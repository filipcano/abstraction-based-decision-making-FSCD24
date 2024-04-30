#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<chrono>
#include<cmath>
#include<cstdlib>
#include<algorithm>
#include<random>
#include<cassert>
#include<cstdlib>
#include<iomanip>
// #include "boost/multi_array.hpp"

using namespace std;

// MACROS //
const int GROUP_A = 0;
const int GROUP_B = 1;
const int ACCEPT = 1;
const int REJECT = 0;
const int INF = 1e8;

typedef vector<int> VI;
typedef vector<float> VD;
typedef vector<VD> V2D;
typedef vector<V2D> V3D;
typedef vector<V3D> V4D;
typedef vector<V4D> V5D;






int sampleNumber(const std::vector<float>& probabilities) {
    // returns an int between 0 and probabilities.size()-1

    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::discrete_distribution<> dist(probabilities.begin(), probabilities.end());

    return dist(gen);
}

int log_a_to_base_b(int a, int b) {
    return log2(a) / log2(b);
}

void PrintMemoryInfo(string extra) {
    std::string line;
    std::ifstream statusFile("/proc/self/status");

    while (getline(statusFile, line)) {
        if (line.substr(0, 6) == "VmSize") {
            std::string memStr = line.substr(7);
            memStr = memStr.substr(0, memStr.find(" kB"));
            int memKb = std::stoi(memStr);
            float memMb = memKb / 1024.0;
            std::cout << extra << " Memory Usage: " << memMb << " MB" << std::endl;
            break;
        }
    }
}

VI convertToBaseBfixedLength(int N, int B, int length) {
    VI digits(length);

    for (int i = 0; i < digits.size(); ++i) {
        digits[i] = N % B;
        N /= B;
    }

    // The digits are stored in reverse order, so we need to reverse them back
    reverse(digits.begin(), digits.end());
    return digits;
}



class MaximallyResponsive {
    V5D VAL;
    // VAL[remaining_decisions][remaining_budget][max streak][current streak][current group]
    // remaining_decisions: 0 ... T
    // remaining_budget: 0 ... B
    // max streak: 0 ... T
    // current streak: 0 ... max streak 
    // current group: {GROUP_A, GROUP_B}

    VD VAL_NAIVE;
    // Val[m] is the value with story m
    // m is an integer in base 4X+1 of up to T digits
    VD get_threshold(int t, int b, int m, int s, int w);
    
public:
    MaximallyResponsive(bool compute_optimized, bool compute_naive);
    int T; // horizon length
    int X; // max price, prices are [1, 2, ..., X]
    int B; // budget
    int defaultVal; //default value to start the VAL table

    V2D Prob; // Prob[G][i] prob that client A gives a price i

    
    float Val(int t, int b, int m, int s, int w);
    float ValNaive(int M);
    void printInputs();
    void prinValNaiveTable();
    V2D make_one_simulation();
    void save_val_to_file(string filename);
    void load_val_from_file(string filename);
    void printCostsAndRewards();
};

MaximallyResponsive::MaximallyResponsive(bool compute_optimized, bool compute_naive) {
    cin >> T >> X >> B;
    defaultVal = -1;
    Prob = V2D(2, VD(X, 0));

    for (int group = 0; group < 2; ++ group) {
        for (int i =0; i < X; ++i) {
            cin >> Prob[group][i];     
        }
    }

    VAL = V5D(T+1, V4D(B+1, V3D(T+1, V2D(T+1, VD(2, defaultVal)))));
    long long int M = 1;
    for (int i =0; i < T;++i) M = M*(2*X*X+1);
    // cout << "M:       " << M << endl;
    // cout << "maxsize: " << VAL_NAIVE.max_size() << endl;
    if (compute_naive) VAL_NAIVE = VD(M+1, defaultVal);
}


void MaximallyResponsive::printInputs() {
    cout << "Printing inputs... \n";
    cout << "T: " << T << ", X: " << X << ", B: " << B <<  endl;
    cout << "Probs Client A: [";
    for (int i = 0; i < X; ++i) cout << Prob[GROUP_A][i] << ", ";
    cout << "]\n";
    cout << "Probs Client B: [";
    for (int i = 0; i < X; ++i) cout << Prob[GROUP_B][i] << ", ";
    cout << "]\n";
    cout << "... finished printing inputs\n";
}

void MaximallyResponsive::save_val_to_file(string filename) {
    // Serialize
    cout << "Saving policy..." << endl;
    std::ofstream file(filename);
    if (file.is_open()) {
        for (int i1 =0; i1 < VAL.size(); ++i1) {
            for (int i2=0; i2 < VAL[0].size(); ++i2) {
                for (int i3=0; i3< VAL[0][0].size(); ++i3) {
                    for (int i4=0; i4< VAL[0][0][0].size(); ++i4) {
                        for (int i5=0; i5< VAL[0][0][0][0].size(); ++i5) {
                            if (VAL[i1][i2][i3][i4][i5] != -1) {
                                file << i1 << " " << i2 << " " << i3 << " " << i4 << " " << i5;
                                file << " " << VAL[i1][i2][i3][i4][i5] << "\n";
                            }
                        }
                    }
                }
            }
        }
        file.close();
    } else {
        std::cerr << "Unable to open file for writing." << std::endl;
    }
}

void MaximallyResponsive::load_val_from_file(string filename) {
    cout << "Loading policy..." << endl;
    std::ifstream file(filename);
    int i1, i2, i3, i4, i5;
    if (file.is_open()) {
        while (file >> i1 >> i2 >> i3 >> i4 >> i5) {
            file >> VAL[i1][i2][i3][i4][i5];
        }
        file.close();
    } else {
        std::cerr << "Unable to open file for reading." << std::endl;
    }
}

float MaximallyResponsive::Val(int t, int b, int m, int s, int w) {
    if (b < 0) return INF;
    if (s + t <= m) return m;
    // cout << "hola1" << endl;
    // cout << t << ", " << b << ", " << m << ", " << s << ", " << w << endl;
    float& res = VAL[t][b][m][s][w];
    // cout << "hola2" << endl;
    if (res != defaultVal) return res;
    
    if (s + t <= m) return res = m;

    res = 0;

    for (int i = 0; i < X; ++i) {
        for (int j = 0; j < X; ++j) {
            float accept;
            float reject;
            if (i >= j) {
                int s_acc = w == GROUP_A ? s+1 : 1;
                accept = Val(t-1, b, max(m, s_acc), s_acc, GROUP_A);
                int s_rej = w == GROUP_B ? s+1 : 1;
                reject = Val(t-1, b-(i-j), max(m, s_rej), s_rej, GROUP_B);
            }
            else {
                int s_acc = w == GROUP_B ? s+1 : 1;
                accept = Val(t-1, b, max(m, s_acc), s_acc, GROUP_B);
                int s_rej = w == GROUP_A ? s+1 : 1;
                reject = Val(t-1, b - (j-i), max(m, s_rej), s_rej, GROUP_A);
            }
            res += Prob[GROUP_A][i]*Prob[GROUP_B][j]*min(accept, reject);
        }
    }
    return res;
}



float MaximallyResponsive::ValNaive(int M) {
    float& res = VAL_NAIVE[M];
    if (res != defaultVal) return res;
    int t, b, m, s, w;
    
    t = T-1;
    if (M != 0) {
        t = T-2 - log_a_to_base_b(M, 4*X+1); 
    }


    if (t == -1) { // this means that the number code corresponds to  a full word of size T, the compute payoff and constraint
        VI word = convertToBaseBfixedLength(M, 2*X*X+1, T);
        b = m = s = w = 0;
        int cost = 0;
        VI As(T, 0); // vector containing values of A
        VI Bs(T, 0); // vector containing values of B
        VD Ds(T, 0); // vector containing decision (GROUP_A or GROUP_B)
        for (int i=0; i < T; ++i) {
            assert(word[i] != 0);
            word[i] = word[i] -1;
            int code = word[i]%4;
            Ds[i] = code%2 == GROUP_A ? GROUP_A : GROUP_B;
            As[i] = (word[i]/2)/X;
            Bs[i] = (word[i]/2)%X;
        }
        m = 0; // max wait time
        s = 0; // current streak
        w = 0; // previous 
        cost = 0; // cost
        for (int i = 0; i < T; ++i) {
            if (Ds[i] == GROUP_A) {
                s = w == GROUP_A ? s+1 : 1;
                w = GROUP_A;
                if (As[i] < Bs[i]) cost += Bs[i] - As[i];
                
            }
            else {
                s = w == GROUP_B ? s+1 : 1;
                w = GROUP_B;
                if (As[i] > Bs[i]) cost += As[i] - Bs[i];
            }
            m = max(m, s);
        }

        if (cost > B) {
            return INF;
        }

        return res = m;
    }

    
    
    res = 0;
    
    for (int price_a = 0; price_a < X; ++price_a) {
        for (int price_b = 0; price_b < X; ++price_b) {
            float takeA = ValNaive(M*(2*X*X+1) + 2*X*price_a + 2*price_b + GROUP_A + 1);
            float takeB = ValNaive(M*(2*X*X+1) + 2*X*price_a + 2*price_b + GROUP_B + 1);
            res += Prob[GROUP_A][price_a]*Prob[GROUP_B][price_b]*min(takeA, takeB);
        }
    }
    
    return res;
}

/*
void MaximallyResponsive::prinValNaiveTable() {
    cout << "Start printing naive table\n";
    for (int m = 0; m < VAL_NAIVE.size(); ++m) {
        VI word = convertToBaseBfixedLength(m, 8*X+1, T);
        int gAseen, gAacc, gBseen, gBacc, cost;        
        int t = word.size();
        if (t == T) {
        gAseen = gAacc = gBseen = gBacc = cost = 0;
        for (int i=0; i < T; ++i) {
            int code = word[i]%8;
            if (code == 0 or code == 4 or code == 1 or code == 5) {
                ++gAseen;
                if (code == 0 or code == 1) ++gAacc; 
            }           
            else {
                ++gBseen;
                if (code == 2 or code == 3) ++gBacc;
            }
            if (code == 1 or code == 3 or code == 4 or code == 6) cost += N[word[i]/8];
        }
        assert(gAseen + gBseen == t);
        float res = abs((float)gAacc/(1+gAseen) - (float)gBacc/(1+gBseen));

        cout << "m: " << m << ", t: " << t << ", gAacc: " <<  gAacc << ", gBacc: " << gBacc;
        cout << ", gAseen: " << gAseen << ", gBseen: " << gBseen << ", res: " << VAL_NAIVE[m] << endl;
        }
    }
    cout << "End printing" << endl;
}

*/



VD MaximallyResponsive::get_threshold(int t, int b, int m, int s, int w) {

    float thresholdA = 0;
    float thresholdB = 0;
    return {thresholdA, thresholdB};
}

V2D MaximallyResponsive::make_one_simulation() {

    V2D result(T, VD(4, 0)); //initialize
    int t = 0; //stage
    int m_greedy = 0; // current max wait time greedy
    int s_greedy = 0; // current streak wait time greedy
    int w_greedy = 0; // previous client greedy
    int reward_greedy = 0;
    int m_optimal = 0;// current max wait time optimal
    int s_optimal = 0; // current streak wait time optimal
    int w_optimal = 0; // previous client optimal
    int reward_optimal = 0;
    int cost_optimal = 0;

    for (t = 0; t < T; ++t) {
        result[t][0] = m_greedy;
        result[t][1] = m_optimal;
        result[t][2] = reward_greedy;
        result[t][3] = reward_optimal;
        int price_A = sampleNumber(Prob[GROUP_A]);
        int price_B = sampleNumber(Prob[GROUP_B]);
        // greedy part
        if (price_A > price_B) { // accept client A
            s_greedy = w_greedy == GROUP_A ? s_greedy+1 : 1;
            w_greedy = GROUP_A;
        }
        else if (price_A < price_B) { // accept client B
            s_greedy = w_greedy == GROUP_B ? s_greedy+1 : 1;
            w_greedy = GROUP_B;
        }
        else { // if they are equal, change client
            w_greedy = w_greedy == GROUP_A ? GROUP_B : GROUP_A;
            s_greedy = 1;
        }
        reward_greedy += w_greedy == GROUP_A ? price_A : price_B;
        m_greedy = max(m_greedy, s_greedy);

        // optimal part
        float acceptA = 0;
        float acceptB = 0;

        if (price_A > price_B) {
            int new_s_A = w_optimal == GROUP_A ? s_optimal + 1 : 1;
            int new_m_A = max(m_optimal, new_s_A);
            acceptA = Val(T-(t)-1, B - cost_optimal, new_m_A, new_s_A, GROUP_A);
            int new_s_B = w_optimal == GROUP_B ? s_optimal + 1 : 1;
            int new_m_B = max(m_optimal, new_s_B);
            acceptB = Val(T-(t)-1, B - cost_optimal - (price_A - price_B), new_m_B, new_s_B, GROUP_B);
            s_optimal = acceptA < acceptB ? new_s_A : new_s_B;
            m_optimal = acceptA < acceptB ? new_m_A : new_m_B;
            cost_optimal += acceptA < acceptB ? 0 : price_A - price_B;
            w_optimal = acceptA < acceptB ? GROUP_A : GROUP_B;
        }
        else {
            int new_s_A = w_optimal == GROUP_A ? s_optimal + 1 : 1;
            int new_m_A = max(m_optimal, new_s_A);
            acceptA = Val(T-(t)-1, B - cost_optimal - (price_B - price_A), new_m_A, new_s_A, GROUP_A);
            int new_s_B = w_optimal == GROUP_B ? s_optimal + 1 : 1;
            int new_m_B = max(m_optimal, new_s_B);
            acceptB = Val(T-(t)-1, B - cost_optimal, new_m_B, new_s_B, GROUP_B);
            s_optimal = acceptA < acceptB ? new_s_A : new_s_B;
            m_optimal = acceptA < acceptB ? new_m_A : new_m_B;
            cost_optimal += acceptA < acceptB ? price_B - price_A : 0;
            w_optimal = acceptA < acceptB ? GROUP_A : GROUP_B;

        }
        reward_optimal += w_optimal == GROUP_A ? price_A : price_B;
    }

    return result;
}




int main(int argc,char **argv) {
    // std::cout << std::fixed;
    // std::cout << std::setprecision(4);
    // srand(10);
    std::mt19937 generator (10);

    bool compute_optimized = false;
    bool compute_naive = false;\
    bool save_policy = false;
    bool load_policy = false;
    bool simulate = false;
    int n_simulations = 1;
    bool input_rewards = false;
    string saved_policy_file = "default";
    for (int i = 0 ; i < argc; ++i) {
        string argument = argv[i];
        if (argument == "--compute_naive") compute_naive = true;
        if (argument == "--compute_optimized") compute_optimized = true;
        if (argument == "--both") compute_naive = compute_optimized = true;
        if (argument == "--save_policy") save_policy = true;
        if (argument == "--load_policy") load_policy = true;
        if (argument == "--simulate") simulate = true; 
        if (argument == "--input_rewards") input_rewards = true; 
        string keyword = "--n_simulations=";
        if (argument.find(keyword) != std::string::npos) 
            n_simulations = stoi(argument.substr(keyword.size(),argument.size()-keyword.size()));
        keyword = "--saved_policy_file=";
        if (argument.find(keyword) != std::string::npos) 
            saved_policy_file = argument.substr(keyword.size(),argument.size()-keyword.size());
    }
    MaximallyResponsive MR = MaximallyResponsive(compute_optimized, compute_naive);
    // MR.printInputs();
    cout << "T: " << MR.T << endl;
    string filename = "saved_policies/" + to_string(MR.T) + "_" + to_string(MR.X) + "_" + to_string(MR.B) + "_enforcer.txt";
    if (saved_policy_file != "default") filename = saved_policy_file;
    if (load_policy) MR.load_val_from_file(filename);
    // MR.load_val_from_file(filename);
    // MR.printCostsAndRewards();
    // return 0;

    if (compute_optimized) {
        auto t_start = std::chrono::high_resolution_clock::now();
        float aux = MR.Val(MR.T, MR.B, 0, 0, 0);
        cout << "Value: " << aux << endl;
        auto t_end = std::chrono::high_resolution_clock::now();
        float elapsed_time_ms = std::chrono::duration<float, std::milli>(t_end-t_start).count();
        cout << "Optimized Elapsed time: " << elapsed_time_ms/1000 << " seconds" << endl;
        PrintMemoryInfo("Optimized");
        if (save_policy) MR.save_val_to_file(filename);
    }
    if (compute_naive) {
        auto t_start = std::chrono::high_resolution_clock::now();
        float auxnaive = MR.ValNaive(0); 
        cout << "Value naive: " << auxnaive << endl;
        auto t_end = std::chrono::high_resolution_clock::now();
        float elapsed_time_ms = std::chrono::duration<float, std::milli>(t_end-t_start).count();
        cout << "Naive Elapsed time: " << elapsed_time_ms/1000 << " seconds" << endl;
        PrintMemoryInfo("Naive");
        // MR.prinValNaiveTable();
    }

    if (simulate) {

        string csvstr = "sim_no,t,";
        csvstr += "m_greedy,m_optimal,cost_greedy,cost_optimal";

        for (int sim_no = 0; sim_no < n_simulations; ++sim_no) {
            float thersholdA = 0;
            float thresholdB =0;
            V2D res = MR.make_one_simulation();
            for (int t = 0; t < res.size(); ++t) {
                string thisline =  "\n" + to_string(sim_no);
                thisline += "," + to_string(t);
                for (int i = 0; i < res[0].size(); ++i)
                    thisline += "," + to_string(res[t][i]);
                csvstr += thisline;
            }
        }

        cout << endl << csvstr << endl;
    }
}