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



class BalacedServer {
    V3D VAL;
    // VAL[remaining_decisions][current imbalance][current group with imbalance favor][reward]
    // remaining_decisions: 0 ... T
    // current imbalance: 0 ... T
    // current group: {GROUP_A, GROUP_B}
    

    V4D VAL_MID;
    // VAL[remaining_decisions][current imbalance][current group with imbalance favor]
    // remaining_decisions: 0 ... T
    // current imbalance: 0 ... T
    // current group: {GROUP_A, GROUP_B}
    // reward: 0 .. X*T

    VD VAL_NAIVE;
    // Val[m] is the value with story m
    // m is an integer in base 4X+1 of up to T digits
    VD get_threshold(int t, int b, int m, int s, int w);
    
public:
    BalacedServer(bool compute_optimized, bool compute_naive, bool compute_mid);
    int T; // horizon length
    int X; // max price, prices are [1, 2, ..., X]
    int B; // budget
    int defaultVal; //default value to start the VAL table
    float pre_balanced_threshold = 0; // threshold for pre-balanced policy

    V2D Prob; // Prob[G][i] prob that client A gives a price i

    
    float Val(int t, int b, int g);
    float ValMid(int t, int b, int g, int reward);
    float ValNaive(int M);
    void printInputs();
    void prinValNaiveTable();
    V2D make_one_simulation(float balanced_threshold);
    void save_val_to_file(string filename);
    void load_val_from_file(string filename);

};

BalacedServer::BalacedServer(bool compute_optimized, bool compute_naive, bool compute_mid) {
    cin >> T >> X >> B;
    defaultVal = -1;
    Prob = V2D(2, VD(X, 0));

    for (int group = 0; group < 2; ++ group) {
        for (int i =0; i < X; ++i) {
            cin >> Prob[group][i];     
        }
    }

    VAL = V3D(T+1, V2D(T+1, VD(2, defaultVal)));
    long long int M = 1;
    for (int i =0; i < T;++i) M = M*(2*X*X+1);
    // cout << "M:       " << M << endl;
    // cout << "maxsize: " << VAL_NAIVE.max_size() << endl;
    if (compute_naive) VAL_NAIVE = VD(M+1, defaultVal);
    if (compute_mid) VAL_MID = V4D(T+1, V3D(T+1, V2D(2, VD(X*T+1, defaultVal))));
    cin >> pre_balanced_threshold;
}


void BalacedServer::printInputs() {
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

void BalacedServer::save_val_to_file(string filename) {
    // Serialize
    cout << "Saving policy..." << endl;
    std::ofstream file(filename);
    if (file.is_open()) {
        for (int i1 =0; i1 < VAL.size(); ++i1) {
            for (int i2=0; i2 < VAL[0].size(); ++i2) {    
                for (int i3=0; i3 < VAL[0][0].size(); ++i3) {   
                    if (VAL[i1][i2][i3] != -1) {
                        file << i1 << " " << i2 << " " << i3;
                        file << " " << VAL[i1][i2][i3] << "\n";
                    }
                }
            }
        }
        file.close();
    } else {
        std::cerr << "Unable to open file for writing." << std::endl;
    }
}

void BalacedServer::load_val_from_file(string filename) {
    cout << "Loading policy..." << endl;
    std::ifstream file(filename);
    int i1, i2, i3;
    if (file.is_open()) {
        while (file >> i1 >> i2 >> i3 ) {
            file >> VAL[i1][i2][i3];
        }
        file.close();
    } else {
        std::cerr << "Unable to open file for reading." << std::endl;
    }
}

float BalacedServer::Val(int t, int b, int g) {
    // cout << "hola1" << endl;
    // cout << t << ", " << b << ", " << m << ", " << s << ", " << w << endl;
    float& res = VAL[t][b][g];
    // cout << "hola2" << endl;
    if (res != defaultVal) return res;
    
    if (t == 0) {
        if (b <= B) return res = 0;
        else return res = -INF;
    }


    res = 0;

    for (int i = 0; i < X; ++i) {
        for (int j = 0; j < X; ++j) {
            float acceptA = 0;
            float acceptB = 0;
            if (b == 0) {
                acceptA = i + Val(t-1, 1, GROUP_A);
                acceptB = j + Val(t-1, 1, GROUP_B);
            }
            else {
                if (g == GROUP_A) {
                    acceptA = i + Val(t-1, b+1, GROUP_A);
                    acceptB = j + Val(t-1, b-1, GROUP_A);
                }
                else {
                    acceptA = i + Val(t-1, b-1, GROUP_B);
                    acceptB = j + Val(t-1, b+1, GROUP_B);
                }
            }
            res += Prob[GROUP_A][i]*Prob[GROUP_B][j]*max(acceptA, acceptB);
        }
    }
    return res;
}


float BalacedServer::ValMid(int t, int b, int g, int reward) {
    // cout << "hola1" << endl;
    // cout << t << ", " << b << ", " << m << ", " << s << ", " << w << endl;
    float& res = VAL_MID[t][b][g][reward];
    // cout << "hola2" << endl;
    if (res != defaultVal) return res;
    
    if (t == 0) {
        if (b <= B) return res = reward;
        else return res = -INF;
    }


    res = 0;

    for (int i = 0; i < X; ++i) {
        for (int j = 0; j < X; ++j) {
            float acceptA = 0;
            float acceptB = 0;
            if (b == 0) {
                acceptA = ValMid(t-1, 1, GROUP_A, reward + i);
                acceptB = ValMid(t-1, 1, GROUP_B, reward + j);
            }
            else {
                if (g == GROUP_A) {
                    acceptA = ValMid(t-1, b+1, GROUP_A, reward + i);
                    acceptB = ValMid(t-1, b-1, GROUP_A, reward + j);
                }
                else {
                    acceptA = ValMid(t-1, b-1, GROUP_B, reward + i);
                    acceptB = ValMid(t-1, b+1, GROUP_B, reward + j);
                }
            }
            res += Prob[GROUP_A][i]*Prob[GROUP_B][j]*max(acceptA, acceptB);
        }
    }
    return res;
}



float BalacedServer::ValNaive(int M) {
    float& res = VAL_NAIVE[M];
    if (res != defaultVal) return res;
    int t = T-1;
    if (M != 0) {
        t = T-2 - log_a_to_base_b(M, 4*X+1); 
    }


    if (t == -1) { // this means that the number code corresponds to  a full word of size T, the compute payoff and constraint
        VI word = convertToBaseBfixedLength(M, 2*X*X+1, T);
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
        int reward = 0;
        int imbalance = 0;

        for (int i = 0; i < T; ++i) {
            if (Ds[i] == GROUP_A) {
                ++imbalance;
                reward += As[i];
            }
            else {
                --imbalance;
                reward += Bs[i];
            }
        }
        if (abs(imbalance) > B) return res = -INF;
        else return res = reward;
    }

    
    
    res = 0;
    
    for (int price_a = 0; price_a < X; ++price_a) {
        for (int price_b = 0; price_b < X; ++price_b) {
            float takeA = ValNaive(M*(2*X*X+1) + 2*X*price_a + 2*price_b + GROUP_A + 1);
            float takeB = ValNaive(M*(2*X*X+1) + 2*X*price_a + 2*price_b + GROUP_B + 1);
            res += Prob[GROUP_A][price_a]*Prob[GROUP_B][price_b]*max(takeA, takeB);
        }
    }
    
    return res;
}




VD BalacedServer::get_threshold(int t, int b, int m, int s, int w) {

    float thresholdA = 0;
    float thresholdB = 0;
    return {thresholdA, thresholdB};
}

V2D BalacedServer::make_one_simulation(float balanced_threshold) {

    V2D result(T, VD(9, 0)); //initialize
    
    int t = 0; //stage

    int b_greedy = 0; // current imbalance greedy
    int b_balanced = 0; // current imbalance pre-balanced
    int b_optimal = 0; // current imbalance optimal

    int g_greedy = 0; // current imbalance group greedy
    int g_balanced = 0; // current imbalance group pre-balanced
    int g_optimal = 0; // current imbalance group optimal

    int reward_greedy = 0;
    int reward_balanced = 0;
    int reward_optimal = 0;

    int cost_optimal = 0;

    for (t = 0; t < T; ++t) {
        result[t][0] = b_greedy;
        result[t][1] = g_greedy;
        result[t][2] = reward_greedy;
        result[t][3] = b_balanced;
        result[t][4] = g_balanced;
        result[t][5] = reward_balanced;
        result[t][6] = b_optimal;
        result[t][7] = g_optimal;
        result[t][8] = reward_optimal;
        
        int price_A = sampleNumber(Prob[GROUP_A]);
        int price_B = sampleNumber(Prob[GROUP_B]);
        // greedy part
        bool takeA = (price_A > price_B) or ((price_A == price_B) and (g_greedy == GROUP_B) );

        if (takeA) {
            reward_greedy += price_A;
            if (g_greedy == GROUP_A) {
                ++b_greedy;
            }
            else if (b_greedy == 0) {
                g_greedy = GROUP_A;
                ++b_greedy;
            }
            else --b_greedy;
        }
        else {
            reward_greedy += price_B;
            if (g_greedy == GROUP_B) {
                ++b_greedy;
            }
            else if (b_greedy == 0) {
                g_greedy = GROUP_B;
                ++b_greedy;
            }
            else --b_greedy;
        }


        // pre-balanced part
        takeA = (price_A - price_B < -balanced_threshold) or ((price_A - price_B == -balanced_threshold) and (g_balanced == GROUP_B) );

        if (takeA) {
            reward_balanced += price_A;
            if (g_balanced == GROUP_A) {
                ++b_balanced;
            }
            else if (b_balanced == 0) {
                g_balanced = GROUP_A;
                ++b_balanced;
            }
            else --b_balanced;
        }
        else {
            reward_balanced += price_B;
            if (g_balanced == GROUP_B) {
                ++b_balanced;
            }
            else if (b_balanced == 0) {
                g_balanced = GROUP_B;
                ++b_balanced;
            }
            else --b_balanced;
        }

        // optimal part
        float acceptA = price_A;
        float acceptB = price_B;

        if (b_optimal > 0) {
            if (g_optimal == GROUP_A) {
                acceptA += Val(T-t-1, b_optimal + 1, GROUP_A);
                acceptB += Val(T-t-1, b_optimal - 1, GROUP_A); 
            }
            else {
                acceptA += Val(T-t-1, b_optimal - 1, GROUP_B);
                acceptB += Val(T-t-1, b_optimal + 1, GROUP_B);
            }
        }
        else {
            acceptA += Val(T-t-1, 1, GROUP_A);
            acceptB += Val(T-t-1, 1, GROUP_B);
        }
        
        takeA = (acceptA > acceptB) or ((acceptA == acceptB) and (g_optimal == GROUP_B));

        if (takeA) {
            reward_optimal += price_A;
            if (g_optimal == GROUP_A) {
                ++b_optimal;
            }
            else if (b_optimal == 0) {
                g_optimal = GROUP_A;
                ++b_optimal;
            }
            else --b_optimal;
        }
        else {
            reward_optimal += price_B;
            if (g_optimal == GROUP_B) {
                ++b_optimal;
            }
            else if (b_optimal == 0) {
                g_optimal = GROUP_B;
                ++b_optimal;
            }
            else --b_optimal;
        }
        
    }

    return result;
    
}



int main(int argc,char **argv) {
    // std::cout << std::fixed;
    // std::cout << std::setprecision(4);
    // srand(10);
    std::mt19937 generator (10);

    bool compute_optimized = false;
    bool compute_naive = false;
    bool compute_mid = false;
    bool save_policy = false;
    bool load_policy = false;
    bool simulate = false;
    int n_simulations = 1;
    bool input_rewards = false;
    string saved_policy_file = "default";
    for (int i = 0 ; i < argc; ++i) {
        string argument = argv[i];
        if (argument == "--compute_naive") compute_naive = true;
        if (argument == "--compute_mid") compute_mid = true;
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
    BalacedServer BS = BalacedServer(compute_optimized, compute_naive, compute_mid);
    // BS.printInputs();
    cout << "T: " << BS.T << endl;
    string filename = "saved_policies/" + to_string(BS.T) + "_" + to_string(BS.X) + "_" + to_string(BS.B) + "_enforcer.txt";
    if (saved_policy_file != "default") filename = saved_policy_file;
    if (load_policy) BS.load_val_from_file(filename);
    // BS.load_val_from_file(filename);
    // BS.printCostsAndRewards();
    // return 0;

    if (compute_optimized) {
        auto t_start = std::chrono::high_resolution_clock::now();
        float aux = BS.Val(BS.T, 0, 0);
        cout << "Value: " << aux << endl;
        auto t_end = std::chrono::high_resolution_clock::now();
        float elapsed_time_ms = std::chrono::duration<float, std::milli>(t_end-t_start).count();
        cout << "Optimized Elapsed time: " << elapsed_time_ms/1000 << " seconds" << endl;
        PrintMemoryInfo("Optimized");
        if (save_policy) BS.save_val_to_file(filename);
    }
    if (compute_mid) {
        auto t_start = std::chrono::high_resolution_clock::now();
        float aux = BS.ValMid(BS.T, 0, 0, 0);
        cout << "Value: " << aux << endl;
        auto t_end = std::chrono::high_resolution_clock::now();
        float elapsed_time_ms = std::chrono::duration<float, std::milli>(t_end-t_start).count();
        cout << "Mid Elapsed time: " << elapsed_time_ms/1000 << " seconds" << endl;
        PrintMemoryInfo("Mid");
        if (save_policy) BS.save_val_to_file(filename);
    }
    if (compute_naive) {
        auto t_start = std::chrono::high_resolution_clock::now();
        float auxnaive = BS.ValNaive(0); 
        cout << "Value naive: " << auxnaive << endl;
        auto t_end = std::chrono::high_resolution_clock::now();
        float elapsed_time_ms = std::chrono::duration<float, std::milli>(t_end-t_start).count();
        cout << "Naive Elapsed time: " << elapsed_time_ms/1000 << " seconds" << endl;
        PrintMemoryInfo("Naive");
        // BS.prinValNaiveTable();
    }

    if (simulate) {

        string csvstr = "sim_no,t,";
        csvstr += "b_greedy,g_greedy,reward_greedy";
        csvstr += ",b_balanced,g_balanced,reward_balanced";
        csvstr += ",b_optimal,g_optimal,reward_optimal";

        for (int sim_no = 0; sim_no < n_simulations; ++sim_no) {
            float thersholdA = 0;
            float thresholdB =0;
            V2D res = BS.make_one_simulation(BS.pre_balanced_threshold);
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