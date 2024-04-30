#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<chrono>
#include<cmath>
#include<cstdlib>
#include<algorithm>
#include<cassert>
#include<random>
// #include "boost/multi_array.hpp"
using namespace std;

const int INF = 1e8;
const int ACCEPT = 1;
const int REJECT = 0;

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

class ValueMaximizer {
    V3D VAL;
    // VAL[remaining_decisions][remaining_budget][appeared]
    // remaining_decisions: 0 ... T
    // remaining_budget: 0 ... B
    // appeared: number of clients that have appeared so far: 0 ... T
    // accepted is not needed (how many clients have been accepted so far), given by t and b

    V4D VAL_MID;
    // VAL_MID[remaining_decisions][remaining_budget][appeared][reward]
    // remaining_decisions: 0 ... T
    // remaining_budget: 0 ... B
    // appeared: number of clients that have appeared so far: 0 ... T
    // Reward observed until now
    // accepted is not needed (how many clients have been accepted so far), given by t and b


    VD VAL_NAIVE;
    // Val[m] is the value with story m
    // m is an integer in base ...

    float get_threshold(int t, int b, int app);
    
public:
    ValueMaximizer(bool compute_optimized, bool compute_naive, bool compute_mid);
    int T; // horizon length
    int X;
    VD N; // sequence of loan values
    int B; // budget
    int defaultVal; //default value to start the VAL table
    VD Prob; //Prob[i] is the probability that a client appears with value i, for i=0...N-1

    
    float Val(int t, int b, int app);
    float ValMid(int t, int b, int app, int reward);
    float ValNaive(int m);
    void printInputs();
    void prinValNaiveTable();
    V2D make_one_simulation(bool use_threshold_policy, float threshold);
    void save_val_to_file(string filename);
    void load_val_from_file(string filename);
};

void ValueMaximizer::printInputs() {
  cout << "T: " << T << ", X: " << X << ", B: " << B << "\nN: [";
  for (int i = 0; i < N.size(); ++i) cout << N[i] << ", ";
  cout << "]\nProb : [";
  for (int i = 0; i < N.size(); ++i) cout << Prob[i] << ", ";
  cout << "]\n";
}

ValueMaximizer::ValueMaximizer(bool compute_optimized, bool compute_naive, bool compute_mid) {
    cin >> T >> X >> B;
    defaultVal = -1;
    VD N_orig(X, 0);
    VD Prob_orig(X, 0);
    
    for (int i=0; i < X; ++i) cin  >> N_orig[i];
    for (int i=0; i < X; ++i) cin  >> Prob_orig[i];
    if (compute_optimized) {
        VAL = V3D(T+1, V2D(B+1, VD(T+1, defaultVal)));
    }
    int maxval = 0;
    for (int i=0; i < X; ++i) maxval = max(maxval, (int)N_orig[i]);
    
    if (compute_mid) {
        VAL_MID = V4D(T+1, V3D(B+1, V2D(T+1, VD(maxval+1, defaultVal))));
    }
    // make N ordered
    vector<pair<int, int>> pairs(0);
    for (int i = 0; i < X; ++i) 
        pairs.push_back({N_orig[i], i});
    sort(pairs.begin(), pairs.end());
    VI indices(0);
    N = VD(0);
    for(const auto &p : pairs) {
        N.push_back(p.first);
        indices.push_back(p.second);
    }
    Prob = VD(X, 0);
    for (int i =0; i < X; ++i) Prob[i] = Prob_orig[indices[i]];




    if (compute_naive) {
        long long int M = 1;
        for (int i = 0; i < T; ++i) M = M*(2*(X+1)+1);
        // cout << "M: " << M << endl;
        VAL_NAIVE = VD(M+1, defaultVal);
    }
}


void ValueMaximizer::save_val_to_file(string filename) {
    // Serialize
    cout << "Saving policy..." << endl;
    std::ofstream file(filename);
    if (file.is_open()) {
        for (int i1 =0; i1 < VAL.size(); ++i1) {
            for (int i2=0; i2 < VAL[0].size(); ++i2) {
                for (int i3=0; i3< VAL[0][0].size(); ++i3) {      
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

void ValueMaximizer::load_val_from_file(string filename) {
    cout << "Loading policy..." << endl;
    std::ifstream file(filename);
    int i1, i2, i3;
    if (file.is_open()) {
        while (file >> i1 >> i2 >> i3) {
            file >> VAL[i1][i2][i3];
        }
        file.close();
    } else {
        std::cerr << "Unable to open file for reading." << std::endl;
    }
}



float ValueMaximizer::get_threshold(int t, int b, int app) {
    int min_it = 0;
    float accept = Val(t-1, b-1, app+1);
    float reject = Val(t-1, b, app+1);
    while ((min_it < N.size()) and (reject > accept + N[min_it])) ++min_it;
    if (reject == accept + N[min_it]) return N[min_it];
    if (min_it == N.size()-1) return N[min_it];
    return ((float)N[min_it] + (float)N[min_it +1])/2;
}

float ValueMaximizer::Val(int t, int b, int app) {
    if (b < 0) return -INF;
    float& res = VAL[t][b][app];
    if (res != defaultVal) return res;
    if (t == 0) {
        return res = 0;
    }
    float acceptance_rate = (float)(1+B-b)/(float)(2+app);
    float accept = 0;
    float reject = 0;
    res = (1-acceptance_rate)*Val(t-1, b, app); // with probability 1-acc_rate, no client appears in the next timestep

    for (int i = 0; i < X; ++i) {
        accept = N[i] + Val(t-1, b-1, app+1);
        reject = Val(t-1, b, app+1);
        // cout << "i: " << i << ", accept: " << accept << ", reject: " << reject;
        // cout << ", acc_rate: " << acceptance_rate << ", Prob: " << Prob[i] << endl; 
        res += acceptance_rate*Prob[i]*max(accept, reject);
    }


    // cout << "t: " << t << ", b: " << b << ", app: " << app << ", res: " << res << endl;
    return res;
}


float ValueMaximizer::ValMid(int t, int b, int app, int reward) {
    if (b < 0) return -INF;
    // cout << "before: " << t << ", "<< b << ", "<< app << ", "<< reward << "\n";
    float& res = VAL_MID[t][b][app][reward];
    // cout << "after: " << t << ", "<< b << ", "<< app << ", "<< reward << "\n";
    if (res != defaultVal) return res;
    if (t == 0) {
        return res = 0;
    }
    float acceptance_rate = (float)(1+B-b)/(float)(2+app);
    float accept = 0;
    float reject = 0;
    res = (1-acceptance_rate)*ValMid(t-1, b, app, reward); // with probability 1-acc_rate, no client appears in the next timestep

    for (int i = 0; i < X; ++i) {
        accept = N[i] + ValMid(t-1, b-1, app+1, reward + N[i]);
        reject = ValMid(t-1, b, app+1, reward);
        // cout << "i: " << i << ", accept: " << accept << ", reject: " << reject;
        // cout << ", acc_rate: " << acceptance_rate << ", Prob: " << Prob[i] << endl; 
        res += acceptance_rate*Prob[i]*max(accept, reject);
    }
    return res;
}

void ValueMaximizer::prinValNaiveTable() {
    cout << "Start printing naive table\n";
    for (int m = 0; m < VAL_NAIVE.size(); ++m) {
        VI word = convertToBaseBfixedLength(m, 2*(X+1)+1, T);
        int app, acc, res;
        app = acc = res = 0;        
        int t = T-1;
        if (m != 0) {
            t = T-2 - log_a_to_base_b(m, 2*(X+1)+1); 
        }
        cout << endl << endl;
        cout << "m: " << m << ", t: " << t;
        cout << ", wordbefore: [";
        for (int i =0; i < word.size(); ++i) {
            cout << word[i] << ", ";
        }
        cout << "]";

        for (int i = 0; i < T; ++i) {
            // assert (word[i] != 0);
            if (word[i] != 0) {
                word[i] = word[i]-1;
                if (word[i]/2 > 0) { // if a customer appeared in stage t
                    ++app;
                    if (word[i]%2 == ACCEPT) { // if custormer is accepeted
                        res += N[word[i]/2-1];
                        ++acc;
                    }
                }
            }        
        }
        
        
        cout << ", word: [";
        for (int i =0; i < word.size(); ++i) {
            cout << word[i] << ", ";
        }
        cout << "]";
        cout << ", app: " << app << ", acc: " << acc;
        // cout << ", t: " << t << ", gAacc: " <<  gAacc << ", gBacc: " << gBacc;
        // cout << ", gAseen: " << gAseen << ", gBseen: " << gBseen;
        cout << ", res: " << VAL_NAIVE[m];
        cout << endl;
    }
    cout << "End printing" << endl;
}


float ValueMaximizer::ValNaive(int m) {
    float& res = VAL_NAIVE[m];
    if (res != defaultVal) return res;

    int t = T-1;
    if (m != 0) {
        t = T-2 - log_a_to_base_b(m, 2*(X+1)+1); 
    }

    int acc = 0; int app = 0;
    res = 0;
    VI word = convertToBaseBfixedLength(m, 2*(X+1)+1, T);   


    for (int i = 0; i < T; ++i) {
        // assert (word[i] != 0);
        if (word[i] != 0) {
            word[i] = word[i]-1;
            if (word[i]/2 > 0) { // if a customer appeared in stage t
                ++app;
                if (word[i]%2 == ACCEPT) { // if custormer is accepeted
                    res += N[word[i]/2-1];
                    ++acc;
                }
            }
        }        
    }

    if (t == -1) { // this means that the number m corresponds to  a full word of size T, the compute payoff and constraint        
        if (acc > B) {
            return res = -INF;
        }
        return res;
    }
    
    float acceptance_rate = (float)(1+acc)/(float)(2+app);

    res = (1-acceptance_rate)*ValNaive(m*(2*(X+1)+1) + 1);



    for (int i = 0; i < X; ++i) {
        float accept = ValNaive(m*(2*(X+1)+1) + 2*(i+1) + ACCEPT + 1);
        float reject = ValNaive(m*(2*(X+1)+1) + 2*(i+1) + REJECT + 1);
        res += acceptance_rate*Prob[i]*max(accept, reject);
    }
    return res;
}

V2D ValueMaximizer::make_one_simulation(bool use_threshold_policy, float threshold) {
    V2D result(T, VD(4, 0)); //initialize
    int t = 0; //stage
    int accepted = 0; //number accepted
    int appeared = 0; //number of clients seen
    int reward = 0; // accumulated reward

    for (t = 0; t < T; ++t) {
        float acceptance_rate = (float)(1+accepted)/(float)(2+appeared);
        result[t][0] = acceptance_rate;
        result[t][1] = B-accepted;
        result[t][2] = reward;
        result[t][3] = use_threshold_policy ? threshold : get_threshold(T-t, B-accepted, appeared);
        //to do threshold

        bool client_appears = sampleNumber({acceptance_rate, (float)1.0 - acceptance_rate}) == 0;
        if (client_appears) {
            ++appeared;
            int client_val = N[sampleNumber(Prob)];
            if (use_threshold_policy) {
                if ((client_val > threshold) and (accepted < B)) {
                    ++accepted;
                    reward += client_val;
                }
            }
            else {
                float accept_val = client_val + Val(T-t, B-(accepted+1), appeared);
                float reject_val = Val(T-t, B-accepted, appeared);
                if (accept_val > reject_val) {
                    ++accepted;
                    reward += client_val;
                }
            }
        }
    }
    return result;
}


int main(int argc,char **argv) {
    bool compute_optimized = false;
    bool compute_naive = false;
    bool compute_mid = false; // compte statistical-DP in its non accumulative form
    bool simulate = false;
    int n_simulations = 1;
    float threshold = 0;
    bool save_policy = false;
    bool load_policy = false;
    bool use_threshold_policy = false;
    string saved_policy_file = "default";
    for (int i = 0 ; i < argc; ++i) {
        string argument = argv[i];
        if (argument == "--compute_naive") compute_naive = true;
        if (argument == "--compute_optimized") compute_optimized = true;
        if (argument == "--compute_mid") compute_mid = true;
        if (argument == "--both") compute_naive = compute_optimized = true;
        if (argument == "--simulate") simulate = true; 
        if (argument == "--save_policy") save_policy = true;
        if (argument == "--load_policy") load_policy = true;
        if (argument == "--use_threshold_policy") use_threshold_policy = true;
        string keyword = "--n_simulations=";
        if (argument.find(keyword) != std::string::npos) 
            n_simulations = stoi(argument.substr(keyword.size(),argument.size()-keyword.size()));
        keyword = "--threshold=";
        if (argument.find(keyword) != std::string::npos) 
            threshold = stof(argument.substr(keyword.size(),argument.size()-keyword.size()));
        keyword = "--saved_policy_file=";
        if (argument.find(keyword) != std::string::npos) 
            saved_policy_file = argument.substr(keyword.size(),argument.size()-keyword.size());
    }
    // compute_optimized = true;
    // compute_naive = true;
    ValueMaximizer VM = ValueMaximizer(compute_optimized, compute_naive, compute_mid);
    cout << "T: " << VM.T << endl;
    string filename = "saved_policies/" + to_string(VM.T) + "_" + to_string(VM.X) + "_" + to_string(VM.B) + "_clientele.txt";
    if (saved_policy_file != "default") filename = saved_policy_file;
    if (load_policy) VM.load_val_from_file(filename);
    if (compute_optimized) {
        auto t_start = std::chrono::high_resolution_clock::now();
        float aux = VM.Val(VM.T, VM.B, 0);
        // VM.printInputs();
        cout << "Value: " << aux << endl;
        auto t_end = std::chrono::high_resolution_clock::now();
        float elapsed_time_ms = std::chrono::duration<float, std::milli>(t_end-t_start).count();
        cout << "Optimized Elapsed time: " << elapsed_time_ms/1000 << " seconds" << endl;
        PrintMemoryInfo("Optimized");
        if (save_policy) VM.save_val_to_file(filename);
    }

    if (compute_mid) {
        auto t_start = std::chrono::high_resolution_clock::now();
        float aux = VM.ValMid(VM.T, VM.B, 0, 0);
        // VM.printInputs();
        cout << "Value: " << aux << endl;
        auto t_end = std::chrono::high_resolution_clock::now();
        float elapsed_time_ms = std::chrono::duration<float, std::milli>(t_end-t_start).count();
        cout << "Mid Elapsed time: " << elapsed_time_ms/1000 << " seconds" << endl;
        PrintMemoryInfo("Mid");
    }

    if (compute_naive) {
        auto t_start = std::chrono::high_resolution_clock::now();
        float auxnaive = VM.ValNaive(0); 
        cout << "Value naive: " << auxnaive << endl;
        auto t_end = std::chrono::high_resolution_clock::now();
        float elapsed_time_ms = std::chrono::duration<float, std::milli>(t_end-t_start).count();
        cout << "Naive Elapsed time: " << elapsed_time_ms/1000 << " seconds" << endl;
        PrintMemoryInfo("Naive");
        // VM.prinValNaiveTable();
    }

    if (simulate) {

        string csvstr = "sim_no, t, acceptance_rate, available_budget, accumulated_reward, local_threshold";

        for (int sim_no = 0; sim_no < n_simulations; ++sim_no) {
            V2D res = VM.make_one_simulation(use_threshold_policy, threshold);
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

    return 0;

}