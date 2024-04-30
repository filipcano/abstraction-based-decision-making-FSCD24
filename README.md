# Abstraction-Based Decision Making for Statistical Properties

Authors: Filip Cano, Thomas A. Henzinger, Bettina Könighofer, Konstantin Kueffner, and Kaushik Mallik

This folder contains the implementation and collected data to reproduce the experiments in the paper Abstraction-Based Decision Making for Statistical Properties presented at the 9th International Conference on Formal Structures for Computation and Deduction.

Experiments run through three examples, as explained in the paper: (i) balanced server, (ii) maximally responsive server, and (iii) clientele-aware server.

## Setup requirements
The scripts should run on any Linux or MacOS machine. In the implementation process, we have used:
- Ubuntu 22.04 as operating system,
- g++ 11.4 as C++ compiler
- python 3.10

## Contents of the folder
- *collected_data/* contains the data used to produce the graphs shown in the paper.
- *inputs/* contains the inputs to the C++ implementation of each model, used to produce the data in the *collected_data/* directory.
- *result_images/* contains the images produced by the notebooks *graphs_exectime.ipynb* and *graphs_simulation.ipynb*. These images correspond to the images in Figures 2 and 3 in the paper.
- *balanced_server.cc* contains the prototype implementation of the balanced server. Builds the policy and simulates it over a period of time. 
- *balanced_server_performance.sh*. This script executes the balanced server example with different inputs to produce the data used to build Figure 3a.
- *balanced_server_resources.sh*. This script executes the balanced server example with different inputs to produce the data used to build Figure 2a.
- *dist_change.cc* contains the prototype implementation of the clientele-aware server. Builds the policy and simulates it over a period of time. 
- *dist_change_performance.sh*. This script executes the clientele-aware server example with different inputs to produce the data used to build Figure 3c.
- *dist_change_resources.sh*. This script executes the clientele-aware server example with different inputs to produce the data used to build Figure 2c.
- *Execution time experiment input generator.ipynb*. Jupyter notebook with scripts to generate input for the different implementations from desired probability distributions and problem parameters.
- *graphs_exectime.ipynb*. Jupyter notebook used to generate images in Figure 2. Requires data found in the *collected_data/* folder.
- *graphs_simulation.ipynb*. Jupyter notebook used to generate images in Figure 3. Requires data found in the *collected_data/* folder.
- *max_response.cc* contains the prototype implementation of the maximally responsive server. Builds the policy and simulates it over a period of time. 
- *max_response_performance.sh*. This script executes the maximally responsive server example with different inputs to produce the data used to build Figure 3b.
- *max_response_resources.sh*. This script executes the maximally responsive server example with different inputs to produce the data used to build Figure 2b.
- *README.md*. This readme file.

## How to reproduce the experiments

### For the experiments in Figure 2:

```
./balanced_server_resources.sh
./dist_change_resources.sh
./max_response_resources.sh
```

This will print the data in the terminal, similar to what can be found in *collected_data/*. To get the figures, save the output of each run to a data file, and run the notebook *graphs_exectime.ipynb*, reading from those data files. 

Warning: Each of these calls will produce processes that may run for up to an hour and take up to 200GB of memory in your system. The whole set of experiments can be done in about 20h in an adequate machine.


### For the experiments in Figure 3:

```
./balanced_server_performance.sh
./dist_change_performance.sh
./max_response_performance.sh
```

This will print the data to concrete files, following the same naming patter as in the *collected_data/* directory. To get the figures, run the notebook *graphs_simulations.ipynb*, making sure to read the data files matching the ones produced in the previous step.

Warning: Since these experiments are not about resource usage, they are much less intense. The whole set of experiments should be done in about 20min, and use less than 15GB or memory.


## How to use the implementation beyond the results in this paper

Each of the three prototype implementations is fairly well documented and can be used to produce and simulate policies for different inputs. In the *inputs/* folder one can find several input examples for the three implementations, as well as three "explainer" inputs policies, where the expected sizes and conditions on the input are documented. These are marked with *explained* in the title. 
