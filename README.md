
# Parallelize Replica Exchange MCMC Sampling 

## Introduction
We designed and implemented an "HPC-version" of Replica Exchange MCMC sampling (Parallel Tempering) based on a hybrid of distributed memory processing and shared-memory processing. In particular, the goal is to simulate an ensemble of Markov chains that exchange states regularly. Each chain is configured to sample a probability distribution (referred to as target distribution) at a certain *temperature*. We deploy two levels of parallelism: (i) each MPI process is used to simulate one or more chains in the ensemble where they exchange states through MPI point-to-point or collective communication; (ii) OpenMP treads are used to simulate next step of each individual chain. We apply the resulting algorithm to two practical applications that require large amount of computation: (i) simulating two-dimensional square-lattice Ising model of interacting magnetic spins; (ii) decipher encrypted text. The experiments are run on Harvard Cannon computing cluster.

## Algorithms and Parallelization Strategy
We will first introduce replica exchange MCMC sampling and why a hybrid of distributed memory processing (MPI) and shared-memory processing (OpenMP) will be a particularly suitable solution for parallelization.

### MCMC Chains
The replica exchange MCMC sampling algorithm run S parallel Markov chains denoted as <img src="https://render.githubusercontent.com/render/math?math=X_t^{(1)}, X_t^{(2)},...,X_t^{(S)}"> at time t=1,2,.... Each step of a chain is referred to as a state and each state could be a number (one-dimensional), an n-array, or a n-by-n matrix. 

Each of these chains are set up to simulate a target distribution <img src="https://render.githubusercontent.com/render/math?math=\pi(x)"> at temperatures <img src="https://render.githubusercontent.com/render/math?math=T_1\ge T_2\ge T_3..."> respectively. To be specific, the i-th chain will simulate <img src="https://render.githubusercontent.com/render/math?math=\pi(x)^{T_i}"> by taking steps as *informed by <img src="https://render.githubusercontent.com/render/math?math=\pi(x)^{T_i}">*. There are two ways for each chain to *take a step* (i.e. moving from <img src="https://render.githubusercontent.com/render/math?math=X_i^{(t)}"> to <img src="https://render.githubusercontent.com/render/math?math=X_i^{(t+1)}">):
1. Metropolis-Hasting kernel: the chain proposes a move to <img src="https://render.githubusercontent.com/render/math?math=X_i^{(t+1)}=y"> uniformly around its current position <img src="https://render.githubusercontent.com/render/math?math=X_i^{(t)}=x"> and accepts with probability <img src="https://render.githubusercontent.com/render/math?math=\min(1, \frac{\pi(y)^{T_i}}{\pi(x)^{T_i}})">. 
2. Gibbs kernel: in the case where each state of the chain is an array or a matrix, we may update each coordinate of the Markov chain state sequentially by the conditional distribution <img src="https://render.githubusercontent.com/render/math?math=\pi(x_i|x_1,...,x_{i-1},x_{i+1},...)">.

### P
This is where we can potentially to speed up 

### Communication of Chains
In last section, we discussed how to 



The replica exchange MCMC sampling algorithm proceeds as follows:
1. 



The instruction to compile and run the project:

1. Clone the repo from github
2. Find the paths to mpicc and mpic++ compiler on your own machine by typing 

`$ which mpicc`

`$ which mpic++`

3. Go to root directory and modify the following lines in the CMakeList.txt file: replace "/usr/local/bin/mpicc" 
and "/usr/local/bin/mpic++" with the paths you find in step 2.

`set(CMAKE_C_COMPILER /usr/local/bin/mpicc)`

`set(CMAKE_CXX_COMPILER /usr/local/bin/mpic++)`

4. In the root directory, type

`$ cmake .`

`$ make`


5. Go to bin directory, type the following to run the code (with 2 tasks)

`$ mpirun -np 2 ./Denigma`

