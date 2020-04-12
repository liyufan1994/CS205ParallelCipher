
# Parallelize Replica Exchange MCMC Sampling 

## Introduction

### Overview
We designed and implemented an "HPC-version" of Replica Exchange MCMC sampling (Parallel Tempering) based on a hybrid of distributed memory processing and shared-memory processing with C++. In particular, the goal is to simulate an ensemble of Markov chains that exchange states regularly. Each chain is configured to sample a probability distribution (referred to as target distribution) at a certain *temperature*. We deploy two levels of parallelism: (i) each MPI process is used to simulate one or more chains in the ensemble where they exchange states through MPI point-to-point or collective communication; (ii) OpenMP threads are used to parallelize simulation of each step of each individual chain. We apply the resulting algorithm to two practical high-dimensional applications that are especially compute-intensive: (i) simulating two-dimensional square-lattice Ising model of interacting magnetic spins; (ii) decipher encrypted text. The experiments are run on Harvard Cannon computing cluster.

### The need for HPC
MCMC algorithm only provides an approximation of the correct sampling. One needs to run the underlying Markov chains sufficiently long to obtain accurate simulation. For high-dimensional problems, each step of the chain could take tens of thousands of float point operation and the chain itself often requires thousands of such steps to converge. For instance, simulating Ising model of moderate size (e.g. 64x64 lattice), each step of the Gibbs sampler requires 64x64x7=28672 floating point operations and convergence could take 10000+ steps. Application in cryptography exhibits similar demands for intensive computation. To exacerbate the problem, one single chain rarely offers optimal performance for such high dimensional problems. Practitioners often need to run an ensemble of S **communicating** Markov chains at different temperatures to attain satisfactory mixing at each temperature level. S ranges from 10 to 50 (or more). The computational burden is S-fold of a single Markov chain for serial implementation. 

## Problem description and comparison with existing work of parallelization
In this section, we give a succinct account of two problems we seek to solve with Replica Exchange MCMC Sampling: (i) Simulating Ising lattice; (ii) decipher encrypted texts. 

### Statistical Mechanics: Simulating Ising Lattice
Consider a NxN matrix of entries of +1 or -1. This matrix is called an Ising lattice and each entry is a site of the lattice and +1, -1 denotes the site's spin. A spin configuration is an assignment of spin value to each lattice site.

<img src="doc/image/IsingSquare.png" width="400" height="400" >

A site has four neighbors: sites that is above, below, right, or left to it. For any two adjacent sites i,j, there is an interaction <img src="https://render.githubusercontent.com/render/math?math=J_{ij}">. The energy of a configuration <img src="https://render.githubusercontent.com/render/math?math=\delta">  is given by the Hamiltonian function:

<img src="https://render.githubusercontent.com/render/math?math=H(\delta)=-\sum_{(i,j)} J_{ij} \sigma_i \sigma_j -\mu \sum_j h_j \sigma_j,">
where the first sum is over pairs of adjacent spins. The magnetic moment is given by <img src="https://render.githubusercontent.com/render/math?math=\mu."> The configuration probability is given by the Boltzmann distribution with inverse temperature <img src="https://render.githubusercontent.com/render/math?math=\beta \ge 0:">

<img src="https://render.githubusercontent.com/render/math?math=P_\beta(\sigma)=\frac{e^{-\beta H(\sigma)}}{Z_\beta},">

where <img src="https://render.githubusercontent.com/render/math?math=Z_\beta"> is normalizing constant. The goal is now to draw sample from target distribution <img src="https://render.githubusercontent.com/render/math?math=P_\beta"> for a range of temperature <img src="https://render.githubusercontent.com/render/math?math=\beta.">

A common approach is to simulate dynamics of the Ising lattice at these temperatures simultaneously using Parallelize Replica Exchange MCMC Sampling where each chain (a Gibbs sampler) targets one temperature level. These chains exchange state frequently depending on the energy level of the entire ensemble so that the chains do not get stuck at local modes. Typically such simulation is computationally expensive for even moderate sized lattice. We consider running chains in the ensemble in parallel through MPI, and parallelizing each individual Gibbs sampler through strip-partition and chess-board-partition of the lattice square.

### Cryptograph: break a substitution cipher using 2-gram MCMC decipher

We consider the problem of decoding a paragraph of text encrypted through a substitution cipher. For example, we have a paragraph of text:

"*Tonight, we will break out of jail. The time we do this will be midnight. Tony and his men will be ready outside.*"

A substitution cipher works by replacing each character in this string with another character. Say, we replace "T" with "n", "o" with "S", "n" with "z", "i" with "&" and so forth. The coded text of this paragraph may look like "gibberish" such as the following: 

"*nSz&Kv&);DS&AKK;z;SKzMKM))0KzMKC)MMK0MLKk)v)U4amsdfUK!0)gSKDz_SKzx&SzU4KDzUK;DSKlz4KDzUK;DSKlz4KDzUK;Dasd*"

What each character in the original text corresponds to in the ciphered text is the key with which a ciphered paragraph may be recovered. The task of a deciphering algorithm is to find this key. For example, in our experiments, we restrict our attention to the 95 printable ASCII character. A ciphering key would be a permutation of these 95 character. To illustrate, consider 4 instead of 95 characters: "H, e, l, o". A permutation is "e, l, o, H". Then we can cipher "Hello" to "elooH". But once this permutation is uncovered, we can easily obtain the original text by reverse the permutation mapping. **There are a total of 95! ( roughly 10 to the power of 148) possible solutions!** So a naitve approach of trying out all possible keys simply do not work.

A strategy to search for the key is to use MCMC algorithm based on 2-gram frequency counts. In particular, we have a long reference text that is preferrably close to the ciphered text in terms of wording an context. For each pair of characters <img src="https://render.githubusercontent.com/render/math?math=\beta_1, \beta_2">(e.g. "h" and "e"), we let <img src="https://render.githubusercontent.com/render/math?math=r(\beta_1, \beta_2)"> denote the number of times that each specific pair (e.g. "he") appear consecutively in the referenced text. Let's use x to denote a putative decryption key to the ciphered text. We then let <img src="https://render.githubusercontent.com/render/math?math=f_x(\beta_1, \beta_2)"> record the number of times that the pair appears when the cipher text is decrypted using the decryption key x. Consider the following score function

<img src="https://render.githubusercontent.com/render/math?math=\pi(x)=\prod_{\beta_1, \beta_2} r(\beta_1, \beta_2)^{f_x(\beta_1, \beta_2)}">

Informally, <img src="https://render.githubusercontent.com/render/math?math=\pi(x)"> here has the interpretation of the likelihood of the decrypted text is correct English (consider r as "true" transition probability from one character to the next). Intuitively, if a pair of character appears frequently in reference text, such pair should also appear frequently in the decrypted text provided that the decryption is correct. In short, we would like to sample from <img src="https://render.githubusercontent.com/render/math?math=\pi(x)"> to obtain keys that are likely to be correct and do some finer processing of the deciphered text. 

Parallelize Replica Exchange MCMC Sampling is often applied to sample <img src="https://render.githubusercontent.com/render/math?math=\pi(x)^\beta"> at different temperature <img src="https://render.githubusercontent.com/render/math?math=\beta">. Each Markov chain follows a Metropolis Hasting algorithm (See below). Roughly speaking, each state of an Markov chain is a candidate decyphering key. It starts from a randomly generated permutation of the 95 ASCII characters and it proposes to take a step forward by permuting two of the characters. It decides if it will accept this proposal by evaluating <img src="https://render.githubusercontent.com/render/math?math=\pi(x)">. As we can see, evaluating <img src="https://render.githubusercontent.com/render/math?math=\pi(x)"> is computationally expensive: in our case, each step of the Markov chain involves 96\*96\*2=18432 floating point operations. We parallelize such evaluation with OpenMP reduction clause within each chain. Another level of parallelism is that we assign one or more chains to one MPI process. Using multiple chains starting at different candidate key is usually necessary since one single chain will likely be stuck at local mode. 



## Algorithms and Parallelization Strategy
We will first introduce replica exchange MCMC sampling and why a hybrid of distributed memory processing (MPI) and shared-memory processing (OpenMP) will be a particularly suitable solution for parallelization.

### Architectural Overview
The overall archetecture is represented in the following diagram:

<img src="doc/image/Achetecture.png" width="600" height="550">

The overall procedure may be summarized as:
1. Assign S Markov chains to N MPI process where each process is hosted on one computing node; Ideally, each node run one Markov chain;
2. Each Markov chain target a specific temperature level of the target distribution (See details below);
3. Each MPI process run a chain through multiple OpenMP threads, parallelizing computation of target density or Gibbs update sites (See details below)
4. The Markov chains are run independently for K steps before they exchange states (See details below); this is repeated for T iterations

### Parallelize Individual MCMC Chains with OpenMP
The replica exchange MCMC sampling algorithm run S parallel Markov chains denoted as <img src="https://render.githubusercontent.com/render/math?math=X_t^{(1)}, X_t^{(2)},...,X_t^{(S)}"> at time t=1,2,.... Each step of a chain is referred to as a state and each state could be a number (one-dimensional), an n-array, or a n-by-n matrix. 

Each of these chains are set up to simulate a target distribution <img src="https://render.githubusercontent.com/render/math?math=\pi(x)"> at temperatures <img src="https://render.githubusercontent.com/render/math?math=T_1\ge T_2\ge T_3..."> respectively. To be specific, the i-th chain will simulate <img src="https://render.githubusercontent.com/render/math?math=\pi(x)^{T_i}"> by taking steps as *informed by <img src="https://render.githubusercontent.com/render/math?math=\pi(x)^{T_i}">*. There are two ways for each chain to *take a step* (i.e. moving from <img src="https://render.githubusercontent.com/render/math?math=X_i^{(t)}"> to <img src="https://render.githubusercontent.com/render/math?math=X_i^{(t+1)}">):
1. Metropolis-Hasting kernel: the chain proposes a move to <img src="https://render.githubusercontent.com/render/math?math=X_i^{(t+1)}=y"> uniformly around its current position <img src="https://render.githubusercontent.com/render/math?math=X_i^{(t)}=x"> and accepts the move with probability <img src="https://render.githubusercontent.com/render/math?math=\min(1, \frac{\pi(y)^{T_i}}{\pi(x)^{T_i}})">. 
2. Gibbs kernel: in the case where each state of the chain is an array or a matrix, we may update each coordinate of the Markov chain state sequentially by the conditional distribution <img src="https://render.githubusercontent.com/render/math?math=\pi(x_i|x_1,...,x_{i-1},x_{i+1},...)">. 

#### Parallelize Gibbs sampler for Ising Model Example
For the Ising lattice application, each state is an N-by-N matrix. We consider two ways to parallelize the Gibbs updates: (i) partition the data matrix by strip; (ii) a checkboard decomposition.

The Gibbs updates for strip partition proceeds as follows:
1. Partition the entire Ising strip into S stripps horizontally and assign each strip to a parallel OpenMP threads;
2. Each OpenMP threads update its strip by sweeping through it row by row **with the exception of the last row**. Each entry is assigned +1 with probability <img src="https://render.githubusercontent.com/render/math?math=\frac{\exp(\beta s_i)}{\exp(\beta s_i) \dagger \exp(-\beta s_i)}"> (otherwise assign -1) where <img src="https://render.githubusercontent.com/render/math?math=s_i"> denotes sum of four neighboring sites of i;
3. An OpenMP barrier is placed until all threads have reached their last row (but has not updated it);
4. All OpenMP threads update the last rows of their strip.

<img src="doc/image/Stripped.png" width="500" height="400">

The reason for this updating schedule is that each sites is dependent on the current value their four neighbors. So it is not possible to update all sites simultaneously at once. The last row of each strip is then left un-updated until the first row of the next strip has been updated by the other OpenMP thread. 

The checkerboard decomposition uses a more efficient partition scheme. The key observation is that each cell's updates depend only on its four neighboring cells (up, down, left and right). For this reason, all the black sites' update are independent of all the white sites and vice versa:

<img src="doc/image/Checkerboard2.png" width="500" height="200">

The updating schedule for checkerboard decomposition is therefore as follows:
1. Update all white sites in parallel using as many OpenMP threads as possible by the similar conditional probability rule: Each entry is assigned +1 with probability <img src="https://render.githubusercontent.com/render/math?math=\frac{\exp(\beta s_i)}{\exp(\beta s_i) \dagger \exp(-\beta s_i)}"> (otherwise assign -1) where <img src="https://render.githubusercontent.com/render/math?math=s_i"> denotes sum of four neighboring sites of i. Note that here each entry 
2. Then update all black sites while keeping white sites unchanged. 

#### Parallelize Metropolis Hasting kernel for cryptograph example
Recall that the goal is to sample from the target distribution <img src="https://render.githubusercontent.com/render/math?math=\pi(x)"> below:

<img src="https://render.githubusercontent.com/render/math?math=\pi(x)=\prod_{\beta_1, \beta_2} r(\beta_1, \beta_2)^{f_x(\beta_1, \beta_2)}">

To illustrate, we suppose that the state space (i.e. all candidate solutions) is permutation of (1,2,3,4). The actual implementation will deal with (1,2,...,95). 

The Metropolis Hasting algorithm proceeds as follows:
1. Each step the chain proposes to permute two entries in the list; for example, if the current state is (1,3,2,4), we may propose to permute 1,2 to obtain (2,3,1,4) as a proposal.
2. Compute <img src="https://render.githubusercontent.com/render/math?math=\pi(x)"> for both current and proposed move and compute acceptance ratio;
  
 The most computation-intensive part is to calculate <img src="https://render.githubusercontent.com/render/math?math=\pi(x)"> which involves a double loop of size 95x95. The parallelization is relatively easy. We simply apply OpenMP reduction clause to the loop.


### Parallelize Replica Exchange with MPI
In the last section, we discussed how to parallelize each individual chain in the ensemble through OpenMP. Now, we will explain how to run the Replica Exchange Sampling ensemble through MPI. 

First consider the following illustration for Replica Exchange MCMC Sampling:

<img src="doc/image/ParTempering.jpg" width="500" height="330">


The replica exchange MCMC sampling algorithm proceeds as follows:
1. There are S chains in the ensemble. Each chain runs either as Gibbs sampler or Metropolis Hasting sampler and will be parallelized through OpenMP. In particular, each chain targets a specific temperature level. This means that the target distribution for the i-th chain will be <img src="https://render.githubusercontent.com/render/math?math=\pi(x)^{T_i}"> where <img src="https://render.githubusercontent.com/render/math?math=T_i"> is the i-th temperature;
2. Suppose that there are N nodes. Each node will be responsible for one or more chains;
3. The chains proceed independently for K steps then do an *exchange of states*. Such exchange of states is achieved by a proposal-acceptance procedure (i.e. Metropolis Hastings step). In particular, we may propose to exchange states of chain i and chain j. Suppose chain i (with temperature <img src="https://render.githubusercontent.com/render/math?math=T_i">) is at <img src="https://render.githubusercontent.com/render/math?math=x_i"> and let <img src="https://render.githubusercontent.com/render/math?math=E_i=\log(\pi(x_i))">, and analogously for chain j. The acceptance probability for the exchange is given by the following:

<img src="https://render.githubusercontent.com/render/math?math=p = \min \left( 1, \frac{ \exp \left( -\frac{E_j}{kT_i} - \frac{E_i}{kT_j} \right) }{ \exp \left( -\frac{E_i}{kT_i} - \frac{E_j}{kT_j} \right) } \right) = \min \left( 1, e^{(E_i - E_j) \left( \frac{1}{kT_i} - \frac{1}{kT_j} \right)} \right)">

We use different exchange schedule for the Ising model application and the deciphering application. The Ising model example will following what we call a "shifting model". In particular, it proposes to exchange chain i and chain i+1 in a sequential fashion. At each iteration, MPI process running chain i+1 will send its state to the process running chain i and that process will compute <img src="https://render.githubusercontent.com/render/math?math=E_i"> and <img src="https://render.githubusercontent.com/render/math?math=E_j"> and determine whether to accept exchange. It will then send back its decision to the MPI process running chain i+1 along with state of chain i if the decision is to accept.

The deciphering application follows a "star model". In particular, the proposal will be to exchange chain 1 and chain i sequentially. The reason for doing this is that we would like chain 1 to be the solution key at the end. As a result, we would like to constantly check if other chains have obtained a more optimal solution.

## Reproduction Instructions and System Specification
This code may be run on local machine or a cluster. Please check the dependencies below. But it was mainly run on CANNON computing cluster at Harvard FAS Research computing. We will provide instructions and slurm commands. 

### Operating System and distribution
Our testing is done on Linux system (CentOS) but it should also work with Mac OS once one has the required compiler and MPI installed.

```
  Operating System: CentOS Linux 7 (Core)
       CPE OS Name: cpe:/o:centos:centos:7
            Kernel: Linux 3.10.0-957.12.1.el7.x86_64
```

### Compiler and MPI
You must have Intel or GCC compiler that supports OpenMP (the "-fopenmp" tag). For this project, we will use GCC of the following version:
```
Thread model: posix
gcc version 9.2.0 (GCC) 
```

We use Open MPI of the following version for this project:
```
(Open MPI) 4.0.2
```

It is also required that cmake is installed on the system:
```
cmake version 2.8.12.2
```

These softwares may be loaded on Cannon Cluster through the following command:
```
module load gcc/9.2.0-fasrc01 openmpi/4.0.2-fasrc01
```

Or if one prefers to use Intel C compiler (better OpenMP threads affinity support):
```
module load intel/19.0.5-fasrc01 openmpi/4.0.1-fasrc01
```

### Hardware Archetecture
The following information is provided at the Harvard FAS research computing website: 

"The new Cannon cluster is primarily comprised of 670 Lenovo SD650 NeXtScale servers, part of their new liquid-cooled Neptune line. Each chassis unit contains two nodes, each containing two Intel 8268 "Cascade Lake" processors and 192GB RAM per node. The nodes are interconnected by HDR 100 Gbps Infiniband (IB) in a single Fat Tree with a 200 Gbps IB core. "

In particular, the CPU information of each node is given below:

```
lscpu
Architecture:          x86_64
CPU op-mode(s):        32-bit, 64-bit
Byte Order:            Little Endian
CPU(s):                48
On-line CPU(s) list:   0-47
Thread(s) per core:    1
Core(s) per socket:    24
Socket(s):             2
NUMA node(s):          2
Vendor ID:             GenuineIntel
CPU family:            6
Model:                 85
Model name:            Intel(R) Xeon(R) Platinum 8268 CPU @ 2.90GHz
Stepping:              7
CPU MHz:               1199.896
CPU max MHz:           3900.0000
CPU min MHz:           1200.0000
BogoMIPS:              5800.00
Virtualization:        VT-x
L1d cache:             32K
L1i cache:             32K
L2 cache:              1024K
L3 cache:              36608K
NUMA node0 CPU(s):     0-23
NUMA node1 CPU(s):     24-47
```

To request nodes, one can either initiate an interative mode where one can interact with the reuqested resource or submit a job through slurm batch file. Our expriment is conducted in the interative mode where the following line is used to request resources:

```
srun -p test -n 8 -N 8 -c 12 --pty --mem 1000 -t 0-30:00 /bin/bash
```

This corresponds to 8 tasks on 8 nodes where each node is equipped with 12 CPUs (i.e. maximum 12 OpenMP threads per node). 



### The instruction to compile and run the project:

1. Clone the repo from github
2. Find the paths to mpicc and mpic++ compiler on your own machine by typing 

```$ which mpicc```

```$ which mpic++```

3. Go to root directory and modify the following lines in the CMakeList.txt file: replace "/usr/local/bin/mpicc" 
and "/usr/local/bin/mpic++" with the paths you find in step 2.

```set(CMAKE_C_COMPILER /usr/local/bin/mpicc)```

```set(CMAKE_CXX_COMPILER /usr/local/bin/mpic++)```

This step tell the cmake which compiler to use

4. In the root directory, type

```$ cmake .```

```$ make```


5. Go to bin directory, you will find the following two executables; Ising is the executable for Ising model application whereas Denigma is the executable for decryption application

```Ising```

```Denigma```

6. Type the following to run the code (with 4 tasks)

```$ mpirun -np 4 ./Denigma```

```$ mpirun -np 4 ./Ising```

7. One may alter number of OpenMP threads by setting OMP_NUM_THREAD environment variable before mpirun;

```$ export OMP_NUM_THREADS=6```


8. One may adjust OpenMP NUMA thread affinity by setting KMP_AFFINITY environment variable (for Intel compiler) or GOMP_CPU_AFFINITY (for GCC). For example,

```$ export GOMP_CPU_AFFINITY="0-5"```

```$ export KMP_AFFINITY=verbose,compact```



## Performance Evaluation 
We apply the paralleled replica exchange MCMC algorithm to two distinct problems: Ising lattice simulation and decryption. The purpose of testing the algorithm of two problems (as opposed to one) is that these two problems tend to reflect different aspects of the performance of the algorithm. In particular, we use different benchmark for these two applications. For the Ising model, we focus on testing raw computing speed and scalability of the algorithm as we increase the problem size (size of the lattice and number of chains in the ensemble). For the decryption problem, we can explicitly evaluate improvement of decryption accuracy under fixed time budget--with MPI-level parallelization we can run more chains to broaden our search and avoid minima trap whereas with OpenMP-level parallelization we can run longer chains to make the search more thorough. So the benchmark now is the accuracy of the algorithm under a fixed time budget. 

### Ising Model

#### Experiment Set-up and Benchmark
The experiment setup is the following: we have a NxN 2D Ising lattice and we run a replica exchange MCMC algorithm of S parallel Markov chains. Temperature level is set to be equal intervals between 0.3 to 0.6 where so that the temperature of the first chain is 0.6, second chain 0.6-(0.6-0.3)/S, the third chain 0.6-2x(0.6-0.3)/S and so forth. We may scale the problem by changing N and S. We assume Note that the code will assign more than one chain to an

#### Scalability with Fixed Problem Size (Strong Scaling)
In this section, we test algorithm performance with a fixed problem size but increasing number of MPI processes and OpenMP threads. Recall that speedup for a fixed problem size with respect to the number of processors is governed by Amdahl's law. 

#### Scalability with Increasing problem size (Weak Scaling)
In this section, we test performance of the algorithm when the problem size is scaled to the number of processors. Recall that this is governed by Gustafson's law.

#### Compare Checkerboard and Strip Decomposition
In this section, we fix both the problem size and number of parallel processors/threads and compare the performance using strip and checkboard decomposition.


### Break Substitution Cipher

#### Accuracy gain with more MPI processes
In this section, we test the improvement in decryption accuracy under fixed time budget when we increase number of parallel MPI processes. This should result in more parallel Markov chains in the ensemble whereas each chain is run for about the same duration (or shorter due to parallelization overhead). This should lead to more consistent result (success in finding global maximum rather than local minimums).

#### Accuracy gain with more OpenMP threads
In this section, we test the improvement in decryption accuracy under fixed time budget when we increase number of parallel OpenMP threads. This should result in longer chains and a more thorough search.

#### Accuracy gain with more MPI processes and more OpenMP threads
In this section, we test the improvement in decryption accuracy under fixed time budget when we increase both the number of parallel MPI processes and OpenMP threads.
