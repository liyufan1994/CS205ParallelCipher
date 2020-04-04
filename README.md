# CS205ParallelCipher
An MCMC cipher based on parallel tempering

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

