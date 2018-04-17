# NKGAclust
This repository contains the source code for the NK Hybrid Genetic Algorithm for Clustering. 

Reference:  Tinós, Zhao, Chicano & Whitley (2018). "NK Hybrid Genetic Algorithm for Clustering". IEEE Transactions on Evolutionary Computation.	

Contact: Renato Tinos <rtinos@ffclrp.usp.br>

Running the code: ./NKGAclust <name_of_instance> <K> <load initial solutions from file>

name_of_instance: an example of the dataset format is given in file Compound.dat (Zahns Compound: https://cs.joensuu.fi/sipu/datasets/ )

K: Parameter of the NK clustering validation 2 (NKCV2) criterion. We recommend to use K=3.

load initial solutions from file: 1 if a file with initial solutions is provided and 0 otherwise. If the parameter is 0, initial random solutions are generated. Examples of files for initial solutions are given in: Compound_1.sol (generated using k-means), Compound_2.sol (generated using DBSCAN), and Compound_3.sol (generated using DP clustering algorithm).

Example for running the code with instance Compound: 

make

./NKGAclust Compound 3 0



Observation 1: Class NK clustering is given in nk.h 

- Function double nk::comp_fitness (int *x) : evaluates a solution x using the NK Clustering Validation 2 (NKCV2) criterion. This function can be used independently from the Genetic Algorithm.
	
- Function double nk::px(int *solution_blue, int *solution_red, int *offspring ): generates an offspring by recombination of two parent solutions using partition crossover for clustering 

- Class nk.h uses two classes for the manipulation of graphs. Those classes were developed in other projects and are reused here. They contain general (and common) functions for the manipulation of graphs. They can be replaced by other classes (codes) for the manipulation of graphs. For each call for one of those functions, a comment was inserted in the code.  
			
Observations 2: file global.cpp contains the parameters of the GA (examples: population size and crossover rate).

Observations 3: NKGAclust generates two files. The first file (sol_<name_of_instance>.dat) contains the best solutions in each run. The second file (sol_<name_of_instance>.dat) contains the fitness of best solutions in each run.
	
