# Overview

This archive contains all data required to reproduce the results in our paper. It also contains detailed information about the systems we performed our evaluation on as well as the used compilers and libraries. 

This archive is structured into three subdirectories:

* cp2k: The git repository of our modified CP2K code. The required commits are given in the individual READMEs in the subdirectories.
* evaluation:
  * software: Software evaluation from Section V of our paper.
  * hardware_accel: Evaluation of GPU and FPGA acceleration from Section VI of our paper.
* exploration: 
  * dim_sparsities: Scripts and tools to evaluate the dimensions of the submatrices and block-based and element-based sparsities discussed in section V.C and shown in Figures 4 and 11.
  * heuristics: Programs to evaluate the k-means-based heuristic and the METIS-based heuristic discussed in section IV.C and shown in Figure 5.

The READMEs in the subdirectories give extensive information on the corresponding topics.
