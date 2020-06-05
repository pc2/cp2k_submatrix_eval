# Overview

This archive contains all data required to reproduce the results in our paper. It also contains detailed information about the systems we performed our evaluation on as well as the used compilers and libraries.

You will also need our modified CP2K code. It can be found in an archived git repository with doi: 10.5281/zenodo.3878669

This archive is structured into two subdirectories:

* evaluation:
  * software: Software evaluation from Section V of our paper (Figures 6-10).
  * hardware_accel: Evaluation of GPU and FPGA acceleration from Section VI of our paper (Figures 12, 13).
* exploration:
  * dim_sparsities: Scripts and tools to evaluate the dimensions of the submatrices and block-based and element-based sparsities discussed in section V.C and shown in Figures 4 and 11.
  * heuristics: Programs to evaluate the k-means-based heuristic and the METIS-based heuristic discussed in section IV.C and shown in Figure 5.
  * truncation_thresholds: Scripts and input files used to obtain Figure 1.

The READMEs in the subdirectories give extensive information on the corresponding topics as well as runtime environments and machine specs.
