# Introduction

This folder contains the required input data and the used scripts to reproduce the software evaluation for Section V of our paper. It also contains the outputs that we have obtained in our evaluation and the scripts used to plot the acquired data.

# Used hardware

The software evaluation has been performed on the Noctua system of the Paderborn Center for Parallel Computing. It has the following specifications:
* CPU-nodes with dual-socket Intel Xeon Gold 6148 (2x20 cores, HT disabled, 6-channel DDR4-2666)
* 192 GB memory per node
* Lustre file system
* CentOS 7.5.1804, kernel 3.10.0-693.2.2.el7.x86_64
* 100 Gbit/s Omni-Path interconnect in fat-tree-topology

The output of the `collect_environment.sh` script can be found in `cpu_node.txt`.

# Evaluated code

The evaluation has been performed using a customized version of CP2K, including experimental features of the submatrix method that have not been included upstream yet. It corresponds to revision c7af285016b4d006bc0ff88bf8720f9b5529cb31 in the provided cp2k repository.

# Used compilers and libraries

We used the Intel C/C++/Fortran compiler in version 19.0.4, including Intel MKL and Intel MPI. Furthermore, we linked CP2K against the following libraries, which as well have been compiled using this environment:

* elpa 2019.05.001
* libint 2.5.0_lmax_6
* libxc 4.3.4
* libxsmm 4c3a692

# Compilation

The ARCH file that we used to compile our modified version of CP2K is included in the code repository under `arch/Linux-x86-64-intel.psmp`. If the environment is setup accordingly (see `env.sh`), compilation can be started by running

```
make -j ARCH="Linux-x86-64-intel" VERSION="psmp"
```

# Overview of folders, each with its own README

* epsilon_filter-sensitivity: inputs, scripts and outputs used to produce Figures 6 and 7
* linear-scaling: inputs, scripts and outputs used to produce Figure 8
* strong-scaling: inputs, scripts and outputs used to produce Figure 9
* weak-scaling: inputs, scripts and outputs used to produce Figure 10

Each of these folders contain their own READMEs.
