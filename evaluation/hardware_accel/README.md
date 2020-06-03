# Introduction

This folder contains the information sources to reproduce the results presented in the chapter "Hardware Acceleration" in the paper.

To make the tests more reproducible, instead of performing a complete calculation in CP2K with the submatrix method, we have exported the matrix for which the submatrix is to be constructed from CP2K into a directory based structure where each file contains an individual block. Due to the different availability we had to use two cluster system of Paderborn Center of Parallel Computing:

## Noctua
### CPU-Nodes
* CPU-nodes with dual-socket Intel Xeon Gold 6148 (2x20 cores, HT disabled, 6-channel DDR4-2666)
* 192 GB memory per node
* Lustre file system
* CentOS 7.5.1804, kernel 3.10.0-693.2.2.el7.x86_64
* 100 Gbit/s Omni-Path interconnect in fat-tree-topology
* These nodes have been used for the CPU-based evaluation of the submatrix method and the export of the matrices for GPU- and FPGA acceleration.

The output of the `collect_environment.sh` script can be found in `cpu_node.txt`.

### FPGA-Nodes
* same as CPU-Nodes except one CPU is a Intel Xeon Gold 6148F (integrated Omni-Path HFI)
* 2x Bittware 520N (Intel Stratix 10 GX2800, 32 GB quad-channel DDR4, PCIe-3.0, electrically x16, working only x8)
* These nodes have been used for the FPGA-based acceleration. Only one FPGA was used.
* CentOS 7.4.1708, kernel 3.10.0-693.2.2.el7.x86_64
* Intel FPGA SDK for OpenCL 19.2, Bittware 520N board support package

The output of the `collect_environment.sh` script can be found in `fpga/fpga_node.txt`.

## OCuLUS
### GPU-Nodes
* GPU-nodes with dual-socket Intel Xeon E5-2670 (2x8 cores, HT disabled, 4-channel DDR3)
* 64 GB memory per node
* one Nvidia RTX 2080 Ti (ZOTAC GAMING GeForce RTX 2080 Ti Blower, 11 GB GDDR6, PCIe-3.0 x16)
* Scientific Linux release 7.2, kernel 3.10.0-693.17.1.el7.x86_64
* CUDA 10.2

The output of the `collect_environment.sh` script can be found in `gpu/gpu_node.txt`.

An additional reason for exporting the matrix from CP2K instead of running a complete calculation in CP2K is that the main memory of the GPU-nodes in Oculus was not sufficient for some of the tests.

# Exporting the Matrices

Thus, the steps are as follows:

1. Compile CP2K from https://github.com/pc2/cp2k/tree/fc11aca067bdf33a249456b22cc6d3f81ad60235 using the architecture file arch/Linux-x86-64-intel.psmp
2. Export the sparse matrix from CP2K.
3. Generate the submatrix from the exported sparse matrix for use in the GPU-based and FPGA-based code.

## Details for Step 1.
Instructions on how to compile CP2K can be found at https://github.com/pc2/cp2k/blob/sc_eval/INSTALL.md.
The export of the matrices has been performed on a CPU-node of Noctua. Thus, we have used CP2K only on the CPU-Nodes of Noctua.

## Details for Step 2.
* In the directory matrix_extraction:
  * Execute job.sh after replacing the path to the compiled cp2k.psmp.
  * This will result in a file "matrix" which holds the orthogonalized Kohn-Sham matrix in a DBCSR-specific binary format.
  * CP2K will also create a directory called "mat" with many subdirectories that contain the individual blocks of the orthogonalized Kohn-Sham matrix.

## Details for Step 3.
* In the directory matrix_extraction:
  * Run the script sm.sh as "bash sm.sh ."
  * The script sm.sh compiles the Fortran program in sm.f90 (with GCC 8.3.0) and constructs the submatrix of the first 32 water molecules. This submatrix is written to ksmat.bin in a binary format and is read by the GPU and FPGA implementations.

# Using the GPU Tensor Core Implementation

* In the directory gpu:
  * Run make_test_interface.sh which will build individual binaries for the four precisions half, halfsingle, single and double. That is, create the files test_interface_half.x, test_interface_halfsingle.x, test_interface_single.x and test_interface_double.x.
  * Copy the created ksmat.bin to the directory where job.sh resides.
  * Run job.sh on a GPU node. This will create the output files 3972_192_[half|halfsingle|single|double]_0.0_[0-20].log.
  * Run "bash analyse.sh". This will create energy.eps (Figure 9) and energy2.eps (Figure 10) but without the FPGA results.

## Measurement of Power Consumption

The power consumption was during calculations was read with the tool "nvidia-smi".

# Using the FPGA-accelerated Implementation
## FPGA Bitstream Synthesis

As discussed in the paper, the matrix-matrix multiplication kernel is based on the OpenCL-based example included in the Intel FPGA SDK for OpenCL 19.2. Our modified version that reaches up to 3.4 TFLOP/s can be found in fpga/synthesis/device/matrix_mult.cl. It was synthesized for the 19.2_hpc BSP for the Bittware 520N board.
The synthesis can be executed, provided the corresponding SDK and BSP in installed with the script fpga/synthesis/job.sh.
The resulting bitstream can be found in fpga/bitstream/matrix_mult.aocx. Please note, that the bitstream was too large to be included in the Github-Repository and is instead hosted on Github's Git Large File Storage (LFS, https://git-lfs.github.com/).

## FPGA-accelerated Implementation

* In the directory fpga:
  * Copy the created bitstream into the directory "bin".
  * Make sure that the common-directory of the examples_aoc-directory supplied with the Intel FPGA SDK for OpenCL 19.2 is found at "../common/" or change the path in build.sh.
  * Run "bash build.sh". This will compile the host-code with the Intel C++ and Fortran compiler and place the binary in bin/hostf.
  * On an FPGA-node run "bash run.sh" to produce the output files 3972_192_single_0.0_[0-20].log.
  * Run "bash analyse.sh" to create energy.eps (Figure 9) and energy2.eps (Figure 10) including the FPGA-results.

## Measurement of Power Consumption

The power consumption was during calculations was measured with a tool "card monitor utility" from Nallatech/Bittware that reads the current and voltage on the PCI-Express bus and the additional 12V power connectors from the FPGA board.
