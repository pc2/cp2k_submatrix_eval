# Dimensions and Sparsity of Submatrices

Please note, the basis sets DVZP-MOLOPT-SR-GTH and SZV-MOLOPT-SR-GTH are abbreviated here as DZVP and SZV respectively.

The obtain the results presented in the paper, several steps have to be executed:

1. Production of raw data:
  * The version of CP2K to export the required data is commit 0448736ad07a8f5103531f3d4751ac53ee533978 of the cp2k-repository at the top level of this archive. Please note that this version of CP2K writes out the complete block structure and all blocks of the orthogonalized Kohn-Sham matrix so that later analysis doesn't require additional CP2K runs. Running this version consumes a significant amount of disk space. It will produce the following files/directories in addition to the conventional output:
     * matrix_ssqrtinv_ks_ssqrtinv_binary: dump of the orthogonalized Kohn-Sham Matrix in the DBCSR-specific binary format
     * The directory matrix_ssqrtinv_ks_ssqrtinv that contains a subdirectory for every column of the matrix. Each of these subdirectories contains files for each non-zero block in this column. Directory and file names denote the indices.
     * The files fort.100, fort.101,... which for every MPI rank contain the list of indices of non-zero blocks.
     * H2O-COORD.XYZ-pos-1.xyz contains the coordinates of the atoms.
  * Raw data for a given system size in terms of three-dimensional repetitions of a 32-molecule basic block (NREP times NREP times NREP basic blocks), the basis set and given matrix thresholds can be generated with the script loop.sh. A generic input for CP2K is H2O-dft-ls.inp and a template for the jobscipt for the Noctua system is in template.sh. The script loop.sh takes H2O-dft-ls.inp and the template jobscript template.sh and creates the actual input files as well as matching job scripts.
  * The naming convention of the directories that will be created by loop.sh is (description)_(basis)_(NREP)_(method)_(matrix threshold)_(number of nodes), e.g. H2OTEST_DZVP_4_SUBMATRIX_0.00001_2.

2. Analysis of raw data:
  * The wrapper script sm.sh takes the name of the data directory of a run, the corresponding NREP-parameter and the blocks dimension (SZV=6, DZVP=23) as commandline arguments and performs the following steps:
     1. Compiling sm.f90 to sm.x.
     2. Concatenation of the fort.* files in the data directory of the calculation to the file blocks.dat which contains a list of the indices of all non-zero blocks in the orthongonalized Kohn-Sham matrix.
     3. Then sm.x performs these steps:
        1. Reading in the indices of non-zero blocks from blocks.dat.
        2. Computing the sparsity of the orthongonalized Kohn-Sham matrix in terms of the number of non-zero blocks. (output "sparsity of full matrix (blocks)")
        3. Calculating the dimension of the first submatrix. (output "single column smdim")
        4. Constructing the first submatrix in terms of non-zero blocks.
        5. Outputing the dimension of the first submatrix in terms of blocks. (output "single column smdim")
        6. Computing the sparsity of the first submatrix in terms of the number of non-zero blocks. (output "sparsity of submatrix (blocks)")
        7. Constructing the first submatrix in terms of elements and computes the number of non-zero elements. (output "non-zero elements")

The output of sm.sh, for example, looks like
```
 $bash sm.sh H2OTEST_SZV_7_SUBMATRIX_0.00001_12 7 6
 running forH2OTEST_SZV_7_SUBMATRIX_0.00001_12                                                                  with NREP=           7  and bs=           6
 reading strcture of matrix
 sparsity of full matrix (blocks)     1027530
 single column smdim=           1          90   729000.00000000000     
 sparsity of submatrix (blocks)        3936
 smdim          90
 reading blocks
 non-zero elements       80510   
```
and, thus, if the procedure outlined above is performed for the different basis sets, number of atoms (represented by NREP) and for a given cutoff of the matrix elements, the information depicted in Figures 4 and 11 can be reproduced.
