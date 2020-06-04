# Heuristics

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

2. Additional Programs
  * The clustering explored here requires two additional programs:
     * METIS 5.1.0 that can be downlaoded from http://glaros.dtc.umn.edu/gkhome/metis/metis/download
     * Scikit-learn 0.23.1 that can be downloaded from https://scikit-learn.org/stable/

3. Analysis of raw data:
  * The wrapper script sm.sh takes the name of the data directory of a run, the corresponding NREP-parameter and the blocks dimension (SZV=6, DZVP=23) as commandline arguments and performs the following steps:
     1. Compiling sm.f90 to sm.x.
     2. Concatenation of the fort.* files in the data directory of the calculation to the file blocks.dat which contains a list of the indices of all non-zero blocks in the orthongonalized Kohn-Sham matrix.
     3. Then sm.x perform these steps:
        1. Reading in the indices of non-zero blocks from blocks.dat.
        2. Writing out of the graph representing the sparsite pattern of the orthogonalized Kohn-Sham matrix to the file graph.
        3. Calling the wrapper metis.sh that runs the commandline utility gpmetis that clusters the graph given in the file graph and writes the resulting partition to the file graph.out.
        4. The file graph.out is then read in the fortran program and the clustering is used to estimate the speedup S defined in the paper. Please note, that not the speedup but the inverse of the speedup 1/S is outputted to the commandline together with number of submatrices.(outputs "metis performance...")
        5. The xyz-file H2O-COORD.XYZ-pos-1.xyz containing the coordinates of atoms in the system is read and the three position of each atom in a molecule are averaged because each molecule corresponds to a block.
        6. A k-means implementation from Scikit-learn is called by executing the wrapper clustering.py which performs a k-means clustering and writes the results to the file, which is automatically read from the fortran program.
        7. The fortran program uses the clusters to compute the estimated speedup S defined in the paper. Please note, that not the speedup but the inverse of the speedup 1/S is outputted to the commandline together with number of submatrices.(outputs "unperiodic kmeans scikit performance...")
        8. A simple fortran implementaion of k-means from https://rosettacode.org/wiki/K-means%2B%2B_clustering#Fortran is used to create comparison results. (output "unperiodic kmeans fortran performance")

The output of sm.sh, for example, looks like
```
$bash sm.sh H2OTEST_SZV_6_SUBMATRIX_0.0000001_2 6 6
 running forH2OTEST_SZV_6_SUBMATRIX_0.0000001_2                                                                 with NREP=           6  and bs=           6
 reading strcture of matrix
 sparsity of full matrix (blocks)     1663992
 single column smdim=           1         246   14886936.000000000     
 sparsity of submatrix (blocks)       27644
 metis performance           2   4.2952247714259286     
 metis performance          71   1.0202767243448181     
 metis performance         140  0.79895402614291267     
 metis performance         209  0.73529494859655153     
 metis performance         278  0.68280818826860579     
 metis performance         347  0.66527191192534363     
 metis performance         416  0.65516258901069568     
 metis performance         485  0.63901485647588929     
 metis performance         554  0.64915197187309570     
 metis performance         623  0.64116118373122954     
 metis performance         692  0.64740111767527897     
 metis performance         761  0.66150522851347360 
...
 unperiodic kmeans fortran performance           2   2.0958989659404095     
 unperiodic kmeans scikit performance           2   4.6130134268960115     
 unperiodic kmeans fortran performance          71   1.0773381218381806     
 unperiodic kmeans scikit performance          71   1.0815809304017776     
 unperiodic kmeans fortran performance         140  0.83463646377449285     
 unperiodic kmeans scikit performance         140  0.84361062719969593     
 unperiodic kmeans fortran performance         209  0.75500002992817117     
 unperiodic kmeans scikit performance         209  0.75474859877375611     
 unperiodic kmeans fortran performance         278  0.70876724164835292     
 unperiodic kmeans scikit performance         278  0.71591606492239768     
 unperiodic kmeans fortran performance         347  0.68806820440945615     
 unperiodic kmeans scikit performance         347  0.68487907867054532     
 unperiodic kmeans fortran performance         416  0.67311221674001931
 unperiodic kmeans scikit performance         416  0.66980859775369894     
 unperiodic kmeans fortran performance         485  0.66563454094479291   
 unperiodic kmeans scikit performance         485  0.66423377515071425     
 unperiodic kmeans fortran performance         554  0.66125491172010398
 unperiodic kmeans scikit performance         554  0.66032088184848692     
 unperiodic kmeans fortran performance         623  0.66241108383184499 
 unperiodic kmeans scikit performance         623  0.65884677052060681
 unperiodic kmeans fortran performance         692  0.65991116035112973
 unperiodic kmeans scikit performance         692  0.65948481924293756     
 unperiodic kmeans fortran performance         761  0.66385885665084576 
 unperiodic kmeans scikit performance         761  0.65892432748237850     
...
```
and, thus, if executed for the example in the paper will create the information depicted in Figure 5.
