# Overview

The `run-*.sh` scripts generate an input file based on a template and run CP2K on different numbers of nodes using slurm. The number used in the file names is the number of nodes used.

The relevant part of the output is the time used to calculate the density matrix, which is the last number in the following line (example):
```
density_matrix_sign_fixed_mu         1  6.0    0.001    0.014   12.353   12.378
```

This folder also contains the aggregated results (`SM.dat`) and the script used for plotting.
