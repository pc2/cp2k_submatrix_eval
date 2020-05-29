# Overview

The `run.sh` script generates input files for differently sized systems based on a template. The number used in the file names is the NREP factor that we use to scale up the system in each dimension.

The relevant part of the output is the time used to calculate the density matrix, which is the last number in the following line (example):
```
density_matrix_sign_fixed_mu         1  6.0    0.000    0.001    0.351    0.352
```

This folder also contains the aggregated results (`SM.dat`) and the script used for plotting.
