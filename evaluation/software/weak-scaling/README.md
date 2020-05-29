# Overview

This part of the evaluation consists of three parts:

* Submatrix: Runs performed with different system sizes and numbers of nodes using the submatrix method
* Newton-Schulz: Runs with different system sizes and numbers of nodes using Newton-Schulz
* plot: Aggregated data and scripts used for plotting

The `run-*.sh` scripts generate input files based on a template and run CP2K with the number of cores using slurm. The number used in the file names is the linear replication factor that we use to scale up the system.

The relevant part of the output is the time used to calculate the density matrix, which is the last number in the following line (example):
```
density_matrix_sign_fixed_mu         1  6.0    0.019    0.022    6.735    6.752
```

# Repition of measurements

For Newton-Schulz we observed that the runtime is relatively sensitive to circumstances that are not under our control like network load in the cluster or the exact placement of processes within the cluster. We therefore repreated all measurements six times and took the minimum run time to draw the graph. For one value, there was still a visible outlier after six runs, so that we ran additional measurements (run7 - run9) at a later point in time.
