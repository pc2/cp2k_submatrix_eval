# Overview

This part of the evaluation consists of four parts:

* precise: Determination of an exactly calculated band structure energy
* Submatrix: Runs performed with different epsilon_filter values using the submatrix method
* Newton-Schulz: Runs with different epsilon_filter values using Newton-Schulz
* plot: Aggregated data and scripts used for plotting

The `run.sh` scripts generate input files based on a template and run CP2K. the `filter*.sh` scripts can be used to extract the relevant data from the output files. The important lines in the output are (examples):
```
Energy Tr(DK)  -12166.978248409592
```
and
```
density_matrix_sign_fixed_mu         1  6.0    0.000    0.001    0.751    0.764
```
where the last value (here: 0.764) divided by the first value (here: 1) denotes the time that was required to calculate the density matrix once.

Note that the energy is output in Hartree. For our plots, we convert the energy differences into meV/atom.
