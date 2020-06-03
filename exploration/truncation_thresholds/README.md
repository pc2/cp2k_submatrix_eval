# Overview

The `run.sh` script generates input files for differently sized systems and different truncation thresholds based on a template. The number at the front of file names is the NREP factor that we use to scale up the system in each dimension. Then there is one reference output `precise` (truncation threshold 1E-12) and multiple outputs `truncated` where the number following corresponds to the exponent of the truncation threshold.

The relevant part of the output is the total energy (example from `4-truncated-4.out`):
```
ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):           -34963.073145484595443
```

This folder also contains the aggregated results (`results.dat`) and the script used for plotting.
