# Simulated 2D Diffusion in presence of a partially reflecting mesh

This source code simulates 2D diffusion (e.g. biological membrane components) in presence of a partially reflective meshgrid.

It was used in the Biophys. Journal publication: Lateral membrane diffusion modulated by a minimal actin cortex, Heinemann F., Vogel SK, Schwille P., 2013 (https://www.cell.com/biophysj/fulltext/S0006-3495(13)00260-9)

Is samples fluorescence traces in 9 Gaussian shaped detection spots and computes fluorescence correlation spectroscopy curves using the multiple tau algorithm.

The simulation was optimized for performance in C++ (with a small inline assembly part for a much faster floor computation).

Dependencies:
- boost (https://www.boost.org/)
(tested with the somewhat old: boost 1_48_0)
- cpp-yaml (https://github.com/jbeder/yaml-cpp)

Environment:
- Windows7/10
- Should in principle also run under Unix

Usage:

>simulated_diffusion.exe simulation_list.yaml


