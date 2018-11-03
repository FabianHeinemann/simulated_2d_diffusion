# Simulated 2D Diffusion in presence of a partially reflecting mesh

This source code simulates 2D diffusion (e.g. biological membrane components) in presence of a partially reflective meshgrid.

- N particles perform a random walk in 2D
- A meshgrid (repesenting a membrane skeleton) reflects particles, which can only pass with a probability p_jump
- 9 Gaussian shaped detection spots are placed on the membrane and simulated fluorescence trances are aquired
- The fluorescence traces are autocorrelated using the multiple tau algorithm, resulting in realistic fluorescence correlation spectroscopy (FCS) curves

The code was used to generate simulations published in the Biophysical Journal paper
<i>Lateral membrane diffusion modulated by a minimal actin cortex</i>, Heinemann F., Vogel SK, Schwille P., 2013 (https://www.cell.com/biophysj/fulltext/S0006-3495(13)00260-9)

The simulation was heavily optimized for performance and is written in C++ with a small inline assembly part for a much faster floor computation.

Dependencies:
- boost (https://www.boost.org/)
(tested with the somewhat old: boost 1_48_0)
- cpp-yaml (https://github.com/jbeder/yaml-cpp)

Environment:
- Tested under Windows 7 and 10
- Should in principle also run under Unix

Usage in console:

>simulated_diffusion.exe simulation_list.yaml

To do:
- describe input files
- Images
