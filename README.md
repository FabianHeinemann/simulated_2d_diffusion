# Simulated 2D Diffusion in presence of a partially reflecting mesh

This source code simulates 2D diffusion (e.g. biological membrane components) in presence of a partially reflective meshgrid.

- N particles perform a random walk in 2D
- A meshgrid (repesenting a membrane skeleton) reflects particles, which can only pass with a probability p_jump
- 9 Gaussian shaped detection spots are placed on the membrane and simulated fluorescence trances are aquired
- The fluorescence traces are autocorrelated using the multiple tau algorithm, resulting in realistic fluorescence correlation spectroscopy (FCS) curves

The code was used to generate simulations published in the Biophysical Journal paper
<i>Lateral membrane diffusion modulated by a minimal actin cortex</i>, Heinemann F., Vogel SK, Schwille P., 2013 (https://www.cell.com/biophysj/fulltext/S0006-3495(13)00260-9)

The simulation was quite heavily optimized for performance using a profiler. It is written in C++ (with a small inline assembly part for a much faster floor computation).

![Simulation_illustration](https://github.com/FabianHeinemann/simulated_2d_diffusion/blob/master/images/Frame_0.png)

Dependencies:
- boost (https://www.boost.org/)
(tested with the somewhat old: boost 1_48_0)
- cpp-yaml (https://github.com/jbeder/yaml-cpp)

Environment:
- Tested under Windows 7 and 10
- Should in principle also run under Unix

Usage in console:

>simulated_diffusion.exe simulation_list.yaml

![Console screenshot](https://github.com/FabianHeinemann/simulated_2d_diffusion/blob/master/images/console.png)

The yaml file contains repeating blocks of 4 parameters. An example is in the repository. Example block:

>voronoimesh: ./mesh/1-2 mean82 sd31.txt
>resultfile: ./1-2 mean82 sd31_p001.txt
>pjump: 0.01
>Tmax: 100

Each block is one simulation and will output one result file. All blocks are computed until the end.

Parameter description:
- voronoimesh: A text file with coordinates of the mesh. See example in repository. I used a custom program performing Voronoi tesselation (https://en.wikipedia.org/wiki/Voronoi_diagram), which I can distrubute on request.
- resultfile: Name of output file to write to. Output will be a simple textfile containing the 9 simulated FCS autocorrelation curves and their average.
- p_jump: Probability to cross a mesh fibre. 
- Tmax time in simulated seconds (I recommend 100-300s, with larger values for dense meshes and / or high pjump).

Example for a simulated FCS curve (without mesh; with mesh the curve will move to right and a second component may appear, see paper https://www.cell.com/biophysj/fulltext/S0006-3495(13)00260-9).

![FCS curve](https://github.com/FabianHeinemann/simulated_2d_diffusion/blob/master/images/fcs_free.png)
