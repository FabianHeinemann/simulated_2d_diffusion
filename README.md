# Simulated 2D Diffusion in presence of a partially reflecting mesh

## Brief description

This source code simulates 2D diffusion (e.g. biological membrane components) under periodic boundary conditions in presence of a partially reflective meshgrid.

- N particles perform a random walk in 2D
- A meshgrid (repesenting a membrane skeleton) reflects particles, which can only pass with a probability p_jump
- 9 Gaussian shaped detection spots are placed on the membrane and simulated fluorescence trances are aquired
- The fluorescence traces are autocorrelated using the multiple tau algorithm, resulting in realistic fluorescence correlation spectroscopy (FCS) curves

The code was used to generate simulations published in the Biophysical Journal paper
_Lateral membrane diffusion modulated by a minimal actin cortex_, Heinemann F., Vogel SK, Schwille P., 2013 (https://www.cell.com/biophysj/fulltext/S0006-3495(13)00260-9)

The simulation was quite heavily optimized for performance using a profiler. It is written in C++ (with a small inline assembly part for a much faster floor computation).

![Simulation_illustration](https://github.com/FabianHeinemann/simulated_2d_diffusion/blob/master/images/Frame_0.png)

## Dependencies

- boost (https://www.boost.org/)
(tested with the somewhat old: boost 1_48_0)
- cpp-yaml (https://github.com/jbeder/yaml-cpp)

## Environment

- Tested under Windows 7 and 10 (should in principle also run under Linux)

## Basic usage

1. Compile the project
2. Create file(s) decribing meshgrid (see mesh/1-2 mean82 sd31.txt)
3. Create / edit yaml file configuring the simulation
4. Open cmd, browse to project folder and run 
>simulated_diffusion.exe simulation_list.yaml
![Console screenshot](https://github.com/FabianHeinemann/simulated_2d_diffusion/blob/master/images/console.png)
5. Drink some coffee and wait...
6. Analyze simulation results written to text files

## Description of configuration yaml

The yaml file contains repeating blocks of 4 parameters. An example is in the repository. Example block:

>voronoimesh: ./mesh/1-2 mean82 sd31.txt

>resultfile: ./1-2 mean82 sd31_p001.txt

>pjump: 0.01

>Tmax: 100

Each block is one simulation and will output one result file. All blocks are computed until the end.

### Parameter description:
_voronoimesh:_ A text file with coordinates of the mesh. See example in repository. I used a custom program performing Voronoi tesselation (https://en.wikipedia.org/wiki/Voronoi_diagram), which I can distrubute on request.
_resultfile:_ Name of output file to write to. Output will be a simple textfile containing the 9 simulated FCS autocorrelation curves and their average.
_p_jump:_ Probability to cross a mesh fibre. 
_Tmax:_ time in simulated seconds (I recommend 100-300s, with larger values for dense meshes and / or high pjump).

## Example result when plotted

A simulated FCS curve (without mesh; with mesh the curve will move to right and a second component may appear, see paper https://www.cell.com/biophysj/fulltext/S0006-3495(13)00260-9).

![FCS curve](https://github.com/FabianHeinemann/simulated_2d_diffusion/blob/master/images/fcs_free.png)

## Note:

CSimulatedDiffusion.cpp contains further physical and simulation parameters, which can be modified, e.g.:
- Diffusion coefficient: D (default: 10 µm²/s)
- Confocal spot size: r0 (default: 250 nm)
- Particle number: nParticles (default: 1000)
