# Implementation of data-driven coupling models at Road Junctions
Codes used in the paper
> Data-driven Models for Traffic Flow at Junctions

by M. Herty & N. Kolbe

A preprint of the paper will be available soon.

Julia codes were written by N. Kolbe

## Contents
The repository contains a Julia project environment, which includes the package `DataDrivenroadjunctions.jl` implementing the coupling models and their parameter estimation, and a directory (`experiments`) with scripts for various use cases. The former implements Riemann solvers, handling of trajectory and macroscopic data, fitting, and a numerical scheme for experiments on 2-to-1 networks. The code relies on other Julia packages, in particular, [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl) used for global optimization in the parameter fitting, [Flux.jl](https://github.com/FluxML/Flux.jl) implementing neural networks and [CentralNetworkScheme.jl](https://github.com/nklb/CentralNetworkScheme) for numerical simulations on the network. The folder `models` contains the best fit models C1, C2, ML1, ML2 and ML3 from the paper, the folder `FCD` is used to store the vehicle trajectory data and the results of the scripts in the `experiments` folder will be stored in the folder `out`.

The scripts `experiments/`

Car trajectory data has not yet been published and is therefore not included in this repository. But it is available upon request. 

## Usage 
Clone the repository on your local machine. To run the scripts activate the project environment, for example, by starting `julia` from the command line within the root folder of the repository, then changing to pkg mode typing `]` and afterwards running `activate .`. Then the experiments can be run by including the corresponding scripts, e.g., `include("experiments/prediction.jl")`.  You can get started by going through the scripts and modifying them as you like. Further documentation will be added in the future. 

In case of questions or problems please contact the authors of the paper or file a GitHub issue in this repository.
