# Confinement and driving effects on continuous and discrete model interfaces ##
Paul Gersberg PhD, Université de Bordeaux & ENS Lyon
July 16th 2020


## Overview

This git records all codes I have developped for my thesis in C++, Python and Bash.

  * Abstract [Thesis](https://www.theses.fr/2020BORD0084)

This thesis examines the properties of the interface between two phases in phase separated systems. We are interested in how finite size effects modify the statistical properties of these interfaces, in particular how the dependence of the free energy on the system size gives rise to long range critical Casimir forces close to the critical point. Often the interfaces in phase separated systems are described by simplified or coarse grained models whose only degrees of freedom are the interface height. We review how the statics and dynamics of these interface models can be derived from microscopic spin models and statistical field theories. We then examine finite size effects for continuous interface models such as the Edwards Wilkinson model and discrete models such as the Solid-On-Solid model and discuss their relevance to the critical Casimir effect.
In the second part of the thesis we examine models of driven interfaces which have nonequilibrium steady states. We develop a model C type model of an interface which shows a nonequlibrium steady state even with constant driving. The resulting nonequlibrium steady state shows properties seen in experiments on sheared colloidal systems, notably the suppression of height fluctuations but an increase in the fluctuations’ correlation length. Finally we propose a new model for one dimensional interfaces which is a modification of the solid-on-solid model and containing an extra entropic term, whose correspondance with physical systems is yet to be found.

## Files

This distribution consists of the following files:

  * Readme.md, the file you are currently reading.

  * Copying.txt, a text file containing the GNU General Public License version 3.

This distribution consists of the following directories:

###   Main codes in C++. Use bash script ./Daedalus to compile the code and launch the corresponding jobs on Slurm, with Python codes for posttreatment
  * Ising/          : Ising simulations in C++ in different contexts, with Glauber and Kawasaki implementations
  * SOS/            : SOS simulations in C++ in different contexts, with Glauber and Kawasaki implementations
  * POP/            : POP model, with still some problems with detailed balance, and GSL matrix diagonalisation

In each directory, if there is a subdirectory 'prog/', all the sources will be found therein.

###   Different uses of main codes
  * Casimir/        : SOS simulations with crossover hamiltonian
  * Equilibrium/    : equilibrium time computation for Ising, SOS, RSOS and POP 
  * CompModels/     : compairson of equilibrium values for Ising, SOS, RSOS and POP 

###   C++ 1D finite element simulation 
  * Langevin/       : 1D finite elements simulations in stochastic equations with stochastic CFL conditions

###   Python codes for scientific computations
  * figuresPyth/    : Some figures of the manuscript with matplotlib, like Fig  2.5 or 4.2
  * TransferMatrix/ : Matrix diagonalisations with numpy.linalg, scipy.optimize.curve_fit 

###   Miscellanous C++ things
  * pRNG/           : a [fast implementation of pseudo-Random Number Generators in C++](https://martin.ankerl.com/2018/12/08/fast-random-bool/)
  * libeingen/      : [Libeigen C++ library](https://gitlab.com/libeigen/eigen)
  * Sauvegardes/    : source files from the [previous PhD](https://www.theses.fr/2015ENSL1025), codes and example of how to use OpenMP and MPI from a formation

###   Python image video processing and finite elements for the miscellanous physics projects
  * IPT/            : some Python simulations made for the [French Physicist's Tournament](https://france.iptnet.info/)



## Support on Linux

The code has been developped for the scientific cluster [Mésocentre de Calcul Intensif d'Aquitaine](https://www.mcia.fr/projects/cluster-curta/wiki)
in a Linux environnement and a Slurm scheduler. All bash scripts and flags have been carefully chosen to match the cluster's configurations, and you might need to adapt them to yours.

## Compiling yourself from the sources

Compilation might be tricky because of the MCIA's configuration, which made linking to libraries a bit hard. For the g++ flags, always refer to compile and Daedalus files.

## Improving or reusing the code

* The latest source can always be found on [GitHub](https://github.com/Bulbille/Curta).

## Terms of use

This code is free, and distributed under the **GNU General Public License version 3**
(GPL v3). Essentially, this means that you are free to do almost exactly
what you want with the program, including distributing it among your
friends, making it available for download from your web site, selling
it (either by itself or as part of some bigger software package), or
using it as the starting point for a software project of your own.

The only real limitation is that whenever you use this code in
some way, you must always include the full source code, or a pointer
to where the source code can be found. If you make any changes to the
source code, these changes must also be made available under the GPL.

If you use this code for your own scientific projects, you must cite the author as a reference in the scientific papers this code was used in.

For full details, read the copy of the GPL v3 found in the file named
*Copying.txt*.
