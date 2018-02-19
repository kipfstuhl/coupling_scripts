Efficient coupling of PDEs based on weak transmission conditions
-------------------------

The MATLAB code in this repository is complementary to the paper *Coupling non-conforming discretizations of PDEs by spectral approximation of the Lagrange multiplier space* (Deparis, S. and Pegolotti, L., 2018). It was used to generate the figures in the *Numerical applications* section. Most of the functionalities we use are implemented in *feamat* (https://github.com/lucapegolotti/feamat), which is a set of routines for the solution of Finite Element problems on structured meshes. Note that the code is not particularly optimized, and therefore *feamat* should be used only for prototyping purposes. The *feamat* repository is hereby included as submodule.

HOW TO USE
-------------------------
The functions of *feamat* must be beforehand included by running the run_addpaths.m script in the root of the project. The repository is structured in multiple folders, each containing the scripts for generating every figure in the paper. In each figure* folder, a run_all.m scrips executes all the scripts that are necessary to generate the data to be plotted.
