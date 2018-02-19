Coupling non-conforming discretizations of PDEs by spectral approximation of the Lagrange multiplier space
-------------------------
The MATLAB code in this repository is complementary to the paper *Coupling non-conforming discretizations of PDEs by spectral approximation of the Lagrange multiplier space* (Deparis, S. and Pegolotti, L., 2018). It was used to generate the figures in the *Numerical results* section. Most of the functionalities we use are implemented in *feamat* (https://github.com/lucapegolotti/feamat), which is a set of MATLAB routines for the solution of Finite Element problems on structured meshes. Note that the code is not particularly optimized, and therefore *feamat* should be used only for prototyping purposes. The *feamat* repository is hereby included as submodule.

Tested version of MATLAB: 2017b. Previous versions might not work properly.

How to clone
-------------------------
Run
```
 git clone --recurse-submodules git@bitbucket.org:pegolotti/coupling_scripts.git
```
to recursively clone the directory and *feamat*.

How to use
-------------------------
The functions of *feamat* must be beforehand included by running the run_addpaths.m script in the root of the project. The repository is structured in multiple folders, each containing the scripts for generating every figure in the paper. In each figure* folder, a run_all.m scrips executes all the scripts that are necessary to generate the data to be plotted.
