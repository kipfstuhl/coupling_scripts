Efficient coupling of PDEs based on weak transmission conditions
---

The MATLAB code in this repository is complementary to the paper *Efficient coupling of PDEs based on weak transmission conditions* (Deparis, S. and Pegolotti, L., 2017). It was used to generate the figures in the *Numerical applications* section. Most of the functionalities we use are implemented in the *feamat* repository (https://github.com/lucapegolotti/feamat), which is a set of routines for the solution of Finite Element problems on structured meshes. Note that the code is not particularly optimized, and therefore *feamat* should be used only for prototyping purposes. The *feamat* repository is hereby included as submodule: see *How to clone* section for details.

HOW TO CLONE
---
...

HOW TO USE
---
The functions of *feamat* must be before hand included by running the run_addpaths.m script in the root of the project. The repository is structured in multiple folders, each containing the scripts for generating the figures of a particular section of the paper. Each script is commented, such that it should be possible to follow the steps without knowledge of the underlying *feamat* functions
