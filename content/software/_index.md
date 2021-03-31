---
title: "Software"

---

## [PeriodicFMM](https://github.com/wenyan4work/PeriodicFMM)
**This package has retired and have been replaced by STKFMM**

PeriodicFMM is a software package that implements the flexibly periodization method for kernel independent FMM described in these two papers:

1. [ Wen Yan, Michael Shelley (2018). Flexibly Imposing Periodicity in Kernel Independent FMM: A Multipole-to-Local Operator Approach. Journal of Computational Physics](https://doi.org/10.1016/j.jcp.2017.11.012)
2. [ Wen Yan, Michael Shelley (2018). Universal Image Systems for Non-Periodic and Periodic Stokes Flows above a No-Slip Wall. Journal of Computational Physics](https://doi.org/10.1016/j.jcp.2018.08.041)

The package provides Stokes single layer (Stokelet) kernel FMM for the following cases:

| Geometry                | Boundary Conditions                         |
| ----------------------- | ------------------------------------------- |
| 3D                      | Non-, Singly-, Doubly-, and Triply-Periodic |
| 3D above a no-slip wall | Non-, Singly-, and Doubly-Periodic          |

This package is based on the 'new_BC' branch of the massively parallel KIFMM package [PVFMM](https://github.com/wenyan4work/pvfmm.git) originally developed by [Dr. Dhairya Malhotra](https://cims.nyu.edu/~malhotra/).
This package provides native interfaces in C++ with full `OpenMP` + `MPI` support and scales beyond thousands of cores. Kernels are hand-written in `AVX/AVX2` SIMD instructions to improve speed. Native C++ interfaces are provided, with C/Fortran/Python wrappers. The wrappers are implemented with the help of [Dr. Florencio Balboa Usabiaga](https://scholar.google.com/citations?user=1lCVNZwAAAAJ&hl=en).

## [STKFMM](https://github.com/wenyan4work/STKFMM)
It computes the classic kernel sum problem: for a given set of single layer sources $s^j$ at points $y_s^j$, double layer sources $d^j$ points $y_d^j$, target points $x_t$, and single layer potential $K_s$, double layer potential $K_d$:
$$p(x_t^i)=\sum_j K_s(x_t^i,y_s^j) s^j +\sum_j K_d (x_t^i,y_d^j)d^j $$

**Note** For some problems the kernels $K_s$ and $K_d$ may not be linear operators applied onto $s^j, d^j$.

This package computes Laplace kernel $L=\dfrac{1}{4\pi r}$ and Stokeslet kernel $G_{ij}=\dfrac{1}{8\pi}\left(\dfrac{1}{r}\delta_{ij}+\dfrac{1}{r^3}r_ir_j\right)$ and their derivatives.

Here is a detailed table, in which the summation $\sum$ is dropped for clarity and the subscript indices $i,j,k,l$ denote the tensor indices.
Einstein summation and comma notation are used to simplify the expressions, for example, $G_{ij,i}f_j=\nabla\cdot (\mathbf{G}\cdot\mathbf{f})$.

In the table:

1. **NA** means input ignored
2. $Q_{ij}$, $D_{ij}$ are 3x3 tensors written as 9-dimension vectors in row-major format
3. $\nabla\nabla p$ is symmetric so it is written as $p_{,xx},p_{,xy},p_{,xz},p_{,yy},p_{,yz},p_{,zz}$.
4. For `RPY`, `StokesRegVel` and `StokesRegVelOmega` kernels, the parameter $b$ and $\epsilon$ can be different for each source point, and the summations are nonlinear functions of $b$ and $\epsilon$. Also $b$ and $\epsilon$ must be much smaller than the lower level leaf box of the adaptive octree, otherwise the convergence property of KIFMM is invalidated.
5. For all kernels, the electrostatic conductivity and fluid viscosity are ignored (set to 1).
6. The regularized Stokeslet is $G_{ij}^\epsilon = \dfrac{1}{8\pi}\dfrac{r^{2}+2 \epsilon^{2}}{\left(r^{2}+\epsilon^{2}\right)^{3 / 2}} \delta_{i j} f_j+\dfrac{1}{\left(r^{2}+\epsilon^{2}\right)^{3 / 2}} r_ir_jf_j$.
7. For Stokes `PVel`, `PVelGrad`, `PVelLaplacian`, and `Traction` kernels, the pressure and velocity fields are:
   $$ p=\frac{1}{4 \pi} \frac{r_{j}}{r^{3}} f_{j} + \frac{1}{4 \pi}\left(-3 \frac{r_{j} r_{k}}{r^{5}}+\frac{\delta_{j k}}{r^{3}}\right) D_{j k}, \quad u_{i}=G_{ij}f_j + \frac{1}{8 \pi \mu}\left(-\frac{r_{i}}{r^{3}} trD\right) + \frac{1}{8 \pi \mu}\left[-\frac{3 r_{i} r_{j} r_{k}}{r^{5}}\right] D_{j k} $$

| Kernel              | Single Layer Source (dim)  | Double Layer Source (dim) | Summation                                       | Target Value (dim)                                                  |
| ------------------- | -------------------------- | ------------------------- | ----------------------------------------------- | ------------------------------------------------------------------- |
| `LapPGrad`          | $q$ (1)                    | $d$ (3)                   | $p=Lq-L_{,j}d_j$                                | $p,\nabla p$ (1+3)                                                  |
| `LapPGradGrad`      | $q$ (1)                    | $d$ (3)                   | $p=Lq-L_{,j}d_j$                                | $p,\nabla p, \nabla\nabla p$ (1+3+6).                               |
| `LapQPGradGrad`     | $Q_{ij}$ (9)               | NA                        | $p=L_{,ij}Q_{ij}$                               | $p,\nabla p, \nabla\nabla p$ (1+3+6).                               |
| `Stokes`            | $f_j$ (3)                  | NA                        | $u_i = G_{ij} f_j$                              | $u_i$ (3)                                                           |
| `RPY`               | $f_j,b$ (3+1)              | NA                        | $u_i = (1+\frac{1}{6}b^2\nabla^2) G_{ij} f_j$   | $u_i,\nabla^2 u_i$ (3+3)                                            |
| `StokesRegVel`      | $f_j,\epsilon$ (3+1)       | NA                        | $u_i = G_{ij}^\epsilon f_j$                     | $u_i$                                                               |
| `StokesRegVelOmega` | $f_k,n_l,\epsilon$ (3+3+1) | NA                        | See Appendix A of doi 10.1016/j.jcp.2012.12.026 | $u_i,w_j$ (3+3)                                                     |
| `PVel`              | $f_j,trD$ (3+1)            | $D_{jk}$ (9)              | see above                                       | $p,u_i$ (1+3)                                                       |
| `PVelGrad`          | $f_j,trD$ (3+1)            | $D_{jk}$ (9)              | see above                                       | $p,u_i,p_{,i},u_{i,j}$ (1+3+3+9)                                    |
| `PVelLapLacian`     | $f_j,trD$ (3+1)            | $D_{jk}$ (9)              | see above                                       | $p,u_i,u_{i,jj}$ (1+3+3)                                            |
| `Traction`          | $f_j,trD$ (3+1)            | $D_{jk}$ (9)              | see above                                       | $\sigma_{ij}=-p \delta_{i j}+\mu\left(u_{i, j}+u_{j, i}\right)$ (9) |

### Features

- All kernels are hand-written with optimized SIMD intrinsic instructions.
- Singly, doubly and triply periodicity in a unified interface.
- Support no-slip boundary condition imposed on a flat wall through image method.
- Single Layer and Double Layer potentials are simultaneously calculated through a single octree.
- M2M, M2L, L2L operations are combined into single layer operations only.
- All PVFMM data structures are wrapped in a single class.
- Multiple kernels can be activated simultaneously.
- Complete MPI and OpenMP support.


This package is based on the 'new_BC' branch of the massively parallel KIFMM package [PVFMM](https://github.com/wenyan4work/pvfmm.git) originally developed by Dr. Dhairya Malhotra. 


## [SimToolbox](https://github.com/wenyan4work/SimToolbox)
SimToolbox is a set of loosely coupled handy tools (a 'toolbox') to simplify the development and maintainance of parallel particle-tracking simulations on HPC

At the lowest level the toolbox relies on `FDPS` and `Trilinos`.
`FDPS` refers to the [Framework for Developing Particle Simulator] (https://github.com/FDPS/FDPS), which provides necessary infrastructures of parallel particle-tracking simulations such as domain decomposition and near neighbor detection.
[`Trilinos`](https://github.com/trilinos/Trilinos) is the massive C++ HPC project for distributed linear algebra and some other handy tool.

At the application level the toolbox implements a stable and efficient collision-resolution algorithm for general smooth-shape particles based on geometric constrained optimization. 
The method is demonstrated to track the system collision stress accurately and is described in detail in the following paper:

1. [Wen Yan, Huan Zhang, Michael J. Shelley (2019). Computing Collision Stress in Assemblies of Active Spherocylinders: Applications of a Fast and Generic Geometric Method. The Journal of Chemical Physics](https://aip.scitation.org/doi/10.1063/1.5080433)

The toolbox also includes some useful code for supportive tasks in simulations, including:

1. Output to binary XML VTK data files
2. Parallel random number generation with full OpenMP and MPI support, with the help of [`TRNG`](https://github.com/rabauke/trng4) library.
3. MPI data directory based on [Zoltan Distributed Directory Utility ](https://cs.sandia.gov/Zoltan/ug_html/ug_util_dd.html)

## [SafeFFT](https://github.com/wenyan4work/SafeFFT)
SafeFFT is a thread-safe c++ wrapper for FFTW and MKL. In FFTW3 (or MKL) the only thread safe functions are `fftw_execute_...()`. This sometimes poses problems on the structure of multithreading code, where each thread may need to perform FFTs with different `fftw_plan`. This simple wrapper around FFTW3 aims at making things easier by maintain a global hash table of `fftw_plan`, where each thread may insert new plans and read already allocated plans. The hash table is locked such that multiple readers can access it but only one thread can insert new entries. This allows multiple threads to reuse already allocated plan simultaneously, without allocating a new plan everytime. I believe this approach has some performance and design advantage because allocating a new plan everytime for every thread requires a mutex lock to allow only one thread to create a plan. This simple wrapper fits a case where a large number of FFTs must be processed, but the total number of different FFT plans are not that large, and it may also be hard to preallocate all possible plans before running any FFTs. 

## [CSDmp](https://github.com/wenyan4work/CSDmp)
This is a pedagogical codebase in Fortran 95 to demonstrate the Conventional Stokesian Dynamics algorithm:

1. Durlofsky, L., Brady, J. F. & Bossis, G. Dynamic Simulation of Hydrodynamically Interacting Particles. Journal of Fluid Mechanics 180, 21–49 (1987).

This package only computes a small number of monodisperse spherical particles in unbounded (Non-periodic boundary condition) Stokes fluid. 
The program relies on standard BLAS and LAPACK routines. Modify the makefile according to your software environment. Mixed precision means the velocity is calculated in single precision and the position is in double precision.

## [DemoSIMD](https://github.com/wenyan4work/DemoSIMD)
This is a collection of short code for beginners to understand SIMD and cache. It includes demos about `memcpy`, `gemm`, fast-inverse-squareroot, and cache blocking. 