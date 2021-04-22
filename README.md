# PeriodicSpectra
Matlab routine to compute the spectrum of a Schrödinger operators with periodic potential defined on the real line. The implementation is based on a method recently developed in https://arxiv.org/abs/2104.09575. A Floquet-Bloch transform is used to replace the real line by a finite domain, after which the operator is discretized by choosing a Fourier basis. As proved in the article, the method is guaranteed to converge.
The code takes the potential function as an input and returns a set of complex numbers, which approximates the spectrum of the Schrödinger operator.

The main building blocks of the routine are
* `Main.m` - Main script that returns approximation of the spectrum;
* `build_lattice(z1, z2, horizontal_resolution)` - computes initial lattice in the complex plane, defined by the two corners `z1, z2`. The number of lattice points in horizontal direction is given by`horizontal_resolution`. ;
* `potential(x)` - returns value of potential function for input x in [0,1];
* `fourier_coefficient(k, n)` - Accepts frequency vector `k`, grid size `n` and returns vector of Fourier coefficients of the potential;
* `compute_potential_matrix(a,N)` - Computes operator matrix of the potential in Fourier representation. Input: Fourier coefficients `a` and size of potential matrix `N`.

Any comments or queries are welcome at https://frank-roesler.github.io/contact/
