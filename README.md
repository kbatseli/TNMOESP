Tensor Network MoESP method for the identification of polynomial state space models in MATLAB&copy;/Octave&copy;
----------------------------------------------------------------------------------------------------------------

This package contains MATLAB/Octave implementations of both the matrix and tensor network implementation of the MOESP algorithm for the identification of specific polynomial state space models.

1. Functions
------------

* demo

Demonstrates the usage of TNmoesp and sim

*  [A,BD,C]=TNmoesp(y,u,d) or [A,BD,C]=TNmoesp(y,u,d,n)

Tensor network implementation of the MOESP algorithm for the identification of a specific polynomial state space model, determined by matrices A,C, a tensor network BD and the degree d of the polynomials. 

* y=simTNss(A,BD,C,x0,u,d)

Simulates the polynomial state space system as specified by the A,C,D matrices, the BD tensor network, the initial state x0, the input vectors u the degree d of the polynomials. 

* [A,B,C,D]=moespd(y,u,d) or [A,B,C,D,e]=moespd(y,u,d,n)

Matrix implementation of the MOESP algorithm for the identification of a specific polynomial state space model, determined by matrices A,B,C,D and the degree d of the polynomials. 

* y=simssd(A,B,C,D,x0,u,d)

Simulates the polynomial state space system as specified by the A,B,C,D matrices, the initial state x0, the input vectors u and the degree d of the polynomials. 


2. Reference
------------

Tensor network subspace identification of polynomial state space models

Authors: Kim Batselier, Ching Yun Ko, Ngai Wong
