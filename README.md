# CSC419 Final Paper Implementation

This project is a libigl-style implementation of the [Appearance Mimicking Surfaces](https://cims.nyu.edu/gcl/papers/mimicking-2014.pdf) paper written by Christian Schuller, Daniele Panozzo, and Olga Sorkine-Hornung.

## Installation, Layout, and Compilation

Ensure that you have installed cmake and a modern c++ compiler on Mac OS X, Linux, or Windows.

In the root project directory (The same directory containing this README file), issue:

    mkdir build
    cd build
    cmake ..
    make 

## Execution

Once built, you can execute the demo from inside the `build` directory using 

    ./ams [path to mesh]

This prject includes two sample meshes in the `data` directory, `bunny.off` and `decimated-knight.off`.

## Implementation

This project implements the main goal of the paper, constraining a 3D object to a specific volume while maintaining its appearance from a specific viewpoint. 

Let $M^0$ be the discrete triangle mesh representing the original surface, let $M$ by the discrete triangle mesh representing the deformed surface, and let $o$ be a viewpoint oriented in the direction of the surface. The paper 
defines $d(M, M^0, o)$ as the "surface similarity" or perceived difference of the meshes $M$ and $M^0$; this is the value the paper seeks to minimize while deforming the mesh. Let $\lambda_i = || v_i - o ||$ where $v_i$ is 
the $i$ th vertex on the deformed mesh. We solve for these $\lambda_i$ values, as we can construct the deformed mesh with them.

Let $A^0_i$ be the Voronoi area of the original mesh associated with the $i$ th vertex, let $L^0$ be the discrete Laplace-Beltrami operator of $M^0$, let $\hat{V}$ be the matrix stacking all normalized vertices of the mesh,
let $D_{\lambda}$ be the diagonal matrix stacking the values of $\lambda_i$ along its diagonal, and define $D_{\lambda^0}$ similarily with values of $\lambda^0$. Finally, define $w_i$ as the weight of the $i$ th vertex of the
mesh. We assign a weight of $1$ to vertices that are visible from the viewpoint and a weight of $0.1$ to vertices that are hidden from the viewpoint. As descrbied in the paper, we can define the surface similarity as follows:

$$
d(M, M^0, o) = \sum_{i \in V} w^2_i A^0_i || (L^0 D_{\lambda} \hat{V})_{i} - (L^0 D_{\lambda^0} \hat{V})_{i} \frac{\lambda_i}{\lambda^0_i} ||^2
$$

This implementation lacks two components of the full scope of paper. Firstly, it does not implement...

## References
1. PAPER REFERENCE HERE
