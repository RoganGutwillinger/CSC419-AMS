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

Let $M<sup>0</sup>$ be the discrete triangle mesh representing the original surface, let M by the discrete triangle mesh representing the deformed surface, and let o be a viewpoint oriented in the direction of the surface. The paper 
defines d(M, M<sup>0</sup>, o) as the "surface similarity" or perceived difference of the meshes M and M<sup>0</sup>; this is the value the paper seeks to minimize while deforming the mesh.

This implementation lacks two components of the full scope of paper. Firstly, it does not implement...

## References
1. PAPER REFERENCE HERE
