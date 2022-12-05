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

This project includes two sample meshes in the `data` directory, `bunny.off` and `decimated-knight.off`.

## Demo

hey

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

To constrain the volume of the deformed mesh, we define $\lambda^{min}\_i$ and $\lambda^{max}_i$ such that:

$$
\lambda^{min}_i \leq \lambda_i \leq \lambda^{max}_i
$$

This gives a set of linear inequalities to account for while solving. To avoid rank deficiency in the optimization, we also constrain the value of $\lambda_i$ for a vertex $i$. This gives a linear equality to account for while solving.

To efficiently solve this problem, the paper outlines how to vectorize this optimization. Let $D_A$ and $D_w$ be diagonal matrices with the areas $A^0_i$ and weights $w_i$ along their respective diagonals. Let
$\tilde{L^0} = L^0 \bigotimes I_3$, let $D_{\hat{V}}$ be the diagonal matrix with row-wise stacked elements of $\hat{V}$, let $S = I_n \bigotimes [1, 1, 1]^T$ where $n$ is the number of vertices, and define $L_{\theta}$ as follows:

$$
L_{\theta} = D^{-1}\_{S \lambda^0} \tilde{L^0} D_{\hat{V}} S \lambda^0
$$

Finally, we have the following vectorized optimization:

$$
argmin\_{\lambda} \hspace{10pt} || D\_A D\_w (\tilde{L^0} D\_{\hat{V}} - D\_{L\_{theta}}) S \lambda ||^2
$$

subject to the constraints we outlined above.

As such, the majority of this project's implementation is constructing these matrices and manipulating them into a form that can be fed into libigl's quadratic program solver, `active_set`.

## Future Improvements

This implementation lacks two main features from the paper description. The first is the lack of z-ordering constraints; as a result, it is possible and likely for the mesh to contain self-intersections. In the future,
these z-ordering constraints can be implemented to preserve the order of vertices with respect to the viewpoint. The second is the lack of support for projecting onto disconnected pieces. Note that, although we did not
include the $mu$ values in the brief derivation above because it does not function correctly with this implementation, much of the implementation has been written as if the feature were present; as a result, this feature
could be added in the future without having to significantly refactor the existing code.

Finally, as previously mentioned, the implementation use's the libigl function `active_set` to solve the quadratic problem. The paper outlines a process to convert this optimization to a linear energy minimization with conic
constraints to improve performance. To my knowledge, libigl does not currently contain a function to solve this kind of optimization, so other tools such as Mosek would need to be used if this improvement were to be implemented.

## References
Christian Schüller, Daniele Panozzo, and Olga Sorkine-Hornung. 2014. Appearance-mimicking surfaces. ACM Trans. Graph. 33, 6, Article 216 (November 2014), 10 pages. https://doi.org/10.1145/2661229.2661267

