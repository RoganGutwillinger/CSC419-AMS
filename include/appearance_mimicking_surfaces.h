#include <Eigen/Core>

// Given a mesh, a viewpoint, and linear constraints:
//      lambdaMin -- combined with lambdaMax this will "squeeze" the mesh into the defined thickness by projecting the vertices
//      lambdaMax
//      mu -- different parts of the shape can be given different thickness by segmenting the vertices into groups
//      fixed vertices
//      weights -- scales the difference of vertex normals so that more visible vertices are given preserved "more" than less visible
//      optional depth ordering
// output a bas-relief deformed mesh whose appearance is preserved 
// from the viewpoint
//
// 
// Inputs:
//    V            #V by 3 list of the vertex positions of the model
//    F            #F by 3 list of triangle indices into V
//    view         3D vector of the coordinates of the viewpoint
//    lambdaMin    #Constraints by 2 list of vertex indices and the lambdaMinValue: [vertexIdx,lambdaMinValue]
//    lambdaMax    #Constraints by 2 list of [vertexIdx,lambdaMaxValue]
//    bf           #fixedVertices vector of indices into V
//    weights      #V length list of vertex weights
//    mu           #V length list of mu indices 
//                      (separates mesh into independent regions each with their own thickness constraint)
//                      (default should be all 1's)
//  Outputs:
//    DV           #V by 3 list of the vertex positions of the deformed model

void appearance_mimicking_surfaces(
        const Eigen::MatrixXd &V, 
        const Eigen::MatrixXi &F, 
        const Eigen::Vector3d &view, 
        const Eigen::MatrixXd &lambdaMin, 
		const Eigen::MatrixXd &lambdaMax, 
        const Eigen::VectorXi &bf, 
        const Eigen::VectorXd &weights, 
        const Eigen::VectorXd &mu, 
        Eigen::MatrixXd &DV);