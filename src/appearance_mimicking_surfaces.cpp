#include "../include/appearance_mimicking_surfaces.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/active_set.h>
#include <igl/cat.h>

void get_voronoi_area_matrix(
  Eigen::SparseMatrix<double>& M,
  Eigen::DiagonalMatrix<double, Eigen::Dynamic>& M_diag)
{
  // Square root each area entry and repeat it three consecutive times along the diagonal
  Eigen::VectorXd diagonal, resized_diagonal;
  diagonal = M.diagonal();
  resized_diagonal.conservativeResize(3 * M.rows());
  for (int i = 0; i < M.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      resized_diagonal(3 * i + j) = sqrt(diagonal(i));
    }
  }
  
  // Return diagonal matrix with repeated area entries along the diagonal
  M_diag.resize(3 * M.rows());
  M_diag.diagonal() = resized_diagonal;
}

void get_weights_matrix(
  const Eigen::VectorXd& W,
  Eigen::DiagonalMatrix<double, Eigen::Dynamic>& W_diag)
{
  // Repeat each area entry three times
  Eigen::VectorXd resized_diagonal;
  resized_diagonal.conservativeResize(3 * W.rows());
  for (int i = 0; i < W.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      resized_diagonal(3 * i + j) = W(i);
    }
  }

  // Return diagonal matrix with repeated weight entries along the diagonal
  W_diag.resize(3 * W.rows());
  W_diag.diagonal() = resized_diagonal;
}

void get_laplace_beltrami_kronecker(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::SparseMatrix<double>& M,
  Eigen::SparseMatrix<double>& LB_kronecker)
{
  // Get cotangent Laplace-Beltrami operator (#V x #V sparse)
  Eigen::SparseMatrix<double> L, LB, Minv;
  igl::cotmatrix(V, F, L);
  igl::invert_diag(M, Minv);
  LB = Minv * L;

  // Return Kronecker product of Laplace-Beltrami matrix and 3-dimensional Identity matrix (3#V x 3#V sparse)
  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve(3 * LB.rows() * LB.rows());
  for (int i = 0; i < LB.rows(); i++) {
    for (int j = 0; j < LB.rows(); j++) {
      double LB_entry = LB.coeff(i, j);
      for (int k = 0; k < 3; k++) {
        tripletList.push_back(Eigen::Triplet<double>(i * 3 + k, j * 3 + k, LB_entry));
      }
    }
  }

  LB_kronecker.conservativeResize(3 * LB.rows(), 3 * LB.rows());
  LB_kronecker.setFromTriplets(tripletList.begin(), tripletList.end());
}

void get_unit_directions(
  const Eigen::MatrixXd& V,
  const Eigen::Vector3d& view,
  Eigen::MatrixXd& V_unit)
{
  // Compute the normalized direction vector from the viewpoint to each vector
  V_unit.conservativeResize(V.rows(), 3);
  for (int i = 0; i < V.rows(); i++) {
    V_unit.row(i) = (V.row(i) - view.transpose()).normalized();
  }
}

void get_diagonal_stacked_unit_directions(
  Eigen::MatrixXd& V_unit,
  Eigen::DiagonalMatrix<double, Eigen::Dynamic>& V_dir_diag)
{
  // Stack the normalized direction vectors along the diagonal
  Eigen::VectorXd resized_diagonal;
  resized_diagonal.conservativeResize(3 * V_unit.rows());
  for (int i = 0; i < V_unit.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      resized_diagonal(3 * i + j) = V_unit(i, j);
    }
  }

  // Return diagonal matrix with unit direction vector components along the diagonal
  V_dir_diag.resize(3 * V_unit.rows());
  V_dir_diag.diagonal() = resized_diagonal;
}

void get_selector_matrix(
  int num_vertices,
  Eigen::SparseMatrix<double>& S)
{
  // The selector matrix is equivalent to the Kronecker product of the n-dimensional Identity matrx and the transpose of [1, 1, 1] (3#V x #V sparse matrix)
  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve(3 * num_vertices);
  for (int i = 0; i < num_vertices; i++) {
    for (int j = 0; j < 3; j++) {
      tripletList.push_back(Eigen::Triplet<double>(i * 3 + j, i, 1));
    }
  }

  S.conservativeResize(3 * num_vertices, num_vertices);
  S.setFromTriplets(tripletList.begin(), tripletList.end());
}

void get_distances_from_view(
  const Eigen::MatrixXd& V,
  const Eigen::Vector3d& view,
  Eigen::VectorXd& lambda
)
{
  lambda.conservativeResize(V.rows());
  for (int i = 0; i < V.rows(); i++) {
    lambda[i] = (V.row(i) - view.transpose()).norm();
  }
}

void get_D_L_theta(
  Eigen::SparseMatrix<double>& S,
  Eigen::VectorXd& lambda_0,
  Eigen::SparseMatrix<double>& L0_tilda,
  Eigen::DiagonalMatrix<double, Eigen::Dynamic>& D_V,
  Eigen::DiagonalMatrix<double, Eigen::Dynamic>& L_theta_diag
)
{
  Eigen::VectorXd Slambda = S * lambda_0;
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> Slambda_diag(Slambda.rows());
  Slambda_diag.diagonal() = Slambda;

  Eigen::VectorXd L_theta = Slambda_diag.inverse() * L0_tilda * D_V * Slambda;

  L_theta_diag.resize(Slambda.rows());
  L_theta_diag.diagonal() = L_theta;
}

void get_fixed_vertex_constraints(
  const Eigen::MatrixXd& V,
  const Eigen::VectorXi& bf,
  Eigen::VectorXd& bc
)
{
  // To eliminate rank deficiency in optimization, we fix lambda for a single vertex in each disconnected group
  // Set these lambdas to the average of lambda_min and lambda_max for that group
  bc.conservativeResize(bf.rows());
  for (int i = 0; i < bf.rows(); i++) {
    bc(i) = V.row(bf(i)).norm();
  }
}

void get_inequality_constraints(
  const Eigen::MatrixXd& lambdaMin,
  const Eigen::MatrixXd& lambdaMax,
  const Eigen::VectorXd& mu,
  Eigen::SparseMatrix<double>& Aieq,
  Eigen::VectorXd& Bieq
)
{
  // Construct inequality constraints
  // We have mu_g lambda_min_i <= lambda_i <= mu_g lambda_max_i. This constraint can be split into the following two constraints:
  // mu_g lambda_min_i - lambda_i <= 0 and lambda_i - mu_g lambda_max_i <= 0
  // The first #lambdaMin rows will cover the former type of constraint, and the last #lambdaMax rows will cover the latter type of constraint
  // Recall that x in this case is [lambda mu]^T, so the matrix will be (#lambdaMin + #lambdaMax) x (#lambda + #mu)
  int vertex_count = mu.rows();
  int mu_count = mu.maxCoeff() + 1;
  int row_count = lambdaMin.rows() + lambdaMax.rows();

  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve(2 * row_count);

  // Upper half of matrix: mu_g lambda_min_i - lambda_i <= 0 constraints
  for (int i = 0; i < vertex_count; i++) {
    tripletList.push_back(Eigen::Triplet<double>(i, i, -1)); // selects -lambda_i
    tripletList.push_back(Eigen::Triplet<double>(i, mu(i) + vertex_count, lambdaMin(i, 1))); // selects mu_g lambda_min_i
  }

  // Lower half of matrix: lambda_i - mu_g lambda_max_i <= 0 constraints
  for (int i = 0; i < vertex_count; i++) {
    tripletList.push_back(Eigen::Triplet<double>(i + vertex_count, i, 1)); // selects lambda_i
    tripletList.push_back(Eigen::Triplet<double>(i + vertex_count, mu(i) + vertex_count, -lambdaMax(i, 1))); // selects -(mu_g lambda_max_i)
  }

  Aieq.conservativeResize(row_count, vertex_count + mu_count);
  Aieq.setFromTriplets(tripletList.begin(), tripletList.end());

  // Right hand side of all inequalities is 0
  Bieq = Eigen::VectorXd::Zero(row_count);
}

void appearance_mimicking_surfaces(
        const Eigen::MatrixXd &V, 
        const Eigen::MatrixXi &F, 
        const Eigen::Vector3d &view, 
        const Eigen::MatrixXd &lambdaMin, 
		    const Eigen::MatrixXd &lambdaMax, 
        const Eigen::VectorXi &bf, 
        const Eigen::VectorXd &weights, 
        const Eigen::VectorXd &mu, 
        Eigen::MatrixXd &DV)
{
  std::cout << "Constructing various matrices for F..." << "\n";

  // Get mass matrix
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_VORONOI, M);

  // Get D_A
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_A;
  get_voronoi_area_matrix(M, D_A);

  // Get D_W - Save as a SparseMatrix for later computation
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_W;
  get_weights_matrix(weights, D_W);
  Eigen::SparseMatrix<double> D_W_sparse = Eigen::SparseMatrix<double>(D_W);

  // Get L0_tilda (3#V x 3#V sparse matrix)
  Eigen::SparseMatrix<double> L0_tilda;
  get_laplace_beltrami_kronecker(V, F, M, L0_tilda);

  // Get V_unit (#V x 3 matrix)
  Eigen::MatrixXd V_unit;
  get_unit_directions(V, view, V_unit);

  // Get D_V (3#V x 3#V diagonal matrix)
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_V;
  get_diagonal_stacked_unit_directions(V_unit, D_V);

  // Get S (3#V x #V sparse matrix)
  Eigen::SparseMatrix<double> S;
  get_selector_matrix(V.rows(), S);

  // Get lambda_0 (#V vector)
  Eigen::VectorXd lambda_0;
  get_distances_from_view(V, view, lambda_0);

  // Get D_L_theta (3#V x 3#V diagonal matrix) - Save as a SparseMatrix for later computation
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_L_theta;
  get_D_L_theta(S, lambda_0, L0_tilda, D_V, D_L_theta);
  Eigen::SparseMatrix<double> D_L_theta_sparse = Eigen::SparseMatrix<double>(D_L_theta);

  // Get number of mu groups; assume mu starts indexing at 0
  int mu_count = mu.maxCoeff() + 1;

  std::cout << "Constructing F and constraint matrices..." << "\n";

  // Construct F
  Eigen::SparseMatrix<double> E;
  E = D_A * D_W_sparse * ((L0_tilda * D_V) - D_L_theta_sparse) * S;
  double alpha = sqrt(pow(10, -7));
  Eigen::SparseMatrix<double> R = (alpha * Eigen::MatrixXd::Identity(mu_count, mu_count)).sparseView();

  Eigen::SparseMatrix<double> Q, Q2, Q_left, Q_right, zeros_1, zeros_2;
  zeros_1 = Eigen::SparseMatrix<double>(mu_count, E.cols());
  zeros_2 = Eigen::SparseMatrix<double>(E.rows(), mu_count);
  igl::cat(1, E, zeros_1, Q_left);
  igl::cat(1, zeros_2, R, Q_right);
  igl::cat(2, Q_left, Q_right, Q);
  Q2 = Q.transpose() * Q;

  // No linear coefficients, set to 0
  Eigen::VectorXd B = Eigen::VectorXd::Zero(V.rows() + mu_count);

  // Fixed vertex constraints
  Eigen::VectorXd bc;
  get_fixed_vertex_constraints(V, bf, bc);

  // No other linear equality constraints, so leave Aeq, Beq empty
  Eigen::SparseMatrix<double> Aeq;
  Eigen::VectorXd Beq;

  // Linear inequality constraints
  Eigen::SparseMatrix<double> Aieq;
  Eigen::VectorXd Bieq;
  get_inequality_constraints(lambdaMin, lambdaMax, mu, Aieq, Bieq);

  std::cout << "Solving for lambda and mu..." << "\n";

  // Solve
  Eigen::VectorXd Z, lx, ux;
  igl::active_set_params as;
  // as.max_iter = 10;
  igl::active_set(Q2, B, bf, bc, Aeq, Beq, Aieq, Bieq, lx, ux, as, Z);

  std::cout << "Finished solving." << "\n";

  // Solution is equivalent to [lambda, mu]^T, so extract lambda and u
  Eigen::VectorXd lambda_sol = Eigen::VectorXd(V.rows());
  for (int i = 0; i < V.rows(); i++) {
    lambda_sol(i) = Z(i);
  }

  Eigen::VectorXd mu_sol = Eigen::VectorXd(mu_count);
  for (int i = 0; i < mu_count; i++) {
    mu_sol(i) = Z(V.rows() + i);
  }

  // Apply lambda and mu to vertices
  // Recall that V_unit has stores a unit vector from the view to each vertex on the mesh
  DV.conservativeResize(V.rows(), 3);

  for (int i = 0; i < V.rows(); i++) {
    // Recall mu(i) gives the index of the mu group that vertex i belongs to. So, mu_sol(mu(i)) is the correct multiplier
    DV.row(i) = (V_unit.row(i) * (lambda_sol(i) * (1 / mu_sol(mu(i))))) + view.transpose();
  }

  return;
}
