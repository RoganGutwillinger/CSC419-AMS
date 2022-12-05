#include "appearance_mimicking_surfaces.h"
#include <igl/cat.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>

void get_AABB_corners(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::RowVector3d& min_corner,
  Eigen::RowVector3d& max_corner) 
{
  min_corner = V.row(0);
  max_corner = V.row(0);
  for (int i = 0; i < V.rows(); i++) {
    Eigen::RowVector3d v = V.row(i);
    min_corner(0) = std::min(min_corner(0), v(0));
    max_corner(0) = std::max(max_corner(0), v(0));
    min_corner(1) = std::min(min_corner(1), v(1));
    max_corner(1) = std::max(max_corner(1), v(1));
    min_corner(2) = std::min(min_corner(2), v(2));
    max_corner(2) = std::max(max_corner(2), v(2));
  }
}

bool plane_intersect(
  Eigen::Vector3d view,
  Eigen::Vector3d unit_dir,
  Eigen::Vector3d plane_point,
  Eigen::Vector3d plane_normal,
  const double min_t, 
  double& t)
{
  double denom = (plane_normal.dot(unit_dir));
  if (denom == 0) { //Checks for 0 in denominator
    return false;
  }

  // Solve for t
  t = plane_normal.dot(plane_point - view) / denom;
  return t >= min_t;
}

bool ray_intersect_triangle(
  Eigen::Vector3d view,
  Eigen::Vector3d unit_dir,
  const Eigen::RowVector3d& A,
  const Eigen::RowVector3d& B,
  const Eigen::RowVector3d& C,
  const double min_t,
  const double max_t,
  double& t)
{
  //2 edges and point to define triangle
  Eigen::RowVector3d edge_1 = B - A;
  Eigen::RowVector3d edge_2 = C - A;

  //Calculate coefficients
  Eigen::Matrix3d m;
  m.col(0) = edge_1.transpose();
  m.col(1) = edge_2.transpose();
  m.col(2) = -unit_dir;
  Eigen::Vector3d coef_sol = m.inverse() * (view - A.transpose());

  double a = coef_sol(0);
  double b = coef_sol(1);
  t = coef_sol(2);

  //Verify constraints
  if (a < 0 || b < 0 || a + b > 1 || t < min_t || t > max_t) {
    return false;
  }

  //Valid intersection
  return true;
}

bool ray_intersect_triangle_mesh(
  Eigen::Vector3d view,
  Eigen::Vector3d unit_dir,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const double min_t,
  const double max_t,
  double& hit_t,
  int& hit_f)
{
  hit_t = std::numeric_limits<double>::infinity();
  int rows = F.rows();
  int cols = F.cols();

  for (int f = 0; f < rows; f++) {
    // Retrieve vertices for this triangle
    Eigen::RowVector3d A = V.row(F(f, 0));
    Eigen::RowVector3d B = V.row(F(f, 1));
    Eigen::RowVector3d C = V.row(F(f, 2));

    double intersect_t;
    if (ray_intersect_triangle(view, unit_dir, A, B, C, min_t, max_t, intersect_t)) {
      //We hit a triange; compare t values to save closest hit
      if (intersect_t < hit_t) {
        hit_t = intersect_t;
        hit_f = f;
      }
    }
  }

  //We initialized hit_t to infinity, so if it is still infinity we have no hit
  if (hit_t == std::numeric_limits<double>::infinity()) {
    return false;
  }

  return true;
}

void get_constraints(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::RowVector3d view_point,
  const Eigen::Vector3d& plane_min_point,
  const Eigen::Vector3d& plane_max_point,
  const Eigen::Vector3d& plane_unit_normal,
  Eigen::MatrixXd& lambdaMin,
  Eigen::MatrixXd& lambdaMax,
  Eigen::VectorXd& weights)
{
  lambdaMin.conservativeResize(V.rows(), 2);
  lambdaMax.conservativeResize(V.rows(), 2);

  // Weights - Set to 1 by default (will be changed to 0.1 if vertex is not visible)
  weights = Eigen::VectorXd::Ones(V.rows());

  for (int i = 0; i < V.rows(); i++) {
    Eigen::RowVector3d dir = V.row(i) - view_point;
    double norm = dir.norm();
    Eigen::RowVector3d unit_dir = dir / norm;

    // Find intersection with min and max plane
    double t_min, t_max;
    plane_intersect(view_point, unit_dir, plane_min_point, plane_unit_normal, 0, t_min);
    plane_intersect(view_point, unit_dir, plane_max_point, plane_unit_normal, 0, t_max);

    lambdaMin(i, 0) = i;
    lambdaMin(i, 1) = t_min;
    lambdaMax(i, 0) = i;
    lambdaMax(i, 1) = t_max;

    // Determine if vertex is visible
    // If there is an intersection with the mesh before reaching the vertex, the vertex is not visible
    double hit_t;
    int hit_f;
    double hit_t_max = norm - 0.001;
    if (ray_intersect_triangle_mesh(view_point, unit_dir, V, F, 0, hit_t_max, hit_t, hit_f)) {
      // Vertex not visible
      weights(i) = 0.1;
    }
  }
}

int main(int argc, char *argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  // Load input meshes
  igl::read_triangle_mesh(
    (argc > 1 ? argv[1] : "../data/decimated-knight.off"), V, F);

  // Only use one continuous group to project to, so only use one group index (0 in this case)
  Eigen::VectorXd mu = Eigen::VectorXd::Zero(V.rows());

  // Since there is only one projection plane, only need to fix one vertex
  Eigen::VectorXi bf = Eigen::VectorXi(1);
  bf(0) = 0;

  std::cout << "Computing mesh intersections with the plane..." << "\n";

  Eigen::RowVector3d min_corner, max_corner;
  get_AABB_corners(V, F, min_corner, max_corner);
  Eigen::RowVector3d mesh_center = (min_corner + max_corner) / 2;
  double depth_x = max_corner(0) - min_corner(0);
  double depth_z = max_corner(2) - min_corner(2);
  double min_dim = std::min(depth_x, depth_z);
  double max_dim = std::max(depth_x, depth_z);

  double radius = 3 * max_dim;
  double deformed_depth = min_dim / 3;

  // We will project onto the plane that corresponds to the back face of the aabb of the mesh
  const Eigen::RowVector3d plane_unit_normal(0, 0, 1);
  const Eigen::RowVector3d plane_min_point = min_corner + plane_unit_normal * deformed_depth;
  const Eigen::RowVector3d plane_max_point = min_corner;

  // Two viewpoints: One in front of the mesh, one offset to the side
  const Eigen::RowVector3d view_point_front = mesh_center + Eigen::RowVector3d(0, 0, 1) * radius;
  const Eigen::RowVector3d view_point_corner = mesh_center + Eigen::RowVector3d(1 / sqrt(2), 0, 1 / sqrt(2)) * radius;

  Eigen::VectorXd weights_front, weights_corner;
  Eigen::MatrixXd lambdaMin_front, lambdaMax_front, lambdaMin_corner, lambdaMax_corner;
  get_constraints(V, F, view_point_front, plane_min_point, plane_max_point, plane_unit_normal, lambdaMin_front, lambdaMax_front, weights_front);
  get_constraints(V, F, view_point_corner, plane_min_point, plane_max_point, plane_unit_normal, lambdaMin_corner, lambdaMax_corner, weights_corner);

  Eigen::MatrixXd DV_front, DV_corner;

  std::cout << "Running AMS for the first viewpoint..." << "\n";

  appearance_mimicking_surfaces(V, F, view_point_front, lambdaMin_front, lambdaMax_front, bf, weights_front, mu, DV_front);

  std::cout << "Running AMS for the second viewpoint..." << "\n";

  appearance_mimicking_surfaces(V, F, view_point_corner, lambdaMin_corner, lambdaMax_corner, bf, weights_corner, mu, DV_corner);

  // Concatenate a plane to the mesh
  // Purposefully invert normal to give contrast against mesh
  Eigen::MatrixXd V0, DV0_front, DV0_corner;
  Eigen::MatrixXi F0;
  double aabb_height = max_corner(1) - min_corner(1);
  double aabb_width = max_corner(0) - min_corner(0);

  Eigen::MatrixXd V_Plane = Eigen::MatrixXd(4, 3);
  V_Plane.row(0) = min_corner - Eigen::RowVector3d(aabb_width, aabb_height, 0);
  V_Plane.row(1) = V_Plane.row(0) + Eigen::RowVector3d(0, 3 * aabb_height, 0);
  V_Plane.row(2) = V_Plane.row(0) + Eigen::RowVector3d(3 * aabb_width, 3 * aabb_height, 0);
  V_Plane.row(3) = V_Plane.row(0) + Eigen::RowVector3d(3 * aabb_width, 0, 0);
  Eigen::MatrixXi F_Plane = Eigen::MatrixXi(2, 3);
  F_Plane.row(0) = Eigen::RowVector3i(V.rows(), V.rows() + 1, V.rows() + 2);
  F_Plane.row(1) = Eigen::RowVector3i(V.rows() + 2, V.rows() + 3, V.rows());

  igl::cat(1, V, V_Plane, V0);
  igl::cat(1, DV_front, V_Plane, DV0_front);
  igl::cat(1, DV_corner, V_Plane, DV0_corner);
  igl::cat(1, F, F_Plane, F0);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;

  Eigen::MatrixXd DV_to_display = V0;

  const auto update = [&]()
  {
    viewer.data().set_mesh(DV_to_display, F0);
  };

  viewer.callback_key_pressed =
    [&](igl::opengl::glfw::Viewer&, unsigned int key, int mod)
  {
    switch (key)
    {
    case '1':
      DV_to_display = V0;
      break;
    case '2':
      DV_to_display = DV0_front;
      break;
    case '3':
      DV_to_display = DV0_corner;
      break;
    default:
      return false;
    }
    update();
    return true;
  };

  viewer.data().set_mesh(DV_to_display, F0);
  viewer.data().set_face_based(true);
     
  // Draw points
  viewer.data().add_points(view_point_front, Eigen::RowVector3d(0, 1, 0));
  viewer.data().add_points(view_point_corner, Eigen::RowVector3d(0, 0, 1));

  update();

  std::cout << R"(
    Viewpoint 1 is the green point, Viewpoint 2 is the blue point
    1  Show base mesh
    2  Show deformed mesh for Viewpoint 1
    3  Show deformed mesh for Viewpoint 2
  )";

  viewer.launch();
}