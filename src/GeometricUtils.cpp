//
// Created by xetql on 10/21/19.
//
#include <GeometricUtils.hpp>
namespace lb{
Vector_3 constraint_force(const Box3& d, const Point_3& p, const Vector_3& f){
    const auto x = p.x();
    const auto y = p.y();
    const auto z = p.z();
    double fx = f.x(), fy = f.y(), fz = f.z();

    if(x <= d.xmin+std::numeric_limits<float>::epsilon()) {
        fx = 0.0;
    }
    if(y <= d.ymin+std::numeric_limits<float>::epsilon()) {
        fy = 0.0;
    }
    if(z <= d.zmin+std::numeric_limits<float>::epsilon()) {
        fz = 0.0;
    }
    if(x >= d.xmax-std::numeric_limits<float>::epsilon()) {
        fx = 0.0;
    }
    if(y >= d.ymax-std::numeric_limits<float>::epsilon()) {
        fy = 0.0;
    }
    if(z >= d.zmax-std::numeric_limits<float>::epsilon()) {
        fz = 0.0;
    }
    return Vector_3(fx, fy, fz);
}
double lb_getTetraederVolumeIndexed(int c1, int c2, int c3, int c4, const std::array<Point_3, 8>& vertices) {
    double dir1_0, dir1_1, dir1_2;
    double dir2_0, dir2_1, dir2_2;
    double dir3_0, dir3_1, dir3_2;

    dir1_0 = vertices[c2].x() - vertices[c1].x();
    dir1_1 = vertices[c2].y() - vertices[c1].y();
    dir1_2 = vertices[c2].z() - vertices[c1].z();

    dir2_0 = vertices[c3].x() - vertices[c1].x();
    dir2_1 = vertices[c3].y() - vertices[c1].y();
    dir2_2 = vertices[c3].z() - vertices[c1].z();

    dir3_0 = vertices[c4].x() - vertices[c1].x();
    dir3_1 = vertices[c4].y() - vertices[c1].y();
    dir3_2 = vertices[c4].z() - vertices[c1].z();

    return (dir1_0 * (dir3_1 * dir2_2 - dir3_2 * dir2_1) +
            dir1_1 * (dir3_2 * dir2_0 - dir3_0 * dir2_2) +
            dir1_2 * (dir3_0 * dir2_1 - dir3_1 * dir2_0));
}
bool in_tetrahedron(const Tetrahedron_3& tet, const Point_3& p){
    auto side = tet.bounded_side(p);
    return side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY;
}
std::array<Plane_3, 12> get_planes(const std::array<Point_3, 8>& vertices){
    return {
            //front face
            Plane_3(vertices[1], vertices[4], vertices[0]),
            Plane_3(vertices[1], vertices[5], vertices[4]),
            //connection to next face
            Plane_3(vertices[2], vertices[1], vertices[0]),
            Plane_3(vertices[2], vertices[3], vertices[1]),
            Plane_3(vertices[3], vertices[7], vertices[5]),
            Plane_3(vertices[3], vertices[5], vertices[1]),

            Plane_3(vertices[4], vertices[2], vertices[0]),
            Plane_3(vertices[4], vertices[6], vertices[2]),
            Plane_3(vertices[5], vertices[7], vertices[6]),
            Plane_3(vertices[5], vertices[6], vertices[4]),
            //back face
            Plane_3(vertices[6], vertices[3], vertices[2]),
            Plane_3(vertices[6], vertices[7], vertices[3]),
    };
}
std::array<Tetrahedron_3, 6> get_tetrahedra(const std::array<Point_3, 8>& vertices) {
    return {
            Tetrahedron_3(vertices[0], vertices[3], vertices[1], vertices[7]),
            Tetrahedron_3(vertices[0], vertices[1], vertices[5], vertices[7]),
            Tetrahedron_3(vertices[0], vertices[2], vertices[3], vertices[7]),
            Tetrahedron_3(vertices[0], vertices[6], vertices[2], vertices[7]),
            Tetrahedron_3(vertices[0], vertices[5], vertices[4], vertices[7]),
            Tetrahedron_3(vertices[0], vertices[4], vertices[6], vertices[7])
    };
}
Point_3 get_center_of_load(const std::vector<double>& weights, const std::vector<Point_3>& elements) {

    const double total_mass = std::accumulate(weights.cbegin(), weights.cend(), 0.0);
    const auto size = elements.size();

    double xcm = 0.0; double ycm = 0.0; double zcm = 0.0;
    for(int i = 0; i < size; ++i) {
        xcm += (weights[i] * elements[i].x());
        ycm += (weights[i] * elements[i].y());
        zcm += (weights[i] * elements[i].z());
    }

    return Point_3(xcm / total_mass, ycm / total_mass, zcm / total_mass);
}
double compute_normalized_load(double my_load, double unit) {
    return my_load / unit;
}
Vector_3 get_vertex_force(const Point_3& vertex, const std::vector<Point_3>& centers_of_load, const std::vector<double>& normalized_loads) {
    Vector_3 f(0,0,0);
    const auto size = centers_of_load.size();
    for(int i = 0; i < size; ++i) {
        f += (normalized_loads[i] - 1.0) * get_direction(vertex, centers_of_load[i]);
    }
    return f;
}

bool lb_isGeometryValid(const std::array<Point_3, 8>& vertices, const std::array<Plane_3, 12>& planes, const double grid_size){
    /* In comparision to the set of rules in the plimpton scheme, these values are less strict
         * Except for the tetrahedron volumes, the checks are not necessary, but prevent the domain
         * from degenerating to severely, which can cause problems in the convergence behavior.
         * The movement of corners is otherwise likely to get stuck in local minima.
         */
    /* Tetrahedral subvolumes in the domain must be positively oriented */
    /* Self-intersecting cubes are bad, very bad*/
    /* Test all four permutations of how the cube can be split into tetrahedrons*/
    double v;
    if (lb_getTetraederVolumeIndexed(0, 5, 4, 7, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(0, 3, 1, 7, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(0, 1, 5, 7, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(0, 4, 6, 7, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(0, 6, 2, 7, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(0, 2, 3, 7, vertices) <= 0) return false;

    if (lb_getTetraederVolumeIndexed(1, 7, 5, 6, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(1, 2, 3, 6, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(1, 3, 7, 6, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(1, 5, 4, 6, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(1, 4, 0, 6, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(1, 0, 2, 6, vertices) <= 0) return false;

    if (lb_getTetraederVolumeIndexed(2, 4, 6, 5, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(2, 1, 0, 5, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(2, 0, 4, 5, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(2, 6, 7, 5, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(2, 7, 3, 5, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(2, 3, 1, 5, vertices) <= 0) return false;

    if (lb_getTetraederVolumeIndexed(3, 6, 7, 4, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(3, 0, 2, 4, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(3, 2, 6, 4, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(3, 7, 5, 4, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(3, 5, 1, 4, vertices) <= 0) return false;
    if (lb_getTetraederVolumeIndexed(3, 1, 0, 4, vertices) <= 0) return false;

    //Additionally prevent the collapse of the corners in the domain
    //This would yield a topological different domain geometry
    if (lb_getTetraederVolumeIndexed(0, 1, 4, 2, vertices) < 0) return false;
    if (lb_getTetraederVolumeIndexed(1, 5, 4, 7, vertices) < 0) return false;
    if (lb_getTetraederVolumeIndexed(2, 4, 6, 7, vertices) < 0) return false;
    if (lb_getTetraederVolumeIndexed(2, 7, 3, 1, vertices) < 0) return false;
    if (lb_getTetraederVolumeIndexed(0, 4, 6, 5, vertices) < 0) return false;
    if (lb_getTetraederVolumeIndexed(0, 1, 5, 3, vertices) < 0) return false;
    if (lb_getTetraederVolumeIndexed(0, 6, 2, 3, vertices) < 0) return false;
    if (lb_getTetraederVolumeIndexed(5, 3, 7, 6, vertices) < 0) return false;

    std::array<std::array<Plane_3, 6>, 8> incident_planes;

    for(int i = 0; i < 8; ++i){
        auto vertex = vertices[i];
        std::copy_if(planes.begin(), planes.end(), incident_planes[i].begin(), [&vertex](auto plane){return plane.has_on(vertex);});
    }

    for(int i = 0; i < 7; ++i){
        // dist from i to incident(7-i)
        for(auto const& plane: incident_planes[7-i]) {
            if( CGAL::squared_distance(vertices[i], plane) <= std::pow(sqrt_3*grid_size, 2)) return false;
        }
    }

    return true;
}
}