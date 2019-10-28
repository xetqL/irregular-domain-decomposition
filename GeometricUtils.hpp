//
// Created by xetql on 10/21/19.
//

#ifndef ADLBIRREG_GEOMETRICUTILS_HPP
#define ADLBIRREG_GEOMETRICUTILS_HPP

#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/squared_distance_3.h>
#include <ostream>

namespace lb {

using Kernel = CGAL::Cartesian<double> ;
using Point_3 = CGAL::Point_3<Kernel>;
using Plane_3 = CGAL::Plane_3<Kernel>;
using Vector_3 = Kernel::Vector_3;
using Tetrahedron_3 = Kernel::Tetrahedron_3;
using Transformation = CGAL::Aff_transformation_3<Kernel>;

struct Box3 {
    double xmin, ymin, zmin, xmax, ymax, zmax;

    friend std::ostream &operator<<(std::ostream &os, const Box3 &box3) {
        os << "xmin: " << box3.xmin << " ymin: " << box3.ymin << " zmin: " << box3.zmin << " xmax: " << box3.xmax
           << " ymax: " << box3.ymax << " zmax: " << box3.zmax;
        return os;
    }
};

Point_3 move_vertex(const Point_3& vertex, const Vector_3& force, double mu);

inline std::pair<int, int> cell_to_global_position(int msx, int msy, long long position){
    return std::make_pair(position % msx, (int) position / msx);
}

template<class A>
inline void linear_to_grid(const A index, const A c, const A r, int& x_idx, int& y_idx, int& z_idx){
    x_idx = (int) (index % (c*r) % c);           // col
    y_idx = (int) std::floor(index % (c*r) / c); // row
    z_idx = (int) std::floor(index / (c*r));     // depth
    assert(x_idx < c);
    assert(y_idx < r);
};

template<class NumericalType>
inline int position_to_cell(Point_3 const& position, const double step, const NumericalType column, const NumericalType row,
        const NumericalType col_shift  = 0.0,
        const NumericalType row_shift  = 0.0,
        const NumericalType depth_shift= 0.0) {
    const std::vector<NumericalType> weight = {1, column, column*row};
    NumericalType idx = 0;

    idx += (NumericalType) std::floor(position.x() / step);
    idx += column * (NumericalType) std::floor(position.y() / step);
    idx += column * row * (NumericalType) std::floor(position.z() / step);

    return idx;
}

/**
 * Force should be nullified along axis of the whole domain iff vertex is on domain border
 * @param d domain
 * @param p vertex
 * @param f force
 * @return force constrained
 */
Vector_3 constraint_force(const Box3& d, const Point_3& p, const Vector_3& f);


template<class InputPointIterator, class InputPointIdIterator>
std::array<int, 4> get_points_on_plane(InputPointIterator beg_points, InputPointIterator end_points,
                                       InputPointIdIterator beg_id, InputPointIdIterator end_id, const Plane_3& p){
    std::array<int, 4> on_plane = {-1, -1, -1, -1};
    int i = 0;
    for(;beg_points != end_points && beg_id != end_id; beg_points++, beg_id++){
        if(p.has_on(*beg_points)){
            on_plane[i] = *beg_id;
            i++;
        }
    }
    return on_plane;
}

double lb_getTetraederVolumeIndexed(int c1, int c2, int c3, int c4, const std::array<Point_3, 8>& vertices);
const double sqrt_3 = std::sqrt(3);

// Domain = U partition;
bool in_tetrahedron(const Tetrahedron_3& tet, const Point_3& p);

std::array<Plane_3, 12> get_planes(const std::array<Point_3, 8>& vertices);

std::array<Tetrahedron_3, 6> get_tetrahedra(const std::array<Point_3, 8>& vertices);

inline Point_3 operator*(const double w, const Point_3& p){
    return Point_3(w*p.x(), w*p.y(), w*p.z());
}

inline Point_3 operator+(const Point_3& p1, const Point_3& p2){
    return Point_3(p1.x()+p2.x(), p1.y()+p2.y(), p1.z()+p2.z());
}

Point_3 get_center_of_load(const std::vector<double>& weights, const std::vector<Point_3>& elements);

/*
std::vector<Point_3> get_centers_of_load_for_vertex(const Point_3& cl, MPI_Comm neighborhood){
    int N;
    MPI_Comm_size(neighborhood, &N);
    std::vector<Point_3> centers_of_load(N);
    std::vector<double>  _centers_of_load(N*3);
    std::array<double,3> my_cl = {cl.x(), cl.y(), cl.z()};
    MPI_Alltoall(my_cl.data(), 3, MPI_DOUBLE, _centers_of_load.data(), 3, MPI_DOUBLE, neighborhood);
    for(int i = 0; i < N; ++i) {
        centers_of_load[i] = Point_3(_centers_of_load[i*3],_centers_of_load[i*3+1],_centers_of_load[i*3+2]);
    }
    return centers_of_load;
}
 */

double compute_normalized_load(double my_load, double unit);
inline Vector_3 get_direction(const Point_3& vertex, const Point_3& center_of_load) {
    return (vertex - center_of_load)/std::sqrt(CGAL::squared_distance(vertex, center_of_load));
}

Vector_3 get_vertex_force(const Point_3& vertex, const std::vector<Point_3>& centers_of_load, const std::vector<double>& normalized_loads);

bool lb_isGeometryValid(const std::array<Point_3, 8>& vertices, const std::array<Plane_3, 12>& planes, const double grid_size);

template<class InputIterator>
std::vector<Point_3> get_points(InputIterator beg, InputIterator end){
    std::vector<Point_3> points;
    unsigned int i = 0;
    while(beg+i != end) {
        points.emplace_back(beg+i, beg+i+1, beg+i+2);
        i+=3;
    }
    return points;
}

}
#endif //ADLBIRREG_GEOMETRICUTILS_HPP
