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

const double sqrt_3 = std::sqrt(3);

struct Box3 {
    double xmin=std::numeric_limits<double>::max(),
           ymin=std::numeric_limits<double>::max(),
           zmin=std::numeric_limits<double>::max(),
           xmax=std::numeric_limits<double>::min(),
           ymax=std::numeric_limits<double>::min(),
           zmax=std::numeric_limits<double>::min();

    double step;

    long long x_idx_min,
              y_idx_min,
              z_idx_min,
              x_idx_max,
              y_idx_max,
              z_idx_max,
              size_x,
              size_y,
              size_z;



    Box3 (const std::array<Point_3, 8>& vertices, double step) : step(step) {
        for(const Point_3& p : vertices) {
            if(p.x() < xmin) xmin = p.x();
            if(p.y() < xmin) ymin = p.y();
            if(p.z() < xmin) zmin = p.z();
            if(p.x() > xmax) xmax = p.x();
            if(p.y() > ymax) ymax = p.y();
            if(p.z() > zmax) zmax = p.z();
        }

        x_idx_min = (xmin / step);
        y_idx_min = (ymin / step);
        z_idx_min = (zmin / step);
        x_idx_max = (xmax / step);
        y_idx_max = (ymax / step);
        z_idx_max = (zmax / step);

        size_x = x_idx_max - x_idx_min;
        size_y = y_idx_max - y_idx_min;
        size_z = z_idx_max - z_idx_min;
    }

    Box3(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double step) :
        xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), zmin(zmin), zmax(zmax), step(step){
        x_idx_min = (xmin / step);
        y_idx_min = (ymin / step);
        z_idx_min = (zmin / step);
        x_idx_max = (xmax / step);
        y_idx_max = (ymax / step);
        z_idx_max = (zmax / step);
        size_x = x_idx_max - x_idx_min;
        size_y = y_idx_max - y_idx_min;
        size_z = z_idx_max - z_idx_min;
    }

    friend std::ostream &operator<<(std::ostream &os, const Box3 &box3) {
        os << "xmin: " << box3.xmin << " ymin: " << box3.ymin << " zmin: " << box3.zmin << " xmax: " << box3.xmax
           << " ymax: " << box3.ymax << " zmax: " << box3.zmax;
        return os;
    }

    long long get_number_of_cells() {
        return (x_idx_max - x_idx_min) * (y_idx_max - y_idx_min) * (z_idx_max - z_idx_min);
    }
};

Point_3 move_vertex(const Point_3& vertex, const Vector_3& force, double mu);

inline std::tuple<int, int, int> cell_to_global_position(int msx, int msy, int msz, long long index){
    auto gidx = index % msx,
         gidy = (long long) std::floor(index % (msx*msy) / msx),
         gidz = (long long) std::floor(index / (msx*msy));
    return std::make_tuple(gidx, gidy, gidz);
}

inline std::tuple<int, int, int> cell_to_local_position(int msx, int msy, int msz, Box3 bbox, long long index){
    auto minx = bbox.x_idx_min,
         maxx = bbox.x_idx_max,
         miny = bbox.y_idx_min,
         maxy = bbox.y_idx_max,
         minz = bbox.z_idx_min,
         maxz = bbox.z_idx_max;

    int gidx =  index % msx,
        gidy = (int) index % (msx*msy) / msx,
        gidz = (int) index / (msx*msy);

    return {gidx - minx,  gidy - miny, gidz - minz};
}

inline void cell_to_local_position(int msx, int msy, int msz, Box3 bbox, long long index, long long* x_idx, long long* y_idx, long long* z_idx){
    auto minx = bbox.x_idx_min,
         maxx = bbox.x_idx_max,
         miny = bbox.y_idx_min,
         maxy = bbox.y_idx_max,
         minz = bbox.z_idx_min,
         maxz = bbox.z_idx_max;

    auto gidx = index % msx,
         gidy = (long long) std::floor(index % (msx*msy) / msx),
         gidz = (long long) std::floor(index / (msx*msy));

    *x_idx = gidx - minx;
    *y_idx = gidy - miny;
    *z_idx = gidz - minz;
}

template<class A>
inline void linear_to_grid(const A index, const A c, const A r, A& x_idx, A& y_idx, A& z_idx){
    x_idx = (A) (index % c);           // col
    y_idx = (A) std::floor(index % (c*r) / c); // row
    z_idx = (A) std::floor(index / (c*r));     // depth
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

template<class NumericalType, class IndexType>
inline int position_to_cell(IndexType x, IndexType y, IndexType z, const double step, const NumericalType column, const NumericalType row,
                            const NumericalType col_shift  = 0.0,
                            const NumericalType row_shift  = 0.0,
                            const NumericalType depth_shift= 0.0) {
    const std::vector<NumericalType> weight = {1, column, column*row};
    NumericalType idx = 0;

    idx += (NumericalType) std::floor(x / step);
    idx += column * (NumericalType) std::floor(y / step);
    idx += column * row * (NumericalType) std::floor(z / step);

    return idx;
}

template<class IndexType>
inline IndexType grid_index_to_cell(IndexType x, IndexType y, IndexType z, const IndexType column, const IndexType row, const IndexType depth) {
    const std::vector<IndexType> weight = {1, column, column*row};
    IndexType idx = 0;

    idx += x;
    idx += column * y;
    idx += column * row * z;

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


inline double lb_getTetraederVolumeIndexed(int c1, int c2, int c3, int c4, const std::array<Point_3, 8>& vertices) {
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

inline bool in_tetrahedron(const Tetrahedron_3& tet, const Point_3& p){
    auto side = tet.bounded_side(p);
    return side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY;
}

std::array<Plane_3, 12> get_planes(const std::array<Point_3, 8>& vertices);

std::array<Tetrahedron_3, 6> get_tetrahedra(const std::array<Point_3, 8>& vertices);

inline Point_3 operator*(const double w, const Point_3& p){
    return Point_3(w*p.x(), w*p.y(), w*p.z());
}

inline Point_3 operator+(const Point_3& p1, const Point_3& p2){
    return Point_3(p1.x()+p2.x(), p1.y()+p2.y(), p1.z()+p2.z());
}

Point_3 get_center_of_load(const std::vector<double>& weights, const std::vector<Point_3>& elements);

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

template<class A>
bool is_there_duplicates(const std::vector<A>& els){
    const auto sz = els.size();
    for(int i = 0; i < sz; ++i){
        for(int j = i+1; j < sz; ++j){
            if(els[i] == els[j]) {
                std::cout << els[i] << " " << els[j] << std::endl;
                return true;
            }
        }
    }
    return false;
}

template<class A>
bool exists(const A& el, const std::vector<A>& els){
    const auto sz = els.size();
    for(int i = 0; i < sz; ++i){
        if(els[i] == el) return true;
    }
    return false;
}

}
#endif //ADLBIRREG_GEOMETRICUTILS_HPP
