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
#include "Types.hpp"
namespace lb {

using namespace type;

using Kernel = CGAL::Cartesian<Real> ;
using Point_3 = CGAL::Point_3<Kernel>;
using Plane_3 = CGAL::Plane_3<Kernel>;
using Vector_3 = Kernel::Vector_3;
using Tetrahedron_3 = Kernel::Tetrahedron_3;
using Transformation = CGAL::Aff_transformation_3<Kernel>;

const double sqrt_3 = std::sqrt(3);

struct Box3 {
    Real xmin=std::numeric_limits<Real>::max(),xmax=std::numeric_limits<Real>::min(),
         ymin=std::numeric_limits<Real>::max(),ymax=std::numeric_limits<Real>::min(),
         zmin=std::numeric_limits<Real>::max(),zmax=std::numeric_limits<Real>::min();

    Real simsize_x, simsize_y, simsize_z;

    Real step;

    DataIndex x_idx_min, x_idx_max,size_x,
              y_idx_min, y_idx_max,size_y,
              z_idx_min, z_idx_max,size_z;

    Box3 (const std::array<Point_3, 8>& vertices, Real step) : step(step) {
        for(const Point_3& p : vertices) {
            if(p.x() < xmin) xmin = p.x();
            if(p.y() < ymin) ymin = p.y();
            if(p.z() < zmin) zmin = p.z();
            if(p.x() > xmax) xmax = p.x();
            if(p.y() > ymax) ymax = p.y();
            if(p.z() > zmax) zmax = p.z();
        }

        //xmin = std::min(xmin, xmin-step);
        //ymin = std::min(ymin, ymin-step);
        //zmin = std::min(zmin, zmin-step);
        xmax += step;
        ymax += step;
        zmax += step;

        hook_to_grid();
        init();
    }

    Box3(Real xmin, Real xmax, Real ymin, Real ymax, Real zmin, Real zmax, Real step) :
            xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), zmin(zmin), zmax(zmax), step(step){
        hook_to_grid();
        init();
    }

    DataIndex get_number_of_cells() {
        return (size_x) * (size_y) * (size_z);
    }

    bool contains(Real x, Real y, Real z) {
        return (xmin <= x && x <= xmax) && (ymin <= y && y <= ymax) && (zmin <= z && z <= zmax);
    }

    friend std::ostream &operator<<(std::ostream &os, const Box3 &box3) {
        os << std::setprecision(15) <<  "xmin: " << box3.xmin << " xmax: " << box3.xmax << " ymin: " << box3.ymin << " ymax: " << box3.ymax
           << " zmin: " << box3.zmin << " zmax: " << box3.zmax << " simsize_x: " << box3.simsize_x << " simsize_y: "
           << box3.simsize_y << " simsize_z: " << box3.simsize_z << " step: " << box3.step << " x_idx_min: "
           << box3.x_idx_min << " x_idx_max: " << box3.x_idx_max  << " y_idx_min: "
           << box3.y_idx_min << " y_idx_max: " << box3.y_idx_max << " z_idx_min: "
           << box3.z_idx_min << " z_idx_max: " << box3.z_idx_max << " size_x: " << box3.size_x << " size_y: " << box3.size_y << " size_z: " << box3.size_z;
        return os;
    }

private:
    void init(){
        simsize_x = xmax-xmin;
        simsize_y = ymax-ymin;
        simsize_z = zmax-zmin;

        x_idx_min = (type::DataIndex) (xmin / step);
        y_idx_min = (type::DataIndex) (ymin / step);
        z_idx_min = (type::DataIndex) (zmin / step);

        x_idx_max = (type::DataIndex) (xmax / step);
        y_idx_max = (type::DataIndex) (ymax / step);
        z_idx_max = (type::DataIndex) (zmax / step);

        size_x    = std::floor(simsize_x / step);
        size_y    = std::floor(simsize_y / step);
        size_z    = std::floor(simsize_z / step);
    }

    void hook_to_grid(){
        xmin = std::round(xmin / step) * step;
        ymin = std::round(ymin / step) * step;
        zmin = std::round(zmin / step) * step;

        xmax = std::round(xmax / step) * step;
        ymax = std::round(ymax / step) * step;
        zmax = std::round(zmax / step) * step;
    }
};

Point_3 move_vertex(const Point_3& vertex, const Vector_3& force, Real mu);

inline std::tuple<DataIndex, DataIndex, DataIndex> cell_to_global_position(DataIndex msx, DataIndex msy, DataIndex msz, DataIndex index){
    auto gidx = index % msx,
         gidy = (DataIndex) std::floor(index % (msx*msy) / msx),
         gidz = (DataIndex) std::floor(index / (msx*msy));
    return std::make_tuple(gidx, gidy, gidz);
}

inline std::tuple<DataIndex, DataIndex, DataIndex> cell_to_local_position(DataIndex msx, DataIndex msy, DataIndex msz, Box3 bbox, DataIndex index){
    auto minx = bbox.x_idx_min,
         miny = bbox.y_idx_min,
         minz = bbox.z_idx_min;

    auto gidx =  index % msx,
         gidy = (DataIndex) index % (msx*msy) / msx,
         gidz = (DataIndex) index / (msx*msy);

    return {gidx - minx,  gidy - miny, gidz - minz};
}

inline void cell_to_local_position(DataIndex msx, DataIndex msy, DataIndex msz, const Box3& bbox, DataIndex index, DataIndex* x_idx, DataIndex* y_idx, DataIndex* z_idx){
    auto gidx = index % msx,
         gidy = (DataIndex) std::floor(index % (msx*msy) / msx),
         gidz = (DataIndex) std::floor(index / (msx*msy));
    *x_idx = gidx - bbox.x_idx_min;
    *y_idx = gidy - bbox.y_idx_min;
    *z_idx = gidz - bbox.z_idx_min;
}

inline void linear_to_grid(const DataIndex index, const DataIndex c, const DataIndex r, DataIndex& x_idx, DataIndex& y_idx, DataIndex& z_idx){
    x_idx = (DataIndex) (index % c);                    // col
    y_idx = (DataIndex) std::floor(index % (c*r) / c);  // row
    z_idx = (DataIndex) std::floor(index / (c*r));      // depth
    if(x_idx >= c) throw std::runtime_error("GROS ENCULEUR1");
    if(y_idx >= r) throw std::runtime_error("GROS_ENCULEUR2");
}

template<class NumericalType>
inline NumericalType position_to_cell(Point_3 const& position, const Real step, const NumericalType column, const NumericalType row,
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

template<class IndexType>
inline IndexType position_to_cell(Real x, Real y, Real z, const Real step, const IndexType column, const IndexType row,
                            const Real col_shift  = 0.0, const Real row_shift  = 0.0, const Real depth_shift= 0.0) {
    const std::vector<IndexType> weight = {1, column, column*row};
    IndexType idx = 0;
    idx += (IndexType) std::floor(x / step);
    idx += column * (IndexType) std::floor(y / step);
    idx += column * row * (IndexType) std::floor(z / step);
    return idx;
}

inline std::tuple<DataIndex, DataIndex, DataIndex> position_to_index(const Real x, const Real y, const Real z, const Real step) {
    DataIndex ix, iy, iz;
    ix = (DataIndex) std::floor(x / step);
    iy = (DataIndex) std::floor(y / step);
    iz = (DataIndex) std::floor(z / step);
    return {ix, iy, iz};
}

template<class IndexType>
inline IndexType grid_index_to_cell(IndexType x, IndexType y, IndexType z, const IndexType column, const IndexType row, const IndexType depth) {
    return x + column * y + column * row * z;
}

/**
 * Force should be nullified along axis of the whole domain iff vertex is on domain border
 * @param d domain
 * @param p vertex
 * @param f force
 * @return force constrained
 */

Vector_3 constraint_force(const Box3& d, const Point_3& p, const Vector_3& f);
Vector_3 constraint_force(type::Real xmin, type::Real ymin, type::Real zmin,
                          type::Real xmax, type::Real ymax, type::Real zmax, const Point_3& p, const Vector_3& f);

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


inline Real lb_getTetraederVolumeIndexed(DataIndex c1, DataIndex c2, DataIndex c3, DataIndex c4, const std::array<Point_3, 8>& vertices) {
    Real dir1_0, dir1_1, dir1_2;
    Real dir2_0, dir2_1, dir2_2;
    Real dir3_0, dir3_1, dir3_2;

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

inline Point_3 operator*(const Real w, const Point_3& p){
    return Point_3(w*p.x(), w*p.y(), w*p.z());
}

inline Point_3 operator+(const Point_3& p1, const Point_3& p2){
    return Point_3(p1.x()+p2.x(), p1.y()+p2.y(), p1.z()+p2.z());
}

Point_3 get_center_of_load(const std::vector<Real>& weights, const std::vector<Point_3>& elements);

Real compute_normalized_load(Real my_load, Real unit);

inline Vector_3 get_direction(const Point_3& vertex, const Point_3& center_of_load) {
    return (vertex - center_of_load)/std::sqrt(CGAL::squared_distance(vertex, center_of_load));
}

Vector_3 get_vertex_force(const Point_3& vertex, const std::vector<Point_3>& centers_of_load, const std::vector<Real>& normalized_loads);
std::array<std::array<VertexIndex, 3>, 12> get_planes_vids(const std::array<VertexIndex, 8>& vids);
bool lb_isGeometryValid(const std::array<Point_3, 8>& vertices, const std::array<Plane_3, 12>& planes, const Real grid_size);

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
