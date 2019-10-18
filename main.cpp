#include <iostream>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/squared_distance_3.h>
#include <mpi.h>
#include <random>
#include <iostream>
#include <chrono>
#include <thread>
#include <set>
#include "Communicator.hpp"

const double sqrt_3 = std::sqrt(3);

int my_rank, world_size;

inline int bitselect(int condition, int truereturnvalue, int falsereturnvalue) {
    return (truereturnvalue & -condition) | (falsereturnvalue & ~(-condition)); //a when TRUE and b when FintLSE
}

using Kernel = CGAL::Cartesian<double> ;
using Point_3 = CGAL::Point_3<Kernel>;
using Plane_3 = CGAL::Plane_3<Kernel>;
using Vector_3 = Kernel::Vector_3;
using Tetrahedron_3 = Kernel::Tetrahedron_3;
using Transformation = CGAL::Aff_transformation_3<Kernel>;

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

/*
inline long long position_to_cell(Point_3 const& position, const double step, const long long column, const long long row) {
    const std::vector<long long> weight = {1, column, column*row};
    long long idx = 0;

    idx += (long long) std::floor(position.x() / step);
    idx += column * (long long) std::floor(position.y() / step);
    idx += column * row * (long long) std::floor(position.z() / step);

    return idx;
}
*/

int get_rank_from_vertices(const std::array<int, 4>& vertices_id, const std::map<int, Communicator>& neighborhoods){
    std::vector<int> ranks;
    std::set<int> unique_ranks;
    for(int vid : vertices_id){
        auto neighbors = neighborhoods.at(vid).get_ranks();
        std::copy(neighbors.begin(), neighbors.end(), std::back_inserter(ranks));
        std::copy(neighbors.begin(), neighbors.end(), std::inserter(unique_ranks, unique_ranks.begin()));
    }
    for(int r : unique_ranks) {
        if(std::count(ranks.cbegin(), ranks.cend(), r) == 4 && r != my_rank) return r;
    }
    throw std::runtime_error("nobody owns the 4 vertices?");
}

inline std::pair<int, int> cell_to_global_position(int msx, int msy, long long position){
    return std::make_pair(position % msx, (int) position / msx);
}

template<class NumericalType>
inline int position_to_cell(Point_3 const& position, const double step, const NumericalType column, const NumericalType row) {
    const std::vector<NumericalType> weight = {1, column, column*row};
    NumericalType idx = 0;

    idx += (NumericalType) std::floor(position.x() / step);
    idx += column * (NumericalType) std::floor(position.y() / step);
    idx += column * row * (NumericalType) std::floor(position.z() / step);

    return idx;
}

inline Point_3 operator*(const double w, const Point_3& p){
    return Point_3(w*p.x(), w*p.y(), w*p.z());
}

inline Point_3 operator+(const Point_3& p1, const Point_3& p2){
    return Point_3(p1.x()+p2.x(), p1.y()+p2.y(), p1.z()+p2.z());
}

Point_3 get_center_of_load(const std::vector<double>& weights, const std::vector<Point_3>& elements) {

    const double total_weight = std::accumulate(weights.cbegin(), weights.cend(), 0.0);
    const auto size = elements.size();

    Point_3 r(0, 0, 0);

    for(int i = 0; i < size; ++i) {
        r = r + (weights[i] * elements[i]);
    }

    return (1.0/total_weight) * r;
}

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

/**
 * To compute the local average load use local communicator, otherwise use MPI_COMM_WORLD
 * @param my_load
 * @param average_load
 * @param neighborhood
 */
inline void compute_average_load(const double my_load, double* average_load, MPI_Comm neighborhood){
    int N;
    MPI_Comm_size(neighborhood, &N);
    MPI_Allreduce(&my_load, average_load, 1, MPI_DOUBLE, MPI_SUM, neighborhood);
    *average_load /= N;
}

inline std::vector<double> get_neighbors_load(double my_load, MPI_Comm neighborhood){
    int N;
    MPI_Comm_size(MPI_COMM_WORLD, &N);
    std::vector<double> all_loads(N);
    MPI_Allgather(&my_load, 1, MPI_DOUBLE, all_loads.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    return all_loads;
}
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

double compute_normalized_load(double my_load, double unit) {
    return my_load / unit;
}

inline Vector_3 get_direction(const Point_3& vertex, const Point_3& center_of_load) {
    return (vertex - center_of_load)/std::sqrt(CGAL::squared_distance(vertex, center_of_load));
}

Vector_3 get_vertex_force(const Point_3& vertex, const std::vector<Point_3>& centers_of_load, const std::vector<double>& normalized_loads) {
    Vector_3 f(0,0,0);
    const auto size = centers_of_load.size();
    for(int i = 0; i < size; ++i) {
        f += (normalized_loads[i] - 1) * get_direction(vertex, centers_of_load[i]);
    }
    return f;
}

// Domain = U partition;
bool in_tetrahedron(const Tetrahedron_3& tet, const Point_3& p){
    auto side = tet.bounded_side(p);
    return side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY;
}

struct Partition;

struct Domain {
    double grid_cell_size = 1.0;

    const double get_grid_cell_size() const {
        return grid_cell_size;
    }

    Point_3 v1, v2, v3, v4, v5, v6, v7, v8;

    inline const double xmin() const {
        return v1.x();
    }

    inline const double ymin() const {
        return v1.y();
    }

    inline const double zmin() const {
        return v1.z();
    }

    inline const double xmax() const {
        return v8.x();
    }

    inline const double ymax() const {
        return v8.y();
    }

    inline const double zmax() const {
        return v8.z();
    }

    inline const double xsize() const {
        return xmax()-xmin();
    }

    inline const double ysize() const {
        return ymax()-ymin();
    }

    inline const double zsize() const {
        return zmax()-zmin();
    }

    std::vector<Partition> partitions;

    Partition& get_my_partition(int rank);

    explicit Domain (std::vector<Partition>& partitions): partitions(partitions) {}

    int num_part;
    Domain (const int x, const int y, const int z):
            v1(0, 0, 0),
            v2(x, 0, 0),
            v3(0, y, 0),
            v4(x, y, 0),
            v5(0, 0, z),
            v6(x, 0, z),
            v7(0, y, z),
            v8(x, y, z)
    {}

    Domain (const Point_3& p1,const Point_3& p2,const Point_3& p3,const Point_3& p4,
            const Point_3& p5,const Point_3& p6,const Point_3& p7,const Point_3& p8):
            v1(p1),v2(p2),v3(p3),v4(p4),
            v5(p5),v6(p6),v7(p7),v8(p8) {
    }

    void bootstrap_partitions(unsigned int nb_partitions) {
        num_part = nb_partitions;
        bootstrap_partitions_cubical(nb_partitions, v1, v2, v3, v4, v5, v6, v7, v8);
    }

    /**
     * number of partitions must be a perfect cube
     */
    void bootstrap_partitions_cubical(unsigned int nb_partitions, const Point_3& p1,const Point_3& p2,const Point_3& p3,const Point_3& p4,
                                      const Point_3& p5,const Point_3& p6,const Point_3& p7,const Point_3& p8){
        const double proc_per_row = std::cbrt(nb_partitions);
        const int    row_size     = (int) proc_per_row;
        const double p_m1         = ((proc_per_row - 1) / proc_per_row);


        Transformation translate_right(  CGAL::TRANSLATION, Vector_3(std::abs(p1.x() - p2.x()) / proc_per_row, 0, 0));
        Transformation translate_fullleft(  CGAL::TRANSLATION, Vector_3(-row_size * std::abs(p1.x() - p2.x()) / proc_per_row, 0, 0));

        Transformation translate_up(     CGAL::TRANSLATION, Vector_3(0, 0, std::abs(p1.z() - p5.z()) / proc_per_row));
        Transformation translate_fulldown(     CGAL::TRANSLATION, Vector_3(0, 0, -row_size * std::abs(p1.z() - p5.z()) * p_m1));

        Transformation translate_forward(CGAL::TRANSLATION, Vector_3(0, std::abs(p1.y() - p3.y()) / proc_per_row,   0));
        Transformation translate_backward(CGAL::TRANSLATION, Vector_3(0, -std::abs(p1.y() - p3.y()) / proc_per_row, 0));

        Transformation translate_fullbackward(CGAL::TRANSLATION, Vector_3(0, -row_size * std::abs(p1.y() - p3.y()) * p_m1, 0));
        Transformation translate_fullforward(CGAL::TRANSLATION, Vector_3(0, row_size * std::abs(p1.y() - p3.y()) * p_m1, 0));

        std::array<Point_3, 8> partition_vertices =
                {p1,
                 translate_right(p1),
                 translate_forward(p1),
                 translate_forward(translate_right(p1)),
                 translate_up(p1),
                 translate_up(translate_right(p1)),
                 translate_forward(translate_up(p1)),
                 translate_forward(translate_up(translate_right(p1)))};

        int id = 0;
        for(int i = 0; i < row_size; ++i) {
            for(int j = 0; j < row_size; ++j) {
                for(int k = 0; k < row_size; ++k) {
                    partitions.emplace_back(Point_3(k ,j, i), this,
                                            partition_vertices[0], partition_vertices[1],
                                            partition_vertices[2], partition_vertices[3],
                                            partition_vertices[4], partition_vertices[5],
                                            partition_vertices[6], partition_vertices[7]);
                    for(auto &pv : partition_vertices) pv = translate_right(pv); //translate forward once
                }
                for(auto &pv : partition_vertices) pv = translate_fullleft(pv);
                for(auto &pv : partition_vertices) pv = translate_forward(pv);
            }
            for(auto &pv : partition_vertices) pv = translate_up(pv);
            for(auto &pv : partition_vertices) pv = translate_fullbackward(pv);
        }
    }

    void print_partitions() {
        std::for_each(partitions.cbegin(), partitions.cend(), [](auto p){std::cout << p << std::endl;});
    }

};

Point_3 move_vertex(const Point_3& vertex, const Vector_3& force, double mu){
    return vertex + mu*force;
}

template<class A>
std::vector<A>& search_in_linear_hashmap(std::array<std::pair<int, std::vector<A>>, 8>& linmap, int key) {
    for(std::pair<int, std::vector<A>>& entry : linmap) {
        if(entry.first == key) return entry.second;
    }
    throw std::runtime_error("doesn't exist");
}

Vector_3 constraint_force(const Domain* d, const Point_3& p, const Vector_3& f);
double compute_mu(const Domain* d, double max_normalized_load);
struct Partition {
    const Domain* d;
    Point_3 coord;

    std::array<Point_3, 8>       vertices;
    std::array<int, 8>           vertices_id;
    std::array<Tetrahedron_3, 6> tetrahedra;
    std::map<int, Communicator>  vertex_neighborhood;
    std::array<Plane_3, 12> planes;


    Partition(Point_3 coord, const Domain* d,
              Point_3 v1, Point_3 v2, Point_3 v3, Point_3 v4,
              Point_3 v5, Point_3 v6, Point_3 v7, Point_3 v8)
              : coord(coord), d(d),
              vertices({std::move(v1), std::move(v2), std::move(v3), std::move(v4),
                        std::move(v5), std::move(v6), std::move(v7), std::move(v8)}){
        for(int j = 0; j < 8; ++j)
            vertices_id[j] = get_vertex_id<int>(vertices[j], d->num_part);
        update();
    }

    void update(){
        construct_tetrahedra();
        construct_planes();
    }

    bool lb_isGeometryValid(const std::array<Point_3, 8>& vertices, const std::array<Plane_3, 12>& planes){
        /* In comparision to the set of rules in the plimpton scheme, these values are less strict
             * Except for the tetrahedron volumes, the checks are not necessary, but prevent the domain
             * from degenerating to severely, which can cause problems in the convergence behavior.
             * The movement of corners is otherwise likely to get stuck in local minima.
             */
        /* Tetrahedral subvolumes in the domain must be positively oriented */
        /* Self-intersecting cubes are bad, very bad*/
        /* Test all four permutations of how the cube can be split into tetrahedrons*/
        double v;
        if (lb_getTetraederVolumeIndexed(0, 5, 4, 7, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(0, 5, 4, 7, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(0, 3, 1, 7, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(0, 3, 1, 7, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(0, 1, 5, 7, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(0, 1, 5, 7, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(0, 4, 6, 7, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(0, 4, 6, 7, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(0, 6, 2, 7, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(0, 6, 2, 7, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(0, 2, 3, 7, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(0, 2, 3, 7, vertices) << std::endl;
            return false;
        }

        if (lb_getTetraederVolumeIndexed(1, 7, 5, 6, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(1, 7, 5, 6, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(1, 2, 3, 6, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(1, 2, 3, 6, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(1, 3, 7, 6, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(1, 3, 7, 6, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(1, 5, 4, 6, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(1, 5, 4, 6, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(1, 4, 0, 6, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(1, 4, 0, 6, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(1, 0, 2, 6, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(1, 0, 2, 6, vertices) << std::endl;
            return false;
        }

        if (lb_getTetraederVolumeIndexed(2, 4, 6, 5, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(2, 4, 6, 5, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(2, 1, 0, 5, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(2, 1, 0, 5, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(2, 0, 4, 5, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(2, 0, 4, 5, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(2, 6, 7, 5, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(2, 6, 7, 5, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(2, 7, 3, 5, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(2, 7, 3, 5, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(2, 3, 1, 5, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(2, 3, 1, 5, vertices) << std::endl;
            return false;
        }

        if (lb_getTetraederVolumeIndexed(3, 6, 7, 4, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(3, 6, 7, 4, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(3, 0, 2, 4, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(3, 0, 2, 4, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(3, 2, 6, 4, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(3, 2, 6, 4, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(3, 7, 5, 4, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(3, 7, 5, 4, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(3, 5, 1, 4, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(3, 5, 1, 4, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(3, 1, 0, 4, vertices) <= 0) {
//            std::cout << lb_getTetraederVolumeIndexed(3, 1, 0, 4, vertices) << std::endl;
            return false;
        }

        //Additionally prevent the collapse of the corners in the domain
        //This would yield a topological different domain geometry
        if (lb_getTetraederVolumeIndexed(0, 1, 4, 2, vertices) < 0) {
//            std::cout << lb_getTetraederVolumeIndexed(0, 1, 4, 2, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(1, 5, 4, 7, vertices) < 0) {
//            std::cout << lb_getTetraederVolumeIndexed(1, 5, 4, 7, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(2, 4, 6, 7, vertices) < 0) {
//            std::cout << lb_getTetraederVolumeIndexed(2, 4, 6, 7, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(2, 7, 3, 1, vertices) < 0) {
//            std::cout << lb_getTetraederVolumeIndexed(2, 7, 3, 1, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(0, 4, 6, 5, vertices) < 0) {
//            std::cout << lb_getTetraederVolumeIndexed(0, 4, 6, 5, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(0, 1, 5, 3, vertices) < 0) {
//            std::cout << lb_getTetraederVolumeIndexed(0, 1, 5, 3, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(0, 6, 2, 3, vertices) < 0) {
//            std::cout << lb_getTetraederVolumeIndexed(0, 6, 2, 3, vertices) << std::endl;
            return false;
        }
        if (lb_getTetraederVolumeIndexed(5, 3, 7, 6, vertices) < 0) {
//            std::cout << lb_getTetraederVolumeIndexed(5, 3, 7, 6, vertices) << std::endl;
            return false;
        }

        std::array<std::array<Plane_3, 6>, 8> incident_planes;

        for(int i = 0; i < 8; ++i){
            auto vertex = vertices[i];
            std::copy_if(planes.begin(), planes.end(), incident_planes[i].begin(), [&vertex](auto plane){return plane.has_on(vertex);});
        }

        for(int i = 0; i < 4; ++i){
            // dist from i to incident(7-i)
            for(auto const& plane: incident_planes[7-i]) {
                if( CGAL::squared_distance(vertices[i], plane) <= std::pow(sqrt_3*d->grid_cell_size, 2)) return false;
            }
        }

        return true;
    }

    bool lb_isGeometryValid() {
        return lb_isGeometryValid(vertices, planes);
    }

    void construct_tetrahedra() {
        tetrahedra = get_tetrahedra(vertices);
    }

/*
 *     6---------7
 *    /|        /|
 *   / |       / |
 *  /  |      /  |
 * 4-- 2 ----5 --3
 * |  /      |  /
 * Y Z       | /
 * |/        |/
 * 0----X----1
 */

    void construct_planes() {
        planes =  get_planes(vertices);
    }



    friend std::ostream &operator<<(std::ostream &os, const Partition &partition) {
        os << " ID(" << partition.vertices_id[0] << "): " << partition.vertices[0]
           << " ID(" << partition.vertices_id[1] << "): " << partition.vertices[1]
           << " ID(" << partition.vertices_id[2] << "): " << partition.vertices[2]
           << " ID(" << partition.vertices_id[3] << "): " << partition.vertices[3]
           << " ID(" << partition.vertices_id[4] << "): " << partition.vertices[4]
           << " ID(" << partition.vertices_id[5] << "): " << partition.vertices[5]
           << " ID(" << partition.vertices_id[6] << "): " << partition.vertices[6]
           << " ID(" << partition.vertices_id[7] << "): " << partition.vertices[7];
        return os;
    }

    double get_volume() {
        return (vertices[1].x()-vertices[0].x()) * (vertices[4].y()-vertices[0].y()) * (vertices[2].z()-vertices[0].z());
    }

    bool is_inside(const Point_3& p) {
         return in_tetrahedron(tetrahedra[0], p) || in_tetrahedron(tetrahedra[1], p) || in_tetrahedron(tetrahedra[2], p) || in_tetrahedron(tetrahedra[3], p) || in_tetrahedron(tetrahedra[4], p) || in_tetrahedron(tetrahedra[5], p);
    }

    template<class CartesianPointTransformer, class A>
    bool is_ghost(const CartesianPointTransformer& c, const A& p) {
        return is_ghost(c.transform(p));
    }

    bool is_ghost(const std::array<Plane_3, 6>& planes, const Point_3& p){
        std::array<double, 12> distances_to_plane;
        std::transform(planes.cbegin(), planes.cend(), distances_to_plane.begin(), [&p](Plane_3 plane){return CGAL::squared_distance(p, plane);});
        auto min_distance = std::sqrt(*std::min_element(distances_to_plane.cbegin(), distances_to_plane.cend()));
        return min_distance <= (std::sqrt(3) * d->grid_cell_size);
    }

    std::vector<Plane_3> find_planes_closer_than(const std::array<Plane_3, 6>& planes, const Point_3& p, double cut_off) {
        std::array<double, 6> distances_to_plane;
        std::transform(planes.cbegin(), planes.cend(), distances_to_plane.begin(), [&p](Plane_3 plane){return std::sqrt(CGAL::squared_distance(p, plane));});
        std::vector<Plane_3> closest_planes;
        for(int i = 0; i < 6; ++i){
            if(distances_to_plane[i] <= cut_off) closest_planes.push_back(planes[i]);
        }
        return closest_planes;
    }

    template<class CartesianPointTransformer, class A>
    void move_data(const std::vector<A>& elements) {
        CartesianPointTransformer c;
        std::array<std::vector<A>, 26> to_move;
        const auto size = elements.size();
        for(int i = 0; i < size; ++i) {
            const auto& e = elements.at(i);
            auto closest_planes = find_planes_closer_than(planes, c.transform(elements[i]), sqrt_3*d->grid_cell_size);
            if(closest_planes.size() > 0)
                ghosts.push_back(std::make_pair(closest_planes, i));
            if(this->is_inside(c.transform(e))) {

            } else {

            }
        }
    }

    template<class CartesianPointTransformer, class A>
    std::vector<std::pair<std::vector<Plane_3>, int> > get_ghosts_and_destination_planes(const std::vector<A>& elements) {
        CartesianPointTransformer c;
        std::vector<std::pair<std::vector<Plane_3>, int>> ghosts;
        const auto size = elements.size();
        for(int i = 0; i < size; ++i) {
            auto closest_planes = find_planes_closer_than(planes, c.transform(elements[i]), sqrt_3*d->grid_cell_size);
            if(closest_planes.size() > 0)
                ghosts.push_back(std::make_pair(closest_planes, i));
        }
        return ghosts;
    }

    std::vector<std::pair<std::vector<int>, int>> compute_ghosts_destination_ranks(const std::vector<std::pair<std::vector<Plane_3>, int>>& ghosts_and_planes) {
        std::vector<std::pair<std::vector<int>, int>> ghost_and_ranks;
        for(auto const& ghost_and_planes : ghosts_and_planes) {
            std::vector<int> neighbors_of;
            for(auto& plane : ghost_and_planes.first) {
                auto vertices_on_plane = get_points_on_plane(vertices.begin(), vertices.end(), vertices_id.begin(), vertices_id.end(), plane);
                neighbors_of.push_back(get_rank_from_vertices(vertices_on_plane, vertex_neighborhood));
            }
            ghost_and_ranks.emplace_back(neighbors_of, ghost_and_planes.second);
        }
        return ghost_and_ranks;
    }

    std::array<std::pair<int, std::vector<Point_3>>, 8> spread_centers_of_load(const Point_3& cl, const std::map<int, Communicator>& neighborhood) {
        int N;
        std::array<std::pair<int, std::vector<Point_3>>, 8> all_cl;
        std::array<std::pair<int, std::vector<double>>,  8> buff_all_cl;
        std::array<MPI_Request, 8> requests;

        std::array<double, 3> my_cl = {cl.x(), cl.y(), cl.z()};
        int idx = 0;
        for(auto&  comm : neighborhood) {
            int i = comm.first;
            N     = comm.second.comm_size;
            all_cl[idx] = std::make_pair(i, std::vector<Point_3>());
            buff_all_cl[idx] = std::make_pair(i, std::vector<double>());
            auto& curr_buff_all_cl = search_in_linear_hashmap(buff_all_cl, i);
            auto& curr_all_cl      = search_in_linear_hashmap(all_cl, i);
            curr_buff_all_cl.resize(3*N);
            curr_all_cl.resize(N);
            comm.second.Allgather(my_cl.data(), 3, MPI_DOUBLE, curr_buff_all_cl.data(), 3, MPI_DOUBLE);
            for(int j = 0; j < N; ++j) {// reconstruct the points
                curr_all_cl[j] = Point_3(curr_buff_all_cl[j*3],curr_buff_all_cl[j*3+1],curr_buff_all_cl[j*3+2]);
            }
            idx++;
        }
        return all_cl;
    }

    template<class LoadComputer, class A>
    double get_load_imbalance(const std::vector<A>& elements, MPI_Comm neighborhood = MPI_COMM_WORLD){
        LoadComputer lc;
        double my_load = lc.compute_load(elements);
        int N; MPI_Comm_size(neighborhood, &N);
        auto all_loads = get_neighbors_load(my_load, neighborhood); //global load balancing with MPI_COMM_WORLD
        auto avg_load  = std::accumulate(all_loads.cbegin(), all_loads.cend(), 0.0) / N;
        return (*std::max_element(all_loads.cbegin(), all_loads.cend()) / avg_load) - 1;
    }

    template<class CartesianPointTransformer, class LoadComputer, class MoveMethod, class A>
    void move_vertices(const std::vector<A>& elements, MPI_Comm neighborhood = MPI_COMM_WORLD) {
        LoadComputer lc;
        CartesianPointTransformer transformer;
        std::vector<Point_3> points;

        double my_load = lc.compute_load(elements);

        const std::vector<double>& weights = lc.get_weights(elements);

        std::transform(elements.cbegin(), elements.cend(), std::back_inserter(points), [&transformer](auto el){return transformer.transform(el);});

        int N; MPI_Comm_size(neighborhood, &N);

        auto all_loads = get_neighbors_load(my_load, neighborhood); //global load balancing with MPI_COMM_WORLD

        auto avg_load  = std::accumulate(all_loads.cbegin(), all_loads.cend(), 0.0) / N;

        std::vector<double> normalized_loads;

        std::transform(all_loads.cbegin(), all_loads.cend(), std::back_inserter(normalized_loads), [&avg_load](auto n){return n/avg_load;});

        auto max_normalized_load = *std::max_element(normalized_loads.cbegin(), normalized_loads.cend());

        double mu = compute_mu(d, max_normalized_load);

        std::transform(all_loads.cbegin(), all_loads.cend(), std::back_inserter(normalized_loads), [&avg_load](auto v){return v / avg_load;});

        auto center_of_load = get_center_of_load(weights, points);

        std::array<std::pair<int, std::vector<Point_3>>, 8> all_cl = spread_centers_of_load(center_of_load, vertex_neighborhood);
        std::vector<int> vertex_local_id_to_move = {0,1,2,3,4,5,6,7,8};

        MoveMethod move_method(vertex_local_id_to_move);
        const auto original_mu = mu;
        const int MAX_TRIAL=1;
        unsigned int iteration = 0;

        while(iteration < 8) {
            auto x = iteration&1; auto y = (iteration&2)>>1; auto z = (iteration&4)>>2;
            int com = (((unsigned int) coord.x()+x) & 1)  +
                      (((unsigned int) coord.y()+y) & 1)*2+
                      (((unsigned int) coord.z()+z) & 1)*4;

            for(int trial = 0; trial < MAX_TRIAL; ++trial) {
                auto new_vertices = vertices;
                auto v   = vertices[com];
                auto vid = vertices_id[com];
                auto cls = search_in_linear_hashmap<Point_3>(all_cl, vid);
                auto f1  = get_vertex_force(v, cls, normalized_loads);
                auto f1_after = constraint_force(d, v, f1);
                new_vertices[com] = move_vertex(v, f1_after, mu);
                int valid = lb_isGeometryValid(new_vertices, get_planes(new_vertices));
                const auto& communicator = vertex_neighborhood[vid];
                std::vector<int> allValid(communicator.comm_size);
                communicator.Allgather(&valid, 1, MPI_INT, allValid.data(), 1, MPI_INT);
                int are_all_valid = std::accumulate(allValid.begin(), allValid.end(), 1, [](auto k, auto v){return k*v;});
                if(are_all_valid) {
                    vertices = new_vertices;
                    update();
                    iteration++;
                    mu = original_mu;
                    break;
                } else {
                    mu /= 2.0;
                    if(trial == MAX_TRIAL-1)
                        iteration++;
                }
            }
        }

    }

    template<class NumericalType>
    const NumericalType get_vertex_id(const Point_3& v, int N) const {
        auto proc_per_row = (NumericalType) std::cbrt(N) + 1;
        const double step = this->d->xsize() / (proc_per_row-1);
        return position_to_cell<NumericalType>(v, step, proc_per_row, proc_per_row);
    }

    void init_communicators(int world_size) {
        int N = world_size, x;

        std::array<int, 8>  sbuff;
        std::vector<int> buff(8*N);
        std::map<int, std::set<int>> neighbors;

        //for(int i = 0; i < 8; ++i)  vertices_id[i] = get_vertex_id<int>(vertices[i], N);
        std::copy(vertices_id.begin(), vertices_id.end(), sbuff.begin());
        std::sort(sbuff.begin(), sbuff.end());

        MPI_Allgather(sbuff.data(), 8, MPI_INT, buff.data(), 8, MPI_INT, MPI_COMM_WORLD);

        for(int j = 0; j < 8; ++j) {
            int vid = vertices_id[j];
            for(int i = 0; i < N; ++i) {
                auto beg = buff.begin() + i * 8;
                auto end = buff.begin() + (i+1) * 8;
                if(std::binary_search(beg, end, vid)) neighbors[vid].insert(i);
            }
            neighbors[vid].insert(my_rank);
            std::vector<int> n_list(neighbors[vid].begin(), neighbors[vid].end());
            Communicator c(n_list);
            vertex_neighborhood[vid] = c;
        }
    }
};

Partition& Domain::get_my_partition(int rank) {
    return partitions.at(rank);
}

/**
 * Force should be nullified along axis of the whole domain iff vertex is on domain border
 * @param d domain
 * @param p vertex
 * @param f force
 * @return force constrained
 */
Vector_3 constraint_force(const Domain* d, const Point_3& p, const Vector_3& f){
    const auto x = p.x();
    const auto y = p.y();
    const auto z = p.z();
    double fx = f.x(), fy = f.y(), fz = f.z();
    if(x == d->xmin()) {
        fx = 0;
    }
    if(y == d->ymin()) {
        fy = 0;
    }
    if(z == d->zmin()) {
        fz = 0;
    }
    if(x == d->xmax()) {
        fx = 0;
    }
    if(y == d->ymax()) {
        fy = 0;
    }
    if(z == d->zmax()) {
        fz = 0;
    }
    return Vector_3(fx,fy,fz);
}

double compute_mu(const Domain* d, double max_normalized_load){
    double sigma_max = max_normalized_load-1;
    return sigma_max == 0 ? 0 : d->get_grid_cell_size()/(sigma_max);
}

template<class A>
inline void linear_to_grid(const A index, const A c, const A r, int& x_idx, int& y_idx, int& z_idx){
    x_idx = (int) (index % (c*r) % c);           // col
    y_idx = (int) std::floor(index % (c*r) / c); // row
    z_idx = (int) std::floor(index / (c*r));     // depth
    assert(x_idx < c);
    assert(y_idx < r);
    assert(z_idx < r);
};

struct CommunicationDatatype {
    MPI_Datatype element_datatype;
    MPI_Datatype minimal_datatype;
    CommunicationDatatype(const MPI_Datatype &el, const MPI_Datatype &min) : element_datatype(el), minimal_datatype(min){}
    void free_datatypes() { MPI_Type_free(&element_datatype); MPI_Type_free(&minimal_datatype);}
};

struct Cell {
    int gid, type; //type = -1:empty, 0:rock, 1:water
    float weight, erosion_probability;
    double average_load = 2000;

    Cell() : gid(0), type(0), weight(1.0), erosion_probability(0) {};
    Cell(int gid, int type, float weight, float erosion_probability)
            : gid(gid), type(type), weight(weight), erosion_probability(erosion_probability) {};
    Cell(int gid, int type, float weight, float erosion_probability, double average_load)
            : gid(gid), type(type), weight(weight), erosion_probability(erosion_probability), average_load(average_load) {};

    template<class NumericType>
    std::array<NumericType, 3> get_position_as_array() const {
        auto position_as_pair = cell_to_global_position(Cell::get_msx(), Cell::get_msy(), gid);
        std::array<NumericType, 2> array = {(NumericType) position_as_pair.first, (NumericType) position_as_pair.second};
        return array;
    }

    std::tuple<int, int, int> get_position_as_pair() const {
        int x,y,z;
        linear_to_grid(gid, Cell::get_msx(), Cell::get_msy(), x, y, z); //cell_to_global_position(Cell::get_msx(), Cell::get_msy(), gid);
        return std::make_tuple(x,y,z);
    }

    std::tuple<double, double, double> get_center() const {
        int x,y,z;
        linear_to_grid(gid, Cell::get_msx(), Cell::get_msy(), x, y, z); //cell_to_global_position(Cell::get_msx(), Cell::get_msy(), gid);
        return std::make_tuple(x+Cell::get_cell_size(), y+Cell::get_cell_size(), z+Cell::get_cell_size());
    }

    void get_center(double* _x, double* _y, double* _z) const {
        int x,y,z;
        linear_to_grid(gid, Cell::get_msx(), Cell::get_msy(), x, y, z); //cell_to_global_position(Cell::get_msx(), Cell::get_msy(), gid);
        *_x=x+Cell::get_cell_size(); *_y=y+Cell::get_cell_size(); *_z =z+Cell::get_cell_size();
    }

    static CommunicationDatatype register_datatype() {

        MPI_Datatype cell_datatype, gid_type_datatype;

        MPI_Aint intex, lb, floatex;

        const int number_of_int_elements    = 2;
        const int number_of_float_elements  = 2;
        const int number_of_double_elements = 1;

        int blockcount_element[3];

        blockcount_element[0] = number_of_int_elements;    // gid, lid, exit, waiting_time
        blockcount_element[1] = number_of_float_elements;  // position <x,y>
        blockcount_element[2] = number_of_double_elements; // position <x,y>
        //int
        MPI_Type_contiguous(number_of_int_elements, MPI_INT, &gid_type_datatype); // position
        MPI_Type_commit(&gid_type_datatype);

        MPI_Datatype blocktypes[3];
        blocktypes[0] = MPI_INT;
        blocktypes[1] = MPI_FLOAT;
        blocktypes[2] = MPI_DOUBLE;

        MPI_Type_get_extent(MPI_INT, &lb, &intex);
        MPI_Type_get_extent(MPI_FLOAT, &lb, &floatex);

        MPI_Aint offset[3];
        offset[0] = static_cast<MPI_Aint>(0);
        offset[1] = 2*intex;
        offset[2] = 2*intex + 2*floatex;

        MPI_Type_create_struct(3, blockcount_element, offset, blocktypes, &cell_datatype);
        MPI_Type_commit(&cell_datatype);

        return {cell_datatype, gid_type_datatype};
    }

    static void set_msx(int _msx){
        static int msx = _msx;
    }
    static void set_msy(int _msy){
        static int msy = _msy;
    }
    static int& get_msx(){
        static int msx;
        return msx;
    }
    static int& get_msy(){
        static int msy;
        return msy;
    }
    static int& get_msz(){
        static int msz;
        return msz;
    }
    static double& get_cell_size(){
        static double size;
        return size;
    }

    friend std::ostream &operator<<(std::ostream &os, const Cell &cell) {
        os << "gid: " << cell.gid << " type: " << cell.type << " weight: " << cell.weight << " erosion_probability: "
           << cell.erosion_probability;
        return os;
    }
};

struct GridElementComputer {
    double compute_load(const std::vector<Cell>& elements){
        return std::accumulate(elements.cbegin(), elements.cend(), 0.0, [](double l, Cell e){return l + e.weight;});
    }
    std::vector<double> get_weights(const std::vector<Cell>& elements) {
        std::vector<double> weights;
        std::transform(elements.cbegin(), elements.cend(), std::back_inserter(weights), [](Cell e) {
            return e.weight;
        });
        return weights;
    }
};

struct GridPointTransformer {
    Point_3 transform(const Cell& element){
        double x,y,z;
        std::tie(x,y,z) = element.get_center();
        return Point_3(x, y, z);
    }
};
struct ColoredMoveMethod {

};
struct GlobalMoveMethod {
    using  LocalVertexId  = int;
    using GlobalVertexId  = int;
    using RejectionCount  = int;

    template<class A>
    using CenterOfLoadsDB = std::array<std::pair<int, std::vector<A>>, 8>;
    std::vector<std::pair<LocalVertexId, RejectionCount>> rejection_count_by_vertex;

    GlobalMoveMethod(const std::vector<LocalVertexId>& vertices_to_move) {
        rejection_count_by_vertex.resize(vertices_to_move.size());
        std::transform(vertices_to_move.cbegin(), vertices_to_move.cend(), rejection_count_by_vertex.begin(), [](auto varrid) {return std::make_pair(varrid, 0);});
    }

    template<class GetNeighborFunc, class GetVertexFunc, class ConstraintForceFunc, class MoveVertexFunc>
    std::array<Point_3, 8> move(const std::array<Point_3, 8>& vertices,
                                    const std::array<int, 8>& vertices_id,
                                    const Domain* d,
                                    const std::vector<double>& normalized_loads,
                                    CenterOfLoadsDB<Point_3> & all_cl,
                                    const double mu,
                                    const std::vector<LocalVertexId>& vertices_to_move,
                                    const GetNeighborFunc& get_neighbors,
                                    const GetVertexFunc& get_vertex_force,
                                    const ConstraintForceFunc& constraint_force,
                                    const MoveVertexFunc& move_vertex ) {
        std::array<Point_3, 8> new_vertices;

        for(int i = 0; i < 8; ++i) {
            auto v  = vertices[i];
            auto vid = vertices_id[i];
            auto cls = get_neighbors(all_cl, vid);
            auto f1 = get_vertex_force(v, cls, normalized_loads);
            auto f1_after = constraint_force(d, v, f1);
            new_vertices[i] = move_vertex(v, f1_after, mu);
        }

        return new_vertices;
    }
};

std::vector<Cell> generate_lattice_single_type( int msx, int msy,
                                                int x_proc_idx, int y_proc_idx,
                                                int cell_in_my_cols, int cell_in_my_rows,
                                                int type, float weight, float erosion_probability) {
    int cell_per_process = cell_in_my_cols * cell_in_my_rows;
    std::vector<Cell> my_cells; my_cells.reserve(cell_per_process);
    for(int j = 0; j < cell_in_my_cols; ++j) {
        for(int i = 0; i < cell_in_my_rows; ++i) {
            int gid = cell_in_my_rows * x_proc_idx + i + msx * (j + (y_proc_idx * cell_in_my_cols));
            my_cells.emplace_back(gid, type, weight, erosion_probability);
        }
    }
    return my_cells;
}

int main() {
    MPI_Init(nullptr, nullptr);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int col_row = (int) std::cbrt(world_size);

    const int DOMAIN_SIZE_X = 10;
    const int DOMAIN_SIZE_Y = 10;
    const int DOMAIN_SIZE_Z = 10;

    Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z);

    d.grid_cell_size = 1;

    Cell::get_msx() = DOMAIN_SIZE_X / d.grid_cell_size;
    Cell::get_msy() = DOMAIN_SIZE_Y / d.grid_cell_size;
    Cell::get_msz() = DOMAIN_SIZE_Z / d.grid_cell_size;
    Cell::get_cell_size() = d.grid_cell_size;

    d.bootstrap_partitions(world_size);

    auto part = d.get_my_partition(my_rank);
    part.init_communicators(world_size);

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> normal_distribution(1.0 + my_rank % 2, 0.2);

    auto bbox = CGAL::bbox_3(part.vertices.begin(), part.vertices.end());
    int cell_in_my_rows = (int) (bbox.xmax() - bbox.xmin()) / d.grid_cell_size;
    int cell_in_my_cols = (int) (bbox.ymax() - bbox.ymin()) / d.grid_cell_size;

    auto cell_per_process = cell_in_my_cols * cell_in_my_rows;

    std::vector<Cell> my_cells; my_cells.reserve(cell_per_process);

    int x_proc_idx, y_proc_idx, z_proc_idx; linear_to_grid(my_rank, col_row, col_row, x_proc_idx, y_proc_idx, z_proc_idx);

    const int xcells = cell_in_my_rows * col_row, ycells = cell_in_my_rows * col_row, zcells = cell_in_my_rows * col_row;

    for(int j = 0; j < cell_in_my_rows; ++j) {
        for(int i = 0; i < cell_in_my_rows; ++i) {
            for(int k = 0; k < cell_in_my_rows; ++k) {
                int gid = cell_in_my_rows * x_proc_idx + i + xcells * (j + (y_proc_idx * cell_in_my_cols)) + ((z_proc_idx * cell_in_my_rows) + k) * xcells * ycells;
                my_cells.emplace_back(gid, 0, normal_distribution(gen), 0.0);
            }
        }
    }
    //std::cout << "BEFORE: " << part.get_load_imbalance<GridElementComputer>(my_cells) << std::endl;

    part.move_vertices<GridPointTransformer, GridElementComputer, GlobalMoveMethod, Cell>(my_cells);

    //std::cout << "AFTER: " << part.get_load_imbalance<GridElementComputer>(my_cells) << std::endl;

    std::cout << my_rank << " " << part << std::endl;

    MPI_Finalize();

    return 0;
}




