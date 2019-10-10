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

int rank, world_size;

inline int bitselect(int condition, int truereturnvalue, int falsereturnvalue) {
    return (truereturnvalue & -condition) | (falsereturnvalue & ~(-condition)); //a when TRUE and b when FintLSE
}

using Kernel = CGAL::Cartesian<double> ;
using Point_3 = CGAL::Point_3<Kernel>;
using Vector_3 = Kernel::Vector_3;
using Tetrahedron_3 = Kernel::Tetrahedron_3;
using Transformation = CGAL::Aff_transformation_3<Kernel>;
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
std::array<MPI_Group, 8> ngr;
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
        //std::cout << (weights[i] * elements[i]) << std::endl;
        r = r + (weights[i] * elements[i]);
    }

    return (1.0/total_weight) * r;
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
    //std::cout << "> "<< my_load << std::endl;
    MPI_Allgather(&my_load, 1, MPI_DOUBLE, all_loads.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    //std::for_each(all_loads.cbegin(), all_loads.cend(), [](auto v){ std::cout << "!!"<<v << std::endl; });
    return all_loads;
}

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
/*
 *     7---------8
 *    /|        /|
 *   / |       / |
 *  /  |      /  |
 * 5-- 3 ----6 --4
 * |  /      |  /
 * Y Z       | /
 * |/        |/
 * 1----X----2
 */
    int num_part;
    Domain (const int x, const int y, const int z):
            v1(0, 0, 0),
            v2(x, 0, 0),
            v3(0, 0, z),
            v4(x, 0, z),
            v5(0, y, 0),
            v6(x, y, 0),
            v7(0, y, z),
            v8(x, y, z){
    }

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
     * What is initial partitioning?
     * @param nb_partitions
     */
    void bootstrap_partitions(int generated_parts, unsigned int nb_partitions, const Point_3& p1,const Point_3& p2,const Point_3& p3,const Point_3& p4,
                              const Point_3& p5,const Point_3& p6,const Point_3& p7,const Point_3& p8) {

        if(generated_parts == nb_partitions) {
            partitions.emplace_back(0, this, p1, p2, p3, p4, p5, p6, p7, p8);
            return;
        }

        int largest_dim = largest_dimension(p1, p2, p3, p4, p5, p6, p7, p8);

        if(largest_dim == 0) { // cut the domain along X
            Transformation translate_right(CGAL::TRANSLATION, Vector_3( std::abs(p1.x() - p2.x()) / 2, 0, 0));
            Transformation translate_left( CGAL::TRANSLATION,  Vector_3(-std::abs(p1.x() - p2.x()) / 2, 0, 0));
            bootstrap_partitions(generated_parts * 2, nb_partitions, p1, translate_left(p2), p3, translate_left(p4), p5, translate_left(p6), p7, translate_left(p8));
            bootstrap_partitions(generated_parts * 2, nb_partitions, translate_right(p1), (p2), translate_right(p3), (p4), translate_right(p5), (p6), translate_right(p7), (p8));
        } else if(largest_dim == 1) { // cut the domain along Y
            Transformation   translate_up(CGAL::TRANSLATION, Vector_3(0,  std::abs(p1.y() - p5.y()) / 2, 0));
            Transformation translate_down(CGAL::TRANSLATION, Vector_3(0, -std::abs(p1.y() - p5.y()) / 2, 0));
            bootstrap_partitions(generated_parts * 2, nb_partitions, translate_up(p1), translate_up(p2), translate_up(p3), translate_up(p4), (p5), (p6), p7, p8);
            bootstrap_partitions(generated_parts * 2, nb_partitions, p1, p2, (p3), (p4), translate_down(p5), translate_down(p6), translate_down(p7), translate_down(p8));
        } else { // cut the domain along Z
            Transformation  translate_forward(CGAL::TRANSLATION, Vector_3(0, 0,  std::abs(p1.z() - p3.z()) / 2));
            Transformation translate_backward(CGAL::TRANSLATION, Vector_3(0, 0, -std::abs(p1.z() - p3.z()) / 2));
            bootstrap_partitions(generated_parts * 2, nb_partitions, translate_forward(p1), translate_forward(p2), p3, p4, translate_forward(p5), translate_forward(p6), p7, p8);
            bootstrap_partitions(generated_parts * 2, nb_partitions, p1, p2, translate_backward(p3), translate_backward(p4), (p5), (p6), translate_backward(p7), translate_backward(p8));
        }
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

        Transformation translate_up(     CGAL::TRANSLATION, Vector_3(0, std::abs(p1.y() - p5.y()) / proc_per_row, 0));
        Transformation translate_fulldown(     CGAL::TRANSLATION, Vector_3(0, -row_size * std::abs(p1.y() - p5.y()) * p_m1, 0));

        Transformation translate_forward(CGAL::TRANSLATION, Vector_3(0, 0, std::abs(p1.z() - p3.z()) / proc_per_row));
        Transformation translate_fullbackward(CGAL::TRANSLATION, Vector_3(0, 0, -row_size * std::abs(p1.z() - p3.z()) * p_m1));

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
                for(int k = 0; k < row_size; ++k, ++id) {
                    partitions.emplace_back(id, this,
                                            partition_vertices[0], partition_vertices[1],
                                            partition_vertices[2], partition_vertices[3],
                                            partition_vertices[4], partition_vertices[5],
                                            partition_vertices[6], partition_vertices[7]);
                    for(auto &pv : partition_vertices) pv = translate_right(pv); //translate forward once
                }
                for(auto &pv : partition_vertices) pv = translate_fullleft(pv);
                for(auto &pv : partition_vertices) pv = translate_up(pv);
            }
            for(auto &pv : partition_vertices) pv = translate_fulldown(pv);
            for(auto &pv : partition_vertices) pv = translate_forward(pv);
        }
    }

    int largest_dimension(const Point_3& p1,const Point_3& p2,const Point_3& p3,const Point_3& p4,
                          const Point_3& p5,const Point_3& p6,const Point_3& p7,const Point_3& p8){
        std::array<double, 3> maxs = {std::abs(p1.x() - p2.x()), std::abs(p1.y() - p5.y()), std::abs(p1.z() - p3.z())};
        int ret = 0; double v = -1.0;
        for(int i = 0; i < 3; ++i) {
            if(maxs[i] > v) {
                ret = i;
                v = maxs[i];
            }
        }
        return ret;
    }

    void print_partitions() {
        std::for_each(partitions.cbegin(), partitions.cend(), [](auto p){std::cout << p << std::endl;});
    }

};

Vector_3 constraint_force(const Domain* d, const Point_3& p, const Vector_3& f);
double compute_mu(const Domain* d, double max_normalized_load);
struct Partition {
    const Domain* d;
    int id;
    std::array<Point_3, 8> vertices;
    std::array<int, 8>     vertices_id;
    std::array<Tetrahedron_3, 6> tetrahedra;
    std::map<int, Communicator> vertex_neighborhood = {
//            {0, MPI_COMM_NULL},
//            {1, MPI_COMM_NULL},
//            {2, MPI_COMM_NULL},
//            {3, MPI_COMM_NULL},
//            {4, MPI_COMM_NULL},
//            {5, MPI_COMM_NULL},
//            {6, MPI_COMM_NULL},
//            {7, MPI_COMM_NULL},
    };

    Partition(int id, const Domain* d, Point_3 v1, Point_3 v2, Point_3 v3, Point_3 v4,
              Point_3 v5, Point_3 v6, Point_3 v7, Point_3 v8)
              : id(id), d(d),
              vertices({std::move(v1), std::move(v2), std::move(v3), std::move(v4),
                        std::move(v5), std::move(v6), std::move(v7), std::move(v8)}){
        for(int j = 0; j < 8; ++j) {
            vertices_id[j] = get_vertex_id<int>(vertices[j], d->num_part);
        }
        construct_tetrahedra();
    }

    void construct_tetrahedra() {
        tetrahedra[0] = Tetrahedron_3(vertices[0], vertices[1], vertices[3], vertices[7]);
        tetrahedra[1] = Tetrahedron_3(vertices[0], vertices[1], vertices[5], vertices[7]);
        tetrahedra[2] = Tetrahedron_3(vertices[0], vertices[2], vertices[3], vertices[7]);
        tetrahedra[3] = Tetrahedron_3(vertices[0], vertices[2], vertices[6], vertices[7]);
        tetrahedra[4] = Tetrahedron_3(vertices[0], vertices[4], vertices[5], vertices[7]);
        tetrahedra[5] = Tetrahedron_3(vertices[0], vertices[4], vertices[6], vertices[7]);

        //Union of tetrahedra must be equal to Partition, hence sum of volumes = volume of partition
        /*auto volume = std::abs(tetrahedra[0].volume()) + std::abs(tetrahedra[1].volume()) +
                      std::abs(tetrahedra[2].volume()) + std::abs(tetrahedra[3].volume()) +
                      std::abs(tetrahedra[4].volume()) + std::abs(tetrahedra[5].volume());
        assert(volume - this->get_volume() <= std::numeric_limits<double>::epsilon());*/
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

    Point_3 move_vertex(const Point_3& vertex, const Vector_3& force, double mu){
        return vertex + mu*force;
    }

    std::array<std::vector<Point_3>, 8> spread_centers_of_load(const Point_3& cl, const std::map<int, Communicator>& neighborhood) {
        int N;
        std::array<std::vector<Point_3>, 8> all_cl;
        std::array<std::vector<double>,  8> buff_all_cl;
        std::array<MPI_Request, 8> requests;

        std::array<double, 3> my_cl = {cl.x(), cl.y(), cl.z()};

        int i = 0;
        for(auto&  comm : neighborhood) {
            //MPI_Comm_size(comm.second, &N);
            N = comm.second.comm_size;
            buff_all_cl[i].resize(3*N);
            all_cl[i].resize(N);
            //MPI_Iallgather(my_cl.data(), 3, MPI_DOUBLE, buff_all_cl[i].data(), 3, MPI_DOUBLE, comm.second, &requests[i]);
            //std::cout<<comm.second.comm_size<<std::endl;
            comm.second.Allgather(my_cl.data(), 3, MPI_DOUBLE, buff_all_cl[i].data(), 3, MPI_DOUBLE);
            ++i;
        }

        //MPI_Waitall(8, requests.data(), MPI_STATUSES_IGNORE);

        for(int i = 0; i < 8; ++i) {
            auto N = all_cl[i].size();
            for(int j = 0; j < N; ++j) {
                all_cl[i][j] = Point_3(buff_all_cl[i][j*3],buff_all_cl[i][j*3+1],buff_all_cl[i][j*3+2]);
            }
        }

        return all_cl;
    }

    template<class CartesianPointTransformer, class LoadComputer, class A>
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

        //std::for_each(all_loads.cbegin(), all_loads.cend(), [](auto v){ std::cout << v << std::endl; });

        decltype(all_loads) normalized_loads;
        std::transform(all_loads.cbegin(), all_loads.cend(), std::back_inserter(normalized_loads), [&avg_load](auto n){return n/avg_load;});

        auto max_normalized_load = *std::max_element(normalized_loads.cbegin(), normalized_loads.cend());

        double mu = compute_mu(d, max_normalized_load);

        std::transform(all_loads.cbegin(), all_loads.cend(), std::back_inserter(normalized_loads), [&avg_load](auto v){return v / avg_load;});

        auto center_of_load = get_center_of_load(weights, points);
        auto all_cl = spread_centers_of_load(center_of_load, vertex_neighborhood);

        for(int i = 0; i < 8; ++i) {
            auto& v = vertices[i];

            auto cls = all_cl[i];

            std::ostringstream vts;

            std::copy(cls.begin(), cls.end(),
                      std::ostream_iterator<Point_3>(vts, ", "));

            //std::cout <<  get_vertex_id<int>(v, 8) << " | " << rank << std::endl;

            auto f1 = get_vertex_force(v, cls, normalized_loads);

            auto f1_after = constraint_force(d, v, f1);

            if(vertices_id[i] == 4){
                auto ranks = vertex_neighborhood[4].get_ranks();
                std::ostringstream vts;
                std::copy(cls.begin(), cls.end(),
                          std::ostream_iterator<Point_3>(vts, "; "));
                std::cout << f1 << " ;" << vts.str()  << std::endl;

            }
            v = move_vertex(v, f1_after, mu);

        }
        construct_tetrahedra();
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

        for(int j = 0; j < 8; ++j){
            int vid = vertices_id[j];
            for(int i = 0; i < N; ++i){
                auto beg = buff.begin() + i * 8;
                auto end = buff.begin() + (i+1) * 8;

                if(std::binary_search(beg, end, vid)) neighbors[vid].insert(i);
            }
            neighbors[vid].insert(rank);
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

    static CommunicationDatatype register_datatype() {

        MPI_Datatype cell_datatype, gid_type_datatype;

        MPI_Aint intex, lb, floatex;

        const int number_of_int_elements    = 2;
        const int number_of_float_elements  = 2;
        const int number_of_double_elements = 1;

        int blockcount_element[3];

        blockcount_element[0] = number_of_int_elements; // gid, lid, exit, waiting_time
        blockcount_element[1] = number_of_float_elements; // position <x,y>
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

        return CommunicationDatatype(cell_datatype, gid_type_datatype);
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
        std::tie(x,y,z) = element.get_position_as_pair();
        return Point_3(x,y,z);
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
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int col_row = (int) std::cbrt(world_size);

    Domain d(10, 10, 10);

    d.grid_cell_size = 1;

    Cell::get_msy() = 10;
    Cell::get_msx() = 10;

    d.bootstrap_partitions(world_size);

    auto part = d.get_my_partition(rank);
    part.init_communicators(world_size);

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> normal_distribution(1.0, 0.2);

    auto bbox = CGAL::bbox_3(part.vertices.begin(), part.vertices.end());

    int cell_in_my_rows = (int) (bbox.xmax() - bbox.xmin()) / d.grid_cell_size;
    int cell_in_my_cols = (int) (bbox.ymax() - bbox.ymin()) / d.grid_cell_size;

    auto cell_per_process = cell_in_my_cols * cell_in_my_rows;

    std::vector<Cell> my_cells; my_cells.reserve(cell_per_process);

    int x_proc_idx, y_proc_idx, z_proc_idx; linear_to_grid(rank, col_row, col_row, x_proc_idx, y_proc_idx, z_proc_idx);

    const int xcells = cell_in_my_rows * col_row, ycells = cell_in_my_rows * col_row, zcells = cell_in_my_rows * col_row;

    for(int j = 0; j < cell_in_my_rows; ++j) {
        for(int i = 0; i < cell_in_my_rows; ++i) {
            for(int k = 0; k < cell_in_my_rows; ++k) {
                int gid = cell_in_my_rows * x_proc_idx + i + xcells * (j + (y_proc_idx * cell_in_my_cols)) + ((z_proc_idx * cell_in_my_rows) + k) * xcells * ycells;
                my_cells.emplace_back(gid, 0, normal_distribution(gen), 0.0);
            }
        }
    }

    part.move_vertices<GridPointTransformer, GridElementComputer, Cell>(my_cells);

    std::cout << part << std::endl;

    MPI_Finalize();

    return 0;
}




