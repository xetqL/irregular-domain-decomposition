#include <iostream>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/squared_distance_3.h>
#include <mpi.h>

inline int bitselect(int condition, int truereturnvalue, int falsereturnvalue)
{
    return (truereturnvalue & -condition) | (falsereturnvalue & ~(-condition)); //a when TRUE and b when FintLSE
}

using Kernel = CGAL::Cartesian<double> ;
using Point_3 = CGAL::Point_3<Kernel>;
using Vector_3 = Kernel::Vector_3;
using Tetrahedron_3 = Kernel::Tetrahedron_3;
using Transformation = CGAL::Aff_transformation_3<Kernel>;

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

    for(int i = 0; i < size; ++i) r = r + (weights[i] * elements[i]);

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
    MPI_Comm_size(neighborhood, &N);
    std::vector<double> all_loads(N);
    MPI_Alltoall(&my_load, 1, MPI_DOUBLE, all_loads.data(), 1, MPI_DOUBLE, neighborhood);
    return all_loads;
}

std::vector<Point_3> get_centers_of_load(const Point_3& cl, MPI_Comm neighborhood){
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

inline Vector_3 get_direction(const Point_3& vertex, const Point_3& center_of_load){
    return (vertex - center_of_load)/std::sqrt(CGAL::squared_distance(vertex, center_of_load));
}

Vector_3 get_vertex_force(const Point_3& vertex, const std::vector<Point_3>& centers_of_load, const std::vector<double>& normalized_loads) {
    Vector_3 f(0,0,0);
    const auto size = centers_of_load.size();
    for(int i = 0; i < size; ++i) f += (normalized_loads[i] - 1) * get_direction(vertex, centers_of_load[i]);
    return f;
}

// Domain = U partition;
bool in_tetrahedron(const Tetrahedron_3& tet, const Point_3& p){
    auto side = tet.bounded_side(p);
    return side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY;
}

struct Domain;
Vector_3 constraint_force(const Domain* d, const Point_3& p, const Vector_3& f);
double compute_mu(const Domain* d, double max_normalized_load);

struct Partition {
    const Domain* d;
    Point_3 v1, v2, v3, v4, v5, v6, v7, v8;
    Tetrahedron_3 t1, t2, t3, t4, t5, t6;

    Partition(const Domain* d, Point_3 v1, Point_3 v2, Point_3 v3, Point_3 v4,
              Point_3 v5, Point_3 v6, Point_3 v7, Point_3 v8)
              : d(d),
                v1(std::move(v1)), v2(std::move(v2)),
                v3(std::move(v3)), v4(std::move(v4)),
                v5(std::move(v5)), v6(std::move(v6)),
                v7(std::move(v7)), v8(std::move(v8)) {
        construct_tetrahedra();
    }

    void construct_tetrahedra() {
        t1 = Tetrahedron_3(v1, v2, v4, v8);
        t2 = Tetrahedron_3(v1, v2, v6, v8);
        t3 = Tetrahedron_3(v1, v3, v4, v8);
        t4 = Tetrahedron_3(v1, v3, v7, v8);
        t5 = Tetrahedron_3(v1, v5, v6, v8);
        t6 = Tetrahedron_3(v1, v5, v7, v8);

        //Union of tetrahedra must be equal to Partition, hence sum of volumes = volume of partition
        auto volume = std::abs(t1.volume()) + std::abs(t2.volume()) + std::abs(t3.volume()) + std::abs(t4.volume()) + std::abs(t5.volume()) + std::abs(t6.volume());
        assert(volume - this->get_volume() <= std::numeric_limits<double>::epsilon());
    }

    friend std::ostream &operator<<(std::ostream &os, const Partition &partition) {
        os << "v1: " << partition.v1 << " v2: " << partition.v2 << " v3: " << partition.v3 << " v4: " << partition.v4
           << " v5: " << partition.v5 << " v6: " << partition.v6 << " v7: " << partition.v7 << " v8: " << partition.v8;
        return os;
    }

    double get_volume() {
        return (v2.x()-v1.x()) * (v5.y()-v1.y()) * (v3.z()-v1.z());
    }

    bool is_inside(const Point_3& p) {
         return in_tetrahedron(t1, p) || in_tetrahedron(t2, p) || in_tetrahedron(t3, p) || in_tetrahedron(t4, p) || in_tetrahedron(t5, p) || in_tetrahedron(t6, p);
    }

    Point_3 move_vertex(const Point_3& vertex, const Vector_3& force, double mu){
        return vertex + mu*force;
    }

    template<class CartesianPointTransformer, class LoadComputer, class A>
    void move_vertices(const std::vector<double>& weights, const std::vector<A>& elements, MPI_Comm neighborhood) {
        LoadComputer lc;
        CartesianPointTransformer transformer;
        std::vector<Point_3> points(elements.size());

        double my_load = lc.compute_load(elements);
        std::transform(elements.cbegin(), elements.cend(), std::back_inserter(points), [&transformer](auto el){return transformer.transform(el);});

        int N; MPI_Comm_size(neighborhood, &N);
        auto all_loads = get_neighbors_load(my_load, neighborhood);
        auto avg_load  = std::accumulate(all_loads.cbegin(), all_loads.cend(), 0.0) / N;
        decltype(all_loads) normalized_loads(N);

        auto max_normalized_load = *std::max_element(normalized_loads.cbegin(), normalized_loads.cend());

        double mu = compute_mu(d, max_normalized_load);

        std::transform(all_loads.cbegin(), all_loads.cend(), std::back_inserter(normalized_loads), [&avg_load](auto v){return v / avg_load;});

        auto center_of_load = get_center_of_load(weights, points);
        auto all_cl = get_centers_of_load(center_of_load, neighborhood);

        auto f1 = get_vertex_force(v1, all_cl, normalized_loads);
        f1 = constraint_force(d, v1, f1);
        v1 = move_vertex(v1, f1, mu);
        auto f2 = get_vertex_force(v2, all_cl, normalized_loads);
        f2 = constraint_force(d, v2, f2);
        v2 = move_vertex(v2, f2, mu);
        auto f3 = get_vertex_force(v3, all_cl, normalized_loads);
        f3 = constraint_force(d, v3, f3);
        v3 = move_vertex(v3, f3, mu);
        auto f4 = get_vertex_force(v4, all_cl, normalized_loads);
        f4 = constraint_force(d, v4, f4);
        v4 = move_vertex(v4, f4, mu);
        auto f5 = get_vertex_force(v5, all_cl, normalized_loads);
        f5 = constraint_force(d, v5, f5);
        v5 = move_vertex(v5, f5, mu);
        auto f6 = get_vertex_force(v6, all_cl, normalized_loads);
        f6 = constraint_force(d, v6, f6);
        v6 = move_vertex(v6, f6, mu);

    }
};

struct Domain {
    double grid_cell_size = 1.0;

    const double get_grid_cell_size() const {
        return grid_cell_size;
    }

    Point_3 v1, v2, v3, v4, v5, v6, v7, v8;

    const double xmin() const {
        return v1.x();
    }
    const double ymin() const {
        return v1.y();
    }
    const double zmin() const {
        return v1.z();
    }
    const double xmax() const {
        return v8.x();
    }
    const double ymax() const {
        return v8.y();
    }
    const double zmax() const {
        return v8.z();
    }

    std::vector<Partition> partitions;

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



    void bootstrap_partitions(unsigned int nb_partitions){
        bootstrap_partitions(1, nb_partitions, v1,v2,v3,v4,v5,v6,v7,v8);
    }

    /**
     * What is initial partitioning?
     * @param nb_partitions
     */

    void bootstrap_partitions(int generated_parts, unsigned int nb_partitions, const Point_3& p1,const Point_3& p2,const Point_3& p3,const Point_3& p4,
                             const Point_3& p5,const Point_3& p6,const Point_3& p7,const Point_3& p8) {

        if(generated_parts == nb_partitions) {
            partitions.emplace_back(this, p1, p2, p3, p4, p5, p6, p7, p8);
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
    return sigma_max == 0 ? 0 : d->get_grid_cell_size()/(max_normalized_load-1);
}
int main() {
    MPI_Init(nullptr, nullptr);

    Domain d(10, 10, 10);
    d.grid_cell_size = 0.1;
    d.bootstrap_partitions(std::pow(2, 4));
    //d.print_partitions();

    return 0;
}




