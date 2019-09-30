#include <iostream>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Aff_transformation_3.h>

inline int bitselect(int condition, int truereturnvalue, int falsereturnvalue)
{
    return (truereturnvalue & -condition) | (falsereturnvalue & ~(-condition)); //a when TRUE and b when FintLSE
}

using Kernel = CGAL::Simple_cartesian<double> ;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;
using Tetrahedron_3 = Kernel::Tetrahedron_3;
using Transformation = CGAL::Aff_transformation_3<Kernel>;

struct Partition {
    Point_3 v1,v2,v3,v4,v5,v6,v7,v8; //vertex of spatial partition, v1 is upper right

    Tetrahedron_3 t1, t2, t3, t4, t5, t6;

    Partition(const Point_3& v1, const Point_3& v2, const Point_3& v3, const Point_3& v4,
              const Point_3& v5, const Point_3& v6, const Point_3& v7, const Point_3& v8)
              : v1(v1), v2(v2), v3(v3), v4(v4), v5(v5), v6(v6), v7(v7), v8(v8) {
        construct_tetrahedra();
    }

    void construct_tetrahedra() {
        t1 = Tetrahedron_3(v1, v2, v4, v8);
        t2 = Tetrahedron_3(v1, v2, v6, v8);
        t3 = Tetrahedron_3(v1, v3, v4, v8);
        t4 = Tetrahedron_3(v1, v3, v7, v8);
        t5 = Tetrahedron_3(v1, v5, v6, v8);
        t6 = Tetrahedron_3(v1, v5, v7, v8);

        auto volume = t1.volume() + t2.volume() + t3.volume() + t4.volume() + t5.volume() + t6.volume();
        std::cout << t1.volume() << std::endl;
        auto ref_volume = v8.x()*v8.y()*v8.z();
        if(volume != ref_volume) {
            std::cout << "ref volume: " << ref_volume << " effective volume: " << volume << std::endl;
        }
    };

    friend std::ostream &operator<<(std::ostream &os, const Partition &partition) {
        os << "v1: " << partition.v1 << " v2: " << partition.v2 << " v3: " << partition.v3 << " v4: " << partition.v4
           << " v5: " << partition.v5 << " v6: " << partition.v6 << " v7: " << partition.v7 << " v8: " << partition.v8;
        return os;
    }

};

struct Domain {
    Point_3 v1,v2,v3,v4,v5,v6,v7,v8;
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
            v5(p5),v6(p6),v7(p7),v8(p8){

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
            partitions.emplace_back(p1, p2, p3, p4, p5, p6, p7, p8);
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

// Domain = U partition;
bool in_tetrahedron(const Tetrahedron_3& tet, const Point_3& p){
    auto side = tet.bounded_side(p);
    return side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY;
}

int main() {
    Domain d(10, 10, 10);
    d.bootstrap_partitions(std::pow(2, 4));
    //d.print_partitions();
    return 0;
}




