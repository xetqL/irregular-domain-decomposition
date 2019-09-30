#include <iostream>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Aff_transformation_3.h>

using Kernel = CGAL::Simple_cartesian<double> ;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;
using Tetrahedron_3 = Kernel::Tetrahedron_3;
using Transformation = CGAL::Aff_transformation_3<Kernel>;
struct Partition {
    Point_3 v1,v2,v3,v4,v5,v6,v7,v8; //vertex of spatial partition, v1 is upper right
    Tetrahedron_3 t1,t2,t3,t4,t5,t6;

    Partition(const Point_3& v1, const Point_3& v2, const Point_3& v3, const Point_3& v4,
              const Point_3& v5, const Point_3& v6, const Point_3& v7, const Point_3& v8)
              : v1(v1), v2(v2), v3(v3), v4(v4), v5(v5), v6(v6), v7(v7), v8(v8) {}

    void construct_tetrahedra(){

    };


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
 * 3-- 5 ----4 --6
 * |  /      |  /
 * Y Z       | /
 * |/        |/
 * 1----X----2
 */
    Domain (const Point_3& p1,const Point_3& p2,const Point_3& p3,const Point_3& p4,
            const Point_3& p5,const Point_3& p6,const Point_3& p7,const Point_3& p8):
            v1(p1),v2(p2),v3(p3),v4(p4),
            v5(p5),v6(p6),v7(p7),v8(p8){

    }



    /**
     * What is initial partitioning?
     * @param nb_partitions
     */
    void bootstrap_partitions(int generated_parts, int nb_partitions, const Point_3& p1,const Point_3& p2,const Point_3& p3,const Point_3& p4,
                             const Point_3& p5,const Point_3& p6,const Point_3& p7,const Point_3& p8) {
        if(generated_parts == nb_partitions/2) {
            partitions.emplace_back(p1, p2, p3, p4, p5, p6, p7, p8);
            return;
        };

        int largest_dim = largest_dimension(p1,p2,p3,p4,p5,p6,p7,p8);

        if(largest_dim == 0) { //cut the domain along X
            Transformation translate_right(CGAL::TRANSLATION, Vector_3(p1.x() - p2.x() / 2, 0, 0));
            Transformation translate_left(CGAL::TRANSLATION, Vector_3(-p1.x() - p2.x() / 2, 0, 0));
            bootstrap_partitions(generated_parts * 2, nb_partitions, p1, translate_left(p2), p3, translate_left(p4), p5, translate_left(p6), p7, translate_left(p8));
            bootstrap_partitions(generated_parts * 2, nb_partitions, translate_right(p1), (p2), translate_right(p3), (p4), translate_right(p5), (p6), translate_right(p7), (p8));
        } else if(largest_dim == 1) { //cut the domain along Y
            Transformation translate_right(CGAL::TRANSLATION, Vector_3(0, 0, 0));
            Transformation translate_left(CGAL::TRANSLATION, Vector_3(0, 0, 0));
            bootstrap_partitions(generated_parts * 2, nb_partitions, p1, translate_left(p2), p3, translate_left(p4), p5, translate_left(p6), p7, translate_left(p8));
            bootstrap_partitions(generated_parts * 2, nb_partitions, translate_right(p1), (p2), translate_right(p3), (p4), translate_right(p5), (p6), translate_right(p7), (p8));
        } else { // cut the domain along Z
            Transformation translate_right(CGAL::TRANSLATION, Vector_3(p1.x() - p2.x() / 2, 0, 0));
            Transformation translate_left(CGAL::TRANSLATION, Vector_3(-p1.x() - p2.x() / 2, 0, 0));
            bootstrap_partitions(generated_parts * 2, nb_partitions, p1, translate_left(p2), p3, translate_left(p4), p5, translate_left(p6), p7, translate_left(p8));
            bootstrap_partitions(generated_parts * 2, nb_partitions, translate_right(p1), (p2), translate_right(p3), (p4), translate_right(p5), (p6), translate_right(p7), (p8));
        }

    }



    int largest_dimension(const Point_3& p1,const Point_3& p2,const Point_3& p3,const Point_3& p4,
                          const Point_3& p5,const Point_3& p6,const Point_3& p7,const Point_3& p8){
        auto x = p1.x() - p2.x();
        auto y = p1.y() - p3.y();
        auto z = p1.z() - p5.z();
        if(x > y && x > z) return 0;
        if(y > x && y > z) return 1;
        if(z > x && z > y) return 2;
        return 0;
    }

};

// Domain = U partition;

bool in_tetrahedron(const Tetrahedron_3& tet, const Point_3& p){
    auto side = tet.bounded_side(p);
    return side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY;
}

int main() {
    Point_3 p = {1,2,3};

    return 0;
}




