//
// Created by xetql on 10/21/19.
//

#ifndef ADLBIRREG_DOMAIN_HPP
#define ADLBIRREG_DOMAIN_HPP

#include <Partition.hpp>

namespace lb{

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

    inline const Box3 get_box() const {
        return {xmin(), ymin(), zmin(), xmax(), ymax(), zmax()};
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
                    partitions.emplace_back(Point_3(k ,j, i), get_box(), grid_cell_size,
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
}
#endif //ADLBIRREG_DOMAIN_HPP
