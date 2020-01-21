//
// Created by xetql on 10/21/19.
//

#ifndef ADLBIRREG_DOMAIN_HPP
#define ADLBIRREG_DOMAIN_HPP

#include <Partition.hpp>

namespace lb {


struct Domain {

    Point_3 v1, v2, v3, v4, v5, v6, v7, v8;
    const type::Real* grid_cell_size;

    std::vector<Partition> partitions;
    int num_part = 0;

    //explicit Domain (std::vector<Partition>& partitions): partitions(partitions) {}
    Domain (const type::Real x, const type::Real y, const type::Real z, const type::Real* cell_size):
            v1(0, 0, 0),
            v2(x, 0, 0),
            v3(0, y, 0),
            v4(x, y, 0),
            v5(0, 0, z),
            v6(x, 0, z),
            v7(0, y, z),
            v8(x, y, z),
            grid_cell_size(cell_size)
    {}

    const type::Real get_grid_cell_size() const {
        return *grid_cell_size;
    }

    inline const type::Real xmin() const {
        return v1.x();
    }

    inline const type::Real ymin() const {
        return v1.y();
    }

    inline const type::Real zmin() const {
        return v1.z();
    }

    inline const type::Real xmax() const {
        return v8.x();
    }

    inline const type::Real ymax() const {
        return v8.y();
    }

    inline const type::Real zmax() const {
        return v8.z();
    }

    inline const type::Real xsize() const {
        return xmax()-xmin();
    }

    inline const type::Real ysize() const {
        return ymax()-ymin();
    }

    inline const type::Real zsize() const {
        return zmax()-zmin();
    }

    inline const Box3 get_bounding_box() const {
        return {xmin(), xmax(), ymin(), ymax(), zmin(), zmax(), *grid_cell_size};
    }

    Partition& get_partition(int rank);

    void bootstrap_partitions(unsigned int nb_partitions) {
        num_part = nb_partitions;
        bootstrap_partitions_cubical(nb_partitions, v1, v2, v3, v4, v5, v6, v7, v8);

    }

    /**
     * number of partitions must be a perfect cube
     */
    void bootstrap_partitions_cubical(unsigned int nb_partitions, const Point_3& p1,const Point_3& p2,const Point_3& p3,const Point_3& p4,
                                      const Point_3& p5,const Point_3& p6,const Point_3& p7,const Point_3& p8) {
        const type::Real proc_per_row = std::cbrt(nb_partitions);
        const int    row_size     = (int) proc_per_row;
        //const type::Real p_m1         = ((proc_per_row - 1) / proc_per_row);

        Transformation translate_right(  CGAL::TRANSLATION, Vector_3(std::abs(p1.x() - p2.x()) / proc_per_row, 0, 0));
        Transformation translate_fullleft(  CGAL::TRANSLATION, Vector_3(-(p2.x()-p1.x()), 0, 0));

        Transformation translate_up(     CGAL::TRANSLATION, Vector_3(0, 0, std::abs(p1.z() - p5.z()) / proc_per_row));
        Transformation translate_fulldown(     CGAL::TRANSLATION, Vector_3(0, 0, -(p5.z() - p1.z())));

        Transformation translate_forward(CGAL::TRANSLATION, Vector_3(0, std::abs(p1.y() - p3.y()) / proc_per_row,   0));
        Transformation translate_backward(CGAL::TRANSLATION, Vector_3(0, -std::abs(p3.y() - p1.y()) / proc_per_row, 0));

        Transformation translate_fullbackward(CGAL::TRANSLATION, Vector_3(0, -std::abs(p1.y() - p3.y()), 0));
        Transformation translate_fullforward(CGAL::TRANSLATION, Vector_3(0, std::abs(p1.y() - p3.y()), 0));

        std::array<Point_3, 8> partition_vertices =
                {p1,
                 translate_right(p1),
                 translate_forward(p1),
                 translate_forward(translate_right(p1)),
                 translate_up(p1),
                 translate_up(translate_right(p1)),
                 translate_forward(translate_up(p1)),
                 translate_forward(translate_up(translate_right(p1)))};

        for(int i = 0; i < row_size; ++i) {
            for(int j = 0; j < row_size; ++j) {
                for(int k = 0; k < row_size; ++k) {
                    partitions.emplace_back(Point_3(k ,j, i), get_bounding_box(), grid_cell_size,
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
