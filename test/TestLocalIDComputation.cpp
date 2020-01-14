//
// Created by xetql on 11/6/19.
//

#include <gtest/gtest.h>
#include "Utils.hpp"
#include "GeometricUtils.hpp"
#include "spatial_elements.hpp"
#include "Cell.hpp"

namespace
{

TEST(CellLocalID, Computation) {
    using Cell = mesh::Cell<elements::Element<3>>;
    Cell::get_cell_size() = 0.5;
    Cell::get_msx()       = 20;
    Cell::get_msy()       = 20;
    Cell::get_msz()       = 20;

    std::array<lb::Point_3, 8> vertices = {
            lb::Point_3(0,  5,  5),
            lb::Point_3(5,  5,  5),
            lb::Point_3(0, 10,  5),
            lb::Point_3(5, 10,  5),
            lb::Point_3(0,  5, 10),
            lb::Point_3(5,  5, 10),
            lb::Point_3(0, 10, 10),
            lb::Point_3(5.26632, 10, 10),
    };

    type::DataIndex gid = 7990;

    lb::Box3 bbox(vertices, Cell::get_cell_size());

    Cell cell(gid, bbox);
    auto cbox = cell.as_box();
    type::DataIndex x,y,z;

    lb::cell_to_local_position(Cell::get_msx(), Cell::get_msy(), Cell::get_msz(), bbox, gid, &x, &y, &z);

    auto lid = mesh::compute_lid(Cell::get_msx(), Cell::get_msy(), Cell::get_msz(), gid, bbox);

    EXPECT_LT(lid, bbox.get_number_of_cells());
}

TEST(ParticleLocalID, Computation) {
    using Cell = mesh::Cell<elements::Element<3>>;
    Cell::get_cell_size() = 0.5;
    Cell::get_msx()       = 20;
    Cell::get_msy()       = 20;
    Cell::get_msz()       = 20;

    std::array<lb::Point_3, 8> vertices = {
            lb::Point_3(0,  4.96129,  0),
            lb::Point_3(5.38698,  5.09485,  0),
            lb::Point_3(0, 10,  0),
            lb::Point_3(5.27085, 10,  0),
            lb::Point_3(0,  5, 5),
            lb::Point_3(5,  5, 5),
            lb::Point_3(0, 10, 5),
            lb::Point_3(5, 10, 5),
    };

    //type::DataIndex gid = 7990;

    lb::Box3 bbox(vertices, Cell::get_cell_size());

    //std::cout << "Bounding box is "<< bbox << std::endl;

    elements::Element<3> el({0,0,0}, {0,0,0});

    auto gid = lb::position_to_cell(el.position[0],el.position[1],el.position[2],
                                Cell::get_cell_size(),
                                Cell::get_msx(),
                                Cell::get_msy());
    std::cout << gid << std::endl;
    std::cout << bbox << std::endl;
    auto lid = mesh::compute_lid(
            Cell::get_msx(),
            Cell::get_msy(),
            Cell::get_msz(),
            gid, bbox);

    std::cout << lid << std::endl;
}

} //end of namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}