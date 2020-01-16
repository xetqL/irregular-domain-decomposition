//
// Created by xetql on 11/6/19.
//

#include <gtest/gtest.h>
#include <Domain.hpp>
#include "Utils.hpp"
#include "GeometricUtils.hpp"
#include "spatial_elements.hpp"
#include "Cell.hpp"
#include "mpi.h"

namespace
{
#define COMPARE(a, b) \
    std::cout << a << " vs. " << b << std::endl


TEST(CellLocalID, Computation) {
    using Cell = mesh::Cell<elements::Element<3>>;

        Cell::get_cell_size()   = 0.5;
        auto step             = Cell::get_cell_size();

        const int DOMAIN_SIZE_X = 10;
        const int DOMAIN_SIZE_Y = 10;
        const int DOMAIN_SIZE_Z = 10;

        Cell::get_msx() = DOMAIN_SIZE_X / Cell::get_cell_size();
        Cell::get_msy() = DOMAIN_SIZE_Y / Cell::get_cell_size();
        Cell::get_msz() = DOMAIN_SIZE_Z / Cell::get_cell_size();
        auto msx              = Cell::get_msx();
        auto msy              = Cell::get_msy();
        auto msz              = Cell::get_msz();

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

    Cell cell(gid, Cell::get_msx(), Cell::get_msy(), Cell::get_msz(), bbox);
    type::DataIndex x,y,z;

    lb::cell_to_local_position(Cell::get_msx(), Cell::get_msy(), Cell::get_msz(), bbox, gid, &x, &y, &z);

    auto lid = mesh::compute_lid(Cell::get_msx(), Cell::get_msy(), Cell::get_msz(), gid, bbox);

    EXPECT_LT(lid, bbox.get_number_of_cells());
}

TEST(ParticleLocalID, Computation) {
    using Cell = mesh::Cell<elements::Element<3>>;

        Cell::get_cell_size()   = 0.5;
        auto step             = Cell::get_cell_size();

        const int DOMAIN_SIZE_X = 10;
        const int DOMAIN_SIZE_Y = 10;
        const int DOMAIN_SIZE_Z = 10;

        Cell::get_msx() = DOMAIN_SIZE_X / Cell::get_cell_size();
        Cell::get_msy() = DOMAIN_SIZE_Y / Cell::get_cell_size();
        Cell::get_msz() = DOMAIN_SIZE_Z / Cell::get_cell_size();
        auto msx              = Cell::get_msx();
        auto msy              = Cell::get_msy();
        auto msz              = Cell::get_msz();

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

    lb::Box3 bbox(vertices, Cell::get_cell_size());

    elements::Element<3> el({0,0,0}, {0,0,0});

    auto gid = lb::position_to_cell(el.position[0],el.position[1],el.position[2],
                                Cell::get_cell_size(),
                                Cell::get_msx(),
                                Cell::get_msy());

    auto lid = mesh::compute_lid(
            Cell::get_msx(),
            Cell::get_msy(),
            Cell::get_msz(),
            gid, bbox);
}

TEST(CellLocalID, areDifferent) {
    using Cell = mesh::Cell<elements::Element<3>>;

        Cell::get_cell_size()   = 0.5;
        auto step             = Cell::get_cell_size();

        const int DOMAIN_SIZE_X = 10;
        const int DOMAIN_SIZE_Y = 10;
        const int DOMAIN_SIZE_Z = 10;

        Cell::get_msx() = DOMAIN_SIZE_X / Cell::get_cell_size();
        Cell::get_msy() = DOMAIN_SIZE_Y / Cell::get_cell_size();
        Cell::get_msz() = DOMAIN_SIZE_Z / Cell::get_cell_size();
        auto msx              = Cell::get_msx();
        auto msy              = Cell::get_msy();
        auto msz              = Cell::get_msz();

    std::array<lb::Point_3, 8> vertices = {
            lb::Point_3(0,  5,  0),
            lb::Point_3(5,  5,  0),
            lb::Point_3(0, 10,  0),
            lb::Point_3(5.29255, 10,  0),
            lb::Point_3(0,  5, 5),
            lb::Point_3(5,  5, 5),
            lb::Point_3(0, 10, 5),
            lb::Point_3(5, 10, 5),
    };

    lb::Box3 bbox(vertices, Cell::get_cell_size());

    elements::Element<3> el({0,0,0}, {0,0,0});

    type::DataIndex gid1 = 390;
    type::DataIndex gid2 = 600;

    auto lid1 = mesh::compute_lid(
            Cell::get_msx(),
            Cell::get_msy(),
            Cell::get_msz(),
            gid1, bbox);

    auto lid2 = mesh::compute_lid(
            Cell::get_msx(),
            Cell::get_msy(),
            Cell::get_msz(),
            gid2, bbox);

    EXPECT_NE(gid1, gid2);
    EXPECT_NE(lid1, lid2);
}


TEST(CellLocalID, areDifferent2) {
    using Cell = mesh::Cell<elements::Element<3>>;

    Cell::get_cell_size()   = 0.5;
    auto step             = Cell::get_cell_size();

    const int DOMAIN_SIZE_X = 10;
    const int DOMAIN_SIZE_Y = 10;
    const int DOMAIN_SIZE_Z = 10;

    Cell::get_msx() = DOMAIN_SIZE_X / Cell::get_cell_size();
    Cell::get_msy() = DOMAIN_SIZE_Y / Cell::get_cell_size();
    Cell::get_msz() = DOMAIN_SIZE_Z / Cell::get_cell_size();
    auto msx              = Cell::get_msx();
    auto msy              = Cell::get_msy();
    auto msz              = Cell::get_msz();

    std::array<lb::Point_3, 8> vertices = {
            lb::Point_3(0,  5,  0),
            lb::Point_3(5,  5,  0),
            lb::Point_3(0, 10,  0),
            lb::Point_3(5.29255, 10,  0),
            lb::Point_3(0,  5, 5),
            lb::Point_3(5,  5, 5),
            lb::Point_3(0, 10, 5),
            lb::Point_3(5, 10, 5),
    };

    lb::Box3 bbox(vertices, Cell::get_cell_size());

    type::DataIndex gid1 = 12;
    type::DataIndex gid2 = 20;
    type::DataIndex x,y,z;
    std::tie(x,y,z) = lb::cell_to_global_position(msx,msy,msz,gid1);
    Cell c1(gid1, Cell::get_msx(), Cell::get_msy(), Cell::get_msz(), bbox);
    Cell c2(gid2, Cell::get_msx(), Cell::get_msy(), Cell::get_msz(), bbox);
    EXPECT_NE(c1.get_lid(), c2.get_lid());
}



TEST(CellLocalID, RearrangingCell) {
    MPI_Init(nullptr, nullptr);
    using Cell = mesh::Cell<elements::Element<3>>;

    Cell::get_cell_size()   = 0.5;

    const type::Real DOMAIN_SIZE_X = 10;
    const type::Real DOMAIN_SIZE_Y = 10;
    const type::Real DOMAIN_SIZE_Z = 10;

    Cell::get_msx() = DOMAIN_SIZE_X / Cell::get_cell_size();
    Cell::get_msy() = DOMAIN_SIZE_Y / Cell::get_cell_size();
    Cell::get_msz() = DOMAIN_SIZE_Z / Cell::get_cell_size();

    int my_rank = 0, worldsize = 8;
    int procs_x = std::cbrt(worldsize), procs_y = std::cbrt(worldsize), procs_z = std::cbrt(worldsize);

    auto cell_per_col     =  Cell::get_msx() / procs_x;

    type::DataIndex x_proc_idx, y_proc_idx, z_proc_idx;

    lb::linear_to_grid(my_rank, procs_x, procs_y, x_proc_idx, y_proc_idx, z_proc_idx);

    auto my_cells = mesh::generate_lattice_single_type<int>(
            Cell::get_msx(), Cell::get_msy(), Cell::get_msz(),
            x_proc_idx, y_proc_idx, z_proc_idx,
            cell_per_col, cell_per_col, cell_per_col,
            mesh::TCellType::REAL_CELL);

    lb::Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z, &Cell::get_cell_size());

    d.bootstrap_partitions(worldsize);
    auto part = d.get_my_partition(my_rank);
    lb::Box3 bbox(part.vertices, Cell::get_cell_size());

    std::cout << bbox << std::endl;
    for(int i = 0; i < my_cells.size(); ++i) {
        auto cell = my_cells[i];
        auto lid = cell.get_lid(Cell::get_msx(), Cell::get_msy(), Cell::get_msz(), bbox);
        EXPECT_EQ(i, lid);
    }
}


TEST(DomainSizeX, Check) {
    using Cell = mesh::Cell<elements::Element<3>>;

    Cell::get_cell_size()   = 0.5;

    const type::Real DOMAIN_SIZE_X = 10;
    const type::Real DOMAIN_SIZE_Y = 10;
    const type::Real DOMAIN_SIZE_Z = 10;

    Cell::get_msx() = DOMAIN_SIZE_X / Cell::get_cell_size();
    Cell::get_msy() = DOMAIN_SIZE_Y / Cell::get_cell_size();
    Cell::get_msz() = DOMAIN_SIZE_Z / Cell::get_cell_size();

    int my_rank = 0, worldsize = 8;
    int procs_x = std::cbrt(worldsize), procs_y = std::cbrt(worldsize), procs_z = std::cbrt(worldsize);

    auto cell_per_col     =  Cell::get_msx() / procs_x;

    type::DataIndex x_proc_idx, y_proc_idx, z_proc_idx;

    lb::linear_to_grid(my_rank, procs_x, procs_y, x_proc_idx, y_proc_idx, z_proc_idx);

    auto my_cells = mesh::generate_lattice_single_type<int>(
            Cell::get_msx(), Cell::get_msy(), Cell::get_msz(),
            x_proc_idx, y_proc_idx, z_proc_idx,
            cell_per_col, cell_per_col, cell_per_col,
            mesh::TCellType::REAL_CELL);

    lb::Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z, &Cell::get_cell_size());

    d.bootstrap_partitions(worldsize);
    auto part = d.get_my_partition(my_rank);
    EXPECT_EQ(d.get_bounding_box().size_x, Cell::get_msx());
    EXPECT_EQ(d.get_bounding_box().size_y, Cell::get_msy());
    EXPECT_EQ(d.get_bounding_box().size_z, Cell::get_msz());
}

TEST(DomainSizeX, ZeroGidZeroLid) {
    using Cell = mesh::Cell<elements::Element<3>>;

    Cell::get_cell_size()   = 0.5;

    const type::Real DOMAIN_SIZE_X = 10;
    const type::Real DOMAIN_SIZE_Y = 10;
    const type::Real DOMAIN_SIZE_Z = 10;

    Cell::get_msx() = DOMAIN_SIZE_X / Cell::get_cell_size();
    Cell::get_msy() = DOMAIN_SIZE_Y / Cell::get_cell_size();
    Cell::get_msz() = DOMAIN_SIZE_Z / Cell::get_cell_size();

    int my_rank = 0, worldsize = 8;
    int procs_x = std::cbrt(worldsize), procs_y = std::cbrt(worldsize), procs_z = std::cbrt(worldsize);

    auto cell_per_col     =  Cell::get_msx() / procs_x;

    type::DataIndex x_proc_idx, y_proc_idx, z_proc_idx;

    lb::linear_to_grid(my_rank, procs_x, procs_y, x_proc_idx, y_proc_idx, z_proc_idx);

    auto my_cells = mesh::generate_lattice_single_type<int>(
            Cell::get_msx(), Cell::get_msy(), Cell::get_msz(),
            x_proc_idx, y_proc_idx, z_proc_idx,
            cell_per_col, cell_per_col, cell_per_col,
            mesh::TCellType::REAL_CELL);

    lb::Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z, &Cell::get_cell_size());

    d.bootstrap_partitions(worldsize);
    auto part = d.get_my_partition(my_rank);

    EXPECT_EQ(mesh::compute_lid(Cell::get_msx(), Cell::get_msy(), Cell::get_msz(), 0, d.get_bounding_box()), 0);
}

} //end of namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}