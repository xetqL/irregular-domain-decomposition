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

    std::cout << lid << std::endl;
}

TEST(CellLocalID, areDifferent) {
    using Cell = mesh::Cell<elements::Element<3>>;
    Cell::get_cell_size() = 0.5;
    Cell::get_msx()       = 20;
    Cell::get_msy()       = 20;
    Cell::get_msz()       = 20;

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
    Cell::get_cell_size() = 0.5;
    Cell::get_msx()       = 20;
    Cell::get_msy()       = 20;
    Cell::get_msz()       = 20;

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

    type::DataIndex gid1 = 390;
    type::DataIndex gid2 = 600;

    Cell c1(gid1, bbox);
    Cell c2(gid2, bbox);
    EXPECT_NE(c1.get_lid(bbox), c2.get_lid(bbox));
}

std::vector<mesh::Cell<elements::Element<3>>> generate_lattice_single_type( long long int msx, long long int msy, long long int msz,
                                                                            long long int x_proc_idx, long long int y_proc_idx, long long int z_proc_idx,
                                                                            long long int cell_in_my_cols, long long int cell_in_my_rows, long long int cell_in_my_depth,
                                                                            mesh::TCellType type) {

    auto cell_per_process = cell_in_my_cols * cell_in_my_rows;
    std::vector<mesh::Cell<elements::Element<3>>> my_cells; my_cells.reserve(cell_per_process);

    auto x_shift = (cell_in_my_rows *  x_proc_idx);
    auto y_shift = (cell_in_my_cols *  y_proc_idx);
    auto z_shift = (cell_in_my_depth * z_proc_idx);
    for(int z = 0; z < cell_in_my_depth; ++z) {
        for(int y = 0; y < cell_in_my_cols; ++y) {
            for(int x = 0; x < cell_in_my_rows; ++x) {
                auto gid = (x_shift + x) + (y_shift + y) * msx + (z_shift + z) * msx * msy;
                my_cells.emplace_back(gid, type);
            }
        }
    }

    return my_cells;
}

TEST(CellLocalID, RearrangingCell) {
    MPI_Init(nullptr, nullptr);
    using Cell = mesh::Cell<elements::Element<3>>;

    Cell::get_cell_size()   = 1.0;

    const int DOMAIN_SIZE_X = 20;
    const int DOMAIN_SIZE_Y = 20;
    const int DOMAIN_SIZE_Z = 20;

    Cell::get_msx() = DOMAIN_SIZE_X / Cell::get_cell_size();
    Cell::get_msy() = DOMAIN_SIZE_Y / Cell::get_cell_size();
    Cell::get_msz() = DOMAIN_SIZE_Z / Cell::get_cell_size();

    int my_rank = 0, worldsize = 8;
    int procs_x = std::cbrt(worldsize), procs_y = std::cbrt(worldsize), procs_z = std::cbrt(worldsize);

    auto cell_per_proc    = (Cell::get_msx() / procs_x) * (Cell::get_msy() / procs_y) * (Cell::get_msz() / procs_z);
    auto cell_per_col     =  Cell::get_msx() / procs_x;
    auto total_cell_count =  Cell::get_msx() * Cell::get_msy() * Cell::get_msz();

    type::DataIndex x_proc_idx, y_proc_idx, z_proc_idx;

    lb::linear_to_grid(my_rank, procs_x, procs_y, x_proc_idx, y_proc_idx, z_proc_idx);
    std::cout << "ProcX: " << x_proc_idx << std::endl;
    std::cout << "ProcY: " << y_proc_idx << std::endl;
    std::cout << "ProcZ: " << z_proc_idx << std::endl;
    std::cout << "CellPerProcs: " << cell_per_proc << std::endl;
    std::cout << "CellX: " << cell_per_col << std::endl;
    std::cout << "CellY: " << cell_per_col << std::endl;
    std::cout << "CellZ: " << cell_per_col << std::endl;
    std::cout << "Total cell count: " << total_cell_count << std::endl;

    auto my_cells = generate_lattice_single_type(Cell::get_msx(), Cell::get_msy(), Cell::get_msz(), x_proc_idx, y_proc_idx, z_proc_idx, cell_per_col, cell_per_col, cell_per_col, mesh::TCellType::REAL_CELL);
    lb::Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z, &Cell::get_cell_size());

    d.bootstrap_partitions(worldsize);
    auto part = d.get_my_partition(my_rank);
    lb::Box3 bbox(part.vertices, Cell::get_cell_size());

    std::cout << bbox << std::endl;
    for(int i = 0; i < my_cells.size(); ++i) {
        auto cell = my_cells[i];
        auto lid = cell.get_lid(bbox);
        EXPECT_EQ(i, lid);
    }
}
} //end of namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}