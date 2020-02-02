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
    MPI_Init(nullptr, nullptr);

    using Cell = mesh::Cell<elements::Element<3>>;
    mesh::GridParams& gp = mesh::GridParams::get_instance();

    gp.set_grid_resolution(0.5);
    auto step             = gp.get_grid_resolution();

    const type::Real DOMAIN_SIZE_X = 10;
    const type::Real DOMAIN_SIZE_Y = 10;
    const type::Real DOMAIN_SIZE_Z = 10;

    gp.set_grid_index_dimensions(
            DOMAIN_SIZE_X / gp.get_grid_resolution(),
            DOMAIN_SIZE_Y / gp.get_grid_resolution(),
            DOMAIN_SIZE_Z / gp.get_grid_resolution());

    auto msx  = gp.get_cell_number_x();
    auto msy  = gp.get_cell_number_y();
    auto msz  = gp.get_cell_number_z();

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

    lb::Box3 bbox(vertices, gp.get_grid_resolution());

    Cell cell(gid, gp.get_cell_number_x(), gp.get_cell_number_y(), gp.get_cell_number_z(), bbox);
    type::DataIndex x,y,z;

    lb::cell_to_local_position(gp.get_cell_number_x(), gp.get_cell_number_y(), gp.get_cell_number_z(), bbox, gid, &x, &y, &z);

    auto lid = mesh::compute_lid(gp.get_cell_number_x(), gp.get_cell_number_y(), gp.get_cell_number_z(), gid, bbox);

    EXPECT_LT(lid, bbox.get_number_of_cells());

    if(lid >= bbox.get_number_of_cells()) {
        type::DataIndex x, y, z;
        lb::cell_to_local_position(gp.get_cell_number_x(), gp.get_cell_number_y(), gp.get_cell_number_z(), bbox, gid, &x, &y, &z);
        std::cout << " >> " << x << " " << y << " " << z << std::endl;
        std::cout << bbox << std::endl;
        std::cout << cell << " " << lid << std::endl;
    }
}

TEST(ParticleLocalID, Computation) {
    using Cell = mesh::Cell<elements::Element<3>>;
    mesh::GridParams& gp = mesh::GridParams::get_instance();
    gp.set_grid_resolution(0.5);
    auto step             = gp.get_grid_resolution();

    const int DOMAIN_SIZE_X = 10;
    const int DOMAIN_SIZE_Y = 10;
    const int DOMAIN_SIZE_Z = 10;

    gp.set_grid_index_dimensions(
            DOMAIN_SIZE_X / gp.get_grid_resolution(),
            DOMAIN_SIZE_Y / gp.get_grid_resolution(),
            DOMAIN_SIZE_Z / gp.get_grid_resolution());

    auto msx              = gp.get_cell_number_x();
    auto msy              = gp.get_cell_number_y();
    auto msz              = gp.get_cell_number_z();

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

    lb::Box3 bbox(vertices, gp.get_grid_resolution());

    elements::Element<3> el({0,0,0}, {0,0,0});

    auto gid = lb::position_to_cell(el.position[0],el.position[1],el.position[2],
                                gp.get_grid_resolution(),
                                gp.get_cell_number_x(),
                                gp.get_cell_number_y());

    auto lid = mesh::compute_lid(
            gp.get_cell_number_x(),
            gp.get_cell_number_y(),
            gp.get_cell_number_z(),
            gid, bbox);
}

TEST(CellLocalID, areDifferent) {
    using Cell = mesh::Cell<elements::Element<3>>;
    mesh::GridParams& gp = mesh::GridParams::get_instance();
    gp.set_grid_resolution(0.5);
    auto step             = gp.get_grid_resolution();

    const int DOMAIN_SIZE_X = 10;
    const int DOMAIN_SIZE_Y = 10;
    const int DOMAIN_SIZE_Z = 10;

    gp.set_grid_index_dimensions(DOMAIN_SIZE_X / gp.get_grid_resolution(), DOMAIN_SIZE_Y / gp.get_grid_resolution(),DOMAIN_SIZE_Z / gp.get_grid_resolution());

    auto msx              = gp.get_cell_number_x();
    auto msy              = gp.get_cell_number_y();
    auto msz              = gp.get_cell_number_z();

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

    lb::Box3 bbox(vertices, gp.get_grid_resolution());

    elements::Element<3> el({0,0,0}, {0,0,0});

    type::DataIndex gid1 = 390;
    type::DataIndex gid2 = 600;

    auto lid1 = mesh::compute_lid(
            gp.get_cell_number_x(),
            gp.get_cell_number_y(),
            gp.get_cell_number_z(),
            gid1, bbox);

    auto lid2 = mesh::compute_lid(
            gp.get_maximum_index_x(),
            gp.get_maximum_index_y(),
            gp.get_maximum_index_z(),
            gid2, bbox);

    EXPECT_NE(gid1, gid2);
    EXPECT_NE(lid1, lid2);
}


TEST(CellLocalID, areDifferent2) {
    using Cell = mesh::Cell<elements::Element<3>>;
    mesh::GridParams& gp = mesh::GridParams::get_instance();
    gp.set_grid_resolution(0.5);
    auto step             = gp.get_grid_resolution();

    const int DOMAIN_SIZE_X = 10;
    const int DOMAIN_SIZE_Y = 10;
    const int DOMAIN_SIZE_Z = 10;

    gp.set_grid_index_dimensions(DOMAIN_SIZE_X / gp.get_grid_resolution(),DOMAIN_SIZE_Y / gp.get_grid_resolution(),DOMAIN_SIZE_Z / gp.get_grid_resolution());
    auto msx              = gp.get_cell_number_x();
    auto msy              = gp.get_cell_number_y();
    auto msz              = gp.get_cell_number_z();

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

    lb::Box3 bbox(vertices, gp.get_grid_resolution());

    type::DataIndex gid1 = 12;
    type::DataIndex gid2 = 20;
    type::DataIndex x,y,z;
    std::tie(x,y,z) = lb::cell_to_global_position(msx,msy,msz,gid1);
    Cell c1(gid1, gp.get_cell_number_x(), gp.get_cell_number_y(), gp.get_cell_number_z(), bbox);
    Cell c2(gid2, gp.get_cell_number_x(), gp.get_cell_number_y(), gp.get_cell_number_z(), bbox);
    EXPECT_NE(c1.get_lid(), c2.get_lid());
}

TEST(Cell, ShouldHaveGidWithinDimensionSize) {

    using Cell = mesh::Cell<elements::Element<3>>;
    mesh::GridParams& gp = mesh::GridParams::get_instance();

    int worldsize = 1,
        procs_x       = std::cbrt(worldsize),
        procs_y       = std::cbrt(worldsize),
        procs_z       = std::cbrt(worldsize);

    gp.set_grid_resolution(0.0625);

    const type::Real _DOMAIN_SIZE_X = 20.625; type::Real DOMAIN_SIZE_X;
    const type::Real _DOMAIN_SIZE_Y = 20.625; type::Real DOMAIN_SIZE_Y;
    const type::Real _DOMAIN_SIZE_Z = 20.625; type::Real DOMAIN_SIZE_Z;

    lb::adjust_simulation_size(_DOMAIN_SIZE_X,_DOMAIN_SIZE_Y,_DOMAIN_SIZE_Z,
                               procs_x, procs_y, procs_z,
                               gp.get_grid_resolution(),
                               &DOMAIN_SIZE_X, &DOMAIN_SIZE_Y, &DOMAIN_SIZE_Z);
    gp.set_grid_index_dimensions(
            DOMAIN_SIZE_X / gp.get_grid_resolution(),
            DOMAIN_SIZE_Y / gp.get_grid_resolution(),
            DOMAIN_SIZE_Z / gp.get_grid_resolution());

    lb::Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z, &gp.get_grid_resolution());

    auto cell_per_col = gp.get_cell_number_x() / procs_x;

    type::DataIndex x_proc_idx, y_proc_idx, z_proc_idx;

    d.bootstrap_partitions(worldsize);

    for(int cpu = 0; cpu < worldsize; ++cpu){
        //std::cout << cpu << std::endl;
        lb::linear_to_grid(cpu, procs_x, procs_y, x_proc_idx, y_proc_idx, z_proc_idx);
        auto my_cells = mesh::generate_lattice_single_type<int>(
                gp.get_cell_number_x(), gp.get_cell_number_y(), gp.get_cell_number_z(),
                x_proc_idx, y_proc_idx, z_proc_idx,
                cell_per_col, cell_per_col, cell_per_col,
                mesh::TCellType::REAL_CELL);
        auto part = d.get_partition(cpu);
        lb::Box3 bbox(part.vertices, gp.get_grid_resolution());
        const auto size = my_cells.size();
        for(int i = 0; i < size; ++i) {
            auto cell = my_cells[i];
            auto gid = cell.gid;
            EXPECT_LT(gid, gp.get_cell_number_x() * gp.get_cell_number_y() * gp.get_cell_number_z());
        }
    }

}

TEST(Cell, ShouldCoverRangeOfAllLid) {

    using Cell = mesh::Cell<elements::Element<3>>;
    mesh::GridParams& gp = mesh::GridParams::get_instance();
    int my_rank = 0, worldsize = 64,
        procs_x       = std::cbrt(worldsize),
        procs_y       = std::cbrt(worldsize),
        procs_z       = std::cbrt(worldsize);
    gp.set_grid_resolution(0.0625);

    const type::Real _DOMAIN_SIZE_X = 20.625; type::Real DOMAIN_SIZE_X;
    const type::Real _DOMAIN_SIZE_Y = 20.625; type::Real DOMAIN_SIZE_Y;
    const type::Real _DOMAIN_SIZE_Z = 20.625; type::Real DOMAIN_SIZE_Z;

    lb::adjust_simulation_size(_DOMAIN_SIZE_X,_DOMAIN_SIZE_Y,_DOMAIN_SIZE_Z,
                               procs_x, procs_y, procs_z,
                               gp.get_grid_resolution(),
                               &DOMAIN_SIZE_X, &DOMAIN_SIZE_Y, &DOMAIN_SIZE_Z);
    gp.set_grid_index_dimensions(
            DOMAIN_SIZE_X / gp.get_grid_resolution(),
            DOMAIN_SIZE_Y / gp.get_grid_resolution(),
            DOMAIN_SIZE_Z / gp.get_grid_resolution());

    lb::Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z, &gp.get_grid_resolution());

    auto cell_per_col = gp.get_cell_number_x() / procs_x;

    type::DataIndex x_proc_idx, y_proc_idx, z_proc_idx;

    d.bootstrap_partitions(worldsize);

    for(int cpu = 0; cpu < worldsize; ++cpu){
        //std::cout << cpu << std::endl;
        lb::linear_to_grid(cpu, procs_x, procs_y, x_proc_idx, y_proc_idx, z_proc_idx);
        auto my_cells = mesh::generate_lattice_single_type<int>(
                gp.get_cell_number_x(), gp.get_cell_number_y(), gp.get_cell_number_z(),
                x_proc_idx, y_proc_idx, z_proc_idx,
                cell_per_col, cell_per_col, cell_per_col,
                mesh::TCellType::REAL_CELL);
        auto part = d.get_partition(cpu);
        lb::Box3 bbox(part.vertices, gp.get_grid_resolution());
        const auto size = my_cells.size();
        for(int i = 0; i < size; ++i) {
            auto cell = my_cells[i];
            auto lid = cell.get_lid(gp.get_cell_number_x(), gp.get_cell_number_y(), gp.get_cell_number_z(), bbox);
            EXPECT_EQ(lid, i);
        }

    }

}

TEST(BoundingBox, isWellDimensionedGlobally){

    using Cell = mesh::Cell<elements::Element<3>>;
    mesh::GridParams& gp = mesh::GridParams::get_instance();

    gp.set_grid_resolution(0.0625);

    int my_rank = 0, worldsize = 1;
    int procs_x = std::cbrt(worldsize),
            procs_y = std::cbrt(worldsize),
            procs_z = std::cbrt(worldsize);

    const type::Real _DOMAIN_SIZE_X = 10.625; type::Real DOMAIN_SIZE_X;
    const type::Real _DOMAIN_SIZE_Y = 10.625; type::Real DOMAIN_SIZE_Y;
    const type::Real _DOMAIN_SIZE_Z = 10.625; type::Real DOMAIN_SIZE_Z;

    lb::adjust_simulation_size(_DOMAIN_SIZE_X,_DOMAIN_SIZE_Y,_DOMAIN_SIZE_Z,
                               procs_x, procs_y, procs_z,
                               gp.get_grid_resolution(),
                               &DOMAIN_SIZE_X, &DOMAIN_SIZE_Y, &DOMAIN_SIZE_Z);
    gp.set_grid_index_dimensions(
            DOMAIN_SIZE_X / gp.get_grid_resolution(),
            DOMAIN_SIZE_Y / gp.get_grid_resolution(),
            DOMAIN_SIZE_Z / gp.get_grid_resolution());

    lb::Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z, &gp.get_grid_resolution());

    d.bootstrap_partitions(worldsize);

    auto part = d.get_partition(0);

    lb::Box3 bbox(part.vertices, gp.get_grid_resolution());

    EXPECT_EQ(bbox.get_number_of_cells(), ((gp.get_cell_number_x())* (gp.get_cell_number_y()) * (gp.get_cell_number_z())));

    std::cout << bbox.get_number_of_cells() << " " << ((gp.get_cell_number_x() )* (gp.get_cell_number_y() ) * (gp.get_cell_number_z() )) << std::endl;
    std::cout << bbox << std::endl;
    std::cout << gp << std::endl;
}
TEST(BoundingBox, isWellDimensionedLocally){

    using Cell = mesh::Cell<elements::Element<3>>;
    mesh::GridParams& gp = mesh::GridParams::get_instance();

    gp.set_grid_resolution(0.0625);

    int my_rank = 0, worldsize = 64;

    int procs_x = std::cbrt(worldsize),
        procs_y = std::cbrt(worldsize),
        procs_z = std::cbrt(worldsize);

    const type::Real _DOMAIN_SIZE_X = 10.625; type::Real DOMAIN_SIZE_X;
    const type::Real _DOMAIN_SIZE_Y = 10.625; type::Real DOMAIN_SIZE_Y;
    const type::Real _DOMAIN_SIZE_Z = 10.625; type::Real DOMAIN_SIZE_Z;

    lb::adjust_simulation_size(_DOMAIN_SIZE_X,_DOMAIN_SIZE_Y,_DOMAIN_SIZE_Z,
                               procs_x, procs_y, procs_z,
                               gp.get_grid_resolution(),
                               &DOMAIN_SIZE_X, &DOMAIN_SIZE_Y, &DOMAIN_SIZE_Z);
    gp.set_grid_index_dimensions(
            DOMAIN_SIZE_X / gp.get_grid_resolution(),
            DOMAIN_SIZE_Y / gp.get_grid_resolution(),
            DOMAIN_SIZE_Z / gp.get_grid_resolution());

    gp.set_simulation_dimension((DOMAIN_SIZE_X), (DOMAIN_SIZE_Y), (DOMAIN_SIZE_Z));

    lb::Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z, &gp.get_grid_resolution());

    d.bootstrap_partitions(worldsize);

    auto part = d.get_partition(0);

    lb::Box3 bbox(part.vertices, gp.get_grid_resolution());

    EXPECT_EQ(bbox.get_number_of_cells(), ((gp.get_cell_number_x())* (gp.get_cell_number_y()) * (gp.get_cell_number_z())) / worldsize);

}

TEST(Cell, uniqueGID) {
    using Cell = mesh::Cell<elements::Element<3>>;

    mesh::GridParams& gp = mesh::GridParams::get_instance();
    gp.set_grid_resolution(0.5);
    auto step             = gp.get_grid_resolution();

    const int DOMAIN_SIZE_X = 10;
    const int DOMAIN_SIZE_Y = 10;
    const int DOMAIN_SIZE_Z = 10;

    gp.set_grid_index_dimensions(DOMAIN_SIZE_X / gp.get_grid_resolution(),DOMAIN_SIZE_Y / gp.get_grid_resolution(),DOMAIN_SIZE_Z / gp.get_grid_resolution());

    int my_rank = 0, worldsize = 1;
    int procs_x = std::cbrt(worldsize), procs_y = std::cbrt(worldsize), procs_z = std::cbrt(worldsize);

    auto cell_per_col     =  gp.get_maximum_index_x() / procs_x;

    type::DataIndex x_proc_idx, y_proc_idx, z_proc_idx;

    lb::linear_to_grid(my_rank, procs_x, procs_y, x_proc_idx, y_proc_idx, z_proc_idx);

    auto my_cells = mesh::generate_lattice_single_type<int>(
            gp.get_maximum_index_x(), gp.get_maximum_index_y(), gp.get_maximum_index_z(),
            x_proc_idx, y_proc_idx, z_proc_idx,
            cell_per_col, cell_per_col, cell_per_col,
            mesh::TCellType::REAL_CELL);

    lb::Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z, &gp.get_grid_resolution());

    d.bootstrap_partitions(worldsize);
    auto part = d.get_partition(my_rank);
    lb::Box3 bbox(part.vertices, gp.get_grid_resolution());

    std::vector<bool> found(my_cells.size(), false);

    for(auto cell : my_cells) {
        EXPECT_EQ(found.at(cell.gid), false);
        found.at(cell.gid) = true;
    }
}

TEST(Domain, ShouldBeTheSameBetweenGridParamsAndBoundingBox) {
    using Cell = mesh::Cell<elements::Element<3>>;
    mesh::GridParams& gp = mesh::GridParams::get_instance();
    gp.set_grid_resolution(0.000625);
    auto step             = gp.get_grid_resolution();

    int my_rank = 0, worldsize = 8;
    int procs_x = std::cbrt(worldsize), procs_y = std::cbrt(worldsize), procs_z = std::cbrt(worldsize);
    const type::Real _DOMAIN_SIZE_X = 0.100625;type::Real DOMAIN_SIZE_X;
    const type::Real _DOMAIN_SIZE_Y = 0.100625;type::Real DOMAIN_SIZE_Y;
    const type::Real _DOMAIN_SIZE_Z = 0.100625;type::Real DOMAIN_SIZE_Z;

    lb::adjust_simulation_size(_DOMAIN_SIZE_X,_DOMAIN_SIZE_Y,_DOMAIN_SIZE_Z,
                               procs_x, procs_y, procs_z,
                               gp.get_grid_resolution(),
                               &DOMAIN_SIZE_X, &DOMAIN_SIZE_Y, &DOMAIN_SIZE_Z);

    gp.set_grid_index_dimensions(
            std::floor(DOMAIN_SIZE_X / gp.get_grid_resolution()),
            std::floor(DOMAIN_SIZE_Y / gp.get_grid_resolution()),
            std::floor(DOMAIN_SIZE_Z / gp.get_grid_resolution()));

    gp.set_simulation_dimension(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z);

    auto cell_per_col  =  gp.get_maximum_index_x() / procs_x;

    type::DataIndex x_proc_idx, y_proc_idx, z_proc_idx;

    lb::linear_to_grid(my_rank, procs_x, procs_y, x_proc_idx, y_proc_idx, z_proc_idx);

    auto my_cells = mesh::generate_lattice_single_type<int>(
            gp.get_maximum_index_x(), gp.get_maximum_index_y(), gp.get_maximum_index_z(),
            x_proc_idx, y_proc_idx, z_proc_idx,
            cell_per_col, cell_per_col, cell_per_col,
            mesh::TCellType::REAL_CELL);

    lb::Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z, &gp.get_grid_resolution());

    d.bootstrap_partitions(worldsize);
    auto part = d.get_partition(my_rank);


    EXPECT_EQ(d.get_bounding_box().size_x, gp.get_cell_number_x());
    EXPECT_EQ(d.get_bounding_box().size_y, gp.get_cell_number_y());
    EXPECT_EQ(d.get_bounding_box().size_z, gp.get_cell_number_z());
}

TEST(Cell, ShouldHaveAZeroLidWithAZeroGid) {
    using Cell = mesh::Cell<elements::Element<3>>;

    mesh::GridParams& gp = mesh::GridParams::get_instance();
    gp.set_grid_resolution(0.000625);
    auto step             = gp.get_grid_resolution();
    int my_rank = 0, worldsize = 8;
    int procs_x = std::cbrt(worldsize), procs_y = std::cbrt(worldsize), procs_z = std::cbrt(worldsize);


    const type::Real _DOMAIN_SIZE_X = 0.100625;type::Real DOMAIN_SIZE_X;
    const type::Real _DOMAIN_SIZE_Y = 0.100625;type::Real DOMAIN_SIZE_Y;
    const type::Real _DOMAIN_SIZE_Z = 0.100625;type::Real DOMAIN_SIZE_Z;

    lb::adjust_simulation_size(_DOMAIN_SIZE_X,_DOMAIN_SIZE_Y,_DOMAIN_SIZE_Z,
                               procs_x, procs_y, procs_z,
                               gp.get_grid_resolution(),
                               &DOMAIN_SIZE_X, &DOMAIN_SIZE_Y, &DOMAIN_SIZE_Z);
    gp.set_grid_index_dimensions(DOMAIN_SIZE_X / gp.get_grid_resolution(),DOMAIN_SIZE_Y / gp.get_grid_resolution(),DOMAIN_SIZE_Z / gp.get_grid_resolution());


    auto cell_per_col     =  gp.get_maximum_index_x() / procs_x;

    type::DataIndex x_proc_idx, y_proc_idx, z_proc_idx;

    lb::linear_to_grid(my_rank, procs_x, procs_y, x_proc_idx, y_proc_idx, z_proc_idx);

    auto my_cells = mesh::generate_lattice_single_type<int>(
            gp.get_maximum_index_x(), gp.get_maximum_index_y(), gp.get_maximum_index_z(),
            x_proc_idx, y_proc_idx, z_proc_idx,
            cell_per_col, cell_per_col, cell_per_col,
            mesh::TCellType::REAL_CELL);

    lb::Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z, &gp.get_grid_resolution());

    d.bootstrap_partitions(worldsize);
    auto part = d.get_partition(my_rank);

    EXPECT_EQ(mesh::compute_lid(gp.get_maximum_index_x(), gp.get_maximum_index_y(), gp.get_maximum_index_z(), 0, d.get_bounding_box()), 0);
}

TEST(Cell, isInside) {
    using Cell = mesh::Cell<elements::Element<3>>;

    mesh::GridParams& gp = mesh::GridParams::get_instance();
    gp.set_grid_resolution(0.000625);
    auto step             = gp.get_grid_resolution();

    const type::Real DOMAIN_SIZE_X = 0.100625;
    const type::Real DOMAIN_SIZE_Y = 0.100625;
    const type::Real DOMAIN_SIZE_Z = 0.100625;


    type::DataIndex x_proc_idx, y_proc_idx, z_proc_idx;
    int my_rank = 0, worldsize = 8;
    /*lb::linear_to_grid(my_rank, procs_x, procs_y, x_proc_idx, y_proc_idx, z_proc_idx);
    lb::adjust_simulation_size(1.0 ,1.0 ,1.0,
                               procs_x, procs_y, procs_z,
                               grid_params.get_grid_resolution(),
                               &DOMAIN_SIZE_X, &DOMAIN_SIZE_Y, &DOMAIN_SIZE_Z);

    gp.set_grid_index_dimensions()(DOMAIN_SIZE_X / gp.get_grid_resolution(),DOMAIN_SIZE_Y / gp.get_grid_resolution(),DOMAIN_SIZE_Z / gp.get_grid_resolution());

    */

    std::array<lb::Point_3, 8> vertices = {
            lb::Point_3(0, 0, 0.524999988265336),
            lb::Point_3(0.524999988265336, 0, 0.524999988265336),
            lb::Point_3(0, 0.524999988265336, 0.524999988265336),
            lb::Point_3(0.524999988265336, 0.524999988265336, 0.524999988265336),
            lb::Point_3(0, 0, 1.04999997653067),
            lb::Point_3(0.377200309550746, 0, 0.909498348366358),
            lb::Point_3(0, 0.524999988265336, 1.04999997653067),
            lb::Point_3(0.524999988265336, 0.524999988265336, 1.04999997653067),
    };

}


} //end of namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}