#include <iostream>
#include <GeometricUtils.hpp>
#include <mpi.h>
#include <random>
#include <iostream>
#include <chrono>
#include <thread>
#include <set>
#include "Communicator.hpp"
#include <Domain.hpp>
#include <Cell.hpp>

int my_rank, world_size;

inline int bitselect(int condition, int truereturnvalue, int falsereturnvalue) {
    return (truereturnvalue & -condition) | (falsereturnvalue & ~(-condition)); //a when TRUE and b when FintLSE
}

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

namespace lb{

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
        std::tie(x,y,z) = element.get_center();
        return Point_3(x, y, z);
    }
};

}

int main() {
    MPI_Init(nullptr, nullptr);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int procs_x = (int) std::cbrt(world_size),
        procs_y = procs_x,
        procs_z = procs_x;

    const int DOMAIN_SIZE_X = 10;
    const int DOMAIN_SIZE_Y = 10;
    const int DOMAIN_SIZE_Z = 10;

    lb::Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z);

    d.grid_cell_size = 1;

    Cell::get_msx() = DOMAIN_SIZE_X / d.grid_cell_size;
    Cell::get_msy() = DOMAIN_SIZE_Y / d.grid_cell_size;
    Cell::get_msz() = DOMAIN_SIZE_Z / d.grid_cell_size;
    Cell::get_cell_size() = d.grid_cell_size;

    d.bootstrap_partitions(world_size);

    auto part = d.get_my_partition(my_rank);
    part.init_communicators(world_size);

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> normal_distribution(1.0 + my_rank % 2, 0.2);

    auto bbox = CGAL::bbox_3(part.vertices.begin(), part.vertices.end());
    int nb_cells_x = (int) (bbox.xmax() - bbox.xmin()) / d.grid_cell_size;
    int nb_cells_y = (int) (bbox.ymax() - bbox.ymin()) / d.grid_cell_size;
    int nb_cells_z = (int) (bbox.zmax() - bbox.zmin()) / d.grid_cell_size;

    auto cell_per_process = nb_cells_x * nb_cells_y * nb_cells_z;

    std::vector<Cell> my_cells; my_cells.reserve(cell_per_process);

    int x_proc_idx, y_proc_idx, z_proc_idx; lb::linear_to_grid(my_rank, procs_x, procs_y, x_proc_idx, y_proc_idx, z_proc_idx);

    const int xcells = nb_cells_x * procs_x, ycells = nb_cells_y * procs_y, zcells = nb_cells_z * procs_z;

    for(int j = 0; j < nb_cells_x; ++j) {
        for(int i = 0; i < nb_cells_y; ++i) {
            for(int k = 0; k < nb_cells_z; ++k) {
                int gid = procs_x * x_proc_idx + i + xcells * (j + (y_proc_idx * procs_y)) +
                        ((z_proc_idx * procs_z) + k) * xcells * ycells;
                my_cells.emplace_back(gid, 0, normal_distribution(gen), 0.0);
            }
        }
    }
    //std::cout << "BEFORE: " << part.get_load_imbalance<GridElementComputer>(my_cells) << std::endl;

    part.move_vertices<lb::GridPointTransformer, lb::GridElementComputer, Cell>(my_cells);

    //std::cout << "AFTER: " << part.get_load_imbalance<GridElementComputer>(my_cells) << std::endl;

    std::cout << my_rank << " " << part << std::endl;

    MPI_Finalize();

    return 0;
}




