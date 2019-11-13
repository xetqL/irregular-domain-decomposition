#include <iostream>
#include <GeometricUtils.hpp>
#include <mpi.h>
#include <random>
#include <iostream>
#include <chrono>
#include <thread>
#include <set>
#include "Communicator.hpp"
#include "DiffusiveOptimizer.hpp"
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
        double x, y, z;
        std::tie(x, y, z) = element.get_center();
        return Point_3(x, y, z);
    }
};

}

void print_load_statitics(lb::Partition::LoadStatistics stats){
    lb::Partition::ElementCount count;
    lb::Partition::Load max_load, avg_load;
    lb::Partition::Imbalance load_imbalance, my_imbalance;
    std::tie(count, max_load, avg_load, load_imbalance, my_imbalance) = stats;
    std::cout << "===============================================\n"
              << "Total number of elements: " << count             << "\n"
              << "            Maximum load: " << max_load          << "\n"
              << "            Average load: " << avg_load          << "\n"
              << "   Global Load Imbalance: " << load_imbalance    << "\n"
              << "       My Load Imbalance: " << my_imbalance      <<
               "\n===============================================" << std::endl;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int procs_x = (int) std::cbrt(world_size),
        procs_y = procs_x,
        procs_z = procs_x;

    const int DOMAIN_SIZE_X = 10;
    const int DOMAIN_SIZE_Y = 10;
    const int DOMAIN_SIZE_Z = 10;

    lb::Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z);

    d.grid_cell_size = std::atof(argv[1]);
    int MAX_TRIAL    = std::atoi(argv[2]);


    Cell::get_msx() = DOMAIN_SIZE_X / d.grid_cell_size;
    Cell::get_msy() = DOMAIN_SIZE_Y / d.grid_cell_size;
    Cell::get_msz() = DOMAIN_SIZE_Z / d.grid_cell_size;
    Cell::get_cell_size() = d.grid_cell_size;
    auto datatype_wrapper = Cell::register_datatype();
    d.bootstrap_partitions(world_size);

    auto part = d.get_my_partition(my_rank);
    part.init_communicators(world_size);

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> normal_distribution(1.0 + my_rank % 2, 0.2);

    int nb_cells_x = (DOMAIN_SIZE_X / (double) procs_x / (double) d.grid_cell_size);
    int nb_cells_y = (DOMAIN_SIZE_Y / (double) procs_y / (double) d.grid_cell_size);
    int nb_cells_z = (DOMAIN_SIZE_Z / (double) procs_z / (double) d.grid_cell_size);

    auto cell_per_process = nb_cells_x * nb_cells_y * nb_cells_z;

    std::vector<Cell> my_cells; my_cells.reserve(cell_per_process);

    int x_proc_idx, y_proc_idx, z_proc_idx; lb::linear_to_grid(my_rank, procs_x, procs_y, x_proc_idx, y_proc_idx, z_proc_idx);

    const int total_cells_x = nb_cells_x * procs_x, total_cells_y = nb_cells_y * procs_y, total_cells_z = nb_cells_z * procs_z;

    auto x_shift = (nb_cells_x * x_proc_idx);
    auto y_shift = (nb_cells_y * y_proc_idx);
    auto z_shift = (nb_cells_z * z_proc_idx);
    lb::GridPointTransformer gpt;

    for(int z = 0; z < nb_cells_z; ++z) {
        for(int y = 0; y < nb_cells_y; ++y) {
            for(int x = 0; x < nb_cells_x; ++x) {
                int gid = (x_shift + x) + (y_shift + y) * total_cells_x + (z_shift + z) * total_cells_x * total_cells_y;
                my_cells.emplace_back(gid, 0, my_rank == 0 ? 10 : 1, 0.0);
            }

        }
    }

    //return 0;
    //std::cout << "BEFORE: " << part.get_load_imbalance<GridElementComputer>(my_cells) << std::endl;
    auto stats = part.get_load_statistics<lb::GridElementComputer>(my_cells);

    if(!my_rank) print_load_statitics(stats);
    /*for(int i = 0; i < world_size; ++i){
        if(my_rank == i) io::scatterplot_3d_output<lb::GridPointTransformer>(i, "debug-domain-decomposition-"+std::to_string(0)+".dat", my_cells);
        MPI_Barrier(MPI_COMM_WORLD);
    }*/
    auto prev_imbalance = stats.global;
    double mu = 1.0;

    auto all_loads = get_neighbors_load(stats.my_load, MPI_COMM_WORLD); //global load balancing with MPI_COMM_WORLD
    auto avg_load  = std::accumulate(all_loads.cbegin(), all_loads.cend(), 0.0) / world_size;

    lb::DiffusiveOptimizer<lb::GridPointTransformer, lb::GridElementComputer, Cell> diff_opt;

    //diff_opt.optimize(part, my_cells, datatype_wrapper.element_datatype, 1.0);

    stats = part.get_load_statistics<lb::GridElementComputer>(my_cells);
    if(!my_rank)
        print_load_statitics(stats);
    std::vector<double> all_mu(8, mu);
    LinearHashMap<int, double, 8> vertices_mu;
    std::zip(part.vertices_id.begin(), part.vertices_id.end(), all_mu.begin(), all_mu.end(), vertices_mu.begin());
    LinearHashMap <int, int,  8> vertices_remaining_trials;
    std::transform(part.vertices_id.begin(), part.vertices_id.end(), vertices_remaining_trials.begin(),
                   [](auto id){return std::make_pair(id, 10);});
    for(int j = 0; j < MAX_TRIAL; ++j){
        auto data = part.move_vertices<lb::GridPointTransformer, lb::GridElementComputer, Cell>(my_cells, datatype_wrapper.element_datatype, avg_load, 1.0, vertices_remaining_trials);
        if(prev_imbalance <= stats.global) mu *= 0.9;
        prev_imbalance = stats.global;

        stats = part.get_load_statistics<lb::GridElementComputer>(my_cells);
        if(!my_rank)
            print_load_statitics(stats);
    }

    MPI_Finalize();

    return 0;
}




