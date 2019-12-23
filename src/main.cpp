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
#include <SlidingWindow.hpp>
#include <Utils.hpp>
#include <SimulationParams.hpp>
#include <CLIParser.hpp>

int my_rank, world_size;

inline int bitselect(int condition, int truereturnvalue, int falsereturnvalue) {
    return (truereturnvalue & -condition) | (falsereturnvalue & ~(-condition)); //a when TRUE and b when FintLSE
}

std::vector<Cell> generate_lattice_single_type( int msx, int msy, int msz,
                                                int x_proc_idx, int y_proc_idx, int z_proc_idx,
                                                int cell_in_my_cols, int cell_in_my_rows, int cell_in_my_depth,
                                                int type, float weight, float erosion_probability) {

    /*for(int j = 0; j < cell_in_my_cols; ++j) {
        for(int i = 0; i < cell_in_my_rows; ++i) {
            for(int k = 0; k < cell_in_my_depth; ++k) {
                int gid = cell_in_my_rows * x_proc_idx + i + msx * (j + (y_proc_idx * cell_in_my_cols));
                my_cells.emplace_back(gid, type, weight, erosion_probability);
            }
        }
    }
    */
    int cell_per_process = cell_in_my_cols * cell_in_my_rows;
    std::vector<Cell> my_cells; my_cells.reserve(cell_per_process);

    auto x_shift = (cell_in_my_rows *  x_proc_idx);
    auto y_shift = (cell_in_my_cols *  y_proc_idx);
    auto z_shift = (cell_in_my_depth * z_proc_idx);

    for(int z = 0; z < cell_in_my_depth; ++z) {
        for(int y = 0; y < cell_in_my_cols; ++y) {
            for(int x = 0; x < cell_in_my_rows; ++x) {
                int gid = (x_shift + x) + (y_shift + y) * msx + (z_shift + z) * msx * msy;
                my_cells.emplace_back(gid, 0, 1, 0.0);
            }
        }
    }

    return my_cells;
}

namespace lb {

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

void generate_lattice_rocks(const int rocks_per_stripe, int msx, int msy,
                            std::vector<Cell>* _cells,
                            float erosion_probability,
                            int begin_stripe, int end_stripe){
    using Radius   = int;
    using XPosition= int;
    using YPosition= int;

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<Cell>& cells = *_cells;
    std::vector<std::tuple<XPosition, YPosition, Radius>> rocks_data(rocks_per_stripe);

    for (int i = 0; i < rocks_per_stripe; ++i) {
//        rocks_data[i] = std::make_tuple((int) std::floor((i+1) * msx / (rocks_per_stripe+1)), (begin_stripe + end_stripe) * (3.0/4.0) : (begin_stripe + end_stripe) / 2, (end_stripe - begin_stripe) / 4);
        rocks_data[i] = std::make_tuple((int) std::floor((i+1) * msx / (rocks_per_stripe+1)),
                //rank == size-1 ? (begin_stripe + end_stripe) / 2 + (end_stripe - begin_stripe) / 4 : (begin_stripe + end_stripe) / 2, (end_stripe - begin_stripe) / 3);
                                        (begin_stripe + end_stripe) / 2, (end_stripe - begin_stripe) / 3);
    }

    for(auto& cell : cells) {
        int cx, cy, cr;
        int gid = cell.gid;
        auto pos = lb::cell_to_global_position(msx, msy, gid);
        for(auto& rock : rocks_data){
            std::tie(cx, cy, cr) = rock;
            if(std::sqrt( std::pow(cx-pos.first, 2) + std::pow(cy - pos.second, 2) ) < cr) {
                cell.type  = Cell::ROCK_TYPE;
                cell.weight= 0.0;
                cell.erosion_probability = erosion_probability;
                break;
            }
        }
    }
}

std::tuple<int, int, int> linear_to_grid(const long long index, const long long c, const long long r){
    int x_idx = (int) (index % (c*r) % c);           // col
    int y_idx = (int) std::floor(index % (c*r) / c); // row
    int z_idx = (int) std::floor(index / (c*r));     // depth
    assert(c==r);
    assert(x_idx < c);
    assert(y_idx < r);
    assert(z_idx < r);
    return std::make_tuple(x_idx, y_idx, z_idx);
};

std::pair<int, int> cell_to_global_position(int msx, int msy, long long position) {
    return std::make_pair(position % msx, (int) position / msx);
}

int main(int argc, char** argv) {
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    rank = my_rank;
    auto cellDatatype = Cell::register_datatype().element_datatype;

    SimulationParams params;
    bool err;

    std::tie(params, err) = parse_cli(argc, argv);

    if (err) {
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_UNKNOWN);
        return EXIT_FAILURE;
    }



    const int DOMAIN_SIZE_X = 10;
    const int DOMAIN_SIZE_Y = 10;
    const int DOMAIN_SIZE_Z = 10;

    lb::Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z);

    //d.grid_cell_size = std::atof(argv[1]);
    //int MAX_ITER     = std::atoi(argv[2]);
    int procs_x = params.xprocs,
            procs_y = params.yprocs,
            procs_z = params.zprocs;
    const unsigned int xprocs = params.xprocs,
                       yprocs = params.yprocs,
                       zprocs = params.zprocs,
                       cell_per_process = params.cell_per_process,
                       MAX_ITER = params.MAX_STEP;

    assert(procs_x * procs_y * procs_z == world_size);
    const int cell_in_my_rows  = (int) std::cbrt(cell_per_process),
              cell_in_my_cols  = cell_in_my_rows,
              cell_in_my_depth = cell_in_my_rows;

    const int xcells = cell_in_my_rows  * xprocs,
              ycells = cell_in_my_cols  * yprocs,
              zcells = cell_in_my_depth * zprocs;

    Cell::get_msx() = xcells;
    Cell::get_msy() = ycells;
    Cell::get_msz() = zcells;

    int& msx = Cell::get_msx();
    int& msy = Cell::get_msy();
    int& msz = Cell::get_msz();

    Cell::get_cell_size() = DOMAIN_SIZE_X / xcells;

    auto datatype_wrapper = Cell::register_datatype();
    d.bootstrap_partitions(world_size);

    auto part = d.get_my_partition(my_rank);
    part.init_communicators(world_size);

    std::random_device rd {};
    std::mt19937 gen{rd()};

    std::normal_distribution<double> normal_distribution(1.0 + my_rank % 2, 0.2);

    int nb_cells_x = (DOMAIN_SIZE_X / (double) procs_x / (double) d.grid_cell_size);
    int nb_cells_y = (DOMAIN_SIZE_Y / (double) procs_y / (double) d.grid_cell_size);
    int nb_cells_z = (DOMAIN_SIZE_Z / (double) procs_z / (double) d.grid_cell_size);

    int x_proc_idx, y_proc_idx, z_proc_idx;
    lb::linear_to_grid(my_rank, procs_x, procs_y, x_proc_idx, y_proc_idx, z_proc_idx);

    std::vector<Cell> my_cells; my_cells.reserve(cell_per_process);

    my_cells = generate_lattice_single_type(msx, msy, msz, x_proc_idx, y_proc_idx, z_proc_idx, cell_in_my_cols, cell_in_my_rows, cell_in_my_depth, Cell::WATER_TYPE, 1.0, 0.0);

    auto stats = part.get_load_statistics<lb::GridElementComputer>(my_cells);

    if(!my_rank) print_load_statitics(stats);

    auto prev_imbalance = stats.global;
    double mu = 1.0;

    auto all_loads = get_neighbors_load(stats.my_load, MPI_COMM_WORLD); //global load balancing with MPI_COMM_WORLD
    auto avg_load  = std::accumulate(all_loads.cbegin(), all_loads.cend(), 0.0) / world_size;

    stats = part.get_load_statistics<lb::GridElementComputer>(my_cells);
    if(!my_rank)
        print_load_statitics(stats);
    LinearHashMap<int, std::vector<double>, 8> com_lb_time;
    std::transform(part.vertices_id.begin(), part.vertices_id.end(), com_lb_time.begin(),
                   [](auto id){return std::make_pair(id, std::vector<double>());});
    std::array<SlidingWindow<double>, 8> arr_window = {
            SlidingWindow<double>(100), SlidingWindow<double>(100),
            SlidingWindow<double>(100), SlidingWindow<double>(100),
            SlidingWindow<double>(100), SlidingWindow<double>(100),
            SlidingWindow<double>(100), SlidingWindow<double>(100)
    };
    LinearHashMap<int, int, 8> com_ncall;
    std::transform(part.vertices_id.begin(), part.vertices_id.end(), com_ncall.begin(),
                   [](auto id){return std::make_pair(id, 0);});
    bool lb_decision = false;
    int pcall = 0;
    for(int j = 0; j < MAX_ITER; ++j) {

        int com = lb::translate_iteration_to_vertex_group(j, part.coord);
        lb_decision = pcall + (*search_in_linear_hashmap<int, int, 8>(com_ncall, com)).second >= j;
        auto vid = part.vertices_id[com];
        const auto& communicator = part.vertex_neighborhood[vid];
        std::vector<int> comm_workloads(communicator.comm_size, 0);
        int workload = my_cells.size();
        communicator.Allgather(&workload, 1, MPI_INT, comm_workloads.data(), 1, MPI_INT, vid);
        SlidingWindow<double>& window = arr_window[com];

        window.add(*std::max_element(comm_workloads.begin(), comm_workloads.end()));

        if(lb_decision) {
            auto start_time = MPI_Wtime();
            auto data = part.move_selected_vertices<lb::GridPointTransformer, lb::GridElementComputer, Cell>(
                    j,                                 // the groups that moves
                    my_cells,                          // the data to LB
                    datatype_wrapper.element_datatype, // the data associated with the data
                    avg_load,                          // the average load
                    mu);                               // the movement factor
            auto lb_time = MPI_Wtime() - start_time;
            //(*search_in_linear_hashmap<int, std::vector<double>, 8>(com_lb_time, com)).second.push_back(lb_time);
            auto C = (*search_in_linear_hashmap<int, std::vector<double>, 8>(com_lb_time, com)).second;
            C.push_back(lb_time);
            std::for_each(C.cbegin(), C.cend(), [&](auto v){std::cout << "->"<< v << std::endl;});
            //std::for_each(C2.cbegin(), C2.cend(), [&](auto v){std::cout << v << std::endl;});
            auto avg_lb_cost = std::accumulate(C.begin(), C.end(), 0) / C.size();
            std::cout << avg_lb_cost << std::endl;
            auto tau = (int) std::sqrt(2*avg_lb_cost/get_slope<double>(window.data_container));
            std::cout << tau << std::endl;


            (*search_in_linear_hashmap<int, int, 8>(com_ncall, com)).second = tau;

            std::cout << "slope: " << (int) get_slope<double>(window.data_container) << std::endl;
            window.data_container.clear();
        }
        stats = part.get_load_statistics<lb::GridElementComputer>(my_cells);
        if(!my_rank) {
            std::cout << "Iteration " << j << std::endl;
            print_load_statitics(stats);
        }
    }

    MPI_Finalize();

    return 0;
}




