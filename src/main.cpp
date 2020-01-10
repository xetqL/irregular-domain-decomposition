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
#include "generate.hpp"
int my_rank, world_size;
constexpr int DIMENSION = 3;
using Particle = elements::Element<DIMENSION>;
using Cell     = mesh::Cell<Particle>;

inline int bitselect(int condition, int truereturnvalue, int falsereturnvalue) {
    return (truereturnvalue & -condition) | (falsereturnvalue & ~(-condition)); //a when TRUE and b when FintLSE
}

std::vector<Cell> generate_lattice_single_type( long long int msx, long long int msy, long long int msz,
                                                long long int x_proc_idx, long long int y_proc_idx, long long int z_proc_idx,
                                                long long int cell_in_my_cols, long long int cell_in_my_rows, long long int cell_in_my_depth,
                                                mesh::TCellType type) {

    int cell_per_process = cell_in_my_cols * cell_in_my_rows;
    std::vector<Cell> my_cells; my_cells.reserve(cell_per_process);

    auto x_shift = (cell_in_my_rows *  x_proc_idx);
    auto y_shift = (cell_in_my_cols *  y_proc_idx);
    auto z_shift = (cell_in_my_depth * z_proc_idx);

    for(int z = 0; z < cell_in_my_depth; ++z) {
        for(int y = 0; y < cell_in_my_cols; ++y) {
            for(int x = 0; x < cell_in_my_rows; ++x) {
                int gid = (x_shift + x) + (y_shift + y) * msx + (z_shift + z) * msx * msy;
                my_cells.emplace_back(gid, type);
            }
        }
    }

    return my_cells;
}

namespace lb {

struct GridElementComputer {
    double compute_load(const std::vector<Cell>& elements){
        return std::accumulate(elements.cbegin(), elements.cend(), 0.0, [](double l, Cell e){return l + 1.0;});
    }

    std::vector<double> get_weights(const std::vector<Cell>& elements) {
        std::vector<double> weights;
        std::transform(elements.cbegin(), elements.cend(), std::back_inserter(weights), [](Cell e) {
            return 1.0;
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

void generate_lattice_rocks(const int rocks_per_stripe, int msx, int msy, int msz,
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
        long long int x,y,z;
        auto pos = lb::cell_to_global_position(msx, msy, msz, gid);
        std::tie(x,y,z) = pos;
        for(auto& rock : rocks_data){
            std::tie(cx, cy, cr) = rock;
            if(std::sqrt( std::pow(cx-x, 2) + std::pow(cy - y, 2) ) < cr) {
                cell.type  = mesh::TCellType::REAL_CELL;
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
}

std::pair<int, int> cell_to_global_position(int msx, int msy, long long position) {
    return std::make_pair(position % msx, (int) position / msx);
}

int main(int argc, char** argv) {



    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    rank = my_rank;
    //auto cellDatatype = Cell::register_datatype().element_datatype;

    SimulationParams params;
    bool err;

    std::tie(params, err) = parse_cli(argc, argv);

    if (err) {
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_UNKNOWN);
        return EXIT_FAILURE;
    }

    const int DOMAIN_SIZE_X = params.simsize_x;
    const int DOMAIN_SIZE_Y = params.simsize_y;
    const int DOMAIN_SIZE_Z = params.simsize_z;

    lb::Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z);

    long long int procs_x = params.xprocs,
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

    auto& msx = Cell::get_msx();
    auto& msy = Cell::get_msy();
    auto& msz = Cell::get_msz();

    Cell::get_cell_size() = (double) DOMAIN_SIZE_X / (double) xcells;

    if(!rank) {
        std::cout << "Cell number (X,Y,Z) = (" << msx<<","<<msy<<","<<msz<<")"<< std::endl;
        std::cout << "Cell size = " << Cell::get_cell_size() << std::endl;
        std::cout << "Domain size X = "<< (cell_in_my_rows  * xprocs * Cell::get_cell_size()) << std::endl;
    }
    //auto datatype_wrapper = Cell::register_datatype();
    d.bootstrap_partitions(world_size);

    auto part = d.get_my_partition(my_rank);
    part.init_communicators(world_size);

    std::random_device rd {};
    std::mt19937 gen{rd()};

    std::normal_distribution<double> normal_distribution(1.0 + my_rank % 2, 0.2);

    //auto nb_cells_x = (long long int) (DOMAIN_SIZE_X / (double) procs_x / d.grid_cell_size);
    //auto nb_cells_y = (long long int) (DOMAIN_SIZE_Y / (double) procs_y / d.grid_cell_size);
    //auto nb_cells_z = (long long int) (DOMAIN_SIZE_Z / (double) procs_z / d.grid_cell_size);

    long long int x_proc_idx, y_proc_idx, z_proc_idx;
    lb::linear_to_grid((long long int) my_rank, procs_x, procs_y, x_proc_idx, y_proc_idx, z_proc_idx);

    std::vector<Cell> my_cells; my_cells.reserve(cell_per_process);

    my_cells = generate_lattice_single_type(msx, msy, msz, x_proc_idx, y_proc_idx, z_proc_idx, cell_in_my_cols, cell_in_my_rows, cell_in_my_depth, mesh::TCellType::REAL_CELL);

    if(!rank){
        std::vector<Particle> particles;
        auto rejection_condition = std::make_shared<initial_condition::lennard_jones::RejectionCondition<DIMENSION>>(
                &particles, params.sig_lj, 6.25 * params.sig_lj * params.sig_lj, params.T0, 0, 0, 0,
                params.simsize_x, params.simsize_y, params.simsize_z
        );
        auto elements_generators = init_generator(rejection_condition, 1, &params);
        int cfg_idx = 0;
        while (!elements_generators.empty()) {
            std::cout << "Starting generation of particles with configuration ("
                      << cfg_idx<<"/"<<elements_generators.size()<<") ..." <<std::endl;
            auto el_gen = elements_generators.front();
            el_gen.first->generate_elements(&particles, el_gen.second, rejection_condition);
            elements_generators.pop();
            std::cout << el_gen.second <<"/"<< params.npart << " particles generated." << std::endl;
            cfg_idx++;
        }
    }

    part.move_data<Particle>(my_cells);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto stats = part.get_load_statistics<lb::GridElementComputer>(my_cells);

    if(!my_rank) print_load_statitics(stats);

    double mu = 1.0;

    auto all_loads = get_neighbors_load(stats.my_load, MPI_COMM_WORLD); //global load balancing with MPI_COMM_WORLD
    //auto avg_load  = std::accumulate(all_loads.cbegin(), all_loads.cend(), 0.0) / world_size;

    std::array<std::vector<double>, 8>  com_lb_time; std::fill(com_lb_time.begin(), com_lb_time.end(), std::vector<double>());
    std::array<double, 8>  com_degradation; std::fill(com_degradation.begin(), com_degradation.end(), 0.0);
    std::array<int, 8>  com_ncall; std::fill(com_ncall.begin(), com_ncall.end(), 0);
    std::array<int, 8>  com_pcall; std::fill(com_pcall.begin(), com_pcall.end(), 0);

    std::array<SlidingWindow<double>, 8> com_max_iteration_time_window = {
            SlidingWindow<double>(100), SlidingWindow<double>(100),
            SlidingWindow<double>(100), SlidingWindow<double>(100),
            SlidingWindow<double>(100), SlidingWindow<double>(100),
            SlidingWindow<double>(100), SlidingWindow<double>(100)
    };

    std::array<SlidingWindow<double>, 8> com_avg_iteration_time_window = {
            SlidingWindow<double>(100), SlidingWindow<double>(100),
            SlidingWindow<double>(100), SlidingWindow<double>(100),
            SlidingWindow<double>(100), SlidingWindow<double>(100),
            SlidingWindow<double>(100), SlidingWindow<double>(100)
    };

    std::vector<double> imbalance_over_time, virtual_times(MAX_ITER, 0.0), cum_vtime(MAX_ITER, 0.0);

    bool lb_decision = false;

    std::sort(my_cells.begin(), my_cells.end(), [](Cell &a, Cell &b) {return a.gid < b.gid;} );

    lb::DiffusiveOptimizer<Particle, lb::GridElementComputer> opt;
#ifdef OUTPUT
    for(int i = 0; i < world_size; ++i){
        if(my_rank == i) {
            io::scatterplot_3d_output<lb::GridPointTransformer>(i, "debug-domain-decomposition-" + std::to_string(0) + ".dat", my_cells);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    int lb_call = 0;
    for(int j = 1; j < MAX_ITER; ++j) {
        int com = lb::translate_iteration_to_vertex_group(j, part.coord);
        std::vector<double> &C = com_lb_time[com];
        int vid = part.vertices_id[com];
        int &pcall = com_pcall[com],
            &ncall = com_ncall[com];
        const Communicator& communicator = part.vertex_neighborhood[vid];
        double workload = 0,
               virtual_time = 0,
               &degradation_since_last_lb = com_degradation[com];

        /* Compute maximum workload among neighborhood */
        for(int com_id = 0; com_id < 8; com_id++) {
            auto vid = part.vertices_id[com];
            const auto& communicator = part.vertex_neighborhood[vid];
            std::vector<double> comm_workloads(communicator.comm_size, 0);
            communicator.Allgather(&workload, 1, MPI_DOUBLE, comm_workloads.data(), 1, MPI_DOUBLE, vid);
            SlidingWindow<double>& max_time_window = com_max_iteration_time_window[com_id];
            max_time_window.add(*std::max_element(comm_workloads.begin(), comm_workloads.end()));
            SlidingWindow<double>& avg_time_window = com_avg_iteration_time_window[com_id];
            avg_time_window.add(mean<double>(comm_workloads.begin(), comm_workloads.end()));
        }

        virtual_time = workload;

        SlidingWindow<double>& max_iteration_time_window = com_max_iteration_time_window[com];
        SlidingWindow<double>& avg_iteration_time_window = com_avg_iteration_time_window[com];
        auto avg_lb_cost = std::accumulate(C.begin(), C.end(), 0.0) / C.size();
        lb_decision = degradation_since_last_lb >= avg_lb_cost;
        std::cout << my_rank << " >> "<< j << "; com = "<< com << " " << (lb_decision ? "True" : "False") << std::endl;

        if(true) {

            auto start_time = MPI_Wtime();
            opt.optimize_neighborhood(com, part, my_cells, mu);
            double lb_time  = MPI_Wtime() - start_time;

            {
                std::vector<double> lbtimes(communicator.comm_size);
                communicator.Allgather(&lb_time, 1, MPI_DOUBLE, lbtimes.data(), 1, MPI_DOUBLE, 122334);
                lb_time = *std::max_element(lbtimes.cbegin(), lbtimes.cend());
            }

            C.push_back(lb_time);
            auto avg_lb_cost = std::accumulate(C.begin(), C.end(), 0.0) / C.size();
            std::cout << avg_lb_cost << std::endl;
            double slope = max_iteration_time_window.data_container.back() - max_iteration_time_window.data_container.front(); //get_slope<double>(iteration_time_window.data_container)
            int tau = (int) std::sqrt(2.0 * avg_lb_cost / slope);
            if(tau < 0)  tau = 40;
            if(tau == 0) tau = 1;
            {
                std::vector<int> com_tau(communicator.comm_size, 0);
                communicator.Allgather(&tau, 1, MPI_INT, com_tau.data(), 1, MPI_INT, vid);
                ncall = *std::max_element(com_tau.cbegin(), com_tau.cend());
            }
            max_iteration_time_window.data_container.clear();
            avg_iteration_time_window.data_container.clear();
            degradation_since_last_lb = 0.0;
            pcall = j+8;
            lb_call++;
            virtual_time += lb_time;
        }

        degradation_since_last_lb +=
                median<double>(max_iteration_time_window.end()-std::min(3, (int) max_iteration_time_window.size()), max_iteration_time_window.end()) -
                median<double>(avg_iteration_time_window.end()-std::min(3, (int) max_iteration_time_window.size()), avg_iteration_time_window.end());

        virtual_times[j] = virtual_time;

        if(!my_rank) {
            /* increase workload */
            cum_vtime[j] = j == 0 ? virtual_time : cum_vtime[j-1] + virtual_time;
        }
#ifdef OUTPUTS
        for(int i = 0; i < world_size; ++i){
            if(my_rank == i)
                io::scatterplot_3d_output<lb::GridPointTransformer>(i, "debug-domain-decomposition-" + std::to_string(j) + ".dat", my_cells);
            MPI_Barrier(MPI_COMM_WORLD);
        }
#endif
    }

    std::vector<int> lb_calls(world_size);
    double total_time = std::accumulate(virtual_times.cbegin(), virtual_times.cend(), 0.0);
    double max_total_time;
    MPI_Reduce(&total_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    gather_elements_on(&lb_call, 1, MPI_INT, lb_calls.data(), 0, MPI_COMM_WORLD);

    if(!my_rank) {
        std::ofstream fimbalance;
        fimbalance.open("imbalance.txt", std::ofstream::out);
        for(auto imbalance : imbalance_over_time){
            fimbalance << imbalance << std::endl;
        }
        fimbalance.close();
        std::ofstream flb_calls;
        flb_calls.open("lb_calls.txt", std::ofstream::out);
        for(auto lb_call : lb_calls){
            flb_calls << lb_call << std::endl;
        }
        flb_calls.close();
        std::ofstream fvtimes;
        fvtimes.open("vtimes.txt", std::ofstream::out);
        for(auto t : virtual_times){
            fvtimes << t << std::endl;
        }
        fvtimes.close();
        std::ofstream fcumvtimes;
        fcumvtimes.open("vcumtimes.txt", std::ofstream::out);
        for(auto t : cum_vtime){
            fcumvtimes << t << std::endl;
        }
        fcumvtimes.close();
        std::cout << "------------------- " << std::endl;
        std::cout << "Total time: "         << max_total_time << std::endl;
        std::cout << "------------------- " << std::endl;
    }

    lb::print_load_statitics(part.get_load_statistics<lb::GridElementComputer>(my_cells));

    MPI_Finalize();

    return 0;
}




