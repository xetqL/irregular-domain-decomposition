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
#include "Types.hpp"

int my_rank, world_size;
constexpr int DIMENSION = 3;
using Particle = elements::Element<DIMENSION>;
using Cell     = mesh::Cell<Particle>;

inline int bitselect(int condition, int truereturnvalue, int falsereturnvalue) {
    return (truereturnvalue & -condition) | (falsereturnvalue & ~(-condition)); //a when TRUE and b when FintLSE
}

namespace lb {

struct GridElementComputer {
    double compute_load(const std::vector<Cell>& elements) {
        if(elements.empty()) return 1.0;
        return std::accumulate(elements.cbegin(), elements.cend(), 1.0, [](double l, const Cell& e){return l + e.number_of_elements();});
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

int main(int argc, char** argv) {
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    rank = my_rank;

    SimulationParams params;
    bool err;

    std::tie(params, err) = parse_cli(argc, argv);

    if (err) {
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_UNKNOWN);
        return EXIT_FAILURE;
    }

    const auto procs_x = params.xprocs,
               procs_y = params.yprocs,
               procs_z = params.zprocs;

    mesh::GridParams& grid_params = mesh::GridParams::get_instance();

    grid_params.set_grid_resolution(2.5 * (params.sig_lj));

    type::Real DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z;

    lb::adjust_simulation_size(params.simsize_x,params.simsize_y,params.simsize_z,
                               procs_x, procs_y, procs_z,
                               grid_params.get_grid_resolution(),
                               &DOMAIN_SIZE_X, &DOMAIN_SIZE_Y, &DOMAIN_SIZE_Z);

    grid_params.set_grid_index_dimensions((type::DataIndex) (DOMAIN_SIZE_X / grid_params.get_grid_resolution()),
                                          (type::DataIndex) (DOMAIN_SIZE_Y / grid_params.get_grid_resolution()),
                                          (type::DataIndex) (DOMAIN_SIZE_Z / grid_params.get_grid_resolution()));

    grid_params.set_simulation_dimension((DOMAIN_SIZE_X), (DOMAIN_SIZE_Y), (DOMAIN_SIZE_Z));

    const auto cell_in_my_rows  = grid_params.get_cell_number_x() / procs_x,
               cell_in_my_cols  = grid_params.get_cell_number_y() / procs_y,
               cell_in_my_depth = grid_params.get_cell_number_z() / procs_z;
    const auto MAX_ITER = (type::DataIndex) params.MAX_STEP;
    const auto msx = grid_params.get_cell_number_x();
    const auto msy = grid_params.get_cell_number_y();
    const auto msz = grid_params.get_cell_number_z();
    const auto total_cell_count =  grid_params.get_cell_number_x() *
                                   grid_params.get_cell_number_y() *
                                   grid_params.get_cell_number_z();

    lb::Domain d(DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z, &grid_params.get_grid_resolution());

    if(!rank) {
        std::cout << "Cell number (X,Y,Z) = (" << msx<<","<<msy<<","<<msz<<")"<< std::endl;
        std::cout << "Total cell count: " << total_cell_count << std::endl;
        std::cout << "Cell size = " << grid_params.get_grid_resolution() << std::endl;
        std::cout << "Domain size X = "<< (cell_in_my_rows  * procs_x * grid_params.get_grid_resolution()) << std::endl;
    }

    d.bootstrap_partitions(world_size);
    auto part = d.get_partition(my_rank);
    part.init_communicators(world_size);
    lb::Box3 bbox(part.vertices, grid_params.get_grid_resolution());
    std::cout << bbox << std::endl;
    type::DataIndex x_proc_idx, y_proc_idx, z_proc_idx;
    lb::linear_to_grid(my_rank, procs_x, procs_y, x_proc_idx, y_proc_idx, z_proc_idx);
    std::vector<Cell> my_cells(bbox.get_number_of_cells(), mesh::EMPTY_CELL);

    std::cout << "Populate lattice " << std::endl;
    mesh::populate_lattice_single_type<Particle>(&my_cells, bbox,
            msx, msy, msz, x_proc_idx, y_proc_idx, z_proc_idx,
            cell_in_my_cols, cell_in_my_rows, cell_in_my_depth,
            mesh::TCellType::REAL_CELL);

    std::vector<Particle> particles;
    if(!rank) {
        auto rejection_condition = std::make_shared<initial_condition::lennard_jones::RejectionCondition<DIMENSION>>(
                &particles, params.sig_lj, 0, params.T0, bbox.xmin, bbox.ymin, bbox.zmin,
                bbox.xmax, bbox.ymax, bbox.zmax
        );
        auto elements_generators = init_generator(rejection_condition, 1, &params, 100000);
        int cfg_idx = 0;
        while (!elements_generators.empty()) {
            std::cout << "Starting generation of particles with configuration ("
                      << cfg_idx<<"/"<<elements_generators.size()<<") ..." <<std::endl;
            auto el_gen = elements_generators.front();
            el_gen.first->generate_elements(&particles, 50, rejection_condition);
            elements_generators.pop();
            std::cout << el_gen.second <<"/"<< params.npart << " particles generated." << std::endl;
            cfg_idx++;
        }
    }

    mesh::insert_or_remove(&my_cells, &particles, bbox);

    /*for(int i = 0; i < my_cells.size(); ++i) {
        auto cell = my_cells[i];
        cell.lid = mesh::compute_lid(
                        grid_params.get_cell_number_x(),
                        grid_params.get_cell_number_y(),
                        grid_params.get_cell_number_z(),
                        cell.gid,
                        bbox);
    }*/

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto stats = part.get_load_statistics<lb::GridElementComputer>(my_cells);
    if(!my_rank) lb::print_load_statitics(stats);

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
            int vid = part.vertices_id[com];
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

        if(true) {

            double start_time = MPI_Wtime();
            auto ghost = opt.optimize_neighborhood(com, part, my_cells, mu);
            double lb_time    = MPI_Wtime() - start_time;

            /* send ghost */
            for(const auto& dest_and_data : ghost) {
                lb::ProcRank               dest = dest_and_data.first;
                std::vector<lb::DataIndex> data = dest_and_data.second;
                std::for_each(data.begin(), data.end(), [&my_cells](auto& lid) { lid = my_cells[lid].gid; });

            }

            /* receive ghost*/

            {
                std::vector<double> lbtimes(communicator.comm_size);
                communicator.Allgather(&lb_time, 1, MPI_DOUBLE, lbtimes.data(), 1, MPI_DOUBLE, 122334);
                lb_time = *std::max_element(lbtimes.cbegin(), lbtimes.cend());
            }

            C.push_back(lb_time);
            auto avg_lb_cost = std::accumulate(C.begin(), C.end(), 0.0) / C.size();
            double slope = max_iteration_time_window.data_container.back() - max_iteration_time_window.data_container.front();
            int tau = (int) std::sqrt(2.0 * avg_lb_cost / slope);

            if(tau < 0)  tau = 40; //
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

        /* increase workload */
        if(!my_rank) {
            cum_vtime[j] = j == 0 ? virtual_time : cum_vtime[j-1] + virtual_time;
        }

        std::vector<Particle > rp(params.npart), sp;
        for(auto& cell : my_cells)
            rp.insert(rp.begin(), cell.begin(), cell.end());

        CommunicationDatatype datatype = Particle::register_datatype();

        gather_elements_on(sp.data(), sp.size(), datatype.element_datatype, rp.data(), 0, MPI_COMM_WORLD);

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
    auto s = part.get_load_statistics<lb::GridElementComputer>(my_cells);

    if(!rank)
        lb::print_load_statitics(s);
    MPI_Finalize();

    return 0;
}




