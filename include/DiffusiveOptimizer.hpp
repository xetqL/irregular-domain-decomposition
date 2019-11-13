//
// Created by xetql on 11/4/19.
//

#ifndef ADLBIRREG_DIFFUSIVEOPTIMIZER_HPP
#define ADLBIRREG_DIFFUSIVEOPTIMIZER_HPP

#include "Partition.hpp"
#include "Utils.hpp"

namespace lb
{
inline void print_load_statitics(Partition::LoadStatistics stats){
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
template<class GridPointTransformer, class GridElementComputer, class Cell>
class DiffusiveOptimizer {

public:
    void optimize(Partition& part, std::vector<Cell>& my_cells, MPI_Datatype datatype, double mu, int overloading = 0) {
        get_MPI_rank(my_rank);
        get_MPI_worldsize(worldsize);
/*
        GridElementComputer lc;
        std::vector<double> loads(worldsize), all_mu(8, mu);
        double my_load = lc.compute_load(my_cells);
        double max_load;

        MPI_Allreduce(&my_load, &max_load, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        double avg_load = max_load / (double) worldsize;
        auto neighbors_load  = get_neighbors_info(my_load,  MPI_DOUBLE, part.neighbor_list, part.vertex_neighborhood);
        std::vector<double> imbalances(8);
        std::transform(neighbors_load.begin(), neighbors_load.end(), imbalances.begin(),
                       [avg_load](auto vals){ return 1000.0; });//(*std::max_element(vals.second.cbegin(), vals.second.cend()) / avg_load) - 1.0; });
        LinearHashMap<int, double, 8> vertices_mu;
        LinearHashMap<int, double, 8> previous_imbalance;
        std::zip(part.vertices_id.begin(), part.vertices_id.end(), imbalances.begin(), imbalances.end(), previous_imbalance.begin());
        auto remaining_it = (int) std::cbrt(worldsize)-2;

        LinearHashMap <int, bool, 8> vertices_status; //active = true, inactive = false
        std::transform(part.vertices_id.begin(), part.vertices_id.end(), vertices_status.begin(),
                       [](auto id){return std::make_pair(id, true);});

        auto neighbors_status = get_neighbors_info(overloading, MPI_INT, part.neighbor_list, part.vertex_neighborhood);
        LinearHashMap<int,  bool, 8> strategy;
        std::transform(neighbors_status.cbegin(), neighbors_status.cend(), strategy.begin(),
           [](std::pair<int, std::vector<int>> statuses) {
               return std::make_pair(statuses.first, std::any_of(statuses.second.cbegin(), statuses.second.cend(), [](int status){return status;}));
           });

        LinearHashMap <int, int, 8> vertices_remaining_trials;
        std::transform(part.vertices_id.begin(), part.vertices_id.end(), vertices_remaining_trials.begin(),
                       [](auto id){return std::make_pair(id, 10);});
        while((remaining_it) > 0){ //|| std::accumulate(vertices_remaining_trials.begin(), vertices_remaining_trials.end(), 0, [](int sum, auto rm) {return sum + rm.second;}) > 0) {
            //auto n_list = filter_active_neighbors(part.vertices_id, vertices_remaining_trials, part.vertex_neighborhood);
#ifdef DEBUG
            std::cout << my_rank << " has "<< remaining_it <<" remaining it and "
            << std::accumulate(vertices_remaining_trials.begin(), vertices_remaining_trials.end(), 0, [](int sum, auto rm) {return sum + rm.second;}) << " active vertices" <<  std::endl;
#endif
            part.move_vertices<GridPointTransformer, GridElementComputer, Cell>(my_cells, datatype, avg_load, 1.0, vertices_remaining_trials);
            my_load = lc.compute_load(my_cells);

            neighbors_load  = get_neighbors_info(my_load,  MPI_DOUBLE, n_list, part.vertex_neighborhood);

            for(int vid : part.vertices_id) {
                //const bool vertex_status   = (*search_in_linear_hashmap<int, bool,   8>(strategy, vid)).second; //ulba or not
                int& vertex_rem_trial = (*search_in_linear_hashmap<int, int,8>(vertices_remaining_trials, vid)).second; // stop or not
                double& prev_imbl = (*search_in_linear_hashmap<int, double, 8>(previous_imbalance, vid)).second;

                if(vertex_rem_trial > 0) { // if continue
                    //compute imbalance
                    auto current_neighborhood_load = (*search_in_linear_hashmap<int, std::vector<double>, 8>(neighbors_load, vid)).second;
                    double imbalance  = (*std::max_element(current_neighborhood_load.cbegin(), current_neighborhood_load.cend()) / avg_load) - 1.0;
                    //new imbalance is good
                    if(imbalance < 0.2) {
                        vertex_rem_trial--;
                    } else {
                        vertex_rem_trial = 10;
                        prev_imbl = imbalance;
                    }
                }

#ifdef DEBUG
                std::cout << my_rank << " with vid "<< vid << " has "<< prev_imbl << " rt: " << vertex_rem_trial << std::endl;
#endif
            }
            remaining_it--;
            //MPI_Barrier(MPI_COMM_WORLD);
        }*/

#ifdef DEBUG
        std::cout << my_rank << " has left"<<std::endl;
#endif

        auto stats = part.get_load_statistics<GridElementComputer>(my_cells);

        if(!my_rank) print_load_statitics(stats);
        /*for(int i = 0; i < world_size; ++i){
            if(my_rank == i) io::scatterplot_3d_output<GridPointTransformer>(i, "debug-domain-decomposition-"+std::to_string(0)+".dat", my_cells);
            MPI_Barrier(MPI_COMM_WORLD);
        }*/
        auto prev_imbalance = stats.global;

        auto all_loads = get_neighbors_load(stats.my_load, MPI_COMM_WORLD); //global load balancing with MPI_COMM_WORLD
        auto avg_load  = std::accumulate(all_loads.cbegin(), all_loads.cend(), 0.0) / worldsize;

        //DiffusiveOptimizer<GridPointTransformer, GridElementComputer, Cell> diff_opt;

        //diff_opt.optimize(part, my_cells, datatype_wrapper.element_datatype, 1.0);

        stats = part.get_load_statistics<GridElementComputer>(my_cells);
        if(!my_rank)
            print_load_statitics(stats);
        std::vector<double> all_mu(8, mu);
        LinearHashMap<int, double, 8> vertices_mu;
        std::zip(part.vertices_id.begin(), part.vertices_id.end(), all_mu.begin(), all_mu.end(), vertices_mu.begin());
        LinearHashMap <int, int, 8> vertices_remaining_trials;
        std::transform(part.vertices_id.begin(), part.vertices_id.end(), vertices_remaining_trials.begin(),
                       [](auto id){return std::make_pair(id, 10);});
        for(int j = 0; j < 10; ++j) {
            auto data = part.move_vertices<GridPointTransformer, GridElementComputer, Cell>(my_cells, datatype, avg_load, 1.0, vertices_remaining_trials);
            if(prev_imbalance <= stats.global) mu *= 0.9;
            prev_imbalance = stats.global;

            stats = part.get_load_statistics<GridElementComputer>(my_cells);
            if(!my_rank)
                print_load_statitics(stats);
        }

    }
};

}

#endif //ADLBIRREG_DIFFUSIVEOPTIMIZER_HPP
