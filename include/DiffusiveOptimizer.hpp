//
// Created by xetql on 11/4/19.
//

#ifndef ADLBIRREG_DIFFUSIVEOPTIMIZER_HPP
#define ADLBIRREG_DIFFUSIVEOPTIMIZER_HPP

#include "Partition.hpp"
#include "Utils.hpp"

namespace lb
{
template<class GridPointTransformer, class GridElementComputer, class Cell>
class DiffusiveOptimizer {

public:
    void optimize(Partition& part, std::vector<Cell>& my_cells, MPI_Datatype datatype, double mu, int overloading = 0) {
        GridElementComputer lc;
        get_MPI_worldsize(worldsize);
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
        auto remaining_it = (int) std::cbrt(worldsize) ;
        LinearHashMap <int, int,  8> vertices_remaining_trials;
        std::transform(part.vertices_id.begin(), part.vertices_id.end(), vertices_remaining_trials.begin(),
                       [&remaining_it](auto id){return std::make_pair(id, 3 + remaining_it);});
        LinearHashMap <int, bool, 8> vertices_status; //active = true, inactive = false
        std::transform(part.vertices_id.begin(), part.vertices_id.end(), vertices_status.begin(),
                       [](auto id){return std::make_pair(id, true);});

        auto neighbors_status = get_neighbors_info(overloading, MPI_INT, part.neighbor_list, part.vertex_neighborhood);
        LinearHashMap<int,  bool, 8> strategy;
        std::transform(neighbors_status.cbegin(), neighbors_status.cend(), strategy.begin(),
           [](std::pair<int, std::vector<int>> statuses) {
               return std::make_pair(statuses.first, std::any_of(statuses.second.cbegin(), statuses.second.cend(), [](int status){return status;}));
           });

        get_MPI_rank(my_rank);
        while((remaining_it) > 0 || std::accumulate(vertices_remaining_trials.begin(), vertices_remaining_trials.end(), 0, [](int sum, auto rm) {return sum + rm.second;}) > 0) {
            std::cout << my_rank << " has "<< remaining_it <<" remaining it and "
            << std::accumulate(vertices_remaining_trials.begin(), vertices_remaining_trials.end(), 0, [](int sum, auto rm) {return sum + rm.second;}) << " active vertices" <<  std::endl;
            part.move_vertices<GridPointTransformer, GridElementComputer, Cell>(my_cells, datatype, avg_load, 1.0, vertices_remaining_trials);
            my_load = lc.compute_load(my_cells);
            neighbors_load  = get_neighbors_info(my_load,  MPI_DOUBLE, part.neighbor_list, part.vertex_neighborhood);
            for(int vid : part.vertices_id) {
                //const bool vertex_status   = (*search_in_linear_hashmap<int, bool,   8>(strategy, vid)).second; //ulba or not
                int& vertex_rem_trial = (*search_in_linear_hashmap<int, int,    8>(vertices_remaining_trials, vid)).second; // stop or not
                double& prev_imbl = (*search_in_linear_hashmap<int, double, 8>(previous_imbalance, vid)).second;

                if(vertex_rem_trial > 0) { // if continue
                    //compute imbalance
                    auto current_neighborhood_load = (*search_in_linear_hashmap<int, std::vector<double>, 8>(neighbors_load, vid)).second;
                    double imbalance  = (*std::max_element(current_neighborhood_load.cbegin(), current_neighborhood_load.cend()) / avg_load) - 1.0;
                    //new imbalance is good
                    if(imbalance < 0.2 || prev_imbl <= imbalance) {
                        vertex_rem_trial--;
                    } else {
                        vertex_rem_trial = 3;
                        prev_imbl = imbalance;
                    }
                }
                std::cout << my_rank << " with vid "<< vid << " has "<< prev_imbl << " rt: " << vertex_rem_trial << std::endl;

            }
            remaining_it--;
            MPI_Barrier(MPI_COMM_WORLD);
        }
        std::cout << my_rank << " has left"<<std::endl;
    }
};

}

#endif //ADLBIRREG_DIFFUSIVEOPTIMIZER_HPP
