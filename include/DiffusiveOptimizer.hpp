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
        std::vector<double> loads(worldsize), all_mu(mu, 8);
        auto my_load = lc.compute_load(my_cells);
        double max_load;
        MPI_Allreduce(&my_load, &max_load, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        auto avg_load = max_load / worldsize;

        LinearHashMap<int, double, 8> vertices_mu;
        LinearHashMap<int, double, 8> previous_imbalance;
        std::transform(part.vertices_id.begin(), part.vertices_id.end(), previous_imbalance.begin(),
                [](auto id){return std::make_pair(id, 10000.0);});
        std::zip(part.vertices_id.begin(), part.vertices_id.end(), all_mu.begin(), all_mu.end(), vertices_mu.begin());

        LinearHashMap <int, int,  8> vertices_remaining_trials;
        std::transform(part.vertices_id.begin(), part.vertices_id.end(), previous_imbalance.begin(),
                       [](auto id){return std::make_pair(id, 3);});
        LinearHashMap <int, bool, 8> vertices_status; //active = true, inactive = false
        std::transform(part.vertices_id.begin(), part.vertices_id.end(), previous_imbalance.begin(),
                       [](auto id){return std::make_pair(id, true);});

        while(std::any_of(vertices_status.begin(), vertices_status.end(), [](auto status) {return status.second;})) {
            part.move_vertices<GridPointTransformer, GridElementComputer, Cell>(my_cells, datatype, avg_load, vertices_mu);
            my_load = lc.compute_load(my_cells);
            LinearHashMap <int, std::vector<int>, 8> neighbors_status = get_neighbors_info(overloading, MPI_INT,    part.neighbor_list, part.vertex_neighborhood);
            auto neighbors_load   = get_neighbors_info(my_load,     MPI_DOUBLE, part.neighbor_list, part.vertex_neighborhood);
            LinearHashMap<int,  bool, 8> strategy;
            /// TODO: populate strategy linear hash map
            // std::transform(neighbors_status.cbegin(), neighbors_status.cend(), strategy.begin(),
            //        [](std::pair<int, std::vector<int>> statuses){ return std::any_of(statuses.second.cbegin(), statuses.second.cend(), [](int status){return (bool) status;});});
            for(int vid : part.vertices_id) {
                bool& vertex_status = (*search_in_linear_hashmap<int, bool, 8>(strategy, vid)).second;
                if(vertex_status){
                    auto current_neighborhood_load = (*search_in_linear_hashmap<int, std::vector<double>, 8>(neighbors_load, vid)).second;
                    double imbalance   = (*std::max_element(loads.cbegin(), loads.cend()) / avg_load) - 1.0; //local for neighborhood
                    double& vertex_mu         = (*search_in_linear_hashmap<int, double, 8>(vertices_mu, vid)).second;
                    vertex_status = imbalance > 0.2;
                    mu = 0.0;
                }
            }
        }
    }
};

}

#endif //ADLBIRREG_DIFFUSIVEOPTIMIZER_HPP
