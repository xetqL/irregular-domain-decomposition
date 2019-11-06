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
    void optimize(const Partition& part, std::vector<Cell>& my_cells, CommunicationDatatype datatype_wrapper, double mu, int overloading = 0) {
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
        int remaining_trial = 3;
        while(remaining_trial > 0){
            part.move_vertices<GridPointTransformer, GridElementComputer, Cell>(my_cells, datatype_wrapper.element_datatype, avg_load, vertices_mu);
            my_load = lc.compute_load(my_cells);
            // adapt mu as function of local balancing
            auto neighbors_status = get_neighbors_info(overloading, MPI_INT,    part.neighbor_list, part.vertex_neighborhood);
            auto neighbors_load   = get_neighbors_info(my_load,     MPI_DOUBLE, part.neighbor_list, part.vertex_neighborhood);
            LinearHashMap<int, bool, 8> strategy;
            std::transform(neighbors_status.cbegin(), neighbors_status.cend(), strategy.begin(),
                    [](auto statuses){ return std::any_of(statuses.cbegin(), statuses.cend());});

            for(int vid : part.vertices_id) {
                bool vertex_status = (*search_in_linear_hashmap<int, bool, 8>(strategy, vid)).second;
                auto loads         = (*search_in_linear_hashmap<int, std::vector<double>, 8>(neighbors_load, vid)).second;
                double imbalance   = *std::max_element(loads.cbegin(), loads.cend()) / avg_load; //local for neighborhood

                //if vertex status is ULBA => check with standard lb metric if we approach convergence w.r.t ULBA
                //  i.e., (1+alphaN(P-N)^-1)mu(i)
                //if not check normally
            }

            //if(prev_imbalance <= stats.global) mu *= 0.9;
            //prev_imbalance = stats.global;

        }
    }
};

}

#endif //ADLBIRREG_DIFFUSIVEOPTIMIZER_HPP
