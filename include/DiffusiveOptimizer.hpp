//
// Created by xetql on 11/4/19.
//

#ifndef ADLBIRREG_DIFFUSIVEOPTIMIZER_HPP
#define ADLBIRREG_DIFFUSIVEOPTIMIZER_HPP

#include "Partition.hpp"
#include "Utils.hpp"

namespace lb
{

template<class A, class GridElementComputer >
class DiffusiveOptimizer {

public:
    void optimize_neighborhood(int com, Partition& part, std::vector<mesh::Cell<A>> & my_cells, Real init_mu) {
        get_MPI_rank(my_rank);
        get_MPI_worldsize(worldsize);
        auto vid = part.vertices_id[com];

#ifdef DEBUG
        std::cout << my_rank << " enters " << __func__ << " with vid "<< vid << " with size " << part.vertex_neighborhood[vid].comm_size << std::endl;
#endif

        auto stats = part.get_neighborhood_load_statistics<GridElementComputer>(com, my_cells);
        auto my_load = stats.my_load;
        auto prev_imbalance = stats.global;
        auto avg_load = stats.avg_load;

        Real delta_load = 0;
        unsigned int remaining_trials = 2;
        Real mu = init_mu;
        while(remaining_trials) {
            part.move_selected_vertices<GridElementComputer, A> (com, my_cells, avg_load, mu, &delta_load);
            my_load += delta_load;
            auto current_stats = part.get_neighborhood_load_statistics(com, my_load, my_cells.size());

            if(prev_imbalance <= current_stats.global) {
                remaining_trials--;
                mu /= 2.0;
            } else {
                remaining_trials = 2;
                prev_imbalance = current_stats.global;
            }
        }
    }
};

}

#endif //ADLBIRREG_DIFFUSIVEOPTIMIZER_HPP
