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
    void optimize_neighborhood(int com, Partition& part, std::vector<Cell>& my_cells, MPI_Datatype datatype, double init_mu) {
        get_MPI_rank(my_rank);
        get_MPI_worldsize(worldsize);
        auto vid = part.vertices_id[com];

#ifdef DEBUG
        std::cout << my_rank << " enters " << __func__ << " with vid "<< vid << " with size " << part.vertex_neighborhood[vid].comm_size << std::endl;
#endif
        GridElementComputer lc;

        auto stats = part.get_neighborhood_load_statistics<GridElementComputer>(com, my_cells);
        auto my_load = stats.my_load;
        auto prev_imbalance = stats.global;
        auto avg_load = stats.avg_load;

        double delta_load = 0;
        unsigned int remaining_trials = 3;
        double mu = init_mu;
        while(remaining_trials) {
            auto data = part.move_selected_vertices<GridPointTransformer, GridElementComputer, Cell>
                    (com, my_cells, datatype, avg_load, mu, &delta_load);
            my_load += delta_load;
            auto current_stats = part.get_neighborhood_load_statistics(com, my_load, my_cells.size());

            if(prev_imbalance <= current_stats.global) {
                remaining_trials--;
                mu /= 2.0;
            } else {
                remaining_trials = 3;
                prev_imbalance = current_stats.global;
            }
        }
    }
};

}

#endif //ADLBIRREG_DIFFUSIVEOPTIMIZER_HPP
