//
// Created by xetql on 10/21/19.
//

#include "Partition.hpp"
namespace lb{
double compute_mu(double grid_size, double max_normalized_load){
    double sigma_max = max_normalized_load-1;
    return sigma_max == 0 ? 0 : grid_size/(sigma_max);
}

int get_rank_from_vertices(const std::array<int, 4>& vertices_id, const std::map<int, Communicator>& neighborhoods){
    get_MPI_rank(my_rank);
    std::vector<int> ranks;
    std::set<int> unique_ranks;
    for(int vid : vertices_id){
        auto neighbors = neighborhoods.at(vid).get_ranks();
        std::copy(neighbors.begin(), neighbors.end(), std::back_inserter(ranks));
        std::copy(neighbors.begin(), neighbors.end(), std::inserter(unique_ranks, unique_ranks.begin()));
    }
    for(int r : unique_ranks) {
        if(std::count(ranks.cbegin(), ranks.cend(), r) == 4 && r != my_rank) return r;
    }
    throw std::runtime_error("nobody owns the 4 vertices?");
}
}