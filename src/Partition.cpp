//
// Created by xetql on 10/21/19.
//

#include "Partition.hpp"
namespace lb {

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

int translate_iteration_to_vertex_group(int physical_iteration, lb::Point_3 coord){
    auto iteration = physical_iteration % 8;
    auto x = iteration&1; auto y = (iteration&2)>>1; auto z = (iteration&4)>>2;
    int com = (((unsigned int) coord.x()+x) & 1)  +
              (((unsigned int) coord.y()+y) & 1)*2+
              (((unsigned int) coord.z()+z) & 1)*4;
    return com;
}

Partition::LoadStatistics Partition::get_neighborhood_load_statistics(int com, double my_load, int count){
    int vid = this->vertices_id[com];
    Communicator communicator = this->vertex_neighborhood[vid];

    unsigned int N = communicator.comm_size;
    std::vector<double> all_loads(N);
    std::vector<int> buf(N);

    communicator.Allgather(&my_load, 1, MPI_DOUBLE, all_loads.data(), 1, MPI_DOUBLE, 87650);
    communicator.Allgather(&count,   1, MPI_INT,    buf.data(),       1, MPI_INT,    87651);

    auto avg_load  = std::accumulate(all_loads.cbegin(), all_loads.cend(), 0.0) / N;
    auto max_load  =*std::max_element(all_loads.cbegin(), all_loads.cend());

    return {std::accumulate(buf.cbegin(), buf.cend(), 0),
            max_load,
            avg_load,
            my_load,
            (max_load / avg_load) - 1.0,
            (my_load / avg_load) - 1.0};
}

}