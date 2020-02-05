//
// Created by xetql on 10/21/19.
//

#include "Partition.hpp"
namespace lb {

Real compute_mu(Real grid_size, Real max_normalized_load){
    Real sigma_max = max_normalized_load-1;
    return sigma_max == 0 ? 0 : grid_size/(sigma_max);
}

int get_rank_from_vertices(const std::array<int, 4>& vertices_id, const std::map<int, Communicator>& neighborhoods){
    get_MPI_rank(my_rank);
    std::vector<int> ranks;
    std::set<int> unique_ranks;
    for(int vid : vertices_id) {
        auto neighbors = neighborhoods.at(vid).get_ranks();
        std::copy(neighbors.begin(), neighbors.end(), std::back_inserter(ranks));
        std::copy(neighbors.begin(), neighbors.end(), std::inserter(unique_ranks, unique_ranks.begin()));
    }
    for(int r : unique_ranks) {
        if(std::count(ranks.cbegin(), ranks.cend(), r) == 4 && r != my_rank) return r;
    }
    throw std::runtime_error("nobody owns the 4 vertices?");
}

std::array<std::set<int>, 12> get_ranks_per_plane(
        const std::array<VertexIndex, 8>& vids,
        const std::map<int, Communicator>& comms) {
    auto planes_vid = get_planes_vids(vids);
    std::array<std::set<int>, 12> ranks_per_plane;
    for(int i = 0; i < 12; ++i) {
        ranks_per_plane[i] = get_ranks_from_vertices<3>(planes_vid[i], comms);
    }
    return ranks_per_plane;
}

int translate_iteration_to_vertex_group(int physical_iteration, lb::Point_3 coord){
    auto iteration = physical_iteration % 8;
    auto x = iteration&1; auto y = (iteration&2)>>1; auto z = (iteration&4)>>2;
    int com = (((unsigned int) coord.x()+x) & 1)  +
              (((unsigned int) coord.y()+y) & 1)*2+
              (((unsigned int) coord.z()+z) & 1)*4;
    return com;
}

LoadStatistics Partition::get_neighborhood_load_statistics(int com, Real my_load, int count){
    int vid = this->vertices_id[com];
    Communicator communicator = this->vertex_neighborhood[vid];

    unsigned int N = communicator.comm_size;
    std::vector<Real> all_loads(N);
    std::vector<int> buf(N);

    communicator.Allgather(&my_load, 1, MPI_TYPE_REAL, all_loads.data(), 1, MPI_TYPE_REAL, 87650);
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