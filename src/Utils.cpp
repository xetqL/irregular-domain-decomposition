//
// Created by xetql on 11/1/19.
//
#include <Utils.hpp>

LinearHashMap<int, std::vector<double>, 8>
get_neighbors_load(double my_load, const std::set<int> &neighbors, const std::map<int, Communicator> &v_neighborhood) {
    const auto nb_neighbors = neighbors.size();
    std::array<std::pair<int, std::vector<double>>, 8> neighborhoods_load;

    std::vector<MPI_Request> srequests(nb_neighbors);
    int cnt = 0;
    for (int neighbor_rank : neighbors){
        MPI_Isend(&my_load, 1, MPI_DOUBLE, neighbor_rank, 4011, MPI_COMM_WORLD, &srequests[cnt]);
        cnt++;
    }

    std::vector<double> neighbors_load(nb_neighbors);
    cnt = 0;
    for (int neighbor_rank : neighbors){
        MPI_Recv(&neighbors_load[cnt], 1, MPI_DOUBLE, neighbor_rank, 4011, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cnt++;
    }

    MPI_Waitall(nb_neighbors, srequests.data(), MPI_STATUSES_IGNORE);

    cnt = 0;
    for (const auto &v_comm : v_neighborhood) {
        const auto &comm_ranks = v_comm.second.get_ranks();
        std::vector<double> loads;
        int id_neighbor = 0;
        for (int neighbor_rank : neighbors) {
            if (std::find(comm_ranks.cbegin(), comm_ranks.cend(), neighbor_rank) != comm_ranks.cend()) { //in here?
                loads.push_back(neighbors_load[id_neighbor]);
            }
            id_neighbor++;
        }
        neighborhoods_load[cnt] = std::make_pair(v_comm.first, loads);
        cnt++;
    }
    return neighborhoods_load;
}

std::pair<int, std::vector<double>> get_neighbors_load(double my_load, int vid, const Communicator &v_comm) {
    auto neighbors = v_comm.get_ranks();
    const auto nb_neighbors = neighbors.size();
    //std::array<std::pair<int, std::vector<double>>, 8> neighborhoods_load;

    std::vector<MPI_Request> srequests(nb_neighbors);
    int cnt = 0;
    for (int neighbor_rank : neighbors){
        MPI_Isend(&my_load, 1, MPI_DOUBLE, neighbor_rank, 4011, MPI_COMM_WORLD, &srequests[cnt]);
        cnt++;
    }

    std::vector<double> neighbors_load(nb_neighbors);
    cnt = 0;
    for (int neighbor_rank : neighbors){
        MPI_Recv(&neighbors_load[cnt], 1, MPI_DOUBLE, neighbor_rank, 4011, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cnt++;
    }

    MPI_Waitall(nb_neighbors, srequests.data(), MPI_STATUSES_IGNORE);

    //for (const auto &v_comm : v_neighborhood) {
    const auto &comm_ranks = v_comm.get_ranks();
    std::vector<double> loads;
    int id_neighbor = 0;
    for (int neighbor_rank : neighbors) {
        if (std::find(comm_ranks.cbegin(), comm_ranks.cend(), neighbor_rank) != comm_ranks.cend()) { //in here?
            loads.push_back(neighbors_load[id_neighbor]);
        }
        id_neighbor++;
    }
    return std::make_pair(vid, loads);
}

std::set<int> filter_active_neighbors(const std::array<type::VertexIndex, 8>& vertices_id,
                                      const LinearHashMap<type::VertexIndex, int, 8>&       vertices_trial,
                                      const std::map<int, Communicator>& vertex_neighborhood) {
    std::set<int> active_neighbors;
    for(int vid : vertices_id) {
        const int vertex_trial = (*search_in_linear_hashmap<type::VertexIndex, int, 8>(vertices_trial, vid)).second;
        if(vertex_trial > 0) {
            auto active_ranks = vertex_neighborhood.at(vid).get_ranks();
            active_neighbors.insert(active_ranks.begin(), active_ranks.end());
        }
    }
    return active_neighbors;
}

bool file_exists(const std::string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}