//
// Created by xetql on 10/21/19.
//

#ifndef ADLBIRREG_UTILS_HPP
#define ADLBIRREG_UTILS_HPP

#include <algorithm>
#include <vector>
#include <array>
#include <mpi.h>
#include <fstream>
#include <Communicator.hpp>
#include <set>
#include "GeometricUtils.hpp"


template<class K, class V, int N>
using LinearHashMap = std::array<std::pair<K, V>, N>;

template<class K, class V, int N>
typename LinearHashMap<K, V, N>::iterator
search_in_linear_hashmap(std::array<std::pair<K, V>, N> &linmap, K key) {
    //for(std::pair<int, std::vector<A>>& entry : linmap) {
    auto entry = linmap.begin();
    for (; entry != linmap.end(); entry++) {
        if ((*entry).first == key) return entry;
    }
    return entry;
}

template<class K, class V, int N>
typename LinearHashMap<K, V, N>::const_iterator
search_in_linear_hashmap(const std::array<std::pair<K, V>, N> &linmap, K key) {
    //for(std::pair<int, std::vector<A>>& entry : linmap) {
    auto entry = linmap.begin();
    for (; entry != linmap.end(); entry++) {
        if ((*entry).first == key) return entry;
    }
    return entry;
}
std::set<int> filter_active_neighbors(const std::array<int, 8>& vertices_id,
                                      const LinearHashMap<int, int, 8>& vertices_trial,
                                      const std::map<int, Communicator>& vertex_neighborhood);
/**
 * To compute the local average load use local communicator, otherwise use MPI_COMM_WORLD
 * @param my_load
 * @param average_load
 * @param neighborhood
 */
inline void compute_average_load(const double my_load, double *average_load, MPI_Comm neighborhood) {
    int N;
    MPI_Comm_size(neighborhood, &N);
    MPI_Allreduce(&my_load, average_load, 1, MPI_DOUBLE, MPI_SUM, neighborhood);
    *average_load /= N;
}

inline std::vector<double> get_neighbors_load(double my_load, MPI_Comm neighborhood) {
    int N;
    MPI_Comm_size(MPI_COMM_WORLD, &N);
    std::vector<double> all_loads(N);
    MPI_Allgather(&my_load, 1, MPI_DOUBLE, all_loads.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    return all_loads;
}

LinearHashMap<int, std::vector<double>, 8>
get_neighbors_load(double my_load, const std::set<int> &neighbors, const std::map<int, Communicator> &v_neighborhood);

template<class A, int N=8>
LinearHashMap<int, std::vector<A>, N>
get_neighbors_info(A my_info, MPI_Datatype datatype, const std::set<int> &neighbors, const std::map<int, Communicator> &v_neighborhood) {
    const auto nb_neighbors = neighbors.size();
    LinearHashMap<int, std::vector<A>, N> neighborhoods_info;

    std::vector<MPI_Request> srequests(nb_neighbors);
    int cnt = 0;
    for (int neighbor_rank : neighbors){
        MPI_Isend(&my_info, 1, datatype, neighbor_rank, 40111, MPI_COMM_WORLD, &srequests[cnt]);
        cnt++;
    }

    std::vector<A> neighbors_load(nb_neighbors);
    cnt = 0;
    for (int neighbor_rank : neighbors){
        MPI_Recv(&neighbors_load[cnt], 1, datatype, neighbor_rank, 40111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cnt++;
    }

    MPI_Waitall(nb_neighbors, srequests.data(), MPI_STATUSES_IGNORE);

    cnt = 0;
    for (const auto &v_comm : v_neighborhood) {
        const auto &comm_ranks = v_comm.second.get_ranks();
        std::vector<A> infos;
        int id_neighbor = 0;
        for (int neighbor_rank : neighbors) {
            if (std::find(comm_ranks.cbegin(), comm_ranks.cend(), neighbor_rank) != comm_ranks.cend()) { //in here?
                infos.push_back(neighbors_load[id_neighbor]);
            }
            id_neighbor++;
        }
        neighborhoods_info[cnt] = std::make_pair(v_comm.first, infos);
        cnt++;
    }
    return neighborhoods_info;
}


namespace std {
template<class InputIt, class BinaryFunction, class T>
T reduce(InputIt begin, InputIt end, BinaryFunction f, T start) {
    T result = start;
    while (begin != end) {
        result = f(result, *begin);
        begin++;
    }
    return result;
}

template<class InputIt1, class InputIt2, class OutputIt>
void zip(InputIt1 begin1, InputIt1 end1, InputIt2 begin2, InputIt2 end2, OutputIt out) {
    while (begin1 != end1 && begin2 != end2) {
        *out = std::make_pair(*begin1, *begin2);
        ++out; ++begin1; ++begin2;
    }
}
}

namespace io {
/**
 * Format of output is:
 * PROCID; x y z
 * */
template<class CartesianPointTransformer, class A>
void scatterplot_3d_output(int rank, std::string filename, const std::vector<A> &elements) {
    CartesianPointTransformer cpf;
    std::ofstream f(filename, std::ofstream::app);
    for (const A &el : elements) {
        lb::Point_3 p = cpf.transform(el);
        f << rank << " " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }
    f.close();
}
}
#endif //ADLBIRREG_UTILS_HPP
