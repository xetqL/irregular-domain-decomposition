//
// Created by xetql on 10/21/19.
//

#ifndef ADLBIRREG_UTILS_HPP
#define ADLBIRREG_UTILS_HPP

#include <algorithm>
#include <vector>
#include <array>
#include <mpi.h>
template<class A>
std::vector<A>& search_in_linear_hashmap(std::array<std::pair<int, std::vector<A>>, 8>& linmap, int key) {
    for(std::pair<int, std::vector<A>>& entry : linmap) {
        if(entry.first == key) return entry.second;
    }
    throw std::runtime_error("doesn't exist");
}

/**
 * To compute the local average load use local communicator, otherwise use MPI_COMM_WORLD
 * @param my_load
 * @param average_load
 * @param neighborhood
 */
inline void compute_average_load(const double my_load, double* average_load, MPI_Comm neighborhood){
    int N;
    MPI_Comm_size(neighborhood, &N);
    MPI_Allreduce(&my_load, average_load, 1, MPI_DOUBLE, MPI_SUM, neighborhood);
    *average_load /= N;
}

inline std::vector<double> get_neighbors_load(double my_load, MPI_Comm neighborhood){
    int N;
    MPI_Comm_size(MPI_COMM_WORLD, &N);
    std::vector<double> all_loads(N);
    MPI_Allgather(&my_load, 1, MPI_DOUBLE, all_loads.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    return all_loads;
}

#endif //ADLBIRREG_UTILS_HPP
