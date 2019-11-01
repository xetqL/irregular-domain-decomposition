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

template<int N, class A>
typename std::array<std::pair<int, std::vector<A>>, N>::iterator
search_in_linear_hashmap(std::array<std::pair<int, std::vector<A>>, N> &linmap, int key) {
    //for(std::pair<int, std::vector<A>>& entry : linmap) {
    auto entry = linmap.begin();
    for (; entry != linmap.end(); entry++) {
        if ((*entry).first == key) return entry;
    }
    return entry;
}

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

std::array<std::pair<int, std::vector<double>>, 8>
get_neighbors_load(double my_load, const std::set<int> &neighbors, const std::map<int, Communicator> &v_neighborhood);

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
