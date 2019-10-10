//
// Created by xetql on 10/9/19.
//
#include <mpi.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <numeric>

#include "Communicator.hpp"

int main(int argc, char** argv){
    MPI_Init(nullptr, nullptr);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    Communicator c({0,  2, 4, 6});
    
    std::vector<long long> x  (3);
    std::iota(x.begin(), x.end(), rank);
    std::vector<long long> x2 (c.comm_size * x.size());

    c.Allgather(x.data(), x.size(), MPI_LONG_LONG, x2.data(), x.size(), MPI_LONG_LONG);

    if(rank == 0) {
        std::for_each(x2.cbegin(), x2.cend(), [&](auto val){std::cout << val << std::endl;});
    }
    
    MPI_Finalize();
}