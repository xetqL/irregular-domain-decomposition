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
    
    Communicator c1({0,  2, 4, 6});
    Communicator c2({1,  3, 5, 7});

    std::vector<long long> x  {rank, rank, rank};

    std::vector<long long> x2 (c1.comm_size * x.size());

    if(rank % 2) {
        c2.Allgather(x.data(), x.size(), MPI_LONG_LONG, x2.data(), x.size(), MPI_LONG_LONG);
    } else {
        c1.Allgather(x.data(), x.size(), MPI_LONG_LONG, x2.data(), x.size(), MPI_LONG_LONG);
    }

    if(rank == 1) {
        std::for_each(x2.cbegin(), x2.cend(), [&](auto val){std::cout << val << std::endl;});
    }
    
    MPI_Finalize();
}