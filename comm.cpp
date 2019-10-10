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
    
    Communicator c({0});
    
    std::vector<long long> x  (20000);
    std::iota(x.begin(), x.end(), rank);
    std::vector<long long> x2 (x.size());

    c.Allgather(x.data(), x.size(), MPI_LONG_LONG, x2.data(), x.size(), MPI_LONG_LONG);

    if(rank == 0) {
        std::for_each(x2.cbegin(), x2.cend(), [&](auto val){std::cout << val << std::endl;});
    }
    
    MPI_Finalize();
}