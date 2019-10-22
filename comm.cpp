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
    
    if(rank == 0) {
        int x;
        MPI_Ssend(nullptr, 0, MPI_INT, 1, 101, MPI_COMM_WORLD);
        std::cout << "0 send!" <<std::endl;

    } else if(rank == 1) {

        MPI_Recv(nullptr, 0, MPI_INT, 0, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::cout << "1 recv!" <<std::endl;
    }
    
    MPI_Finalize();
}