//
// Created by xetql on 10/10/19.
//

#ifndef ADLBIRREG_COMMUNICATOR_HPP
#define ADLBIRREG_COMMUNICATOR_HPP

#include <algorithm>
#include <mpi.h>
#include <vector>
#include <cmath>

class Communicator {
    std::vector<int> comm_to_world;
    std::map<int, int> world_to_comm;

public:
    int world_rank, comm_size;

    Communicator(): comm_to_world({}){}
    Communicator(std::vector<int> neighbors) : comm_to_world(std::move(neighbors)) {
        std::sort(comm_to_world.begin(), comm_to_world.end());
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        for(int i = 0; i < comm_to_world.size(); ++i) {
            world_to_comm[comm_to_world[i]] = i;
        }
        comm_size = comm_to_world.size();
    };
    Communicator(std::initializer_list<int> neighbors) : comm_to_world((neighbors)) {
        std::sort(comm_to_world.begin(), comm_to_world.end());
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        for(int i = 0; i < comm_to_world.size(); ++i) {
            world_to_comm[comm_to_world[i]] = i;
        }
        comm_size = comm_to_world.size();
    };

    void Bcast_hypercube(void* buffer, int count, MPI_Datatype type, int root)  const {
        if(world_to_comm.find(world_rank) == world_to_comm.end()) {
            return;
        }

        if(comm_to_world.size() == 1) {
            return;
        }

        int virtual_rank = world_to_comm.at(world_rank) ^ world_to_comm.at(root);
        int dim   = std::ceil(std::log2(comm_to_world.size()));
        int mask = (int) std::pow(2, dim) - 1;

        for(int i = 0; i < dim; ++i) {
            mask = mask ^ ((int) std::pow(2, i));
            if((virtual_rank & mask) == 0) {
                if((virtual_rank & ((int) std::pow(2, i))) == 0) {
                    int virtual_dest   = virtual_rank ^ ((int) std::pow(2, i));
                    MPI_Send(buffer, count, type, comm_to_world.at(virtual_dest ^ world_to_comm.at(root)), 0, MPI_COMM_WORLD);
                } else {
                    int virtual_src   = virtual_rank ^ ((int) std::pow(2, i));
                    MPI_Recv(buffer, count, type, comm_to_world.at( virtual_src ^ world_to_comm.at(root)), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
    }

    void Bcast(void* buffer, int count, MPI_Datatype type, int root) const {
        if (world_to_comm.find(world_rank) == world_to_comm.end()) {
            return;
        }
        if (comm_to_world.size() == 1) {
            return;
        }
        std::vector<MPI_Request> reqs(comm_size);
        int i = 0;
        if(root == world_rank)
            for(auto const n : comm_to_world) {
                MPI_Isend(buffer, count, type, n, 123456, MPI_COMM_WORLD, &reqs[i]);
                ++i;
            }
        else MPI_Recv(buffer, count, type, root, 123456, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Waitall(i, reqs.data(), MPI_STATUSES_IGNORE);
    }

    void Allgather(void* sendbuffer, int sendcount, MPI_Datatype sendtype, void* recvbuffer, int recvcount, MPI_Datatype recvtype) const {
        if (world_to_comm.find(world_rank) == world_to_comm.end()) {
            return;
        }
        if (comm_to_world.size() == 1) {
            return;
        }
        int i = 0;
        short size;
        if(sendtype == MPI_INT)size         = sizeof(int);
        if(sendtype == MPI_DOUBLE)size      = sizeof(double);
        if(sendtype == MPI_FLOAT)size       = sizeof(float);
        if(sendtype == MPI_SHORT)size       = sizeof(short);
        if(sendtype == MPI_LONG)size        = sizeof(long);
        if(sendtype == MPI_LONG_LONG)size   = sizeof(long long);

        int dim = (int)(std::log2(comm_to_world.size()));
        if(std::pow(2, dim) == comm_to_world.size())
            for(auto root : comm_to_world) {
                Bcast_hypercube(root == world_rank ? sendbuffer : ((char*) recvbuffer + i*recvcount*size), sendcount, sendtype, root);
                if(root == world_rank){
                    std::copy((char*) sendbuffer, (char*)sendbuffer + sendcount * size, ((char*) recvbuffer + i*recvcount*size));
                }
                i++;
            }
        else
            for(auto root : comm_to_world) {
                Bcast(root == world_rank ? sendbuffer : ((char*) recvbuffer + i*recvcount*size), sendcount, sendtype, root);
                if(root == world_rank){
                    std::copy((char*) sendbuffer, (char*)sendbuffer + sendcount * size, ((char*) recvbuffer + i*recvcount*size));
                }
                i++;
            }
    }

    std::vector<int> get_ranks() { return comm_to_world; }
};


#endif //ADLBIRREG_COMMUNICATOR_HPP
