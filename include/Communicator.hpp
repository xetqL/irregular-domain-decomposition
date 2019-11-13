//
// Created by xetql on 10/10/19.
//

#ifndef ADLBIRREG_COMMUNICATOR_HPP
#define ADLBIRREG_COMMUNICATOR_HPP

#include <algorithm>
#include <mpi.h>
#include <vector>
#include <cmath>
#include <assert.h>
#include <map>

#define get_MPI_rank(rank_var)\
    int rank_var;\
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_var)

#define get_MPI_worldsize(worldsize_var)\
    int worldsize_var;\
    MPI_Comm_size(MPI_COMM_WORLD, &worldsize_var)

inline short get_type_size(MPI_Datatype type) {
    if(type == MPI_INT)       return sizeof(int);
    if(type == MPI_DOUBLE)    return sizeof(double);
    if(type == MPI_FLOAT)     return sizeof(float);
    if(type == MPI_SHORT)     return sizeof(short);
    if(type == MPI_LONG)      return sizeof(long);
    if(type == MPI_LONG_LONG) return sizeof(long long);
    return -1;
}

struct CommunicationDatatype {
    MPI_Datatype element_datatype;
    MPI_Datatype minimal_datatype;
    CommunicationDatatype(const MPI_Datatype &el, const MPI_Datatype &min) : element_datatype(el), minimal_datatype(min){}
    void free_datatypes() { MPI_Type_free(&element_datatype); MPI_Type_free(&minimal_datatype);}
};

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


    void Bcast_hypercube(void* buffer, int count, MPI_Datatype type, int root, int tag)  const {
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
                    MPI_Send(buffer, count, type, comm_to_world.at(virtual_dest ^ world_to_comm.at(root)), tag, MPI_COMM_WORLD);
                } else {
                    int virtual_src   = virtual_rank ^ ((int) std::pow(2, i));
                    MPI_Recv(buffer, count, type, comm_to_world.at( virtual_src ^ world_to_comm.at(root)), tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
    }

    void Bcast(void* buffer, int count, MPI_Datatype type, int root, int tag) const {
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
                MPI_Isend(buffer, count, type, n, tag, MPI_COMM_WORLD, &reqs[i]);
                ++i;
            }
        else MPI_Recv(buffer, count, type, root, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Waitall(i, reqs.data(), MPI_STATUSES_IGNORE);
    }

    void Allgather(void* sendbuffer, int sendcount, MPI_Datatype sendtype, void* recvbuffer, int recvcount, MPI_Datatype recvtype, int tag) const {
        if (world_to_comm.find(world_rank) == world_to_comm.end()) {
            return;
        }
        if (comm_to_world.size() == 1) {
            return;
        }
        int i = 0;
        short size = get_type_size(sendtype);

        int dim = (int)(std::log2(comm_to_world.size()));
        if(std::pow(2, dim) == comm_to_world.size())
            for(auto root : comm_to_world) {
                Bcast_hypercube(root == world_rank ? sendbuffer : ((char*) recvbuffer + i*recvcount*size), sendcount, sendtype, root, tag);
                if(root == world_rank){
                    std::copy((char*) sendbuffer, (char*)sendbuffer + sendcount * size, ((char*) recvbuffer + i*recvcount*size));
                }
                i++;
            }
        else
            for(auto root : comm_to_world) {
                Bcast(root == world_rank ? sendbuffer : ((char*) recvbuffer + i*recvcount*size), sendcount, sendtype, root, tag);
                if(root == world_rank){
                    std::copy((char*) sendbuffer, (char*)sendbuffer + sendcount * size, ((char*) recvbuffer + i*recvcount*size));
                }
                i++;
            }
    }

    // Peer-to-peer communications are just wrapper around MPI directives
    int Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag) const {
        assert(std::find(comm_to_world.cbegin(), comm_to_world.cend(), dest) != comm_to_world.cend());
        return MPI_Send(buf, count ,datatype, dest, tag, MPI_COMM_WORLD);
    }

    int Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Request *request) const {
        assert(std::find(comm_to_world.cbegin(), comm_to_world.cend(), dest) != comm_to_world.cend());
        return MPI_Isend(buf, count ,datatype, dest, tag, MPI_COMM_WORLD, request);
    }

    int Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag) const {
        assert(std::find(comm_to_world.cbegin(), comm_to_world.cend(), source) != comm_to_world.cend());
        return MPI_Recv(buf, count, datatype, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    int Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Request *request) const {
        assert(std::find(comm_to_world.cbegin(), comm_to_world.cend(), source) != comm_to_world.cend());
        return MPI_Irecv(buf, count, datatype, source, tag, MPI_COMM_WORLD, request);
    }

    const std::vector<int> get_ranks() const { return comm_to_world; }
    const int translate_rank_to_local(int world_rank) const { return world_to_comm.at(world_rank);}
    const int translate_rank_to_global(int local_rank) const { return comm_to_world[local_rank];}

    //const std::vector<int> get_neighbors() const { return comm_to_world; }

};


#endif //ADLBIRREG_COMMUNICATOR_HPP
