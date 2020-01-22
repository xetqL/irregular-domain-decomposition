//
// Created by xetql on 10/21/19.
//

#include "Domain.hpp"
#include "cmath"

namespace lb {

Partition& Domain::get_partition(int rank) {
    return partitions.at(rank);
}

Point_3 move_vertex(const Point_3& vertex, const Vector_3& force, double mu){
    return vertex + mu*force;
}

void adjust_simulation_size(type::Real _DOMAIN_SIZE_X, type::Real _DOMAIN_SIZE_Y, type::Real _DOMAIN_SIZE_Z,
                            int procs_x, int procs_y, int procs_z,
                            type::Real grid_resolution,
                            type::Real* DOMAIN_SIZE_X, type::Real* DOMAIN_SIZE_Y, type::Real* DOMAIN_SIZE_Z){

    *DOMAIN_SIZE_X = std::ceil(_DOMAIN_SIZE_X / (procs_x * grid_resolution)) * (procs_x * grid_resolution);
    *DOMAIN_SIZE_Y = std::ceil(_DOMAIN_SIZE_Y / (procs_y * grid_resolution)) * (procs_y * grid_resolution);
    *DOMAIN_SIZE_Z = std::ceil(_DOMAIN_SIZE_Z / (procs_z * grid_resolution)) * (procs_z * grid_resolution);
}

}
