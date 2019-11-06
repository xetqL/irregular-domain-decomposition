//
// Created by xetql on 10/21/19.
//

#include "Domain.hpp"

namespace lb {

Partition& Domain::get_my_partition(int rank) {
    return partitions.at(rank);
}

Point_3 move_vertex(const Point_3& vertex, const Vector_3& force, double mu){
    return vertex + mu*force;
}

}
