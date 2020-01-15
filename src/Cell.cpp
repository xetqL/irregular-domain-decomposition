//
// Created by xetql on 10/21/19.
//

#include "Cell.hpp"

type::DataIndex mesh::compute_lid(type::DataIndex msx, type::DataIndex msy, type::DataIndex msz, type::DataIndex gid, lb::Box3 bbox){
    type::DataIndex x, y, z;
    lb::cell_to_local_position(msx, msy, msz, bbox, gid, &x, &y, &z);
    return lb::grid_index_to_cell(x, y, z, bbox.size_x+1, bbox.size_y+1, bbox.size_z+1);
}

type::DataIndex mesh::get_lid(type::DataIndex msx, type::DataIndex msy, type::DataIndex msz, type::DataIndex gid, lb::Box3 bbox){
    type::DataIndex x, y, z;
    lb::cell_to_local_position(msx, msy, msz, bbox, gid, &x, &y, &z);
    return lb::grid_index_to_cell(x, y, z, bbox.size_x, bbox.size_y, bbox.size_z);
}