//
// Created by xetql on 10/21/19.
//

#ifndef ADLBIRREG_CELL_HPP
#define ADLBIRREG_CELL_HPP

#include <GeometricUtils.hpp>
#include <Communicator.hpp>
#include "Types.hpp"
#include <zupply.hpp>

namespace mesh {
    enum TCellType {REAL_CELL = 1, EMPTY_CELL=2, GHOST_CELL=3};

template<class ContainedElement>
struct Cell {
    using IndexType = type::DataIndex;
    using value_type= ContainedElement;
    using Real      = type::Real;

    IndexType gid, lid;
    TCellType type;
private:
    Cell() : gid(0), lid(0), type(TCellType::REAL_CELL) {};
    std::vector<ContainedElement> elements;
public:

    Cell(TCellType type) : gid(0), lid(0), type(type) {};

    Cell(IndexType gid, lb::Box3 bbox) : gid(gid), type(TCellType::REAL_CELL) {
        update_lid(bbox);
    };

    Cell(IndexType gid, TCellType type) : gid(gid), lid(0), type(type) {};

    Cell(IndexType gid, lb::Box3 bbox, TCellType type) : gid(gid), type(type) {
        update_lid(bbox);
    };

    inline size_t get_number_of_elements(){
        return elements.size();
    }

    static Cell get_empty_cell(lb::Box3 bbox, IndexType lid){
        Cell c;
        c.lid  = lid;
        c.type = EMPTY_CELL;
        IndexType x,y,z;
        lb::linear_to_grid(lid, bbox.size_x, bbox.size_y, x, y, z);
        x+=bbox.x_idx_min;
        y+=bbox.y_idx_min;
        z+=bbox.z_idx_min;
        c.gid = lb::grid_index_to_cell(x, y, z, Cell::get_msx(), Cell::get_msy(), Cell::get_msz());
        return c;
    }

    void update_lid(lb::Box3 bbox){
        IndexType x, y, z;
        lb::cell_to_local_position(Cell::get_msx(), Cell::get_msy(), Cell::get_msz(), bbox, gid, &x, &y, &z);
        lid  = lb::grid_index_to_cell(x, y, z, bbox.size_x, bbox.size_y, bbox.size_z);
    }

    std::tuple<Real, Real, Real> get_center() const {
        IndexType x, y, z;
        lb::linear_to_grid(gid, Cell::get_msx(), Cell::get_msy(), x, y, z); //cell_to_global_position(Cell::get_msx(), Cell::get_msy(), gid);
        return std::make_tuple(x*Cell::get_cell_size()+Cell::get_cell_size()/2.0, y*Cell::get_cell_size()+Cell::get_cell_size()/2.0, z*Cell::get_cell_size()+Cell::get_cell_size()/2.0);
    }

    lb::Point_3 get_center_point() const {
        IndexType x, y, z;
        lb::linear_to_grid(gid, Cell::get_msx(), Cell::get_msy(), x, y, z); //cell_to_global_position(Cell::get_msx(), Cell::get_msy(), gid);
        return {x*Cell::get_cell_size()+Cell::get_cell_size()/2.0, y*Cell::get_cell_size()+Cell::get_cell_size()/2.0, z*Cell::get_cell_size()+Cell::get_cell_size()/2.0};
    }

    std::tuple<Real, Real, Real> get_coordinates() const {
        IndexType x, y, z;
        lb::linear_to_grid(gid, Cell<ContainedElement>::get_msx(), Cell<ContainedElement>::get_msy(), x, y, z); //cell_to_global_position(Cell::get_msx(), Cell::get_msy(), gid);
        return std::make_tuple(x*Cell<ContainedElement>::get_cell_size(), y*Cell<ContainedElement>::get_cell_size(), z*Cell<ContainedElement>::get_cell_size());
    }

    void add(ContainedElement e) {
        elements.push_back(e);
    }

    size_t number_of_elements() const {
        return elements.size();
    }

    typename std::vector<ContainedElement>::iterator begin(){
        return elements.begin();
    }

    typename std::vector<ContainedElement>::iterator end(){
        return elements.end();
    }

    static IndexType& get_msx(){
        static IndexType msx;
        return msx;
    }

    static IndexType& get_msy(){
        static IndexType msy;
        return msy;
    }

    static IndexType& get_msz(){
        static IndexType msz;
        return msz;
    }

    static Real& get_cell_size(){
        static Real size;
        return size;
    }

    friend bool operator==(const Cell &lhs, const Cell &rhs) {
        return lhs.gid == rhs.gid;
    }

    friend bool operator!=(const Cell &lhs, const Cell &rhs) {
        return !(rhs == lhs);
    }

    friend std::ostream &operator<<(std::ostream &os, const Cell &cell) {
        Real x,y,z, cx, cy,cz;
        std::tie(cx,cy,cz) = cell.get_center();
        std::tie(x,y,z)    = cell.get_coordinates();
        os << "gid: " << cell.gid << " CENTER("<<cx<<";"<<cy<<";"<<cz <<")" << " COORD("<<x<<";"<<y<<";"<<z <<")" << " type: " << cell.type;
        return os;
    }

    lb::Box3 as_box(){
        Real x, y, z;
        std::tie(x, y, z) = get_coordinates();
        return {x, x+Cell<ContainedElement>::get_cell_size(),
                y, y+Cell<ContainedElement>::get_cell_size(),
                z, z+Cell<ContainedElement>::get_cell_size(),
                     Cell<ContainedElement>::get_cell_size()};
    }

};

type::DataIndex compute_lid(type::DataIndex msx, type::DataIndex msy, type::DataIndex msz, type::DataIndex gid, lb::Box3 bbox);


template<class T>
void insert_or_remove(std::vector<Cell<T>>* _cells, std::vector<T>* _elements, lb::Box3 bbox) {
    std::vector<Cell<T>>& cells = *_cells;
    std::vector<T>& elements = *_elements;
    for(T& el : elements) {
        if(bbox.contains(el.position[0], el.position[1], el.position[2])) {
            //migrate the data
            auto indexes = lb::position_to_index(el.position[0],el.position[1],el.position[2], mesh::Cell<T>::get_cell_size());
            type::DataIndex ix = std::get<0>(indexes), iy = std::get<1>(indexes), iz = std::get<2>(indexes);
            type::DataIndex lid = (ix - bbox.x_idx_min) + bbox.size_x * (iy - bbox.y_idx_min) + bbox.size_y * bbox.size_x * (iz - bbox.z_idx_min);
            lb::Box3 b = cells.at(lid).as_box();
            cells.at(lid).add(el);
        }
    }
    elements.clear();
}


}
#endif //ADLBIRREG_CELL_HPP
