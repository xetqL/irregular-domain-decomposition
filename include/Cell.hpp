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
    enum TCellType {REAL_CELL = 1, EMPTY_CELL=2, GHOST_CELL=3, UNITIALIZED_CELL=-1};
    type::DataIndex compute_lid(type::DataIndex msx, type::DataIndex msy, type::DataIndex msz, type::DataIndex gid, lb::Box3 bbox);

template<class ContainedElement>
struct Cell {
    using IndexType = type::DataIndex;
    using value_type= ContainedElement;
    using Real      = type::Real;

    IndexType gid, lid;
    TCellType type;
private:
    std::vector<ContainedElement> elements;
public:
    Cell() : gid(0), lid(0), type(TCellType::UNITIALIZED_CELL) {};

    Cell(TCellType type) : gid(0), lid(0), type(type) {};

    Cell(IndexType gid, IndexType msx, IndexType msy, IndexType msz, lb::Box3 bbox) : gid(gid), type(TCellType::REAL_CELL) {
        update_lid(msx, msy, msz, bbox);
    };

    Cell(IndexType gid, TCellType type) : gid(gid), lid(0), type(type) {};

    Cell(IndexType gid, IndexType msx, IndexType msy, IndexType msz, lb::Box3 bbox, TCellType type) : gid(gid), type(type) {
        update_lid(msx, msy, msz, bbox);
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


    IndexType get_lid(IndexType msx, IndexType msy, IndexType msz, lb::Box3 bbox){
        return mesh::compute_lid(msx, msy, msz, gid, bbox);
    }

    IndexType get_lid(){
        return lid;
    }

    void update_lid(IndexType msx, IndexType msy, IndexType msz, lb::Box3 bbox){
        lid = get_lid(msx, msy, msz, bbox);
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
        return std::make_tuple(x, y, z);
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


template<class T>
void insert_or_remove(std::vector<Cell<T>>* _cells, std::vector<T>* _elements, lb::Box3 bbox) {
    std::vector<Cell<T>>& cells = *_cells;
    std::vector<T>& elements = *_elements;
    int nbp = 0;
    for(T& el : elements) {
        if(bbox.contains(el.position[0], el.position[1], el.position[2])) {
            //migrate the data
            auto indexes = lb::position_to_index(el.position[0],el.position[1],el.position[2], mesh::Cell<T>::get_cell_size());
            type::DataIndex ix = std::get<0>(indexes), iy = std::get<1>(indexes), iz = std::get<2>(indexes);
            type::DataIndex lid = (ix - bbox.x_idx_min) + bbox.size_x * (iy - bbox.y_idx_min) + bbox.size_y * bbox.size_x * (iz - bbox.z_idx_min);
            cells.at(lid).add(el);
            nbp++;
        }
    }
    std::cout << "Number of particle added: " << nbp << "/" << elements.size()<<  std::endl;
    elements.clear();
}
template<class T>
std::vector<Cell<T>> generate_lattice_single_type(type::DataIndex msx, type::DataIndex msy, type::DataIndex msz,
                                                  type::DataIndex x_proc_idx, type::DataIndex y_proc_idx, type::DataIndex z_proc_idx,
                                                  type::DataIndex cell_in_my_cols, type::DataIndex cell_in_my_rows, type::DataIndex cell_in_my_depth,
                                                  mesh::TCellType type) {
    using Cell = Cell<T>;
    auto cell_per_process = cell_in_my_cols * cell_in_my_rows;
    std::vector<Cell> my_cells; my_cells.reserve(cell_per_process);

    auto x_shift = (cell_in_my_rows *  x_proc_idx);
    auto y_shift = (cell_in_my_cols *  y_proc_idx);
    auto z_shift = (cell_in_my_depth * z_proc_idx);

    for(int z = 0; z < cell_in_my_depth; ++z) {
        for(int y = 0; y < cell_in_my_cols; ++y) {
            for(int x = 0; x < cell_in_my_rows; ++x) {
                auto gid = (x_shift + x) + (y_shift + y) * msx + (z_shift + z) * msx * msy;
                my_cells.emplace_back(gid, type);
            }
        }
    }

    return my_cells;
}

}
#endif //ADLBIRREG_CELL_HPP
