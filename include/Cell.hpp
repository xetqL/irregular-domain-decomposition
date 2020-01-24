//
// Created by xetql on 10/21/19.
//

#ifndef ADLBIRREG_CELL_HPP
#define ADLBIRREG_CELL_HPP

#include <GeometricUtils.hpp>
#include <Communicator.hpp>
#include "Types.hpp"
#include <zupply.hpp>
#include <ostream>

namespace mesh {

enum TCellType {REAL_CELL = 1, EMPTY_CELL=2, GHOST_CELL=3, UNITIALIZED_CELL=-1};
inline type::DataIndex compute_lid(type::DataIndex msx, type::DataIndex msy, type::DataIndex msz, type::DataIndex gid, const lb::Box3& bbox){
    type::DataIndex x, y, z;
    lb::cell_to_local_position(msx, msy, msz, bbox, gid, &x, &y, &z);
    return lb::grid_index_to_cell(x, y, z, bbox.size_x, bbox.size_y, bbox.size_z);
}

class GridParams {
    type::DataIndex max_cell_index_x, max_cell_index_y, max_cell_index_z,
                    cell_number_x, cell_number_y, cell_number_z;
    type::Real      grid_resolution;
    type::Real      simulation_size_x, simulation_size_y, simulation_size_z;
    GridParams(){};
public:
    static GridParams& get_instance() {
        static GridParams gp;
        return gp;
    }

    GridParams(GridParams const&)      = delete;
    void operator=(GridParams const&)  = delete;

    void set_grid_resolution(type::Real _grid_resolution){
        grid_resolution = _grid_resolution;
    }

    void set_grid_index_dimensions(type::DataIndex _msx, type::DataIndex _msy, type::DataIndex _msz){
        cell_number_x = _msx;
        cell_number_y = _msy;
        cell_number_z = _msz;
        max_cell_index_x = _msx-1;
        max_cell_index_y = _msy-1;
        max_cell_index_z = _msz-1;
    }

    void set_simulation_dimension(type::Real msx, type::Real msy, type::Real msz){
        simulation_size_x = msx;
        simulation_size_y = msy;
        simulation_size_z = msz;
    }

    const type::DataIndex get_cell_number_x() const {
        return cell_number_x;
    }

    const type::DataIndex get_cell_number_y() const {
        return cell_number_y;
    }

    const type::DataIndex get_cell_number_z() const {
        return cell_number_z;
    }

    const type::Real& get_grid_resolution() const {
        return grid_resolution;
    }
    const type::DataIndex get_maximum_index_x() const {
        return max_cell_index_x;
    }
    const type::DataIndex get_maximum_index_y() const {
        return max_cell_index_y;
    }
    const type::DataIndex get_maximum_index_z() const {
        return max_cell_index_z;
    }
    const type::Real get_maximum_simulation_size_x() const {
        return simulation_size_x;
    }
    const type::Real get_maximum_simulation_size_y() const {
        return simulation_size_y;
    }
    const type::Real get_maximum_simulation_size_z() const {
        return simulation_size_z;
    }

    friend std::ostream &operator<<(std::ostream &os, const GridParams &params) {
        os << "max_cell_index_x: " << params.max_cell_index_x << " max_cell_index_y: " << params.max_cell_index_y
           << " max_cell_index_z: " << params.max_cell_index_z << " grid_resolution: " << params.grid_resolution
           << " simulation_size_x: " << params.simulation_size_x << " simulation_size_y: " << params.simulation_size_y
           << " simulation_size_z: " << params.simulation_size_z;
        return os;
    }
};

static GridParams& grid_params = GridParams::get_instance();

template<class ContainedElement>
struct Cell {
    using IndexType  = type::DataIndex;
    using value_type = ContainedElement;
    using Real       = type::Real;

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
        c.gid = lb::grid_index_to_cell(x, y, z, grid_params.get_cell_number_x(), grid_params.get_cell_number_y(), grid_params.get_cell_number_z());
        return c;
    }

    inline IndexType get_lid(IndexType msx, IndexType msy, IndexType msz, const lb::Box3& bbox) {
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
        lb::linear_to_grid(gid, grid_params.get_cell_number_x(), grid_params.get_cell_number_y(), x, y, z); //cell_to_global_position(msx, msy, gid);
        return std::make_tuple(x*grid_params.get_grid_resolution()+grid_params.get_grid_resolution()/2.0, y*grid_params.get_grid_resolution()+grid_params.get_grid_resolution()/2.0, z*grid_params.get_grid_resolution()+grid_params.get_grid_resolution()/2.0);
    }

    lb::Point_3 get_center_point() const {
        IndexType x, y, z;
        lb::linear_to_grid(gid, grid_params.get_cell_number_x(), grid_params.get_cell_number_y(), x, y, z); //cell_to_global_position(msx, msy, gid);
        return {x*grid_params.get_grid_resolution()+grid_params.get_grid_resolution()/2.0, y*grid_params.get_grid_resolution()+grid_params.get_grid_resolution()/2.0, z*grid_params.get_grid_resolution()+grid_params.get_grid_resolution()/2.0};
    }

    std::tuple<Real, Real, Real> get_coordinates() const {
        IndexType x, y, z;
        lb::linear_to_grid(gid, grid_params.get_cell_number_x(), grid_params.get_cell_number_y(), x, y, z); //cell_to_global_position(msx, msy, gid);
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
        os <<std::setprecision(15)<< "gid: " << cell.gid << " CENTER("<<cx<<";"<<cy<<";"<<cz <<")" << " COORD("<<x<<";"<<y<<";"<<z <<")" << " type: " << cell.type;
        return os;
    }
/*
    lb::Box3 as_box(){
        Real x, y, z;
        std::tie(x, y, z) = get_coordinates();
        return {x, x+grid_params.get_grid_resolution(),
                y, y+grid_params.get_grid_resolution(),
                z, z+grid_params.get_grid_resolution(),
                     grid_params.get_grid_resolution()};
    }*/

};


template<class T>
void insert_or_remove(std::vector<Cell<T>>* _cells, std::vector<T>* _elements, lb::Box3 bbox) {
    std::vector<Cell<T>>& cells = *_cells;
    std::vector<T>& elements = *_elements;
    int nbp = 0;
    for(T& el : elements) {
        if(bbox.contains(el.position[0], el.position[1], el.position[2])) {
            //migrate the data
            auto indexes = lb::position_to_index(el.position[0], el.position[1], el.position[2], grid_params.get_grid_resolution());
            type::DataIndex ix  = std::get<0>(indexes), iy = std::get<1>(indexes), iz = std::get<2>(indexes);
            type::DataIndex lid = (ix - bbox.x_idx_min) + bbox.size_x * (iy - bbox.y_idx_min) + bbox.size_y * bbox.size_x * (iz - bbox.z_idx_min);
            try {
                cells.at(lid).add(el);

            } catch(...) {
                std::cout << ix << " "<< iy<<" " <<iz << std::endl;
                std::cout << bbox << std::endl;
                std::cout << el << std::endl;
                throw;
            }
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
    auto cell_per_process = cell_in_my_cols * cell_in_my_rows * cell_in_my_depth;
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
