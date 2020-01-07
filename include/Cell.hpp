//
// Created by xetql on 10/21/19.
//

#ifndef ADLBIRREG_CELL_HPP
#define ADLBIRREG_CELL_HPP

#include <GeometricUtils.hpp>
#include <Communicator.hpp>
namespace mesh {
    enum TCellType {REAL_CELL = 1, EMPTY_CELL=2, GHOST_CELL=3};

template<class ContainedElement>
struct Cell {

    std::vector<ContainedElement> elements;

    int gid;
    TCellType type;

    Cell() : gid(0), type(TCellType::REAL_CELL) {};

    Cell(int gid, TCellType type) : gid(gid), type(type){};

    template<class NumericType>
    std::array<NumericType, 3> get_position_as_array() const {
        auto position_as_pair = lb::cell_to_global_position(Cell::get_msx(), Cell::get_msy(), gid);
        std::array<NumericType, 2> array = {(NumericType) position_as_pair.first, (NumericType) position_as_pair.second};
        return array;
    }

    std::tuple<int, int, int> get_position_as_pair() const {
        int x,y,z;
        lb::linear_to_grid(gid, Cell::get_msx(), Cell::get_msy(), x, y, z); //cell_to_global_position(Cell::get_msx(), Cell::get_msy(), gid);
        return std::make_tuple(x, y, z);
    }

    std::tuple<double, double, double> get_center() const {
        int x, y, z;
        lb::linear_to_grid(gid, Cell::get_msx(), Cell::get_msy(), x, y, z); //cell_to_global_position(Cell::get_msx(), Cell::get_msy(), gid);
        return std::make_tuple(x*Cell::get_cell_size()+Cell::get_cell_size()/2.0, y*Cell::get_cell_size()+Cell::get_cell_size()/2.0, z*Cell::get_cell_size()+Cell::get_cell_size()/2.0);
    }

    static CommunicationDatatype register_datatype() {

        MPI_Datatype cell_datatype, gid_type_datatype;

        MPI_Aint intex, lb, floatex;

        const int number_of_int_elements    = 1;

        int blockcount_element[1];

        blockcount_element[0] = number_of_int_elements;    // gid, lid, exit, waiting_time

        //int
        MPI_Type_contiguous(number_of_int_elements, MPI_INT, &gid_type_datatype);
        MPI_Type_commit(&gid_type_datatype);

        return {gid_type_datatype, gid_type_datatype};
    }

    void add(ContainedElement e){
        elements.push_back(e);
    }

    static void set_msx(int _msx){
        static int msx = _msx;
    }

    static void set_msy(int _msy){
        static int msy = _msy;
    }

    static int& get_msx(){
        static int msx;
        return msx;
    }

    static int& get_msy(){
        static int msy;
        return msy;
    }

    static int& get_msz(){
        static int msz;
        return msz;
    }

    static double& get_cell_size(){
        static double size;
        return size;
    }

    friend bool operator==(const Cell &lhs, const Cell &rhs) {
        return lhs.gid == rhs.gid;
    }

    friend bool operator!=(const Cell &lhs, const Cell &rhs) {
        return !(rhs == lhs);
    }

    friend std::ostream &operator<<(std::ostream &os, const Cell &cell) {
        double x,y,z;
        std::tie(x,y,z) = cell.get_center();
        os << "gid: " << cell.gid << " ("<<x<<";"<<y<<";"<<z <<") type: " << cell.type;
        return os;
    }

};

template<class T>
void insert_in_proper_cell(std::vector<Cell<T>>* cells, T&& element) {
#ifdef DEBUG
    auto dim = std::cbrt(cells->size());
    assert(dim == (int) dim); //dim is a perfect cubical root required for searching in linear 3D-array
#endif


}

}
#endif //ADLBIRREG_CELL_HPP
