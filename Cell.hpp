//
// Created by xetql on 10/21/19.
//

#ifndef ADLBIRREG_CELL_HPP
#define ADLBIRREG_CELL_HPP

#include <GeometricUtils.hpp>
#include <Communicator.hpp>

struct Cell {
    int gid, type; // type = -1:empty, 0:rock, 1:water
    float weight, erosion_probability;
    double average_load = 2000;

    Cell() : gid(0), type(0), weight(1.0), erosion_probability(0) {};
    Cell(int gid, int type, float weight, float erosion_probability)
            : gid(gid), type(type), weight(weight), erosion_probability(erosion_probability) {};
    Cell(int gid, int type, float weight, float erosion_probability, double average_load)
            : gid(gid), type(type), weight(weight), erosion_probability(erosion_probability), average_load(average_load) {};

    template<class NumericType>
    std::array<NumericType, 3> get_position_as_array() const {
        auto position_as_pair = lb::cell_to_global_position(Cell::get_msx(), Cell::get_msy(), gid);
        std::array<NumericType, 2> array = {(NumericType) position_as_pair.first, (NumericType) position_as_pair.second};
        return array;
    }

    std::tuple<int, int, int> get_position_as_pair() const {
        int x,y,z;
        lb::linear_to_grid(gid, Cell::get_msx(), Cell::get_msy(), x, y, z); //cell_to_global_position(Cell::get_msx(), Cell::get_msy(), gid);
        return std::make_tuple(x,y,z);
    }

    std::tuple<double, double, double> get_center() const {
        int x, y, z;
        lb::linear_to_grid(gid, Cell::get_msx(), Cell::get_msy(), x, y, z); //cell_to_global_position(Cell::get_msx(), Cell::get_msy(), gid);
        return std::make_tuple(x*Cell::get_cell_size()+Cell::get_cell_size()/2.0, y*Cell::get_cell_size()+Cell::get_cell_size()/2.0, z*Cell::get_cell_size()+Cell::get_cell_size()/2.0);
    }

    void get_center(double* _x, double* _y, double* _z) const {
        int x,y,z;
        lb::linear_to_grid(gid, Cell::get_msx(), Cell::get_msy(), x, y, z); //cell_to_global_position(Cell::get_msx(), Cell::get_msy(), gid);
        *_x=x+Cell::get_cell_size(); *_y=y+Cell::get_cell_size(); *_z =z+Cell::get_cell_size();
    }

    static CommunicationDatatype register_datatype() {

        MPI_Datatype cell_datatype, gid_type_datatype;

        MPI_Aint intex, lb, floatex;

        const int number_of_int_elements    = 2;
        const int number_of_float_elements  = 2;
        const int number_of_double_elements = 1;

        int blockcount_element[3];

        blockcount_element[0] = number_of_int_elements;    // gid, lid, exit, waiting_time
        blockcount_element[1] = number_of_float_elements;  // position <x,y>
        blockcount_element[2] = number_of_double_elements; // position <x,y>
        //int
        MPI_Type_contiguous(number_of_int_elements, MPI_INT, &gid_type_datatype); // position
        MPI_Type_commit(&gid_type_datatype);

        MPI_Datatype blocktypes[3];
        blocktypes[0] = MPI_INT;
        blocktypes[1] = MPI_FLOAT;
        blocktypes[2] = MPI_DOUBLE;

        MPI_Type_get_extent(MPI_INT, &lb, &intex);
        MPI_Type_get_extent(MPI_FLOAT, &lb, &floatex);

        MPI_Aint offset[3];
        offset[0] = static_cast<MPI_Aint>(0);
        offset[1] = 2*intex;
        offset[2] = 2*intex + 2*floatex;

        MPI_Type_create_struct(3, blockcount_element, offset, blocktypes, &cell_datatype);
        MPI_Type_commit(&cell_datatype);

        return {cell_datatype, gid_type_datatype};
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
        os << "gid: " << cell.gid << " type: " << cell.type << " weight: " << cell.weight << " erosion_probability: "
           << cell.erosion_probability;
        return os;
    }
};

#endif //ADLBIRREG_CELL_HPP
