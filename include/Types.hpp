//
// Created by xetql on 1/10/20.
//

#ifndef ADLBIRREG_TYPES_HPP
#define ADLBIRREG_TYPES_HPP

#include <mpi.h>
namespace type
{
    /* Type for indexing the data */
    using DataIndex   = unsigned long;
    extern MPI_Datatype MPI_TYPE_DATA_INDEX;
    /* Type for enumerating Processor Rank (MPI REQUIREMENT) */
    using ProcRank    = int;
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /* Shortcut for real datatype. this should be used throughout the project for type consistency */
    using Real        = double;
    extern MPI_Datatype MPI_TYPE_REAL;
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /* Type for indexing vertices */
    using VertexIndex = unsigned int;
    extern MPI_Datatype MPI_TYPE_VERTEX_INDEX;
}

#endif //ADLBIRREG_TYPES_HPP
