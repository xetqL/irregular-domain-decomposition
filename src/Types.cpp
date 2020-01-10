//
// Created by xetql on 1/10/20.
//

#ifndef ADLBIRREG_TYPES_HPP
#define ADLBIRREG_TYPES_HPP

#include <mpi.h>
namespace type
{
    /* Type for indexing the data */
    MPI_Datatype MPI_TYPE_DATA_INDEX = MPI_UNSIGNED_LONG;
    /* Type for enumerating Processor Rank (MPI REQUIREMENT) */
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /* Shortcut for real datatype. this should be used throughout the project for type consistency */
    MPI_Datatype MPI_TYPE_REAL= MPI_DOUBLE;
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /* Type for indexing vertices */
    MPI_Datatype MPI_TYPE_VERTEX_INDEX = MPI_UNSIGNED;
}
#endif //ADLBIRREG_TYPES_HPP
