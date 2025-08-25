#ifndef GRID_H
#define GRID_H


#include "StateVector.h"
#include "Vector3.h"
#include "Array2D.h"
#include <vector>


//Grid structure: 
//Takes num_x and num_y as constructor parameters
//Has all persistent and intermediate arrays as members
//Initialises arrays on construction 
//Contains a template function for converting between 2D indices and 1D arrays
//Contains various template functions for boundary conditions

struct Grid{
    
    size_t ghost_cells; //Won't realistically be changing 
    size_t num_xcells;
    size_t num_ycells;

    //Constructor

    Grid(size_t nx, size_t ny,size_t g_=2);

    //These values will be calculated in initialise function
    std::vector<double> x;
    std::vector<double> y;

    double dx;
    double dy;
    double c_h;
    double dt;
    double B_timestep;

    // RESISTIVITY

    Array2D<double> eta;
    Array2D<Vector3> grad_eta;

    
    //Persistent arrays
    Array2D<ConservedStateVector> U;
    Array2D<PrimitiveStateVector> Prim;
    
    //Named intermediate arrays (memory buffers)

    //  - FOR MUSCL HANCOCK
    Array2D<ConservedStateVector> uBarL;
    Array2D<ConservedStateVector> uBarR;

    Array2D<PrimitiveStateVector> primL;
    Array2D<PrimitiveStateVector> primR;

    //Buffer arrays
    Array2D<Vector3> B;
    Array2D<double> E;
    Array2D<Vector3> LaplacianB;
    Array2D<Vector3> J;





    void reset_intermediate_arrays();
};


















#endif