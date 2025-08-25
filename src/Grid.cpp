#include "Grid.h"


Grid::Grid(size_t nx, size_t ny,size_t g_)
    :ghost_cells(g_),
    num_xcells(nx),
    num_ycells(ny),
    

    x(nx),
    y(ny),

    eta(nx,ny,g_),
    grad_eta(nx,ny,g_),
 
    
    U(nx,ny,g_),
    Prim(nx,ny,g_),
    uBarL(nx,ny,g_),
    uBarR(nx,ny,g_),
    primL(nx,ny,g_),
    primR(nx,ny,g_),
    
    B(nx,ny,g_),
    E(nx,ny,g_),
    LaplacianB(nx,ny,g_),
    J(nx,ny,g_)
   

    {};

void Grid::reset_intermediate_arrays(){
    //only acts on the keynamed arrays
 
    for(size_t i=0;i<num_xcells+2*ghost_cells;i++){
        for(size_t j=0;j<num_ycells+2*ghost_cells;j++){
        uBarL(i,j)*=0;
        uBarR(i,j)*=0;
        primL(i,j)*=0;
        primR(i,j)*=0; 
}}
}

