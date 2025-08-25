#ifndef SLIC_H
#define SLIC_H

#include "Grid.h"
#include "SimulationConfig.h"

using PSV=PrimitiveStateVector;
using CSV=ConservedStateVector;

//Each step has its own dedicated function and then a single function collects all 
void set_w_x(Grid& grid);
void set_w_y(Grid& grid);

StateVector get_Delta(const PSV& prim_plus, const PSV& prim_minus);
StateVector get_r(const PSV& prim0,const PSV& prim_plus, const PSV& prim_minus,const double TOL=1e-12);

void do_VL_limiting(StateVector& r);

void calc_ubar(bool x,const Grid& grid, const SimulationConfig& cfg,CSV& uBarL, CSV& uBarR,PSV& primL, PSV& primR,const PSV& Prim,const StateVector& Delta, const StateVector& Chi);

void calc_ubar_plus(bool x,const Grid& grid, CSV& uBarL, CSV& uBarR,PSV& primL, PSV& primR);

void do_SLIC_xupdate(Grid& grid, const SimulationConfig& cfg);

void do_SLIC_yupdate(Grid& grid, const SimulationConfig& cfg);







#endif