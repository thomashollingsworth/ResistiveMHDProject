#ifndef HLLD_H
#define HLLD_H

#include "Grid.h"
#include "SimulationConfig.h"
#include "c_f.h"

using PSV=PrimitiveStateVector;
using CSV=ConservedStateVector;
//used for deciding what flux each point should take
enum class Region{
    L,
    L_star,
    L_doublestar,
    R_doublestar,
    R_star,
    R
};

//Update the uBarLR and PrimLR arrays in place
void redefine_LRx(Grid& grid);
void redefine_LRy(Grid& grid);


void set_xtilde_vals(Grid& grid);
void set_ytilde_vals(Grid& grid);


//Single cell level operations (storing doubles is cheap and easy)

std::tuple<double,double> calc_S_LR_x(const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,const double gamma);
//(use calc_cfx for primL and primR)
std::tuple<double,double> calc_S_LR_y(const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,const double gamma);
//(use calc_cfx for primL and primR)

double calc_S_M_x(const double S_L,const double S_R,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR);
double calc_S_M_y(const double S_L,const double S_R,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR);

std::tuple<double,double> calc_S_stars_x(const double S_L,const double S_R,const double S_M,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR);
//Requires calc of density star states 
std::tuple<double,double> calc_S_stars_y(const double S_L,const double S_R,const double S_M,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR);


void calc_u_starLR_x(bool L,ConservedStateVector& output,const double S_M,const double S_L,const double S_R,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,const ConservedStateVector& uL,const ConservedStateVector& uR);
void calc_u_starLR_y(bool L,ConservedStateVector& output,const double S_M,const double S_L,const double S_R,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,const ConservedStateVector& uL,const ConservedStateVector& uR);
//Full in place star state calc, can be either L or R state
//Must do L or R for star regions, L and R for double star regions

void calc_u_doublestarLR_x(bool L,ConservedStateVector& output,ConservedStateVector& u_starL,ConservedStateVector& u_starR);
void calc_u_doublestarLR_y(bool L,ConservedStateVector& output,ConservedStateVector& u_starL,ConservedStateVector& u_starR);

void calc_HLLD_flux_x(const Grid& grid,StateVector& flux_out, const PSV& pL,const PSV& pR,const CSV& uL,const CSV& uR,const double gamma);
void calc_HLLD_flux_y(const Grid& grid,StateVector& flux_out,const PSV& pL,const PSV& pR,const CSV& uL,const CSV& uR,const double gamma);


void do_HLLD_x_update(Grid& grid, const SimulationConfig& cfg);
void do_HLLD_y_update(Grid& grid, const SimulationConfig& cfg);

void UpdatePrim(Grid& grid, const SimulationConfig& cfg);





#endif