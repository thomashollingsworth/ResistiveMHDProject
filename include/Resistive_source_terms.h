#ifndef RESISTIVESOURCETERM_H
#define RESISTIVESOURCETERM_H

#include "Array2D.h"
#include "VectorCalc.h"
#include "SimulationConfig.h"
#include "Grid.h"
#include "BoundaryConditions.h"


//Handling source terms that can't be solved analytically
//Want to use an explicit scheme and optionally an IMEX scheme
void do_ResistiveRK2_subcycle(Grid& grid,const SimulationConfig& cfg);
void do_RK_step(Grid& grid,double dt);
void update_J_buffers(Grid& grid);
void update_LaplacianB_buffers(Grid& grid);
void buffers_to_U(Grid& grid);
void assign_Resistive_buffers(Grid& grid);
double get_E_timestep(Grid& grid);

#endif