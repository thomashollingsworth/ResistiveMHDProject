#include "Grid.h"
#include "SimulationConfig.h"
#include "Initialise.h"
#include "SavingRoutine.h"
#include "MHD_flux.h"
#include "c_f.h"
#include "TimeStep.h"
#include "SLIC.h"
#include "HLLD.h"
#include "psi_source_term.h"
#include <iostream>
#include "BoundaryConditions.h"
#include "Resistive_source_terms.h"



int main() {
//Grid will hold all 2D arrays and data for intermediate values
Grid grid(128,256);
//Simulation config is used to control parameters of the test
SimulationConfig cfg;

cfg.C=0.8;

cfg.save_directory="KelvinHelmholtzData";

cfg.initialcondition= SimulationConfig::InitialCondition::KelvinHelmholtz;

cfg.save_interval=10;
cfg.plot_interval=20;



//Initialise function sets up grid using the simulation config settings

initialise(grid,cfg);

double time=0;
int iteration=0;

do{
    
    if(iteration%cfg.save_interval==0){
        save_to_file(grid,cfg,iteration,time);
        std::cout<<"Iteration "<<iteration<< " saved"<<std::endl;
    }
    
    get_time_step(grid,cfg); //requires an updated Prim and U array
    time+=grid.dt;
    iteration+=1;

    do_half_psi_update(grid);
    update_bcs(grid,cfg,grid.U);//ALL U cells are now updated
    
    //Resistive term updates only update U
    //do_ResistiveRK2_subcycle(grid,cfg);
    UpdatePrim(grid,cfg);
    
    do_SLIC_xupdate(grid,cfg);//Requires up to date U and Prim

    do_HLLD_x_update(grid,cfg);//HLLD only updates U

    grid.reset_intermediate_arrays();
    
    UpdatePrim(grid,cfg);
    
    do_SLIC_yupdate(grid,cfg);//Requires up to date U and Prim

    do_HLLD_y_update(grid,cfg);//HLLD only updates U

    grid.reset_intermediate_arrays();

    do_half_psi_update(grid);
    //Resistive term updates only update U
    //do_ResistiveRK2_subcycle(grid,cfg);
    UpdatePrim(grid,cfg);

    std::cout<<"Iteration " <<iteration << " completed t= "<<time<<std::endl;
}while(time<cfg.t_end);

}
