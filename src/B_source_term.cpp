#include "Resistive_source_terms.h"
#include <algorithm>
#include <iostream>


double get_E_timestep(Grid& grid){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    
    double max_Bsquared=0;
    for(size_t i=g;i<nx+g;i++){
        for(size_t j=g;j<ny+g;j++){
            double B_squared =grid.B(i,j).magsquared();
            if(B_squared > max_Bsquared){max_Bsquared=B_squared;}
        }}
    return grid.B_timestep/std::sqrt(max_Bsquared);
}

//Updating buffer arrays for reused intermediates 

void assign_Resistive_buffers(Grid& grid){
    //Buffer Array2D<double> store Bx,By,Bz for better cache access
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=0;i<nx+2*g;i++){
        for(size_t j=0;j<ny+2*g;j++){
            grid.B(i,j)= Vector3(grid.U(i,j).B().x(),grid.U(i,j).B().y(),grid.U(i,j).B().z());
            grid.E(i,j)=grid.U(i,j).energy();
}}}

void buffers_to_U(Grid& grid){
    //Buffer Array2D<double> store Bx,By,Bz for better cache access
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=g;i<nx+g;i++){
        for(size_t j=0;j<ny+g;j++){
            grid.U(i,j).B().x()= grid.B(i,j).x();
            grid.U(i,j).B().y()= grid.B(i,j).y();
            grid.U(i,j).B().z() = grid.B(i,j).z();
            grid.U(i,j).energy() = grid.E(i,j);
}}}

void update_LaplacianB_buffers(Grid& grid){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;

for(size_t i=g;i<nx+g;i++){
        for(size_t j=g;j<ny+g;j++){
            grid.LaplacianB(i,j)=Laplacian2D<Array2D<Vector3>, Vector3>(i,j,grid.B,grid.dx,grid.dy);}}
    }

void update_J_buffers(Grid& grid){
    //Buffer Array2D<double> store Jx,Jy,Jz for better cache access
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
            grid.J(i,j)=curl2D(grid.B(i+1,j),grid.B(i,j+1),grid.B(i-1,j),grid.B(i,j-1),grid.dx,grid.dy);         
}}
}



void do_RK_step(Grid& grid,double dt){

    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    
    
    for(size_t i=g;i<nx+g;i++){
        for(size_t j=g;j<ny+g;j++){
            Vector3 JcrossGradEta =cross(grid.J(i,j),grid.grad_eta(i,j));
            //Update E first as it depends on B
            grid.E(i,j)+=dt*(dot(grid.B(i,j),JcrossGradEta) +grid.eta(i,j)*(dot(grid.B(i,j),grid.LaplacianB(i,j))+grid.J(i,j).magsquared()));
            //Now safe to update B array
            grid.B(i,j)+=dt*(JcrossGradEta+grid.eta(i,j)*grid.LaplacianB(i,j));
            
}}
}

void do_ResistiveRK2_subcycle(Grid& grid,const SimulationConfig& cfg){
    double subtime=0;
    size_t subcycle_count=0;
    
    double tf= grid.dt/2.; //Do a half step either side of main HLLD loop
    assign_Resistive_buffers(grid);

    while(subtime<tf){
        
        double E_time=get_E_timestep(grid);
        double timestep=std::min({E_time,grid.B_timestep,(tf-subtime)});
        subtime+=timestep;
        //Perform RK1 loop
        update_LaplacianB_buffers(grid);
        update_J_buffers(grid);
        do_RK_step(grid,timestep/2);
        
        // update_LaplacianB_buffers(grid);
        // update_J_buffers(grid);
        // do_RK_step(grid,timestep/2);
        subcycle_count+=1;
    }
    std::cout<<"Completed "<< subcycle_count << " subcycles" << std::endl;
    buffers_to_U(grid);
    update_bcs(grid,cfg, grid.U);

}