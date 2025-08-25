#include "SLIC.h"
#include "StateVector.h"
#include "MHD_flux.h"
#include "BoundaryConditions.h"

using PSV=PrimitiveStateVector;
using CSV=ConservedStateVector;

void set_w_x(Grid& grid){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
            //B_n and psi are temporarily replaced with their characteristic variables for slope limiting process
            double w_plus=grid.Prim(i,j).psi()+grid.c_h*grid.Prim(i,j).B().x();
            double w_minus=grid.Prim(i,j).psi()-grid.c_h*grid.Prim(i,j).B().x();
            grid.Prim(i,j).B().x()=w_plus;
            grid.Prim(i,j).psi()=w_minus;
        }}}

void set_w_y(Grid& grid){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
            //B_n and psi are temporarily replaced with their characteristic variables for slope limiting process
            double w_plus=grid.Prim(i,j).psi()+grid.c_h*grid.Prim(i,j).B().y();
            double w_minus=grid.Prim(i,j).psi()-grid.c_h*grid.Prim(i,j).B().y();
            grid.Prim(i,j).B().y()=w_plus;
            grid.Prim(i,j).psi()=w_minus;
        }}}


StateVector get_Delta(const PSV& prim_plus, const PSV& prim_minus){
    StateVector output = 0.5*(prim_plus - prim_minus);      
    return output; }


StateVector get_r(const PSV& prim0,const PSV& prim_plus, const PSV& prim_minus,const double TOL){
    //Slope limiting using primitive variables 
    StateVector output;
    for(size_t k=0;k<9;k++){
            double denom=(prim_plus[k]- prim0[k]);
            if(std::fabs(denom)<TOL){denom=((denom>=0) ? TOL  : -TOL);};
        
            output[k]=(prim0[k]- prim_minus[k])/denom;      }
    return output;
}


void do_VL_limiting(StateVector& r){
    //Do Van Leer limiting in place
    for(size_t k=0;k<9;k++){
            
        if(r[k]<=0){r[k]=0;} 
        else{r[k]=std::min(2*r[k]/(1+r[k]),2/(1+r[k]));}
}}


void calc_ubar(bool x,const Grid& grid, const SimulationConfig& cfg,CSV& uBarL, CSV& uBarR,PSV& primL, PSV& primR,const PSV& Prim,const StateVector& Delta, const StateVector& Chi){
    //updates uBar L&R
    primL= Prim-0.5*Chi*Delta;
    primR= Prim+0.5*Chi*Delta;

    
    //Convert back from characteristic variables to B_n and psi
    //________________________________________________________________
    size_t index = (x ? 5 : 6);
    double B_nL= (primL[index]-primL[8])/(2*grid.c_h);
    double B_nR= (primR[index]-primR[8])/(2*grid.c_h);
    double psi_L= (primL[index]+primL[8])/2.;
    double psi_R= (primR[index]+primR[8])/2.;

    primL[index]=B_nL;
    primR[index]=B_nR;

    primL.psi()=psi_L;
    primR.psi()=psi_R;
    //________________________________________________________________
            
    uBarL=primL.prim_to_con(cfg.gamma);
    uBarR=primR.prim_to_con(cfg.gamma);
}


void calc_ubar_plus(bool x,const Grid& grid, CSV& uBarL, CSV& uBarR,PSV& primL, PSV& primR){
    StateVector fluxL;
    StateVector fluxR;
    if(x){
        MHD_xflux(fluxL,uBarL,primL,grid.c_h);
        MHD_xflux(fluxR,uBarR,primR,grid.c_h);
    }else{
        MHD_yflux(fluxL,uBarL,primL,grid.c_h);
        MHD_yflux(fluxR,uBarR,primR,grid.c_h);
    }
    
    uBarL-=0.5*grid.dt/grid.dx*(fluxR-fluxL);
    uBarR-=0.5*grid.dt/grid.dx*(fluxR-fluxL);
}


void do_SLIC_xupdate(Grid& grid, const SimulationConfig& cfg){
    //Performs entire SLIC process starting from U and ending with 
    //putting Ubarplus arrays in place
    set_w_x(grid);

    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
            StateVector Delta=get_Delta(grid.Prim(i+1,j),grid.Prim(i-1,j));
            StateVector slopes=get_r(grid.Prim(i,j),grid.Prim(i+1,j),grid.Prim(i-1,j));
            do_VL_limiting(slopes);
            calc_ubar(true,grid,cfg,grid.uBarL(i,j),grid.uBarR(i,j),grid.primL(i,j),grid.primR(i,j),grid.Prim(i,j),Delta,slopes);
            calc_ubar_plus(true,grid,grid.uBarL(i,j),grid.uBarR(i,j),grid.primL(i,j),grid.primR(i,j));
        }}
}
void do_SLIC_yupdate(Grid& grid, const SimulationConfig& cfg){
    //Performs entire SLIC process starting from U and ending with 
    //putting Ubarplus arrays in place
    set_w_y(grid);

    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
            StateVector Delta=get_Delta(grid.Prim(i,j+1),grid.Prim(i,j-1));
            StateVector slopes=get_r(grid.Prim(i,j),grid.Prim(i,j+1),grid.Prim(i,j-1));
            do_VL_limiting(slopes);
            calc_ubar(false,grid,cfg,grid.uBarL(i,j),grid.uBarR(i,j),grid.primL(i,j),grid.primR(i,j),grid.Prim(i,j),Delta,slopes);
            calc_ubar_plus(false,grid,grid.uBarL(i,j),grid.uBarR(i,j),grid.primL(i,j),grid.primR(i,j));
        }}
}

//Ubar vals are up to date, prim are not, boundaries haven't been updated