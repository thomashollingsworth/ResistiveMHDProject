#include "HLLD.h"
#include <cmath>
#include "MHD_flux.h"
#include "BoundaryConditions.h"
#include <iostream>
#include <stdexcept> // for std::runtime_error


using PSV=PrimitiveStateVector;
using CSV=ConservedStateVector;



//L and R states in HLLD are defined differently to those in SLIC
void redefine_LRx(Grid& grid){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+g;i++){
        for(size_t j=1;j<ny+g;j++){
            auto tempL=grid.uBarR(i,j);
            auto tempR=grid.uBarL(i+1,j);

            grid.uBarL(i,j)=tempL;
            grid.uBarR(i,j)=tempR;
}}  
}

void redefine_LRy(Grid& grid){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+g;i++){
        for(size_t j=1;j<ny+g;j++){
            auto tempL=grid.uBarR(i,j);
            auto tempR=grid.uBarL(i,j+1);

            grid.uBarL(i,j)=tempL;
            grid.uBarR(i,j)=tempR;;
}}
}

void set_xtilde_vals(Grid& grid){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    //precompute constants
    const double factorB=0.5*1./grid.c_h;
    const double factorPsi=0.5*grid.c_h;
    
    for(size_t i=1;i<nx+g;i++){
        for(size_t j=1;j<ny+g;j++){
            double B_xtilde=0.5*(grid.uBarL(i,j).B().x()+grid.uBarR(i,j).B().x())-factorB*(grid.uBarR(i,j).psi()-grid.uBarL(i,j).psi());
            double psi_tilde=0.5*(grid.uBarR(i,j).psi()+grid.uBarL(i,j).psi())-factorPsi*(grid.uBarR(i,j).B().x()-grid.uBarL(i,j).B().x());
            
            grid.uBarR(i,j).B().x()=B_xtilde;
            grid.uBarL(i,j).B().x()=B_xtilde;
            grid.uBarR(i,j).psi()=psi_tilde;
            grid.uBarL(i,j).psi()=psi_tilde;
}}
}
void set_ytilde_vals(Grid& grid){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    //precompute constants
    const double factorB=0.5*1./grid.c_h;
    const double factorPsi=0.5*grid.c_h;

    for(size_t i=1;i<nx+g;i++){
        for(size_t j=1;j<ny+g;j++){
            double B_ytilde=0.5*(grid.uBarL(i,j).B().y()+grid.uBarR(i,j).B().y())-factorB*(grid.uBarR(i,j).psi()-grid.uBarL(i,j).psi());
            double psi_tilde=0.5*(grid.uBarR(i,j).psi()+grid.uBarL(i,j).psi())-factorPsi*(grid.uBarR(i,j).B().y()-grid.uBarL(i,j).B().y());
            
            grid.uBarR(i,j).B().y()=B_ytilde;
            grid.uBarL(i,j).B().y()=B_ytilde;
            grid.uBarR(i,j).psi()=psi_tilde;
            grid.uBarL(i,j).psi()=psi_tilde;
}}
}

//All following operations can be defined at single cell level to try and best utilise the register

std::tuple<double,double> calc_S_LR_x(const PSV& pL,const PSV& pR,const double gamma){
    double c_fL;
    double c_fR;
    calc_cf_x(gamma, c_fL, pL);
    calc_cf_x(gamma, c_fR, pR);
    double SL=std::min(pL.velocity().x()-c_fL,pR.velocity().x()-c_fR);
    double SR=std::max(pL.velocity().x()+c_fL,pR.velocity().x()+c_fR);
    return{SL,SR};
}

std::tuple<double,double> calc_S_LR_y(const PSV& pL,const PSV& pR,const double gamma){
    double c_fL;
    double c_fR;
    calc_cf_y(gamma, c_fL, pL);
    calc_cf_y(gamma, c_fR, pR);
    double SL=std::min(pL.velocity().y()-c_fL,pR.velocity().y()-c_fR);
    double SR=std::max(pL.velocity().y()+c_fL,pR.velocity().y()+c_fR);
    return{SL,SR};
}


double calc_S_M_x(const double S_L,const double S_R,const PSV& pL,const PSV& pR){
    double S_M;
    S_M= ((S_R-pR.velocity().x())*pR.density()*pR.velocity().x()-(S_L-pL.velocity().x())*pL.density()*pL.velocity().x() -pR.pressure_T()+pL.pressure_T())/((S_R-pR.velocity().x())*pR.density()-(S_L-pL.velocity().x())*pL.density());
    return S_M;
}

double calc_S_M_y(const double S_L,const double S_R,const PSV& pL,const PSV& pR){
    double S_M;
    S_M= ((S_R-pR.velocity().y())*pR.density()*pR.velocity().y()-(S_L-pL.velocity().y())*pL.density()*pL.velocity().y() -pR.pressure_T()+pL.pressure_T())/((S_R-pR.velocity().y())*pR.density()-(S_L-pL.velocity().y())*pL.density());
    return S_M;
}

std::tuple<double,double> calc_S_stars_x(const double S_L,const double S_R,const double S_M,const PSV& pL,const PSV& pR){
//Requires calc of density star states
double density_starL=pL.density()*(S_L-pL.velocity().x())/(S_L-S_M);
double density_starR=pR.density()*(S_R-pR.velocity().x())/(S_R-S_M);

double S_starL= S_M-std::fabs(pL.B().x())/std::sqrt(density_starL);
double S_starR= S_M+std::fabs(pR.B().x())/std::sqrt(density_starR);

return {S_starL,S_starR};
}

std::tuple<double,double> calc_S_stars_y(const double S_L,const double S_R,const double S_M,const PSV& pL,const PSV& pR){
//Requires calc of density star states
double density_starL=pL.density()*(S_L-pL.velocity().y())/(S_L-S_M);
double density_starR=pR.density()*(S_R-pR.velocity().y())/(S_R-S_M);

double S_starL= S_M-std::fabs(pL.B().y())/std::sqrt(density_starL);
double S_starR= S_M+std::fabs(pR.B().y())/std::sqrt(density_starR);

return {S_starL,S_starR};
}

void calc_u_starLR_x(bool L,CSV& output,const double S_M,const double S_L,const double S_R,const PSV& pL,const PSV& pR,const CSV& uL,const CSV& uR){
    double TOL=1e-8; //for detecting 0 denominator
    const PSV& p = (L ? pL : pR);
    const CSV& u = (L ? uL : uR);
    double S = (L ? S_L : S_R);

    output.density()=p.density()*(S-p.velocity().x())/(S-S_M);
    output.momentum().x()=output.density()*S_M;
    //Check denominators
    double vel_denom= (p.density()*(S-p.velocity().x())*(S-S_M)-p.B().x()*p.B().x());
    if (std::fabs(vel_denom)<TOL){
        output.momentum().y()=output.density()*p.velocity().y();
        output.momentum().z()=output.density()*p.velocity().z();
    }
    else{
    output.momentum().y()= output.density()*(p.velocity().y()-((S_M-p.velocity().x())*p.B().x()*p.B().y())/vel_denom);
    output.momentum().z()= output.density()*(p.velocity().z()-((S_M-p.velocity().x())*p.B().x()*p.B().z())/vel_denom);
    }
    output.B().x()=p.B().x();

    double B_denom=(p.density()*(S-p.velocity().x())*(S-S_M)-p.B().x()*p.B().x());
    if (std::fabs(B_denom)<TOL){
        output.B().y()=0.;
        output.B().z()=0.;
    }
    else{
    output.B().y()=p.B().y()*(p.density()*(S-p.velocity().x())*(S-p.velocity().x())-p.B().x()*p.B().x())/B_denom;
    output.B().z()=p.B().z()*(p.density()*(S-p.velocity().x())*(S-p.velocity().x())-p.B().x()*p.B().x())/B_denom;
    }
    
    double p_T_star=((S_R-pR.velocity().x())*pR.density()*pL.pressure_T()-(S_L-pL.velocity().x())*pL.density()*pR.pressure_T()+pL.density()*pR.density()*(S_R-pR.velocity().x())*(S_L-pL.velocity().x())*(pR.velocity().x()-pL.velocity().x()))/((S_R-pR.velocity().x())*pR.density()-(S_L-pL.velocity().x())*pL.density());



    output.energy()=((S-p.velocity().x())*u.energy()-p.pressure_T()*p.velocity().x()+p_T_star*S_M + p.B().x()*(dot(p.velocity(),p.B())-1/output.density()*dot(output.momentum(),output.B())))/(S-S_M);
    output.psi()=p.psi();

}

void calc_u_starLR_y(bool L,CSV& output,const double S_M,const double S_L,const double S_R,const PSV& pL,const PSV& pR,const CSV& uL,const CSV& uR){
    double TOL=1e-8; //for detecting 0 denominator
    const PSV& p = (L ? pL : pR);
    const CSV& u = (L ? uL : uR);
    double S = (L ? S_L : S_R);

    output.density()=p.density()*(S-p.velocity().y())/(S-S_M);
    output.momentum().y()=output.density()*S_M;
    
    double vel_denom= (p.density()*(S-p.velocity().y())*(S-S_M)-p.B().y()*p.B().y());
    if (std::fabs(vel_denom)<TOL){
        output.momentum().x()=output.density()*p.velocity().x();
        output.momentum().z()=output.density()*p.velocity().z();
    }
    else{
    output.momentum().x()= output.density()*(p.velocity().x()-((S_M-p.velocity().y())*p.B().y()*p.B().x())/vel_denom);
    output.momentum().z()= output.density()*(p.velocity().z()-((S_M-p.velocity().y())*p.B().y()*p.B().z())/vel_denom);
    }
    output.B().y()=p.B().y();

    double B_denom=(p.density()*(S-p.velocity().y())*(S-S_M)-p.B().y()*p.B().y());
    if (std::fabs(B_denom)<TOL){
        output.B().x()=0.;
        output.B().z()=0.;
    }
    else{
    output.B().x()=p.B().x()*(p.density()*(S-p.velocity().y())*(S-p.velocity().y())-p.B().y()*p.B().y())/B_denom;
    output.B().z()=p.B().z()*(p.density()*(S-p.velocity().y())*(S-p.velocity().y())-p.B().y()*p.B().y())/B_denom;
    }

    double p_T_star=((S_R-pR.velocity().y())*pR.density()*pL.pressure_T()-(S_L-pL.velocity().y())*pL.density()*pR.pressure_T()+pL.density()*pR.density()*(S_R-pR.velocity().y())*(S_L-pL.velocity().y())*(pR.velocity().y()-pL.velocity().y()))/((S_R-pR.velocity().y())*pR.density()-(S_L-pL.velocity().y())*pL.density());
 
    output.energy()=((S-p.velocity().y())*u.energy()-p.pressure_T()*p.velocity().y()+p_T_star*S_M + p.B().y()*(dot(p.velocity(),p.B())-1/output.density()*dot(output.momentum(),output.B())))/(S-S_M);
    output.psi()=p.psi();

}

void calc_u_doublestarLR_x(bool L,CSV& output,CSV& u_starL,CSV& u_starR){
    const CSV& u = (L ? u_starL : u_starR);
    int pm =(L? -1.:1.);

    output.density()=u.density();
    output.momentum().x()=u.momentum().x();
    
    int sgnBx=(0.0 <= u.B().x()) - (u.B().x() < 0.0);
    double sqrt_ustarL_density=std::sqrt(u_starL.density());
    double sqrt_ustarR_density=std::sqrt(u_starR.density());
    
    output.momentum().y()=output.density()*(u_starL.momentum().y()/sqrt_ustarL_density+u_starR.momentum().y()/sqrt_ustarR_density+(u_starR.B().y()-u_starL.B().y())*sgnBx)/(sqrt_ustarR_density+sqrt_ustarL_density);
    output.momentum().z()=output.density()*(u_starL.momentum().z()/sqrt_ustarL_density+u_starR.momentum().z()/sqrt_ustarR_density+(u_starR.B().z()-u_starL.B().z())*sgnBx)/(sqrt_ustarR_density+sqrt_ustarL_density);

    output.B().x()=u.B().x();

    output.B().y()=(sqrt_ustarL_density*u_starR.B().y()+sqrt_ustarR_density*u_starL.B().y()+sqrt_ustarL_density*sqrt_ustarR_density*(u_starR.momentum().y()/u_starR.density()-u_starL.momentum().y()/u_starL.density())*sgnBx)/(sqrt_ustarR_density+sqrt_ustarL_density);
    
    output.B().z()=(sqrt_ustarL_density*u_starR.B().z()+sqrt_ustarR_density*u_starL.B().z()+sqrt_ustarL_density*sqrt_ustarR_density*(u_starR.momentum().z()/u_starR.density()-u_starL.momentum().z()/u_starL.density())*sgnBx)/(sqrt_ustarR_density+sqrt_ustarL_density);

    output.energy()=u.energy()+pm*1/std::sqrt(u.density())*(dot(u.momentum(),u.B())-dot(output.momentum(),output.B()))*sgnBx;

    output.psi()=u.psi();
}

void calc_u_doublestarLR_y(bool L,CSV& output,CSV& u_starL,CSV& u_starR){
    const CSV& u = (L ? u_starL : u_starR);
    int pm =(L? -1:1);

    output.density()=u.density();
    output.momentum().y()=u.momentum().y();
    
    int sgnBy=(0.0 <= u.B().y()) - (u.B().y() < 0.0);
    double sqrt_ustarL_density=std::sqrt(u_starL.density());
    double sqrt_ustarR_density=std::sqrt(u_starR.density());
    
    output.momentum().x()=output.density()*(u_starL.momentum().x()/sqrt_ustarL_density+u_starR.momentum().x()/sqrt_ustarR_density+(u_starR.B().x()-u_starL.B().x())*sgnBy)/(sqrt_ustarR_density+sqrt_ustarL_density);
    output.momentum().z()=output.density()*(u_starL.momentum().z()/sqrt_ustarL_density+u_starR.momentum().z()/sqrt_ustarR_density+(u_starR.B().z()-u_starL.B().z())*sgnBy)/(sqrt_ustarR_density+sqrt_ustarL_density);

    output.B().y()=u.B().y();

    output.B().x()=(sqrt_ustarL_density*u_starR.B().x()+sqrt_ustarR_density*u_starL.B().x()+sqrt_ustarL_density*sqrt_ustarR_density*(u_starR.momentum().x()/u_starR.density()-u_starL.momentum().x()/u_starL.density())*sgnBy)/(sqrt_ustarR_density+sqrt_ustarL_density);
    output.B().z()=(sqrt_ustarL_density*u_starR.B().z()+sqrt_ustarR_density*u_starL.B().z()+sqrt_ustarL_density*sqrt_ustarR_density*(u_starR.momentum().z()/u_starR.density()-u_starL.momentum().z()/u_starL.density())*sgnBy)/(sqrt_ustarR_density+sqrt_ustarL_density);

    output.energy()=u.energy()+pm*1/std::sqrt(u.density())*(dot(u.momentum(),u.B())-dot(output.momentum(),output.B()))*sgnBy;

    output.psi()=u.psi();
}

void calc_HLLD_flux_x(const Grid& grid,StateVector& flux_out, const PSV& pL,const PSV& pR,const CSV& uL,const CSV& uR,const double gamma){
    Region region;
    //Temporary objects that can be handled in register or L1
    StateVector flux_L;
    StateVector flux_R;
    CSV u_starL;
    CSV u_starR;
    CSV u_doublestar;
    
    MHD_xflux(flux_L,uL,pL,grid.c_h);
    MHD_xflux(flux_R,uR,pR,grid.c_h);

    auto [S_L,S_R]= calc_S_LR_x(pL,pR,gamma);
    double S_M = 0.0, S_starL = 0.0, S_starR = 0.0;
    
    if(S_L>0){region=Region::L;}
    else if(S_R<0){region=Region::R;}
    else{
        S_M = calc_S_M_x(S_L,S_R,pL,pR);
        std::tie(S_starL,S_starR) = calc_S_stars_x(S_L,S_R,S_M,pL,pR);
        if(S_starL>0){
            region=Region::L_star;
            calc_u_starLR_x(true,u_starL,S_M,S_L,S_R,pL,pR,uL,uR);
        }else if(S_M>0){
            region=Region::L_doublestar;
            calc_u_starLR_x(true,u_starL,S_M,S_L,S_R,pL,pR,uL,uR);//left star states
            calc_u_starLR_x(false,u_starR,S_M,S_L,S_R,pL,pR,uL,uR);//right star states

            calc_u_doublestarLR_x(true,u_doublestar,u_starL,u_starR); //Calculate left version of double star state
        }else if(S_starR>0){
            region=Region::R_doublestar;
            calc_u_starLR_x(true,u_starL,S_M,S_L,S_R,pL,pR,uL,uR);//left star states
            calc_u_starLR_x(false,u_starR,S_M,S_L,S_R,pL,pR,uL,uR);//right star states

            calc_u_doublestarLR_x(false,u_doublestar,u_starL,u_starR);

        }else{
            region=Region::R_star;
            calc_u_starLR_x(false,u_starR,S_M,S_L,S_R,pL,pR,uL,uR);

        }
    }
    
    switch(region){
        case Region::L:flux_out=flux_L;                                                                break;
        case Region::L_star:flux_out=(flux_L+S_L*(u_starL-uL));                                        break;
        case Region::R_star:flux_out=(flux_R+S_R*(u_starR-uR));                                        break;
        case Region::L_doublestar:flux_out=flux_L+S_L*(u_starL-uL)+S_starL*(u_doublestar-u_starL);     break;
        case Region::R_doublestar:flux_out=(flux_R+S_R*(u_starR-uR)+S_starR*(u_doublestar-u_starR));   break;
        case Region::R:flux_out=flux_R;                                                                break;
    
    }  
}


void calc_HLLD_flux_y(const Grid& grid,StateVector& flux_out,const PSV& pL,const PSV& pR,const CSV& uL,const CSV& uR,const double gamma){
    Region region;
    StateVector flux_L;
    StateVector flux_R;
    CSV u_starL;
    CSV u_starR;
    CSV u_doublestar;

    MHD_yflux(flux_L,uL,pL,grid.c_h);
    MHD_yflux(flux_R,uR,pR,grid.c_h);

    auto [S_L,S_R]= calc_S_LR_y(pL,pR,gamma);
    double S_M = 0.0, S_starL = 0.0, S_starR = 0.0;
    
    if(S_L>0){region=Region::L;}
    else if(S_R<0){region=Region::R;}
    else{
        S_M = calc_S_M_y(S_L,S_R,pL,pR);
        std::tie(S_starL,S_starR) = calc_S_stars_y(S_L,S_R,S_M,pL,pR);
        
        if(S_starL>0){region=Region::L_star;
            
            calc_u_starLR_y(true,u_starL,S_M,S_L,S_R,pL,pR,uL,uR);
        }
        else if(S_M>0){region=Region::L_doublestar;
            
            calc_u_starLR_y(true,u_starL,S_M,S_L,S_R,pL,pR,uL,uR);//left star states
            calc_u_starLR_y(false,u_starR,S_M,S_L,S_R,pL,pR,uL,uR);//right star states
            calc_u_doublestarLR_y(true,u_doublestar,u_starL,u_starR); //Calculate left version of double star state
        }
        else if(S_starR>0){region=Region::R_doublestar;
            
            calc_u_starLR_y(true,u_starL,S_M,S_L,S_R,pL,pR,uL,uR);//left star states
            calc_u_starLR_y(false,u_starR,S_M,S_L,S_R,pL,pR,uL,uR);//right star states
            calc_u_doublestarLR_y(false,u_doublestar,u_starL,u_starR);
        }
        else{region=Region::R_star;
            calc_u_starLR_y(false,u_starR,S_M,S_L,S_R,pL,pR,uL,uR);
        }
    }
    switch(region){
        case Region::L:flux_out=flux_L;                                                                break;
        case Region::L_star:flux_out=(flux_L+S_L*(u_starL-uL));                                        break;
        case Region::R_star:flux_out=(flux_R+S_R*(u_starR-uR));                                        break;
        case Region::L_doublestar:flux_out=flux_L+S_L*(u_starL-uL)+S_starL*(u_doublestar-u_starL);     break;
        case Region::R_doublestar:flux_out=(flux_R+S_R*(u_starR-uR)+S_starR*(u_doublestar-u_starR));   break;
        case Region::R:flux_out=flux_R;                                                                break;
    
    }  
}


//Composite functions to perform the entire HLLD update

void do_HLLD_x_update(Grid& grid, const SimulationConfig& cfg){
    //Starting with ubarplus states after SLIC method
   
    redefine_LRx(grid);
    set_xtilde_vals(grid);

    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    StateVector flux_temp;
    for(size_t i=1;i<nx+2*g-2;i++){
        for(size_t j=1;j<ny+2*g-2;j++){
            //Fast cell level operations in register/L1 cache
            grid.uBarL(i,j).con_to_prim(grid.primL(i,j),cfg.gamma);
            grid.uBarR(i,j).con_to_prim(grid.primR(i,j),cfg.gamma);
            
            calc_HLLD_flux_x(grid,flux_temp,grid.primL(i,j),grid.primR(i,j),grid.uBarL(i,j),grid.uBarR(i,j),cfg.gamma);

            grid.U(i,j)-=grid.dt/grid.dx*flux_temp;
            grid.U(i+1,j)+=grid.dt/grid.dx*flux_temp;
        }}

update_bcs(grid,cfg,grid.U);
}


void do_HLLD_y_update(Grid& grid, const SimulationConfig& cfg){
    //Starting with ubarplus states after SLIC method
    redefine_LRy(grid);
    set_ytilde_vals(grid);
    
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    StateVector flux_temp;
    
    for(size_t i=1;i<nx+2*g-2;i++){
        for(size_t j=1;j<ny+2*g-2;j++){
            //Fast cell level operations in register/L1 cache
            grid.uBarL(i,j).con_to_prim(grid.primL(i,j),cfg.gamma);
            grid.uBarR(i,j).con_to_prim(grid.primR(i,j),cfg.gamma);
            calc_HLLD_flux_y(grid,flux_temp,grid.primL(i,j),grid.primR(i,j),grid.uBarL(i,j),grid.uBarR(i,j),cfg.gamma);
            
            grid.U(i,j)-=grid.dt/grid.dy*flux_temp;
            grid.U(i,j+1)+=grid.dt/grid.dy*flux_temp;

}}
update_bcs(grid,cfg,grid.U);
}

void UpdatePrim(Grid& grid, const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    
    for(size_t i=0;i<nx+2*g;i++){
        for(size_t j=0;j<ny+2*g;j++){
            grid.U(i,j).con_to_prim(grid.Prim(i,j),cfg.gamma);//Updates primitive data
}}
}
