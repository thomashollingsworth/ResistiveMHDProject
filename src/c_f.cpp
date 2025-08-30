#include "c_f.h"
#include <cmath>



void calc_cf_x(double gamma, double& output, const PSV& p){
    double factor= (gamma*p.pressure()+dot(p.B(),p.B()))/p.density();
    output=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*p.pressure()*p.B().x()*p.B().x()/(p.density()*p.density()))));
}


void calc_cf_y(double gamma, double& output, const PSV& p){
    double factor= (gamma*p.pressure()+dot(p.B(),p.B()))/p.density();
    output=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*p.pressure()*p.B().y()*p.B().y()/(p.density()*p.density()))));
}


void calc_cf_z(double gamma, double& output, const PSV& p){
    double factor= (gamma*p.pressure()+dot(p.B(),p.B()))/p.density();
    output=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*p.pressure()*p.B().z()*p.B().z()/(p.density()*p.density()))));
}


void set_c_h(Grid& grid, const SimulationConfig& cfg,const Array2D<PSV>& prim_array){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;

    double c_fx;
    double c_fy;
    double c_fz;
    double max_ch=0;
  
    for(size_t i=g;i<nx+g;i++){
        for(size_t j=g;j<ny+g;j++){
            calc_cf_x(cfg.gamma,c_fx,prim_array(i,j));
            calc_cf_y(cfg.gamma,c_fy,prim_array(i,j));
            calc_cf_z(cfg.gamma,c_fz,prim_array(i,j));
        double new_max=std::max({
            std::fabs(prim_array(i,j).velocity().x())+c_fx,
            std::fabs(prim_array(i,j).velocity().y())+c_fy,
            std::fabs(prim_array(i,j).velocity().z())+c_fz});

            if(new_max> max_ch){
            max_ch=new_max;
        } }}
  grid.c_h= max_ch;
}

