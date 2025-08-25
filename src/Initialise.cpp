#include "Initialise.h"
#include "BoundaryConditions.h"
#include <iostream>
#include <cmath>

void set_xy(Grid& grid,const double x0,const double xf,const double y0,const double yf,const size_t nx,const size_t ny){
    //Initialising the x and y arrays
        grid.dx=(xf-x0)/nx;
        for(size_t i=0;i<nx;i++){
            grid.x[i]=x0+(0.5+i)*grid.dx;}
        
        grid.dy=(yf-y0)/nx;
        for(size_t i=0;i<ny;i++){
            grid.y[i]=y0+(0.5+i)*grid.dy;}
}

void initialise(Grid& grid, SimulationConfig& cfg){
    using test=SimulationConfig::InitialCondition;

    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;

    //Declaring parameters
    double x0;
    double xf;
    double y0;
    double yf;
    double tf;
    double gamma;


    //constants
    double pi = 4*std::atan(1);

    //Assigns boundary conditions and intialises U,Prim,x,y arrays 
    switch(cfg.initialcondition){
    //______________________________________________________________________
    case test::BrioWu_x:{
        //Boundary Conditions
        cfg.bcs_x0=SimulationConfig::BoundaryCondition::Transmissive;
        cfg.bcs_xf=SimulationConfig::BoundaryCondition::Transmissive;
        cfg.bcs_y0=SimulationConfig::BoundaryCondition::Transmissive;
        cfg.bcs_yf=SimulationConfig::BoundaryCondition::Transmissive;
        //Domain conditions
         x0=0.;
         xf=800.;
         y0=0.;
         yf=800.;
         tf=80.;
        
        //Other parameters
        gamma=2.;

        //Setting values
        cfg.t_end=tf;
        cfg.gamma=gamma;
        
        //Initialising the x and y arrays
        set_xy(grid,x0,xf,y0,yf,nx,ny);

        //Move from large x to small x with field
        for(size_t i=g; i<nx+g;i++){
        for(size_t j=g;j<ny+g;j++){
        if(grid.x[i-g]<=400 && grid.y[j-g]<=400){
            //bottom left corner
            grid.Prim(i,j)[0]=0.125;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=0.1;
            grid.Prim(i,j)[5]=0.75;
            grid.Prim(i,j)[6]=-1.;
            grid.Prim(i,j)[7]=0.0;
            grid.Prim(i,j)[8]=0.;
        }
        else if(grid.x[i-g]<=400 && grid.y[j-g]>400){
            //top left corner
            grid.Prim(i,j)[0]=0.125;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=0.1;
            grid.Prim(i,j)[5]=0.75;
            grid.Prim(i,j)[6]=-1.;
            grid.Prim(i,j)[7]=0.0;
            grid.Prim(i,j)[8]=0.;
        }
        else if(grid.x[i-g]>400 && grid.y[j-g]<400){
            //bottom right corner
            grid.Prim(i,j)[0]=1.;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=1.;
            grid.Prim(i,j)[5]=0.75;
            grid.Prim(i,j)[6]=1.;
            grid.Prim(i,j)[7]=0.;
            grid.Prim(i,j)[8]=0.;
            
        }else if(grid.x[i-g]>400 && grid.y[j-g]>400){
            //top right corner
            grid.Prim(i,j)[0]=1.;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=1.;
            grid.Prim(i,j)[5]=0.75;
            grid.Prim(i,j)[6]=1.;
            grid.Prim(i,j)[7]=0.;
            grid.Prim(i,j)[8]=0.;
        }
        grid.U(i,j)=grid.Prim(i,j).prim_to_con(cfg.gamma);}}
        update_bcs(grid,cfg,grid.U);
        update_bcs(grid,cfg,grid.Prim);
        break;}
    //_____________________________________________________________________________________
    
    case test::BrioWu_y:{
        //Boundary Conditions
        cfg.bcs_x0=SimulationConfig::BoundaryCondition::Transmissive;
        cfg.bcs_xf=SimulationConfig::BoundaryCondition::Transmissive;
        cfg.bcs_y0=SimulationConfig::BoundaryCondition::Transmissive;
        cfg.bcs_yf=SimulationConfig::BoundaryCondition::Transmissive;
        //Domain conditions
         x0=0.;
         xf=800.;
         y0=0.;
         yf=800.;
         tf=80.;
        
        //Other parameters
         gamma=2.;

        //Setting values
        cfg.t_end=tf;
        cfg.gamma=gamma;
        
        //Initialising the x and y arrays
        set_xy(grid,x0,xf,y0,yf,nx,ny);

        //Move from large x to small x with field
        for(size_t i=g; i<nx+g;i++){
        for(size_t j=g;j<ny+g;j++){
        //Move from small y to large y with B field
        if(grid.x[i-g]<=400 && grid.y[j-g]<=400){
            //bottom left corner
            grid.Prim(i,j)[0]=1.;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=1.;
            grid.Prim(i,j)[6]=0.75;
            grid.Prim(i,j)[5]=1.;
            grid.Prim(i,j)[7]=0.;
            grid.Prim(i,j)[8]=0.;
        }
        else if(grid.x[i-g]<=400 && grid.y[j-g]>400){
            //top left corner
            grid.Prim(i,j)[0]=0.125;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=0.1;
            grid.Prim(i,j)[6]=0.75;
            grid.Prim(i,j)[5]=-1.;
            grid.Prim(i,j)[7]=0.0;
            grid.Prim(i,j)[8]=0.;
        }
        else if(grid.x[i-g]>400 && grid.y[j-g]<400){
            //bottom right corner
            grid.Prim(i,j)[0]=1.;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=1.;
            grid.Prim(i,j)[6]=0.75;
            grid.Prim(i,j)[5]=1.;
            grid.Prim(i,j)[7]=0.;
            grid.Prim(i,j)[8]=0.;
        }else if(grid.x[i-g]>400 && grid.y[j-g]>400){
            //top right corner
            grid.Prim(i,j)[0]=0.125;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=0.1;
            grid.Prim(i,j)[6]=0.75;
            grid.Prim(i,j)[5]=-1.;
            grid.Prim(i,j)[7]=0.0;
            grid.Prim(i,j)[8]=0.;
        }
        grid.U(i,j)=grid.Prim(i,j).prim_to_con(cfg.gamma);}}
        update_bcs(grid,cfg,grid.U);
        update_bcs(grid,cfg,grid.Prim);
        break;}
    //_____________________________________________________________________________________
    
    case test::OrzsagTang:{
        //Boundary Conditions
        cfg.bcs_x0=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_xf=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_y0=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_yf=SimulationConfig::BoundaryCondition::Periodic;
        //Domain conditions
         x0=0.;
         xf=1.;
         y0=0.;
         yf=1.;
         tf=1.1;
        
        //Other parameters
         gamma=5./3.;

        //Setting values
        cfg.t_end=tf;
        cfg.gamma=gamma;
        
        //Initialising the x and y arrays
        set_xy(grid,x0,xf,y0,yf,nx,ny);

        //Move from large x to small x with field
        
        for(size_t i=g; i<nx+g;i++){
        for(size_t j=g;j<ny+g;j++){
            grid.Prim(i,j)[0]=cfg.gamma*cfg.gamma;
            grid.Prim(i,j)[1]=-std::sin(2*pi*grid.y[j-g]);
            grid.Prim(i,j)[2]=std::sin(2*pi*grid.x[i-g]);
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=cfg.gamma;
            grid.Prim(i,j)[5]=-std::sin(2*pi*grid.y[j-g]);
            grid.Prim(i,j)[6]=std::sin(4*pi*grid.x[i-g]);
            grid.Prim(i,j)[7]=0.;
            grid.Prim(i,j)[8]=0.;
   
        grid.U(i,j)=grid.Prim(i,j).prim_to_con(cfg.gamma);}}
        
        update_bcs(grid,cfg,grid.U);
        update_bcs(grid,cfg,grid.Prim);
        break;}
    //_____________________________________________________________________________________
    
    case test::KelvinHelmholtz:{
        //Boundary Conditions
        cfg.bcs_x0=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_xf=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_y0=SimulationConfig::BoundaryCondition::Reflective;
        cfg.bcs_yf=SimulationConfig::BoundaryCondition::Reflective;
        //Domain conditions
         x0=0.;
         xf=1.;
         y0=-1.;
         yf=1.;
         tf=20.;
        
        //Other parameters
         gamma=5./3.;

        //Setting values
        cfg.t_end=tf;
        cfg.gamma=gamma;
        
        //Initialising the x and y arrays
        set_xy(grid,x0,xf,y0,yf,nx,ny);

        //Move from large x to small x with field
         pi = 4*std::atan(1);
        
        for(size_t i=g; i<nx+g;i++){
        for(size_t j=g;j<ny+g;j++){
            grid.Prim(i,j)[0]=cfg.gamma*cfg.gamma;
            grid.Prim(i,j)[1]=-std::sin(2*pi*grid.y[j-g]);
            grid.Prim(i,j)[2]=std::sin(2*pi*grid.x[i-g]);
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=cfg.gamma;
            grid.Prim(i,j)[5]=-std::sin(2*pi*grid.y[j-g]);
            grid.Prim(i,j)[6]=std::sin(4*pi*grid.x[i-g]);
            grid.Prim(i,j)[7]=0;
            grid.Prim(i,j)[8]=0;
        grid.U(i,j)=grid.Prim(i,j).prim_to_con(cfg.gamma);}}
        
        update_bcs(grid,cfg,grid.U);
        update_bcs(grid,cfg,grid.Prim);
        break;}
    //______________________________________________________________________
    
    case test::DoubleCurrentSheet:{
        //Boundary Conditions
        cfg.bcs_x0=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_xf=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_y0=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_yf=SimulationConfig::BoundaryCondition::Periodic;

        //Initial Parameters:
        double B_0=1.;
        double u_0=0.1;
        double eta_0=1e-5;
    

        gamma=5./3.; 
        
        //Domain conditions
         x0=-0.5;
         xf=0.5;
         y0=-0.5;
         yf=0.5;
         tf=10.;

        //Setting values
        cfg.t_end=tf;
        cfg.gamma=gamma;
        
        //Setting resistivity
        cfg.UseResistivity=true;
        cfg.resistivity_profile= std::make_unique<SimulationConfig::Uniform>(eta_0);

        //Initialising the x and y arrays
        set_xy(grid,x0,xf,y0,yf,nx,ny);

        for(size_t i=g; i<nx+g;i++){
        for(size_t j=g;j<ny+g;j++){
            grid.Prim(i,j)[0]=1.;
            grid.Prim(i,j)[1]=u_0*std::sin(2*pi*grid.y[j-g]);
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=0.05;
            grid.Prim(i,j)[5]=0.;
            grid.Prim(i,j)[7]=0.;
            grid.Prim(i,j)[8]=0.;
        if(std::fabs(grid.x[i-g])>0.25){grid.Prim(i,j)[6]=B_0;}
        else{grid.Prim(i,j)[6]=-B_0;}
        
        grid.U(i,j)=grid.Prim(i,j).prim_to_con(cfg.gamma);}}
        
        update_bcs(grid,cfg,grid.U);
        update_bcs(grid,cfg,grid.Prim);
        break;}
    
    //______________________________________________________________________
    
    case test::GEM:{
        //Boundary Conditions
        cfg.bcs_x0=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_xf=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_y0=SimulationConfig::BoundaryCondition::Reflective;
        cfg.bcs_yf=SimulationConfig::BoundaryCondition::Reflective;

        //Initial Parameters:
        double pi = 4*std::atan(1);
        double ion_inertial_length=1.;
        double L_y=12.8*ion_inertial_length;
        double L_x=2*L_y;
        
        double rho_0=1.;
        double rho_inf=0.2;
        double lambda=0.5;
        double B_0=1.;
        double psi_0=0.1;

        gamma=5./3.; 
        
        //Domain conditions
         x0=-L_x/2;
         xf=L_x/2;
         y0=-L_y/2;
         yf=L_y/2;
         tf=200.;

        //Setting values
        cfg.t_end=tf;
        cfg.gamma=gamma;
        
        //Setting resistivity
        cfg.UseResistivity=true;
        cfg.resistivity_profile= std::make_unique<SimulationConfig::Uniform>(0.001);

        //Initialising the x and y arrays
        set_xy(grid,x0,xf,y0,yf,nx,ny);

        for(size_t i=g; i<nx+g;i++){
        for(size_t j=g;j<ny+g;j++){
            grid.Prim(i,j)[0]=rho_0+(1/std::cosh(grid.y[j-g]/lambda))*(1/std::cosh(grid.y[j-g]/lambda))+ rho_inf;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=B_0*B_0/(2*rho_0)*grid.Prim(i,j).density();
            grid.Prim(i,j)[5]=B_0*std::tanh(grid.y[j-g]/lambda)+ pi/L_y*psi_0*std::cos(2*pi*grid.x[i-g]/L_x)*std::sin(pi*grid.y[j-g]/L_y);
            grid.Prim(i,j)[6]=2*pi/L_x*psi_0*std::sin(2*pi*grid.x[i-g]/L_x)*std::cos(pi*grid.y[j-g]/L_y);
            grid.Prim(i,j)[7]=0;
            grid.Prim(i,j)[8]=0;
    
        grid.U(i,j)=grid.Prim(i,j).prim_to_con(cfg.gamma);}}
        
        update_bcs(grid,cfg,grid.U);
        update_bcs(grid,cfg,grid.Prim);
        break;}
    //______________________________________________________________________
    
    case test::SingleCurrentSheet:{
        //Boundary Conditions
        cfg.bcs_x0=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_xf=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_y0=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_yf=SimulationConfig::BoundaryCondition::Periodic;

        //Initial Parameters:
        double B_0=1.;
        double eta_0=1e-5;

        gamma=5./3.; 
        
        //Domain conditions
         x0=-0.5;
         xf=0.5;
         y0=-0.5;
         yf=0.5;
         tf=5.0;

        //Setting values
        cfg.t_end=tf;
        cfg.gamma=gamma;
        
        //Setting resistivity
        cfg.UseResistivity=true;
        cfg.resistivity_profile= std::make_unique<SimulationConfig::Uniform>(eta_0);

        //Initialising the x and y arrays
        set_xy(grid,x0,xf,y0,yf,nx,ny);

        for(size_t i=g; i<nx+g;i++){
        for(size_t j=g;j<ny+g;j++){
            grid.Prim(i,j)[0]=1.;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=0.1;
            grid.Prim(i,j)[5]=0.;
            grid.Prim(i,j)[7]=0.;
            grid.Prim(i,j)[8]=0.;
        if(grid.x[i-g]<=0){grid.Prim(i,j)[6]=B_0;}
        else{grid.Prim(i,j)[6]=-B_0;}

        grid.U(i,j)=grid.Prim(i,j).prim_to_con(cfg.gamma);}}
        
        update_bcs(grid,cfg,grid.U);
        update_bcs(grid,cfg,grid.Prim);
        break;}

}

    //Setting up resistivity and grad_resistivity arrays
    if(cfg.UseResistivity){
    
    double min_dX2 = std::min(grid.dx*grid.dx, grid.dy*grid.dy);
    double max_eta=0;
    
    for(size_t i=g; i<nx+g;i++){
        for(size_t j=g;j<ny+g;j++){
            double new_eta=(*cfg.resistivity_profile)(i-g,j-g);//Virtual functions can give overhead but this is only called for intialisation
            grid.eta(i,j)=new_eta;
            if(std::fabs(new_eta)> max_eta){
                max_eta=std::fabs(new_eta);
        }}
    update_bcs(grid,cfg,grid.eta);
    grid.B_timestep= min_dX2/(4*max_eta);

    for(size_t i=1; i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
            grid.grad_eta(i,j).x()=1./(2*grid.dx)*(grid.eta(i+1,j)-grid.eta(i-1,j));
            grid.grad_eta(i,j).y()=1./(2*grid.dy)*(grid.eta(i,j+1)-grid.eta(i,j-1));
            grid.grad_eta(i,j).z()=0.;
        }}
}
    }
    grid.reset_intermediate_arrays();
}