#ifndef SIMULATIONCONFIG_H
#define SIMULATIONCONFIG_H

#include <string>
//Structure that is used to choose the initial conditions, boundary conditions etc for the simulation

struct SimulationConfig{

    SimulationConfig(){};
    
    enum class InitialCondition{
        BrioWu_x,
        BrioWu_y,
        OrzsagTang,
        KelvinHelmholtz,
        GEM,
        DoubleCurrentSheet,
        SingleCurrentSheet,
    };

    enum class BoundaryCondition{
        Transmissive,
        Periodic,
        Reflective,
        Conducting
    };

    //Test parameter members initialised with default values

    InitialCondition initialcondition=InitialCondition::OrzsagTang;
    double t_end=1.0;
    double C=0.8;
    double gamma=5./3.;
    BoundaryCondition bcs_x0=BoundaryCondition::Periodic;
    BoundaryCondition bcs_xf=BoundaryCondition::Periodic;
    BoundaryCondition bcs_y0=BoundaryCondition::Periodic;
    BoundaryCondition bcs_yf=BoundaryCondition::Periodic;

    //Saving Parameters set to a default
    bool save=true;
    int save_interval=10;
    bool plot=true;
    int plot_interval=20;

    std::string save_directory="Results";

    bool UseResistivity=true; //controls whether to bother updating the resistive source terms

    //Runtime polymorphism structure for Resistivity
    struct ResProfile{
        virtual ~ResProfile() = default;
        enum class Type {
            Uniform,
        };
        Type profile_type;
        // Return the resistivity at a given point (x,y) 
        virtual double operator()(double x, double y) const = 0;
    };
    struct Uniform : ResProfile{
        Type profile_type= Type::Uniform;
        double const_value=1e-4;
        Uniform(double v) : const_value(v) {};
        double operator()(double /*x*/, double /*y*/) const override{
            return const_value; //Each x,y point is set to const_value
        }
        

    };
    //Config has a member 'resistivity_profile' which is a smart pointer to a specific ResProfile structure, defaults to a uniform profile
    std::unique_ptr<ResProfile> resistivity_profile = std::make_unique<Uniform>(1e-4);
    
    
    
  
    
    
    
    
    

    
    


};











#endif