#ifndef VECTORCALC_H
#define VECTORCALC_H

#include <tuple>
#include "Vector3.h"

//generic templated functions that handle Array2D<double> or ScalarFieldViews of Array2Ds
//template functions need to be defined inside header file!


template<typename Vector>//Can be any derivative structure from Vector3Template
Vector3 grad2D(const Vector& u_xplus,const Vector& u_xminus,const Vector& u_yplus,const Vector& u_yminus, const double dx,const double dy){
    //Returns {grad_x(u),grad_y(u),grad_z(u)} 
    Vector3 output;
    output.x()=1/(2*dx)*(u_xplus-u_xminus);
    output.y()=1/(2*dy)*(u_yplus-u_yminus); 
    output.z()=0.;  
    return output;       
}

template<typename Vector>
double div2D(const Vector& F_xplus,const Vector& F_yplus,const Vector& F_xminus,const Vector& F_yminus, const double dx,const double dy){
    double output;
    output=1/(2*dx)*(F_xplus.x()-F_xminus.x());
    output+=1/(2*dy)*(F_yplus.y()-F_yminus.y());       
}
template<typename Vector>
Vector3 curl2D(const Vector& F_xplus,const Vector& F_yplus,const Vector& F_xminus,const Vector& F_yminus, const double dx,const double dy){
    Vector3 output;

    output.x()=1/(2*dy)*(F_yplus.z()-F_yminus.z());
    output.y()=-1/(2*dx)*(F_xplus.z()-F_xminus.z());
    output.z()=1/(2*dx)*(F_xplus.y()-F_xminus.y())-1/(2*dy)*(F_yplus.x()-F_yminus.x());
    return output;}

template<typename Array,typename Out>
Out Laplacian2D(const size_t i,const size_t j, const Array& u, const double dx,const double dy){
    //Returns 9-point laplacian of a scalar field u at the point i,j
    // u must be accessible by indices (i,j)
    Out output;
    output=1/6*(4*((u(i+1,j)-2*u(i,j)+u(i-1,j))/(dx*dx)+(u(i,j+1)-2*u(i,j)+u(i,j-1))/(dy*dy))+(u(i+1,j+1)+u(i+1,j-1)+u(i-1,j+1)+u(i-1,j-1)-4*u(i,j))/(dx*dy));
    return output; 
}
 

#endif