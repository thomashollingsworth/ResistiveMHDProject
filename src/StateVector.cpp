#include "StateVector.h"



//Local aliases for shorthand

using CSV=ConservedStateVector;
using PSV=PrimitiveStateVector;
using SV=StateVector;


//Accessing scalar quantities by name
double& PSV::density()   {return this->data[0];}
double& CSV::density()   {return this->data[0];}

const double& PSV::density() const   {return this->data[0];}
const double& CSV::density() const    {return this->data[0];}

double& PSV::pressure()   {return this->data[4];}
double& CSV::energy()   {return this->data[4];}

const double& PSV::pressure() const   {return this->data[4];}
const double& CSV::energy() const   {return this->data[4];}

double& PSV::psi()   {return this->data[8];}
double& CSV::psi()   {return this->data[8];}

const double& PSV::psi() const   {return this->data[8];}
const double& CSV::psi() const   {return this->data[8];}


double PSV::pressure_T(){return this->data[4] + 0.5*dot(B(),B()); }
double PSV::pressure_T() const {return this->data[4] + 0.5*dot(B(),B()); }


//Accessing vector quantities by name
Vector3View PSV::B() {return Vector3View(&this->data[5],&this->data[6],&this->data[7]);}
Vector3View CSV::B() {return Vector3View(&this->data[5],&this->data[6],&this->data[7]);}
const Vector3ConstView PSV::B() const {return Vector3ConstView(&this->data[5],&this->data[6],&this->data[7]);}
const Vector3ConstView CSV::B() const {return Vector3ConstView(&this->data[5],&this->data[6],&this->data[7]);}

Vector3View PSV::velocity() {return Vector3View(&this->data[1],&this->data[2],&this->data[3]);}
Vector3View CSV::momentum() {return Vector3View(&this->data[1],&this->data[2],&this->data[3]);}
const Vector3ConstView PSV::velocity() const {return Vector3ConstView(&this->data[1],&this->data[2],&this->data[3]);}
const Vector3ConstView CSV::momentum() const {return Vector3ConstView(&this->data[1],&this->data[2],&this->data[3]);}




//Converting between primitive and conservative

CSV PSV::prim_to_con(const double gamma) const{
    CSV conserved(
        density(),
        density() * velocity().x(),
        density() * velocity().y(),
        density() * velocity().z(),
        pressure()/((gamma-1)) + 0.5* density() * dot(velocity(),velocity()) + 0.5* dot(B(),B()),
        B().x(),
        B().y(),
        B().z(),
        psi()
    );
    return conserved;
}
void PSV::prim_to_con(ConservedStateVector& out,double gamma) const{
        out[0]=density();
        out[1]=density() * velocity().x();
        out[2]=density() * velocity().y();
        out[3]=density() * velocity().z();
        out[4]=pressure()/((gamma-1)) + 0.5* density() * dot(velocity(),velocity()) + 0.5* dot(B(),B());
        out[5]=B().x();
        out[6]=B().y();
        out[7]=B().z();
        out[8]=psi();
}


PSV CSV::con_to_prim(const double gamma) const{
    PSV primitive(
        density(),
        momentum().x()/density(),
        momentum().y()/density(),
        momentum().z()/density(),
        (gamma-1)*(energy()-0.5*dot(momentum(),momentum())/density()-0.5*dot(B(),B())),
        B().x(),
        B().y(),
        B().z(),
        psi()
    );
    return primitive;
}
void CSV::con_to_prim(PrimitiveStateVector& out,double gamma) const{
        out[0]=density();
        out[1]=momentum().x()/density();
        out[2]=momentum().y()/density();
        out[3]=momentum().z()/density();
        out[4]=(gamma-1)*(energy()-0.5*dot(momentum(),momentum())/density()-0.5*dot(B(),B()));
        out[5]=B().x();
        out[6]=B().y();
        out[7]=B().z();
        out[8]=psi();
}




