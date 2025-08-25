#ifndef ARRAY2D_H
#define ARRAY2D_H
#include <vector>

//This structure will be a general template for a std::vector that represents a 2D array
//The Array2D can have any type of component e.g. double StateVector double*

template<typename Component>
struct Array2D{
   
    size_t nx;
    size_t ny;
    size_t nghost;
    std::vector<Component> data;
   
    //Constructor
    Array2D(size_t nx,size_t ny,size_t nghost=2)
    : nx(nx), ny(ny), nghost(nghost), data((nx+2*nghost)*(ny+2*nghost)) {}

    //set_default_vals, updates all components to their default constructor value
    void set_default_vals() {
    data.assign((nx + 2*nghost) * (ny + 2*nghost), Component{});
}

    
    
    //Access indices by 2D indexing
    Component& operator()(size_t i, size_t j){
    return data[(ny+2*nghost)*i+j];}
    
    const Component& operator()(size_t i, size_t j) const {
    return data[(ny+2*nghost)*i+j];}

};
//Add additional class that views a 2D scalar field by indexing a single parameter from a larger Array2D
//For component "K"; ScalarFieldofK(i,j) = Array2D(i,j)[K]
//Element of the Array2D must have a well defined [] operator, i.e. can be a StateVector, Vector3, std::array but NOT a double 

template<typename Array>
struct ScalarFieldView {
    size_t nx;
    size_t ny;
    size_t nghost;
    Array& data; //Takes a view/reference of an array and an index for the relevant component
    size_t index;

    // Constructor
    ScalarFieldView(Array& a, size_t idx)
        : data(a), index(idx),
          nx(a.nx), ny(a.ny), nghost(a.ghost_cells) {}

    //Access indices by 2D indexing
    auto& operator()(size_t i, size_t j){
    return data(i,j)[index];}
    
    const auto& operator()(size_t i, size_t j) const {
    return data(i,j)[index];}

};

template<typename Array>
struct Vector3FieldView{
    size_t nx;
    size_t ny;
    size_t nghost;
    Array& data; //Takes a view/reference of an array and an index for the relevant component
    size_t index_x;
    size_t index_y;
    size_t index_z;

    // Constructor
    Vector3FieldView(Array& a, size_t idx_x,size_t idx_y,size_t idx_z)
        : data(a), index_x(idx_x),index_y(idx_y),index_z(idx_z),
          nx(a.nx), ny(a.ny), nghost(a.ghost_cells) {}

    //Access indices by 2D indexing
    auto& operator()(size_t i, size_t j){
    return Vector3View(&this->data(i,j)[index_x],&this->data(i,j)[index_y],&this->data(i,j)[index_z]);}
   
    const auto& operator()(size_t i, size_t j) const {
    return Vector3ConstView(&this->data(i,j)[index_x],&this->data(i,j)[index_y],&this->data(i,j)[index_z]);}
};




#endif