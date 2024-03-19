#ifndef BUILD_GEOMETRY_HPP
#define BUILD_GEOMETRY_HPP

#include "MyDataStruct.hpp"
#include "Point.hpp"
#include <vector>
#include <iostream>
#include <fstream> 
#include <string>
#include <cmath>

class Build_Geometry{

   private:

     MyDataStruct my_data;

   public:
     
     //CONSTRUCTOR
     Build_Geometry(MyDataStruct d): my_data(d) {};

     //METHODS TO COMPUTE
     std::vector<Point> compute_profile() const;   //this method compute the points that made the airfoil profile
     std::vector<Point> compute_emitter() const;   //this method compute the points that made the emitter
     std::vector<Point> compute_domain()  const;   //this method compute the points taht made the rectangular domain
     
     //METHODS TO WRITE IN THE OUTPUT FILE
     void write_head(std::ofstream & ofs) const;               //this method writes in the output file the head of the .geo file
     void write_profile(std::ofstream & ofs) const;            //this method writes in the output file the points that compose the airfoil and the emitter
     void write_emitter(std::ofstream & ofs) const;            //this method writes in the output file the points that compose the emitter geometry
     void write_domain(std::ofstream & ofs) const;             //this method writes in the output file the points that compose the rectanguklar domain
     void write_loops(std::ofstream & ofs) const;              //this method writes in the output file the loops that define the airfoil, the emitter and the domain
     void write_surface(std::ofstream & ofs) const;            //this method writes in the output file the surface that will be meshed
     void write_physical_groups(std::ofstream & ofs) const;    //this method writes in the output file the physical groups that describe the regions of the domain
     void write_boundary_layer(std::ofstream & ofs) const;     //this method writes in the output thr boundary layer that wraps the airfoil

     
};


#endif //BUILD_GEOMETRY_HPP