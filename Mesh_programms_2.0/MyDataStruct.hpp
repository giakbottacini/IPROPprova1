#ifndef MY_DATA_STRUCT_HPP
#define MY_DATA_STRUCT_HPP

#include <string>

struct MyDataStruct{               //this struct contains all the data of the problem

   std::string airfoil_type;                    //string that stores the name of the type of the airfoil  eg: "NACA0012"
   int last_two_digit;                          //int that stores the last two digit of the NACA type, it's used to compute the thickness of the airfoil
   double chord_length;                         //double that stores the length of the chord
   int NACA_points;                             //int that stores the number of the points that we want to use to draw half NACA profile (serve pari! vedi funzione profilo airfoil)
   double radius_emitter;                       //double that stores the radius length of the emitters
   double distance_emitter_collector;           //double that stores the distance between emitter(circle) and collector(airfoil)
   double distance_Tedge_outlet;                //double that represent the distance between the trailing edge and the end of the rectangular domain
   double distance_emitter_inlet;               //double that stores the distance between the emitter and the inlet of the domain
   double distance_emitter_up_bottom;           //double that stores the distance between the emitter and the bottom/up part of the domain
   double mesh_ref_1;                           //double that stores the mesh refinemet
};


//
//      |----------------------------------------------------------------------------------------|
//      |                                                                                        |
//      |                                                                                        |
//      |                                                                                        |
//      |                                                                                        |
//      |                                                                                        |
//      |                                                                                        |
//      |                                                                                        |
//      |                                                                                        |
//      |                                                                                        |
//      |----------------------------------------------------------------------------------------|
//
//


#endif //MY_DATA_STRUCT_HPP