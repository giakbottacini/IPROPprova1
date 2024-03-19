#include "Build_Geometry.hpp"

// This method compute the points (Point instances identified by an unique tag) that wrap the x,y coordinates of the airfoil profile, it stores them in a vector

std::vector<Point> Build_Geometry::compute_profile() const{

    //you want to retrive the coordinates of the airfoil, the front part of the airfoil in in (0,0)
    // create vectors of double to store the coordinates, one for the x coordinates, the other two for the upper
    // and the lower profile of the NACA respectively 
    std::vector<double> x_coord;
    std::vector<double> up_coord;
    std::vector<double> low_coord;

    // reserve memory for the vectors
    x_coord.reserve(this->my_data.NACA_points);
    up_coord.reserve(this->my_data.NACA_points);
    low_coord.reserve(this->my_data.NACA_points);

    
    // we compute the x coordinates supposing an equal spacing between the points, the interval is (first,second). 
    // first is zero and second is equal to the length of the chord
    double first = 0;                            
    double second = this->my_data.chord_length;               
    
    int n_sub_intervals = this->my_data.NACA_points - 1;   //number of sub intervals in the chord

    double h = (second - first)/(n_sub_intervals);  //length of the sub interval

    for(size_t i = 0 ; i < this->my_data.NACA_points; ++i){

        x_coord[i] = i*h; //just this since we start from zero

    }


    // we compute now the the coordinates of both the profiles
    double t = this->my_data.last_two_digit; //in the struct we save the last two digit, but we need the ratio of the thickness
    t = t/100;


    //lambda function that specify the profile of the NACA airfoil
    auto F_y = [t](double x) -> double {

        return 5*t*(0.2969*std::sqrt(x) - 0.1260*x - 0.3516*std::pow(x,2) + 0.2843*std::pow(x,3) - 0.1015*std::pow(x,4));

    };

    
    //we now compute the y coordinates of the upper and lower profile
    for(size_t i = 0; i< this->my_data.NACA_points ; ++i){

        up_coord[i] = F_y( x_coord[i]);
        low_coord[i] = -1*up_coord[i];

    }

    //finally we create a vector of Points to return, firstly we add the points of the upper profile and then the points of the lower part
    //NB: in order to have less problem with the creation of the lines, we insert the lower points starting from the last one in vector "low_coord"

    std::vector<Point> Points;
    Points.reserve(this->my_data.NACA_points -2);

    for(size_t i = 0; i<this->my_data.NACA_points; ++i){

        Point temp( x_coord[i] ,  up_coord[i] ,0.0, this->my_data.mesh_ref_1);
        Points.push_back(temp);

    }

    for(size_t i = this->my_data.NACA_points -2; i>0; --i){

        Point temp( x_coord[i] , low_coord[i] ,0.0, this->my_data.mesh_ref_1);
        Points.push_back(temp);

    }

    return Points;

}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method compute the 5 points that are used to define a circular arc in gmsh enviroment

std::vector<Point> Build_Geometry::compute_emitter() const{

    //we just need 5 points in order to create the emitter geometry
    std::vector<Point> Points;
    Points.reserve(5);
    
    //retrive some useful data
    double r = this->my_data.radius_emitter;
    double d = this->my_data.distance_emitter_collector;
    double ref = this->my_data.mesh_ref_1;

    double c_x = -d-r; // x coordinate of the center of the circunference
    double c_y = 0.0; //y coordinate of the center of the circunference

    Point p1(c_x,c_y,0.0,ref); //center Point
    Point p2(-d,c_y,0.0,ref);
    Point p3(c_x,r,0.0,ref);
    Point p4(c_x-r,c_y,0.0,ref);
    Point p5(c_x,c_y-r,0.0,ref);
    

    //then load the Points
    Points.push_back(p1);
    Points.push_back(p2);
    Points.push_back(p3);
    Points.push_back(p4);
    Points.push_back(p5);

    return Points;

}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------
//this method compute the position of the rectangular domain where the simulation takes place
std::vector<Point> Build_Geometry::compute_domain()  const{
    
   double ref = this->my_data.mesh_ref_1;

   //half of the heigh of the rectangular domain
   double r = this->my_data.radius_emitter;
   double dist_up = this->my_data.distance_emitter_up_bottom;
   double H = r + dist_up;

   //x position of the inlet and the outlet
   double dist_e_c = this->my_data.distance_emitter_collector;
   double dist_e_i = this->my_data.distance_emitter_inlet;
   double chord = this->my_data.chord_length ;
   double dist_Tedge =  this->my_data.distance_Tedge_outlet;

   double inlet = -dist_e_c -2*r - dist_e_i;
   double outlet = chord + dist_Tedge;

   //now we define the four points and we store them in a vector
   Point p1(inlet,H,0.0,ref);
   Point p2(inlet,-H,0.0,ref);
   Point p3(outlet,-H,0.0,ref);
   Point p4(outlet,H,0.0,ref);

   std::vector<Point> Points;
   Points.reserve(4);
   Points.push_back(p1);
   Points.push_back(p2);
   Points.push_back(p3);
   Points.push_back(p4);

   return Points;


}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
//this method simply write in the outfile the head of the .geo file

void Build_Geometry::write_head(std::ofstream & ofs) const{

   ofs << "// ===========================================" <<std::endl;
   ofs << "// ==================================MESH FILE" <<std::endl;
   ofs << "// ===========================================" <<std::endl;

   return;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
//this method writes in the output file the Points and the lines that compose the airfoil profile

void Build_Geometry::write_profile(std::ofstream & ofs) const{
    
    //first we retrive the points that we need
    std::vector<Point> airfoil_points = this->compute_profile();
    

    //we start with a comment in order to have a more readable .geo file
    ofs << std::endl;
    ofs << "//AIRFOIL GEOMETRY"<<std::endl;
    

    //we then write all the points with the following sintax
    for(size_t i = 0; i<airfoil_points.size(); ++i){
        
        ofs << "Point("<<airfoil_points[i].get_tag()<<") = {"<<airfoil_points[i].get_x()<<", "<<airfoil_points[i].get_y()<<", "<<airfoil_points[i].get_z()<<", "<<airfoil_points[i].get_local_mesh_ref()<<"};"<<std::endl;

    }
    
    //we pass to the lines
    ofs << std::endl;
    ofs << "//AIRFOIL CURVE"<<std::endl;

    //we write in the output file the four SPline that define the profile of the airfoil
    int h = this->my_data.NACA_points;
    h = h/2;


    ofs << "Spline(1) = {"<<airfoil_points[0].get_tag()<<":"<<airfoil_points[h-1].get_tag()<<"};"<<std::endl;
    ofs << "Spline(2) = {"<<airfoil_points[h-1].get_tag()<<":"<<airfoil_points[2*h-1].get_tag()<<"};"<<std::endl;
    ofs << "Spline(3) = {"<<airfoil_points[2*h-1].get_tag()<<":"<<airfoil_points[3*h-1].get_tag()<<"};"<<std::endl;
    ofs << "Spline(4) = {"<<airfoil_points[3*h-1].get_tag()<<":"<<airfoil_points[4*h -3].get_tag()<<", "<<airfoil_points[0].get_tag()<<"};"<<std::endl;


    return;
  
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
//this method writes in the output file the Points and the lines that compose the emitter geometry

void Build_Geometry::write_emitter(std::ofstream & ofs) const{
    
    //we start by retriving the points
    std::vector<Point> emitter_points = this->compute_emitter();
    

    //comment for the readability of the .geo file
    ofs <<std::endl;
    ofs << "//EMITTER GEOMETRY"<<std::endl;
    
    //write the Points in the output file
    for(size_t i = 0; i<emitter_points.size(); ++i){
        
        ofs << "Point("<<emitter_points[i].get_tag()<<") = {"<<emitter_points[i].get_x()<<", "<<emitter_points[i].get_y()<<", "<<emitter_points[i].get_z()<<", "<<emitter_points[i].get_local_mesh_ref()<<"};"<<std::endl;

    }

    //now we write the istructions to build a circunference exploiting gmsh commands
    ofs << std::endl;
    ofs << "//CIRCULAR ARCS"<<std::endl;

    ofs << "Circle(5) = {" << emitter_points[1].get_tag() <<"   , "<<emitter_points[0].get_tag()<<"   , "<<emitter_points[2].get_tag()<<"};"<<std::endl;
    ofs << "Circle(6) = {" << emitter_points[2].get_tag() <<"   , "<<emitter_points[0].get_tag()<<"   , "<<emitter_points[3].get_tag()<<"};"<<std::endl;
    ofs << "Circle(7) = {" << emitter_points[3].get_tag() <<"   , "<<emitter_points[0].get_tag()<<"   , "<<emitter_points[4].get_tag()<<"};"<<std::endl;
    ofs << "Circle(8) = {" << emitter_points[4].get_tag() <<"   , "<<emitter_points[0].get_tag()<<"   , "<<emitter_points[1].get_tag()<<"};"<<std::endl;

    return;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------
//this method writes in the output file the points and the lines that build the rectangular domain
void Build_Geometry::write_domain(std::ofstream & ofs) const{

    // first we retrive the domain points
    std::vector<Point> domain_points = this->compute_domain();


    //comment for the .geo file
    ofs <<std::endl;
    ofs << "//RECTANGULAR DOMAIN"<<std::endl;

    //we write all the points
    for(size_t i = 0; i<domain_points.size(); ++i){
        
        ofs << "Point("<<domain_points[i].get_tag()<<") = {"<<domain_points[i].get_x()<<", "<<domain_points[i].get_y()<<", "<<domain_points[i].get_z()<<", "<<domain_points[i].get_local_mesh_ref()<<"};"<<std::endl;

    }
    
    //and then the lines
    ofs << std::endl;
    ofs << "//DOMAIN LINES"<<std::endl;

    ofs << "Line(9) = {"<<domain_points[0].get_tag()<<", "<<domain_points[1].get_tag()<<"};"<<std::endl;
    ofs << "Line(10) = {"<<domain_points[1].get_tag()<<", "<<domain_points[2].get_tag()<<"};"<<std::endl;
    ofs << "Line(11) = {"<<domain_points[2].get_tag()<<", "<<domain_points[3].get_tag()<<"};"<<std::endl;
    ofs << "Line(12) = {"<<domain_points[3].get_tag()<<", "<<domain_points[0].get_tag()<<"};"<<std::endl;
 
    return;

}
//---------------------------------------------------------------------------------------------------------------------------------------------------------
//this method write the loops
void Build_Geometry::write_loops(std::ofstream & ofs) const{

    ofs << std::endl;
    ofs << "//LOOPS"<<std::endl;

    ofs <<"Curve Loop(1) = {9, 10, 11, 12};"<<std::endl; //curve loop of the domain counter-clokwise

    ofs <<"Curve Loop(2) = {5, 6, 7, 8};"<<std::endl; //curve loop of the circular emitter counter-clockwise

    ofs <<"Curve Loop(3) = {-1, -4, -3, -2};"<<std::endl; //curve loop of the airfoil counter-clockwise

    return;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method writes the surface that will be meshed
void Build_Geometry::write_surface(std::ofstream & ofs) const{

    ofs << std::endl;
    ofs << "//SURFACES"<<std::endl;

    ofs << "Plane Surface(1) = {1, 2, 3};"<<std::endl; //surface that will be meshed

    return;

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//This method define all the physical quantities that will be used by gmsh to understand the domain
void Build_Geometry::write_physical_groups(std::ofstream & ofs) const{

    ofs << std::endl;
    ofs << "//PHYSICAL GROUPS"<<std::endl;
    ofs << std::endl;
    ofs << "//PHYSICAL LINES"<<std::endl;

    //ofs << "Physical Line(""FARFIELD"") = {10, 12};"<<std::endl;             // these are the lines that define the farfield, namely up and bottom endges of the rectangular domain
    //ofs << "Physical Line(""INLET"") = {9};"<<std::endl;                      // this is the edge of the inlet
    //ofs << "Physical Line(""OUTLET"") = {11};"<<std::endl;                    // this is the edge of the outlet
    //ofs << "Physical Line(""AIRFOIL"") = {-1, -4, -3, -2};"<<std::endl;      // these are the lines that define the airfoil profile
    //ofs << "Physical Line(""EMITTER"") = {5, 6, 7, 8};"<<std::endl;           // these are the curves that define th emitter

    ofs << std::endl;
    ofs << "//PHYSICAL SURFACE"<<std::endl;

    ofs << "Physical Surface(1) = {1};"<<std::endl;                           //this is the physical surface, it takes the only surface that we have

    return;

}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// andranno parametrizzati i parametri, see gmsh sito per documentazione
void Build_Geometry::write_boundary_layer(std::ofstream & ofs) const{

    ofs << std::endl;
    ofs << "//BOUNDARY LAYER"<<std::endl;

    ofs << "Field[1]=BoundaryLayer;"<<std::endl;
    ofs << "Field[1].CurvesList={1,2,3,4};"<<std::endl;          //Tags of curves in the geometric model for which a boundary layer is needed
    ofs << "Field[1].Quads=1;"<<std::endl;                       //Generate recombined elements in the boundary layer
    ofs << "Field[1].Ratio=1.2;"<<std::endl;                     //Size ratio between two successive layers
    ofs << "Field[1].Size= 0.002336;"<<std::endl;                //Mesh size normal to the curve
    ofs << "Field[1].Thickness=0.06;"<<std::endl;                //Maximal thickness of the boundary layer
    ofs << "Field[1].FanPointsList={100};"<<std::endl;           //Tags of points in the geometric model for which a fan is created
    ofs << "Field[1].FanPointsSizesList={11};"<<std::endl;       //Number of elements in the fan for each fan point. If not present default value Mesh.BoundaryLayerFanElements
    ofs << "BoundaryLayer Field = 1;"<<std::endl;

    return;

}