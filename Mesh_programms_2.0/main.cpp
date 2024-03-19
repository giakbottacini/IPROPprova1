#include "json.hpp"
#include "Build_Geometry.hpp"

using json = nlohmann::json;


int main(int argc, char** argv){


//##################### RETRIVE DATA FROM .JSON FILE ###############################################################################################################



// Open a file stream for reading
std::ifstream inFile("Data.json");

// Check if the file stream is open
if (!inFile.is_open()) {
  std::cerr << "Failed to open the file for reading." << std::endl;
  return 1;
}

// Read JSON data from the file
json json_data;       //object type json
inFile >> json_data;

// Close the file stream
inFile.close();

// Access the data from the JSON object and store them in MyDataStruct object

MyDataStruct s_data;  //structured data

s_data.airfoil_type = json_data["airfoil_type"];       
s_data.last_two_digit = json_data["last_2_digit_NACA"];           
s_data.chord_length = json_data["chord_length"];            
s_data.NACA_points = json_data["NACA_points"];
s_data.radius_emitter = json_data["radius_emitter"];
s_data.distance_emitter_collector=json_data["distance_emitter_collector"];
s_data.distance_Tedge_outlet=json_data["distance_Tedge_outlet"];
s_data.distance_emitter_inlet=json_data["distance_emitter_inlet"];
s_data.distance_emitter_up_bottom=json_data["distance_emitter_up_bottom"];
s_data.mesh_ref_1=json_data["mesh_ref_1"];        



//##################### WRITE THE COORDINATES IN A .GEO OUTPUT FILE ##################################################################################

// Create an ofstream object to write to a file
std::ofstream outFile("example.geo");  


// Recall that if the "example.geo" file already exists, its content will be erased! (check if the file already exists)

// Check if the file is opened successfully
if (!outFile.is_open()) {
  std::cerr << "Error opening the file." << std::endl;
  return 1; // return an error code
}


// Create a Build_Geometry object

Build_Geometry my_geometry(s_data);

// Write data to the file exploiting the methods

my_geometry.write_head(outFile);
my_geometry.write_profile(outFile);
my_geometry.write_emitter(outFile);
my_geometry.write_domain(outFile);
my_geometry.write_loops(outFile);
my_geometry.write_surface(outFile);
my_geometry.write_physical_groups(outFile);
my_geometry.write_boundary_layer(outFile);


// Close the file
outFile.close();

std::cout << "Data has been written to the file successfully." << std::endl;

      
return 0;
}