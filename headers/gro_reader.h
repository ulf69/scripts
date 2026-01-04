
#ifndef gro_reader
#define gro_reader

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <tuple>
#include <iomanip>

std::vector<std::string> infile_to_vector(const std::string& in_gro);
std::tuple<double, double, double> get_minimum_coordinates(const std::vector<std::string>&vec);
void write_gro_file(const std::vector<std::string>& vec, std::string out_gro, int startmol, int particles, double yshift,double zshift);
void make_fitting_box(std::string in_gro);
void sort_molecules(std::string top_file, std::string gro_file, std::string out_file); 
#endif
