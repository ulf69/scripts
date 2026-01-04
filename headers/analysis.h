#ifndef	analysis
#define analysis

#include <iostream>
#include <string>
#include <fstream>
#include <sstream> 
#include <utility>
#include <vector>

//bool skip_lines(char first_letter_of_line);
void timeframe_min_observable(const std::string& in_file);
void timeframe_below_value(const std::string& in_file, double threshold_value);
std::pair<double, double>  get_range_histogramm(const std::string& filename);
void make_histogramm(const std::string& in_file);
//  get_range_histogramm(const std::string& filename);
std::vector<double> sort_values(const std::string& filename);
double get_coordination_number_from_rdf(const std::string& filename, const std::string& cutoff);

#endif
