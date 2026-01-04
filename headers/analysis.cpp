#ifndef analysis
#define analysis

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream> 
#include <utility> 
#include <cmath> 



bool skip_lines(char first_letter_of_line){
        return (first_letter_of_line=='#'||first_letter_of_line=='@');  //skip first lines

}


	



double integrator(std::string file_name, double startx, double endx){
double sum = 0;
auto [x,y] = dh::get_x_y(file_name);
double binwidth=x[1]-x[0];
for(int i = 0;i<x.size();i++){
	if(x[i]<startx){
		continue;
	}
	else if(x[i]>endx){
		break;
	}
	else{
	if(i>0){
		binwidth=x[i]-x[i-1];
	}
		sum +=y[i]*binwidth;	
		//dh::print(i);
	}
}
return sum;
	
}

std::pair<double, double> ads_keg(std::string file_name, double middle, int keg){
double normalization =	integrator(file_name, 0, 100);

double keg_layer_one =	integrator(file_name, 0, middle)/normalization*keg;
double keg_layer_two =	integrator(file_name, middle,100)/normalization*keg;
return std::make_pair(keg_layer_one,keg_layer_two);
}



std::tuple<double, double, double> ads_keg2(std::string file_name, double st_min, double nd_min, int keg){
double normalization =	integrator(file_name, 0, 100);

double keg_layer_one =	integrator(file_name, 0, st_min)/normalization*keg;
double keg_bulk =	integrator(file_name,  st_min,nd_min)/normalization*keg;
double keg_layer_two =	integrator(file_name, nd_min,100)/normalization*keg;
return {keg_layer_one,keg_bulk,keg_layer_two};
}



double get_coordination_number_from_rdf(const std::string& filename, const std::string& cutoff){


std::string line;
double radius=std::stod(cutoff);

double V_sphere;
double x;
double dx;
double coord_number = 0;
double bin_width    = 0.002;
std::ifstream observable(filename);

while(std::getline(observable, line)){
	if(skip_lines(line[0])) continue;
	int i =0;
	std::istringstream iss(line);
	std::string col;
	while(iss>>col){
		i++;

		if(i%2==1){
		if(std::stod(col)>=radius)
		{
			break;	
		}
		
			std::cout<<col<<'\n';
			x  = std::stod(col)+0.001;
			dx = x + bin_width;
			V_sphere = 4*M_PI*(dx*dx*dx-x*x*x);

		}	
		else{
			std::cout<<col<<'\n';
			coord_number+=V_sphere*std::stod(col);
		}
	}

}
return coord_number;
}
#endif
