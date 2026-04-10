#ifndef filewriter 
#define filewriter 

#include <vector>
#include <string>
#include <fstream>


namespace fw{

void write_y_x_file(std::vector<double> y, std::vector<double> x, std::string filename){
    std::ofstream outfile(filename, std::ios::app);	
	int vec_size=y.size();
	for(int i = 0; i < vec_size; i++){
		outfile<<x[i]<<'\t'<<y[i]<<'\n';
	}
    outfile.close();
}






void append_y_x(double fx, double x, std::string filename){
    std::ofstream outfile(filename, std::ios::app);	
	outfile<<x<<'\t'<<fx<<'\n';
    outfile.close();
}
void append_string(std::string str1,std::string str2,std::string filename){
    std::ofstream outfile(filename, std::ios::app);	
	outfile<<str1<<'\t'<<str2<<'\n';
    outfile.close();
}
void append_y_x_errors(double fx, double dfx, double x, double dx,  std::string filename){
		std::ofstream outfile(filename, std::ios::app);	
		outfile<<x<<'\t'<<dx<<'\t'<<fx<<'\t'<<dfx<<'\n';
	    outfile.close();
}
void append_y_x_dy(double fx, double dfx, double x,   std::string filename){
		if (std::isnan(fx)){return;}
		std::ofstream outfile(filename, std::ios::app);	
		outfile<<x<<'\t'<<fx<<'\t'<<dfx<<'\n';
	    outfile.close();
}


void write_subvectors(std::vector<double> sorted_gy, std::vector<double> permuted_fy, std::string out){

	std::vector<double> subvector_gy;
	std::vector<double> subvector_fy;
	int vector_size = sorted_gy.size();
	for (int i = 1; i < vector_size;i++){
			if(i==vector_size-1){
				std::string filename = out+std::to_string(subvector_gy[0])+".xvg";
				write_y_x_file( subvector_fy, subvector_gy, filename);	
				}
			else if (sorted_gy[i] == sorted_gy[i-1]){
				subvector_gy.push_back(sorted_gy[i]);
				subvector_fy.push_back(permuted_fy[i]);
			}	
			else{
			std::string filename = out+std::to_string(sorted_gy[i-1])+".xvg";
			write_y_x_file( subvector_fy, subvector_gy, filename);
			subvector_fy.clear();
			subvector_gy.clear();
			
			subvector_gy.push_back(sorted_gy[i]);
			subvector_fy.push_back(permuted_fy[i]);
			}
		}	
}

}
#endif
