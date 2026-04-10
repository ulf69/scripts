#ifndef math
#define math

#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <utility>
#include <algorithm>

#include "datahandle.cpp"
#include "abbrevations.cpp"

class stat 
{

public:
	std::vector<std::vector<double>> data;
	stat(const std::string& filename)
	  :data(dh::read_multi_col_file(filename))
	{}	
	
	 std::vector<double> add_cols(){
		int siz = data.size();
		int siz2 = data[0].size();
		std::vector<double> vec;
		double sum;
		for(int j =0;j<siz2;j++){
		 	 for(int i =1; i<siz;i++){
				sum +=data[i][j];
		 	 }	
	 	 	vec.push_back(sum);
	 	 	sum=0;
	 	 	
	 	 	}
		return vec;
	 }

std::pair<std::vector<double>, std::vector<double>> make_bins(double startx, int data_points){

		std::vector<double> vecy;
		std::vector<double> vecx;
		int siz =data[0].size(); 
		int i =0;
		double sum;
		while(data[0][i] < startx){
			vecy.push_back(data[1][i]);
			vecx.push_back(data[0][i]);
			i++;
		}
		int bin_size = (siz -i )/data_points;
		for(int j = i+1; j<siz;j++){
			sum+=data[1][j];
			if((j-i)%bin_size == 0){
				vecy.push_back(sum/bin_size);
				sum=0;
				vecx.push_back(data[0][j]/2-data[0][j-bin_size]/2+data[0][j-bin_size]);
			}
		}
		std::vector<double> vec_out=scale_to_null(vecy);
	
		return std::make_pair(vecx, vec_out);
	 }
		
	 

 std::vector<double>	  scale_to_null(std::vector<double> vec){
	  	int siz=vec.size();
	 	int start=siz-10;
	 	
	 	double sum;
	 	
	 	for(int i = start;i<siz-1;i++){
	 		sum += vec[i];
	 	}
	 	sum = sum/9;
		std::vector<double> vec_out;
	 	for(int i = 0;i<siz;i++){
	 		vec_out.push_back(vec[i]-sum);
	 	}
	 	return vec_out;
	 }

	 static std::vector<double> add_vec(std::vector<double> vec1, std::vector<double> vec2){
		int  vector_size = vec1.size();
		int vector2_size = vec2.size();

		if(vector_size-vector2_size!=0){
			std::cout<<"we have got a problem cap"<<'\n';
			std::cout<<"vectors have differnt sizes"<<'\n';	
		}
		std::vector<double> vec;
	 	for(int i = 0 ; i < vector_size; i++){
			vec.push_back(vec1[i]+vec2[i]);
	 	}
	 	return vec;
	 }


	
	static double average(const std::vector<double> vec){
		int n = vec.size();
		double sum =0;
		for(int i =0; i < n; i++){
			sum+=vec[i];
		}
			return sum/n;

	}


	static double sd(const std::vector<double> vec){
		
		double av = average(vec);
		double st =0;
		int siz = vec.size();
		
		for(int i = 0; i < siz; i++){
			st += (vec[i]- av)*(vec[i]- av);
		}
		
		return std::sqrt(st/siz);
	}
	static double sem(const std::vector<double> vec ){
		double err = sd(vec);
		return err/sqrt(vec.size());
	}

	static double block_av(const std::vector<double> vec, int blocks){
		int block_size = vec.size()/blocks;
		std::vector<double> block(block_size);
		std::vector<double> block_sd;
		
		for(int b =0; b< blocks; b++){
			for(int i = 0 ;i < block_size; i++){
				block[i]=vec[i+b*block_size];
			}
			block_sd.push_back(sd(block));
		}

		
		return average(block_sd)/sqrt(block_sd.size());	
	}

	static double integrator(std::vector<double> x, std::vector<double> y){
		double sum = 0;
		double binwidth;
		int siz = x.size();
		for(int i = 0;i<siz;i++){
			if(i == siz){
				binwidth = x[siz]-x[siz-1];
			}
			else{		
				binwidth = x[i+1]-x[i];
			}
				sum +=y[i]*binwidth;	
			}
		return sum;
	}
	
	static double sum( std::vector<double> y){
			double sum = 0;
			
			int siz = y.size();
			for(int i = 0;i<siz;i++){
					sum +=y[i];	
				}
			return sum;
		}

	static 
std::pair<std::vector<double>, std::vector<double>> function_composition(std::vector<double> fy, std::vector<double> gy ){
		int vector_size  = gy.size();
		int vector1_size = fy.size();		
		if(vector_size-vector1_size!=0){
			std::cout<<"we have got a problem cap"<<'\n';
			std::cout<<"vectors have differnt sizes"<<'\n';	
		}
		// Create indices 0..n-1
		std::vector<int> idx(vector_size);
		for (int i = 0; i< vector_size; i++) idx[i]=i;

		//Sort indices according to gy values
		std::sort(idx.begin(), idx.end(), [&gy](int i1, int i2){return gy[i1]<gy[i2];});
		
		std::vector<double> sorted_gy(vector_size);
		std::vector<double> permuted_fy(vector_size);
		for (int i = 0; i < vector_size; ++i) {
		       sorted_gy[i]   = gy[idx[i]];
		       permuted_fy[i] = fy[idx[i]];
		}	
    return std::make_pair(sorted_gy, permuted_fy);


		
}


static double  overlap_hist(std::string hist1, std::string hist2){

	auto [x_hist1,y_hist1] 	  = dh::get_x_y(hist1);
	auto [x_hist2,y_hist2] 	  = dh::get_x_y(hist2);	
	int starting_position  	  = dh::get_common_min_val(x_hist1,x_hist2);
	double normalization_area = stat::sum(y_hist1)+stat::sum(y_hist2);
	double overlap_area    	  = 0;
	if(x_hist1.empty() || x_hist2.empty() || y_hist1.empty() || y_hist1.empty()){
		abb::print("-1000");
		return -1000;
	}
	else if(x_hist1[0]>x_hist2[0]){
		abb::print("-1000");
		return -1000;
	}
	for(int i =starting_position; i<x_hist1.size(); i++){
		
		if(y_hist2[i-starting_position]<y_hist1[i]){
			overlap_area+=y_hist2[i-starting_position];
		}
		else{
			overlap_area+=y_hist1[i];
		}
	}
	
	abb::print(static_cast<int>(2*overlap_area/(normalization_area)*100));
	return 2*overlap_area/(normalization_area);
}

static std::vector<double> vec_sub(std::vector<double> vec1, std::vector<double> vec2){
	int vector_size=vec1.size();
	int vector1_size=vec2.size();
	std::vector<double> vec_out;
	if(vector_size-vector1_size!=0){
				std::cout<<"we have got a problem cap"<<'\n';
				std::cout<<"vectors have differnt sizes"<<'\n';	
			}
			for(int i = 0; i < vector_size;i++){
				vec_out.push_back(vec1[i]-vec2[i]);
			}
		return vec_out;
}


static std::pair<std::vector<double>, std::vector<double>> sort_func(std::vector<double> vecx,std::vector<double> vecy ){
		int sizx = vecx.size();		
		int sizy = vecy.size();		
	if(sizx-sizy!=0){
					std::cout<<"we have got a problem cap"<<'\n';
					std::cout<<"vectors have differnt sizes"<<'\n';	
				}	
		
		// Create indices 0..n-1
		std::vector<int> idx(sizx);
		for (int i = 0; i< sizx; i++) idx[i]=i;

		//Sort indices according to gy values
		std::sort(idx.begin(), idx.end(), [&vecx](int i1, int i2){return vecx[i1]<vecx[i2];});
		
		std::vector<double> sorted_x(sizx);
		std::vector<double> permuted_y(sizx);
		for (int i = 0; i < sizx; ++i) {
		       sorted_x[i]   = vecx[idx[i]];
		       permuted_y[i] = vecy[idx[i]];
		}	
    return std::make_pair(sorted_x, permuted_y);
}

};
#endif
