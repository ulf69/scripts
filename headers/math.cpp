#ifndef math
#define math

#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <utility>
#include <algorithm>

#include "datahandle.cpp"

class stat 
{
private:
	std::vector<double> vecx;
	std::vector<double> vecy;
public:
	stat(const std::string& filename)
	{
		auto result=dh::get_x_y(filename);
		vecx = std::move(result.first);	
		vecy = std::move(result.second);	
	}
	const std::vector<double>& get_vecx() const { return vecx; }
	const std::vector<double>& get_vecy() const { return vecy; }

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





};
#endif
