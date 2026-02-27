#ifndef abbrevations
#define abbrevations

#include <iostream>

namespace abb{
template <typename T>
void print(const T&var){
	std::cout<<var<<'\n';
}

template <typename T>
void printt(const T&var){
	std::cout<<var<<'\t';
}

void show(std::vector<std::string> array){
        for(int i =0;i<array.size();i++){
                std::cout<<array[i]<<'\n';
        }
}


void show_int(std::vector<int> array){
        for(int i =0;i<array.size();i++){
                std::cout<<array[i]<<'\n';
        }
}

void show_dou(std::vector<double> array){
        for(int i =0;i<array.size();i++){
                std::cout<<array[i]<<'\n';
        }
}
void show_xy(std::vector<double> array,std::vector<double> array2){
        for(int i =0;i<array.size();i++){
                printt(array[i]);
                print(array2[i]);
        }
}

	
} 
#endif
