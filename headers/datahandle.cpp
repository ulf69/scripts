#ifndef datahandle 
#define datahandle 

#include <sstream>
#include <fstream>
#include "filewriter.cpp"
#include <vector>
#include <iostream>
namespace dh{

std::vector<std::string> iss(std::string str){
	std::istringstream iss(str);
	std::string col;
	std::vector<std::string> vec;
	while(iss>>col){
		vec.push_back(col);
	}
	return vec;
}



std::vector<std::string> whole_file_to_vector(const std::string& in_gro){
           std::ifstream infile(in_gro);           //opens inputfile
           std::vector<std::string> written_file;                  //vector containing whole file
           std::string line;
                   while (std::getline(infile, line)){ //fill vector line per line
                           written_file.push_back(line);
                   }
                   infile.close();
                   return written_file;
}

std::vector<std::string> infile_to_vector(const std::string& in_gro){
                std::ifstream infile(in_gro);           //opens inputfile
                std::vector<std::string> written_file;                  //vector containing whole file
                std::string line;

                        while (std::getline(infile, line)){ //fill vector line per line
                                if(line.substr(0,1) == "#" || line.substr(0,1) == "@" ){
                                        continue;
                                }
                                else{
                                        written_file.push_back(line);
                                }
                        }
                        infile.close();
                        return written_file;
        }

        
std::pair<std::vector<double>, std::vector<double>> get_x_y(std::string in_file){
	std::ifstream infile(in_file);          			 //opens inputfile
    std::string line;
    std::vector<double> x_val, y_val;

    while (std::getline(infile, line)){ //fill vector line per line
            if(line.substr(0,1) == "#" || line.substr(0,1) == "@" ){
                    continue;
            }
            else{
				            
                    std::vector<std::string> vec = iss(line); 
                    x_val.push_back(stod(vec[0]));
                    y_val.push_back(stod(vec[1]));
            }
    }
    infile.close();
    return std::make_pair(x_val, y_val);
    }



std::vector<std::vector<double>> read_multi_col_file(std::string in_file){ 
	std::ifstream infile(in_file);          			 //opens inputfile
    std::string line;
    
    std::vector<std::vector<double>> file; 
    int j = 0;   
    while (std::getline(infile, line)){ //fill vector line per line
    	
       if(line.substr(0,1) == "#" || line.substr(0,1) == "@" || line.substr(0,1) == "&"){
               continue;
       }
       
       else{
     		
               std::vector<std::string> vec = iss(line);
               
              for(int i = 0; i<vec.size();i++){
              	
              	if(j==0){
       		 		file.resize(i+1);      
       		 		} 
       		 		try
       		 		{
              			file [i].push_back(std::stod(vec[i]));
              		} catch(const std::exception& e)
              		{
              			std::cout<<vec[i]<<'\n';
              			file [i].push_back(0);
              			fw::append_string(in_file,vec[i],"wtf.xvg");
              		}
              			
              }
              j=1;
      }
      }
      return file;
}





std::vector<std::string> merge_vecs(std::vector<std::string> vec1, std::vector<std::string> vec2){
	for(int i=0;i<vec2.size();i++){
		vec1.push_back(vec2[i]);
	}
	return vec1;
}

int get_common_min_val(std::vector<double> vec1, std::vector<double> vec2){
for(int i =0; i < vec1.size(); i++){
	if(vec1[i]==vec2[0]){
		return i;
	}
}	
	return -1;
}










}//end datahandle
#endif
