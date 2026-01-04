#ifndef datahandle 
#define datahandle 


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

std::vector<std::string> grofile_to_vector(const std::string& in_gro){
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







std::vector<std::string> merge_vecs(std::vector<std::string> vec1, std::vector<std::string> vec2){
	for(int i=0;i<vec2.size();i++){
		vec1.push_back(vec2[i]);
	}
	return vec1;
}
double sum(std::vector<double> vec){
	double sum=0;
	for(int i=0;i<vec.size();i++){
		sum+=vec[i];
	}
	return sum;
}

int get_common_min_val(std::vector<double> vec1, std::vector<double> vec2){
for(int i =0; i < vec1.size(); i++){
	if(vec1[i]==vec2[0]){
		return i;
	}
}	
	return -1;
}

}


#endif
