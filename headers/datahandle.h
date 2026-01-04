#ifndef datahandle 
#define datahandle 


namespace datahandle  {


//makes a vector out of a string e. g. :
// string = "1	2	3	4" --> vec(4)={1,2,3,4} 
std::vector<std::string> iss(std::string str);

//	std::cout<<var<<'\n';
void printn(const T&var);

//	std::cout<<var<<'\t';
void printt(const T&var);


void show(std::vector<std::string> array);
void show_int(std::vector<int> array);
void show_dou(std::vector<double> array);


//this method inserts the whole gro file line for line in an vector
std::vector<std::string> grofile_to_vector(const std::string& in_gro);

//this one scips lines starting with # or @
std::vector<std::string> infile_to_vector(const std::string& in_gro);

//merge two vectors
std::vector<std::string> merge_vecs(std::vector<std::string> *vec1, std::vector<std::string> *vec2);


}
#endif
