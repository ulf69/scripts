#ifndef gro_reader
#define gro_reader

#include "datahandle.cpp"
#include "abbrevations.cpp"

#include <fstream>
#include <vector>
#include <algorithm>
 
class gro_manipulator 
{
private:
	std::vector<std::string> topfile;
	std::vector<std::string> grofile;
	
public:
	gro_manipulator(const std::string& top, const std::string& gro)
		{
			topfile = dh::whole_file_to_vector(top);
			grofile = dh::whole_file_to_vector(gro);
		}
		
std::pair<std::vector<std::string>, std::vector<int>> get_topfile_molecules(){


	std::vector<std::string> lines; 	
	std::string line;	
	int siz = topfile.size();

	for(int i = siz-1;i>0;i--){
		if(topfile[i].empty() || topfile[i][0]==';'||topfile[i][0]=='['){
			break;
		}
		lines.push_back(topfile[i]);
		
	}
	std::vector<std::string>	top_residues;
	std::vector<int>			top_res_counts;
		for (int i = lines.size() -1 ; i >= 0; i-- ){			
				//dh::merge_vecs(top_residues, dh::iss(lines[i])); 
				
				std::vector <std::string> vec= dh::iss(lines[i]);
				top_residues.insert(top_residues.end(), vec.begin(), vec.end());

				top_res_counts.push_back(std::stoi(top_residues.back())); //problem here

				top_residues.pop_back();
		}
	return std::make_pair(top_residues, top_res_counts);
}





void sort_molecules(){
	auto [top_residues, top_res_counts] =get_topfile_molecules();

	//change top.residue to .gro format
	
	for(int i = 0; i<top_residues.size();i++){
		top_residues[i] = top_residues[i]+ std::string(5-top_residues[i].size(), ' ');
	}

////#########################################from here on ################################################
////#########################################gro file		################################################
////#########################################manipulation ################################################

std::vector<int> residue_size;
std::vector<std::vector<std::string>> coordinates(top_residues.size());
std::vector<std::vector<std::string>> atoms(top_residues.size());




	for (int i=2; i <grofile.size()-1 ;i++){
	//int particle_count;
	
		std::string residue_name;
		residue_name    = grofile[i].substr(5,5);
	//	particle_count  = std::stoi(grofile[i].substr(15,5)); 
		for(int j = 0; j<top_residues.size();j++){
			if(residue_name==top_residues[j]){
			//	std::cout<<residue_name<<"1"<<'\n';
			//	std::cout<<top_residues[j]<<"1"<<'\n';
				//abb::print(residue_name);
				//abb::print(top_residues[j]);
				atoms[j].push_back(grofile[i].substr(5,10));
				coordinates[j].push_back(grofile[i].substr(20));
			}
		}	
	}

for(int i= 0; i  < atoms.size();i++){
	abb::printt(atoms[i].size());
	abb::print(top_res_counts[i]);
	residue_size.push_back((atoms[i].size())/top_res_counts[i]);
}	
//abb::show(top_residues);
//abb::show_int(residue_size);

	int residue_count = 1;
	int atom_count = 0;
	std::ofstream outfile("sorted.gro");
	outfile<<grofile[0]<<'\n';
	outfile<<grofile[1]<<'\n';
for(int i= 0; i  < atoms.size();i++){
	for(int j =0 ; j<atoms[i].size();j++){
		atom_count++;
		outfile<<std::string(5-std::to_string(residue_count).size(),' ')<<residue_count<<atoms[i][j]<<std::string(5-std::to_string(atom_count).size(),' ')<<atom_count<<coordinates[i][j]<<'\n';
		if((j+1)%residue_size[i]==0){
			residue_count++;
		}
	}
}
	outfile<<grofile[grofile.size()-1]<<'\n';
outfile.close();
}

};


	
////#########################################from here on ################################################
////#########################################gro file		################################################
////#########################################manipulation ################################################
//std::vector<int> residue_size;
//std::vector<std::vector<std::string>> coordinates(top_residues.size());
//std::vector<std::vector<std::string>> atoms(top_residues.size());
//
//
//
//
//	for (int i=2; i <grofile.size()-1 ;i++){
//	int particle_count;
//	
//	std::string residue_name;
//		residue_name    = grofile[i].substr(5,5);
//		particle_count  = std::stoi(grofile[i].substr(15,5)); 
//		for(int j = 0; j<top_residues.size();j++){
//			if(residue_name==top_residues[j]){
//				atoms[j].push_back(grofile[i].substr(5,10));
//				coordinates[j].push_back(grofile[i].substr(20));
//			}
//		}	
//	}
//
//for(int i= 0; i  < atoms.size();i++){
//	residue_size.push_back((atoms[i].size())/top_res_counts[i]);
//}	
//
//
//	int residue_count = 1;
//	int atom_count = 0;
//	std::ofstream outfile(out_file);
//	outfile<<grofile[0]<<'\n';
//	outfile<<grofile[1]<<'\n';
//for(int i= 0; i  < atoms.size();i++){
//	for(int j =0 ; j<atoms[i].size();j++){
//		atom_count++;
//		outfile<<std::string(5-std::to_string(residue_count).size(),' ')<<residue_count<<atoms[i][j]<<std::string(5-std::to_string(atom_count).size(),' ')<<atom_count<<coordinates[i][j]<<'\n';
//		if((j+1)%residue_size[i]==0){
//			residue_count++;
//		}
//	}
//}
//	outfile<<grofile[grofile.size()-1]<<'\n';
//outfile.close();
//}
//
//
//
//
//
//
//
//double overlap_hist(std::string hist1, std::string hist2){
//
//	auto [x_hist1,y_hist1] 	  = dh::get_x_y(hist1);
//	auto [x_hist2,y_hist2] 	  = dh::get_x_y(hist2);	
//	int starting_position  	  = dh::get_common_min_val(x_hist1,x_hist2);
//	double normalization_area = stat::sum(y_hist1)+stat::sum(y_hist2);
//	double overlap_area    	  = 0;
//	if(x_hist1.empty() || x_hist2.empty() || y_hist1.empty() || y_hist1.empty()){
//		abb::print("-1000");
//		return -1000;
//	}
//	else if(x_hist1[0]>x_hist2[0]){
//		abb::print("-1000");
//		return -1000;
//	}
//	for(int i =starting_position; i<x_hist1.size(); i++){
//		
//		if(y_hist2[i-starting_position]<y_hist1[i]){
//			overlap_area+=y_hist2[i-starting_position];
//		}
//		else{
//			overlap_area+=y_hist1[i];
//		}
//	}
//	
//	abb::print(static_cast<int>(2*overlap_area/(normalization_area)*100));
//	return 2*overlap_area/(normalization_area);
//}
//
//
//
//
//
//
//
//
//
//
//std::vector<double> get_monolayer_z_coordinates(std::string obs_ervable){
//
//std::vector<std::string> observable =      dh::infile_to_vector(obs_ervable);
//double x_value=0;
//double maxdens=0;
//int        count =0;
//bool start_analysis=false;
//std::vector<double> monolayer_z_coordinates(2);
//        for(int i =0; i<observable.size(); i++){
//                std::istringstream iss(observable[i]);
//                std::string col;
//                int j =0;
//
//                while(iss>>col){
//                        j++;
//						
//                        if(j%2==1){
//                                        x_value=stod(col);
//                                }
//                        else if(maxdens<stod(col)){
//                        		if(stod(col)>1){
//                        			start_analysis=true;
//            					}
//                                maxdens=stod(col);
//
//                                if(stod(col)>monolayer_z_coordinates[count]){
//                                        monolayer_z_coordinates[count]=x_value;
//                                }
//                        }
//                                else if(stod(col)==0&&start_analysis==true){
//                                        if(count==0){
//                                                maxdens=monolayer_z_coordinates[0]/2;
//                                        }
//                                        else{
//                                                continue;
//                                        }
//                                        count = 1;
//
//                                }
//
//                        }
//                }
//		abb::show_dou(monolayer_z_coordinates);
//        return monolayer_z_coordinates;
//}
//
//
//
//void function_composition(std::string f_x,std::string g_x){
//
//auto [fx,fy] 	  = dh::get_x_y(f_x);
//auto [gx,gy] 	  = dh::get_x_y(g_x);
//
//fx.clear();
//gx.clear();
//
//
//
//int vector_size = gy.size();
//
//// Create indices 0..n-1
//std::vector<int> idx(vector_size);
//for (int i = 0; i< vector_size; i++) idx[i]=i;
//
//
////Sort indices according to gy values
//std::sort(idx.begin(), idx.end(), [&gy](int i1, int i2){return gy[i1]<gy[i2];});
//
//
//std::vector<double> sorted_gy(vector_size);
//std::vector<double> permuted_fy(vector_size);
//
//for (int i = 0; i < vector_size; ++i) {
//       sorted_gy[i]   = gy[idx[i]];
//       permuted_fy[i] = fy[idx[i]];
//}
////hier gehts weiter
//
//
//std::vector<double> subvector_gy;
//std::vector<double> subvector_fy;
//
//
//for (int i = 1; i < vector_size;i++){
//	if(i==vector_size-1){
//		std::string filename = "../U_d/U_"+std::to_string(subvector_gy[0])+".xvg";
//		dh::write_y_x_file( subvector_fy, subvector_gy, filename);	
//		}
//	else if (sorted_gy[i] == sorted_gy[i-1]){
//		subvector_gy.push_back(sorted_gy[i]);
//		subvector_fy.push_back(permuted_fy[i]);
//	}	
//	else{
//	std::string filename = "../U_d/U_"+std::to_string(sorted_gy[i-1])+".xvg";
//	dh::write_y_x_file( subvector_fy, subvector_gy, filename);
//	subvector_fy.clear();
//	subvector_gy.clear();
//	
//	subvector_gy.push_back(sorted_gy[i]);
//	subvector_fy.push_back(permuted_fy[i]);
//	}
//}
//}
//
#endif
