#ifndef gro_reader
#define gro_reader

#include "datahandle.cpp"
#include "abbrevations.cpp"

#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>
 
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
			
				std::vector <std::string> vec= dh::iss(lines[i]);
				top_residues.insert(top_residues.end(), vec.begin(), vec.end());

				top_res_counts.push_back(std::stoi(top_residues.back())); //problem here

				top_residues.pop_back();
		}
	return std::make_pair(top_residues, top_res_counts);
}





void sort_molecules(){
	auto [top_residues, top_res_counts] = get_topfile_molecules();

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
	
		std::string residue_name;
		residue_name    = grofile[i].substr(5,5);
		for(int j = 0; j<top_residues.size();j++){
			if(residue_name==top_residues[j]){
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

std::vector<std::string> get_particles_in_box(double x_min, double x_max,double y_min, double y_max,double z_min, double z_max){
	std::vector<std::string> particles_in_box;
	std::string particle_ph;
	for(int i=2; i < grofile.size() - 1; i++){

		double x_coord = std::stod(grofile[i].substr(20,8));
		double y_coord = std::stod(grofile[i].substr(28,8));
		double z_coord = std::stod(grofile[i].substr(36,8));
		if(x_coord < x_max && x_coord > x_min && y_coord < y_max && y_coord > y_min && z_coord < z_max && z_coord > z_min ){
			std::string particle = grofile[i].substr(0,5);
			if(particle == particle_ph){
				continue;
			}
			else{
				particles_in_box.push_back(particle);
			}
			particle_ph = particle;
		}
	
	}
	
		return particles_in_box;
}

void translate_box(std::string box){

	double x_min = 100000000 ;
	double y_min = 100000000 ;
	double z_min = 100000000 ;
	double x_max = 0 ;
	double y_max = 0 ;
	double z_max = 0 ;

	std::vector<std::string>	boxfile = dh::whole_file_to_vector(box);
	
	int siz = boxfile.size();
	
	for(int i = 2; i < siz-1;i++ ){
		double x_coord = std::stod(boxfile[i].substr(20,8));
		double y_coord = std::stod(boxfile[i].substr(28,8));
		double z_coord = std::stod(boxfile[i].substr(36,8));
		if(x_coord>x_max){
			x_max = x_coord;
		}
		else if(x_coord < x_min){
			x_min = x_coord;
		}
		if(y_coord>y_max){
			y_max = y_coord;
		}
		else if(y_coord < y_min){
			y_min = y_coord;
		}
		if(z_coord > z_max){
			z_max  = z_coord;
		}
		else if(z_coord < z_min){
			z_min = z_coord;
		}
	}
	std::ofstream outfile("cut_trans_box.gro");
		outfile<<boxfile[0]<<'\n';
		outfile<<siz -3<<'\n';


	
	for(int i = 2; i < siz-1;i++ ){
	
		double x_coord = std::stod(boxfile[i].substr(20,8)) - x_min;
		double y_coord = std::stod(boxfile[i].substr(28,8)) - y_min;
		double z_coord = std::stod(boxfile[i].substr(36,8)) - z_min;
		
		
		outfile<<std::fixed << std::setprecision(3) <<boxfile[i].substr(0,20)<<std::string(11-std::to_string(x_coord).size(),' ')<<x_coord<<std::string(11-std::to_string(y_coord).size(),' ')<<y_coord<<std::string(11-std::to_string(z_coord).size(),' ')<<z_coord<<'\n';	
	}
	outfile << "   "<< std::fixed << std::setprecision(5) << x_max-x_min<< "   "<< y_max-y_min<< "   "<< z_max-z_min << "\n";
	outfile.close();

}


void cut_box(double x_min, double x_max,double y_min, double y_max,double z_min, double z_max){

	std::vector<std::string> particles = get_particles_in_box( x_min,  x_max, y_min,  y_max, z_min,  z_max);
	std::ofstream outfile("cut_box.gro");
	int j = 0;
	int count = 0;
	outfile<<grofile[0]<<'\n';
	outfile<<particles.size()<<'\n';
	
	for(int i = 0; i < grofile.size();i++){
		if(grofile[i].substr(0,5)==particles[j]){
		
			count++;
			std::string spaces =  std::string(5-std::to_string(j+1).size(),' ');		
			outfile<<j+1<<spaces<<grofile[i].substr(5,10)<< std::string(5-std::to_string(count).size(),' ') <<count << grofile[i].substr(20,24)<<'\n';			
		}
		if(j+1<particles.size() &&grofile[i].substr(0,5)==particles[j+1]){
			j++;
			count++;
				std::string spaces =  std::string(5-std::to_string(j+1).size(),' ');		
				outfile<<j+1<<spaces<<grofile[i].substr(5,10)<< std::string(5-std::to_string(count).size(),' ') <<count << grofile[i].substr(20,24)<<'\n';			
		} 		
	}
	outfile<<grofile[grofile.size()-1]<<'\n';
	outfile.close();
	translate_box("cut_box.gro");

}

void multiply_box(std::string box){
	std::vector<std::string>	boxfile = dh::whole_file_to_vector(box);

	std::vector<int> x_shift = { 1,  2, 0, 2, 1,  1,  2};
	std::vector<int> y_shift = { 0,  0, 1, 2, 1,  2,  1};
	std::ofstream outfile("multiplied_box.gro");
	outfile<<boxfile[0]<<'\n';
	outfile<<8*stod(boxfile[1])<<'\n';
	int siz = boxfile.size();

	for(int i = 2; i < siz-1; i++){
		outfile	<< boxfile[i]<<'\n';	
	}
	int particles_in_starting_box = std::stoi(boxfile[siz-2].substr(0,5));
	int atom_in_starting_box = std::stoi(boxfile[siz-2].substr(15,5));
	double box_x = std::stod(boxfile[siz-1].substr(0,10));
	double box_y = std::stod(boxfile[siz-1].substr(10,10));
	double box_z = std::stod(boxfile[siz-1].substr(20,10));
	
	int step = 0;
	for(int i = 0; i <x_shift.size();i++ ){
			step += 1;
			for(int k = 2; k < siz-1; k++){
				int particle = step*particles_in_starting_box+std::stoi(boxfile[k].substr(0,5));
				std::string residue_name = boxfile[k].substr(5,5);
				std::string atom_name = boxfile[k].substr(10,5);
				int atom_number = stoi(boxfile[k].substr(15,5))+step*atom_in_starting_box;
				double x_coord = std::stod(boxfile[k].substr(20,8)) + x_shift[i]*box_x ;
				double y_coord = std::stod(boxfile[k].substr(28,8)) + y_shift[i]*box_y;
				double z_coord = std::stod(boxfile[k].substr(36,8));

				
				outfile<<particle<< std::string(5-std::to_string(particle).size(),' ')<<residue_name<<atom_name<<std::string(5-std::to_string(atom_number).size(),' ')<< atom_number;
				 
				outfile<<std::fixed << std::setprecision(3) <<std::string(11-std::to_string(x_coord).size(),' ')<<x_coord<<std::string(11-std::to_string(y_coord).size(),' ')<<y_coord<<std::string(11-std::to_string(z_coord).size(),' ')<<z_coord<<'\n';					
			}		
			
	}
	
	outfile << "   " << std::fixed << std::setprecision(5) << box_x+2*box_x<< "   "<< box_y+box_y*2<< "   "<< box_z << "\n";
	
	outfile.close();
}




void mirror_box_z(std::string box){
	std::vector<std::string>	boxfile = dh::whole_file_to_vector(box);
	std::ofstream outfile("mirrored_box.gro");
	
	outfile<<boxfile[0]<<'\n';
	outfile<<boxfile[1]<<'\n';
	int siz = boxfile.size();
	double box_z = std::stod(boxfile[siz-1].substr(24,10));
	for(int i = 2; i < siz-1; i++){
		double z_coord = stod(boxfile[i].substr(36,8))*(-1) + box_z ;
		outfile	<<std::fixed << std::setprecision(3)<< boxfile[i].substr(0,36)<<std::string(11-std::to_string(z_coord).size(),' ') <<z_coord<<'\n';	
	}
    	
	outfile <<boxfile[siz-1]<<'\n';
	outfile.close();
}


void merge_boxes(std::string box1, std::string box2){
	std::vector<std::string>	boxfile1 = dh::whole_file_to_vector(box1);
	std::vector<std::string>	boxfile2 = dh::whole_file_to_vector(box2);

	std::ofstream outfile("merged_box.gro");		
	outfile<<boxfile1[0]<<'\n';
	int all_particles= stoi(boxfile1[1])+stoi(boxfile2[1]);
	outfile<<all_particles<<'\n';

	int siz1 = boxfile1.size();

	for(int i = 2; i < siz1-1; i++){	
		outfile	<< boxfile1[i]<<'\n';	
	}
	
	int particles = stoi(boxfile1[siz1-2].substr(0,5));
	int atoms     = siz1-3;
	
	
	int siz2 = boxfile2.size();
	for(int i = 2; i < siz2-1; i++){	
		
		int particles_box_2 =  stoi(boxfile2[i].substr(0,5))+particles;
		int atoms_box_2     =  stoi(boxfile2[i].substr(15,5))+atoms;
		double new_z		=  stod(boxfile2[i].substr(36,8))+stod(boxfile1[siz1-1].substr(24,7));
		abb::print(new_z);
		outfile	<< particles_box_2<<std::string(5-std::to_string(particles_box_2).size(),' ') << boxfile2[i].substr(5,10) <<atoms_box_2<< std::string(5-std::to_string(atoms_box_2).size(),' ') <<boxfile2[i].substr(20,16)<<std::string(11-std::to_string(new_z).size(),' ') <<new_z<< '\n';	
	}

	double new_box_z = stod(boxfile1[siz1-1].substr(24,7)) + stod(boxfile2[siz2-1].substr(24,7));
	
	outfile<<std::fixed << std::setprecision(5) << boxfile1[siz1-1].substr(0,24)<< new_box_z;
} 
};




#endif
