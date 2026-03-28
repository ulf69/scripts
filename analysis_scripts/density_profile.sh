
#!/bin/bash

source /usr/local/gromacs/2025.2_cuda_ompi/bin/GMXRC


#######################################################################################################
##############################declare your input data for the residue analysis#########################
#######################################################################################################

exe="gmx_mpi"									#change here the gmx version you need
list=(KEG KEG) 									#change hier the name of the residues 
axis=z 				 							#declare here the axis you have a look on
ptf=/home/scw27692/Masterthesis/POPC-KEG/10KEG	#change here the path to your file
file=($ptf/pr/10Keg $ptf/pr2/10Keg)				#declare here your files you want to analyzes
sl=1000 										#declare amount of bins you want to use in the histogramm file
b=499000										#equillibration time
o=(test/pr test/pr2)							#output directory
mkdir test										#create dir
 
#######################################################################################################
##############################density profile of all residues in list##################################
#######################################################################################################
i=0
for res in "${list[@]}";do
rm ${o[i]}/${res}.xvg
echo $res  | $exe density -f ${file[i]}.xtc -b $b -s  ${file[i]}.tpr  -d $axis -sl $sl -o ${o[i]}${res}.xvg
((i++))

done


#######################################################################################################
##############################declare your input data for the atom analysis###########################
#######################################################################################################

atoms=(P N) 									#change hier the name of the atoms 
#axis=z 				 						#declare here the axis you have a look on
#ptf=/home/scw27692/Masterthesis/POPC-KEG/10KEG	#change here the path to your file
#file=($ptf/pr/10Keg $ptf/pr2/10Keg)			#declare here your files you want to analyzes
#sl=1000 										#declare amount of bins you want to use in the histogramm file
#b=499000										#equillibration time
#o=(test/pr test/pr2)							#output directory
#mkdir test
i=0
for a in "${atoms[@]}";do
rm ${o[i]}${a}.xvg
printf "a $a \n q\n" | $exe make_ndx    -f  ${file[i]}.tpr 
printf "a $a \n q\n"

echo $a  | $exe density -f ${file[i]}.xtc -b $b -s  ${file[i]}.tpr -n index -d $axis -sl $sl -o ${o[i]}${a}.xvg
rm index.ndx
((i++))
done
