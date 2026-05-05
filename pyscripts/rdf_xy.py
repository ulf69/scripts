import sys
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from MDAnalysis.lib.distances import distance_array



#load topology,  trajectory 
outfilename = sys.argv[5]
topol = sys.argv[1]
if topol == "--help":
    print("usage : python3 rdf_xy.py topol.tpr traj.xtc z_max z_min outfilename")
    print("whereas z_max and z_min refer to the borders of the monolayer")
    sys.exit()


    
traj = sys.argv[2]
#define positions of monolayers
z_max =float(sys.argv[3]) 
z_min =float(sys.argv[4])

u     = mda.Universe(topol, traj)#loads trajectory
group = u.select_atoms("resname KEG")#choose all atoms of all keggin ions
#Histogramm-paramter-monolayer-1
nbins = 100                                 #define amount of bins
if u.dimensions[0] < u.dimensions[1]:
    rmax = u.dimensions[0]/20
else:
    rmax = u.dimensions[1]/20
                                            #define max radius
rdf_hist = np.zeros(nbins)                  #returns an array with size nbins
bin_edges = np.linspace(0, rmax,nbins+1)    #creates an array with evenly spaced numbers 
n_total  = 0 

#Histogramm-paramter-monolayer-2

rdf_hist_2 = np.zeros(nbins)
n2_total = 0 

frame    = 0
for ts in u.trajectory: 
    #print current frame
    frame += 1
   # if frame > 10000:
   #     break
    print(frame)
    
    kegs   = group.residues                                               #groups all atoms of keggin ions into molecules
    coords = np.array([res.atoms.center_of_mass() for res in kegs])       #takes center of mass of all keggin ions   
    coords = coords / 10.0                                                 #scale from Angstrom to nanometers
    box_xy = u.dimensions
    box_3d = np.array([box_xy[0], box_xy[1], 1.0, 90.0 , 90.0, 90.0])     #dummy z box
    
    #upper leaflet
    coords_xy = coords[coords[:, 2 ] > z_min][:, :2]                     #only use x y coordinates
    n_total  += len(coords_xy) - 1
    coords_xy_3d = np.column_stack([coords_xy,np.zeros(len(coords_xy))]) #add z coordinates = 0 "dummz z-coordinates"
    dists     = distance_array(coords_xy_3d, coords_xy_3d, box=box_3d)      #creates an N * N matrix with diagonal elements = 0 as that are the distance from particle i to itself
                                                                         #in addition it calculates every distance twice (N*N-N)/2 = N*(N-1)/2
    dists_flat = dists[np.triu_indices(len(coords_xy), k = 1)]        
    counts, _= np.histogram(dists_flat, bins=bin_edges)
    rdf_hist += counts



    #lower leaflet
    coords_xy_2 = coords[coords[:, 2 ] < z_max][:, :2]
    n2_total  += len(coords_xy_2)-1
    
    coords_xy_3d_2  = np.column_stack([coords_xy_2,np.zeros(len(coords_xy_2))])       #add z coordinates = 0 "dummz z-coordinates"
    dists_2         = distance_array(coords_xy_3d_2, coords_xy_3d_2, box=box_3d)      #creates an N * N matrix with diagonal elements = 0 as that are the distance from particle i to itself
    dists_flat_2    = dists_2[np.triu_indices(len(coords_xy_2), k = 1)]        
    counts_2, _= np.histogram(dists_flat_2, bins=bin_edges)
    rdf_hist_2 += counts_2
    
    
    

    
radii = 0.5*(bin_edges[1:]+ bin_edges[:-1])
areas = np.pi * (bin_edges[1:]**2-bin_edges[:-1]**2)


n_avg = n_total / len(u.trajectory)
density = n_avg / (rmax*rmax)
rdf =  rdf_hist / (areas*density*len(u.trajectory))





n2_avg    = n2_total/len(u.trajectory)
density_2 = n2_avg / (rmax*rmax)
rdf_2     = rdf_hist_2 / (areas*density_2*len(u.trajectory))

n2_avg += 1
n_avg  += 1
file="rdf/22upper_leaflet" + outfilename + ".xvg"
with open(file, "w") as f:
    f.write(f"#average amount of keggin ions: {n_avg}\n")
    f.write(f"#keggin ions have z value >:  {z_min}\n")
    np.savetxt(f, np.column_stack([radii ,rdf])  , header="r_xy[nm] g_xy(r)")

    
file1="rdf/22lower_leaflet" + outfilename + ".xvg"
with open(file1, "w") as f1:
    f1.write(f"#average amount of keggin ions:{n2_avg}\n")
    f1.write(f"#keggin ions have z value <:{z_max}\n" )
    np.savetxt(f1, np.column_stack([radii ,rdf_2])  , header="r_xy[nm] g_xy(r)")


#np.savetxt("rdf/22monolayer1" + outfilename + ".xvg", np.column_stack([radii ,rdf])  , header="r_xy[nm] g_xy(r)")
#np.savetxt("rdf/22monolayer1" + outfilename + ".xvg", np.column_stack([radii ,rdf_2]), header="r_xy[nm] g_xy(r)")
