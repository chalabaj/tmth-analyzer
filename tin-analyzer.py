#!/usr/bin/env python

''' ANALYZE MOVIES
working in anaconda WITH python 3.6.5
PROGRAM DO GEOMETRY ANALYSIS OF MOVIE-LIKE (in XYZ format) FILES

FOR ANY MOLECULAR SYSTEM ONLY /natoms & molecule/ VARABLES NEED TO BE MODIFIED IN FUNCTION input_check() TO GET DISTANCE MATRICES

1/ EACH FILE IS CHECKED IF IT IS COMPLETE (NUMBER OF LINE MUST MATCH NATOMS+2 CRITERIA) :: def geoms_check
2/ GO THROUGH EACH MOVIE ONE BY ONE GEMEOTRY -> CREATE DISTANCE MATRIX // should be fast according to measuring of walltime
3/ PROCESS DISTANCE MATRIX ACCORDING TO PARTICULAR SYSTEM (ANALYSIS IS ALWAYS DIFFERENT AND DEPENDS ON POSSIBLE REACTION CHENNEL) 

- TO DO - CREATE  CLUSTERING ALGORITHM INDEPENDENT OF SYSTEM WHICH WILL FIND MOLECULAR FRAGMENT/SEPARATE MOLECULES 
- IF DISTANCE MATRICES ARE STORED FOR FUTURE USE, JAGGED MATRICES MIGHT BE BETTER  AND SAVE MEMORY INSTEAD OF USING THE FULL NATOMS X NATOMS DISTANCE MATRIX
'''
import math
import sys
import numpy as np
import random
import time
import os
import subprocess
import string
import itertools
import csv
##############################################
##############################################
##############################################

np.set_printoptions(linewidth  = 150)  # avoid text wrapping in console when printing np.array for checks

def input_check():
    global molecule,natoms,results_file   # same number and molecule for all movies and geoms
    if len(sys.argv) < 3:
     print("Error: not enought parameters.\nUsage: python ",sys.argv[0]," th/tm/molecule movie.xyz movie2.xyz....")
     sys.exit(1)

    movies  = []
    lines   = []
    geoms   = []
    molecule     = sys.argv[1]   # TMTH or THMS
    movs         = sys.argv[3:]  # list of movies to process
    results_file = sys.argv[2]   # file used for saving final data
    
    ##### MODIFIE HERE FOR SPECIFIC MOLECULE: #####
    if molecule == "tm":
       lines_per_mol = 17
       natoms = 15
    elif molecule == "th":
       lines_per_mol = 13
       natoms = 11
    else: sys.exit("Wrong system (tm/th)")
       
    print('Number of movie files requested:', len(sys.argv)-2)

    for mov in movs:
       cwd = os.getcwd()
       movie_path = os.path.join(cwd,mov)

       if os.path.isfile(movie_path):
         movies.append(movie_path)
         gc = geoms_check(movies[-1],lines_per_mol) #no need to specify order
         #print(gc)
         lines.append(gc[0])
         geoms.append(gc[1])
         print(mov,' OK, Nlines: ', lines[-1], 'Ngeoms: ', geoms[-1])
       else:
         print(mov ,' File NOT EXISTS.')    
    print('File check finished.',"\n",'##########FILES:############################')
    return molecule,movies,geoms
#end input check

def geoms_check(mov,lines_per_mol):   # fast number_of_lines_ reader exploiting limited buffer                 
    lines = 0
    geoms = 0
    buf_size = 1024 * 1024
    with open(mov,'r') as f:
      read_f = f.read
      buf = read_f(buf_size)
      while buf:
         lines += buf.count('\n')  
         #\n) is left at the end of the string, and is only omitted on the last line of the file
         buf = read_f(buf_size)
      f.close()
    if not (lines % lines_per_mol): geoms = lines / lines_per_mol  # 0 FALSE
    else:
       print('Nlines is not divisible by l_p_m: ',lines, lines_per_mol, mov)
       print('Check if there is empty line at the end of the file \n')
       sys.exit(1)
    return lines,geoms

# MAIN ROUTINE TO GO THROUGH EACH MOVIE - READ XYZ, CALCULATE DISTANCE, ANALYZE GEOMETRY
def process_movies(movies,geoms):
    """
     Expecting .xyz file 
     first line = natoms
     second line = comment + time/timestep information, might require change in timestep assingment split index []
    """
    analyzed_geoms = np.array([[0, 0]])                                          # main array with time and reaction channel for each geometry
    for m,mov in enumerate(movies):                                              # iterate over movies
         print("Processing ",m+1,"movie: ",mov)
         with open(mov,'r') as f:
                                       
          for g in range(1,int(geoms[m])+1):                                     # iterate over geoms in each mov file, first index is inclusive, last exclusive!
              #natoms = int(f.readline())
              atoms = f.readline()                                               # atoms
              timestep = f.readline().split()[2]                                 # comment + time 
              if g == 1:
                 xyz = np.zeros(shape=(natoms,3))
                 if os.path.isfile(os.path.join(os.getcwd(),'dist_mat.dat')): os.remove('dist_mat.dat')     
              #print('geometry: ',g)

              for at in range(0,natoms):                                         # iterate over atoms in each geometry
                  line = f.readline().split()
                  xyz[at]=[float(line[1]),float(line[2]),float(line[3])]
              #print(time,"\n",xyz)
              
              dist_mat = distance_matrix(xyz)                                    #cal dist matrix
              
              ##### MODIFIE HERE FOR SPECIFIC MOLECULE: #####
              if   molecule == "tm"  : channel  = analyze_tm(dist_mat)[0]        # analyze geometry
              elif molecule == "th"  : channel  = analyze_th(dist_mat)[0] 
                                            
              analyzed_geoms = np.append(analyzed_geoms, [[int(timestep),int(channel)]], axis = 0)  # save analyzed data for statistics
          f.close()
          
    return(analyzed_geoms)
     
# DISTANCE MATRIX
def distance_matrix(xyz):    
# all combinations of pairs: list(itertools.combinations(range(natoms),2)) - yet still need to loop over two-indices to call dist func brute force number of combinations len(<-)
# with open('dist_mat.dat','a') as file_dist_save:  # save dist_mat in file for check if needed
# np.savetxt(file_dist_save, dist_mat, newline='\n', fmt='%.8e',footer =" ")
    dist_mat = np.zeros(shape=(natoms-1,natoms)) # create empty dist matrix - matrix is not stored for future
    for k in range(0,natoms):
          for l in range(k+1,natoms):
                v1, v2 = np.array(xyz[k]), np.array(xyz[l])
                dist = [(a - b)**2 for a, b in zip(v1, v2)]
                dist_mat[k][l] = math.sqrt(sum(dist))
                # print(v1,v2,l,k) # combination check

    return dist_mat

###############################################
#ConstantS - CRITERIA FOR GEOMETRY ANALYSIS 
H_diss_dist  =  3.000 
OO_bond_dist  = 4.500
CH_bond_dist  = 3.000
HH_bond_dist  = 1.500 
SnH_bond_dist = 3.000
SnX_bond_dist = 5.000   # x = C or O 
###############################################

# GEOMETRY ANALYSIS  
def analyze_th(dist_mat):
#print(molecule,natoms)
    channel = 0
    
    return channel

# GEOMETRY ANALYSIS  
def analyze_tm(dist_mat):
  """
  atom order:
  0    Sn      
  1    C       Sn-C  = 1,2
  2    C       Sn-C  = 1,3
  3    C       Sn-C  = 1,4
  4    O       Sn-O  = 1,5
  5-14 H
  """
#1) WHERE ARE HYDROGEN ATOMS

  channel = 9
  me_diss = 0
  oh_diss = 0
  h_diss  = 0                      # number of dissciated hydrogen atom
  h_diss_index = []                # which H atoms are dissociated
  h_bonds = []                     # list of X - H bonds to test for shortest distance 
  h_ats_on_heavies = [0,0,0,0,0]   # how many H atoms are on each heavy atom
  
  for hydrogen_atom in range(5,natoms):
     for heavy_atom in range(0,5):                                    #  last index excluded, upper diagonal l matrix, first index < second one
         h_bonds.append(dist_mat[heavy_atom][hydrogen_atom])
         #print(hydrogen_atom,heavy_atom," : ",dist_mat[heavy_atom][hydrogen_atom])
    
     shortest_bond = min((j,i) for i,j in enumerate(h_bonds))         #  find the smallest bod and print heavy atom related to it, enumerate over heavy atoms 0 - 5
     h_bonds.clear()  # dont need anymore 
     
     #1a) how many hydrogens are on each heavy atom
     if shortest_bond[0] < H_diss_dist: 
       h_ats_on_heavies[shortest_bond[1]] = h_ats_on_heavies[shortest_bond[1]] + 1
     else : 
       h_diss = h_diss + 1 
       h_diss_index.append(shortest_bond[1])
       if h_diss >= 2: 
         print("2 diss H CAREFULL")
  #print("Sn,C,C,C,O: ",h_ats_on_heavies) 
   
#2) Where are the heavy atoms 
  oh_diss = 0   # OH group diss 0/1
  me_diss = 0   # Methyl group diss / if h_diss = 0 otherwise CH2, CH1 possible
  
  # Sn-O
  if dist_mat[0][4] > SnX_bond_dist: oh_diss = oh_diss + 1
     
  # Sn-C
  for heavy_atom in range(1,4):  
     if dist_mat[0][heavy_atom] > SnX_bond_dist:  
        me_diss = me_diss + 1
        if (oh_diss == 0 and h_ats_on_heavies[heavy_atom] == 2) : channel = 8


  """
  Channels:
  0 nothing happened
  1 1 Methyl diss
  2 2 Methyl diss  
  3 3 Methyl diss  
  4 OH + Methyl diss
  5 OH + 2 or more Methyl diss
  6 H diss + komplex
  7 H diss (komplex + O + H)
  8 H diss (komplex + CH2 + H)
  9 OH diss
  """       
  if h_diss == 0:
     if   oh_diss == 0:   
          if   me_diss == 1: channel = 1 
          elif me_diss == 2: channel = 2     
          elif me_diss == 3: channel = 3  
     elif oh_diss == 1:
          if   me_diss == 0: channel = 9 
          elif me_diss == 1: channel = 4 
          elif me_diss == 2: channel = 5     
          elif me_diss == 3: channel = 5
  elif h_diss == 1:
     if (me_diss == 0 and oh_diss == 0): channel = 6
     if (me_diss == 0 and oh_diss == 1 and h_ats_on_heavies[4] == 0) : channel = 7    #H from O group  
  #print(' channel,h_diss,me_diss,oh_diss,sum(h_ats_on_heavies:',channel,h_diss,me_diss,oh_diss,sum(h_ats_on_heavies))  
  #print("----------------------------------")
  
  #if channel == 9:  print("unknown geom or nothing happened")                                     
  return channel,h_diss         
           
def channel_statistics(analyze_geoms):
    """
    MODIFIE PARAMETERS FOR EACH TYPE OF MOLECULE (n_channels)
    nstep, timestep depends on simulation number of steps (e.g. nsteps in input.in)
    """
    AU_TO_FS   = 0.024189
    n_channels = 10
    n_steps    = 2100 + 1  # number of simulation steps, +1 since upper limit index is exluded
    timestep   = 10        # 
    procentual = 1         # 0 - 1 or 0-100
 
    channel_pop = np.zeros(shape=(n_steps,n_channels))   # 2D array, 0 column time, rest {1,n_channel} are channels
    #totpop     = np.zeros(shape=(n_steps))              # no need to store totpop in each step
    print("Total number of geoms: ",len(analyze_geoms)-1)
    for rec in range(1,len(analyze_geoms)):              # first row is 0,0 entry from array init
      channel = int(analyze_geoms[rec][1])
      step    = int(analyze_geoms[rec][0])
      time    = (step * timestep) * AU_TO_FS
      channel_pop[step][channel] = channel_pop[step][channel] + 1  
    for step in range(1,n_steps):    
        totpop = sum(channel_pop[step])
        time    = (step * timestep) * AU_TO_FS
        for chan in range(0,n_channels):
            channel_pop[step][chan] = (channel_pop[step][chan]/totpop) * procentual         
        line = ( str('%.4f ' %time) + ("  ".join("%.3f" %n for n in channel_pop[step])))
        
        #WRITE EACH LINE
        with open(results_file, 'a') as res_file:
         res_file.write(line)
         res_file.close()
          
##############################################
     ##########  MAIN   ##########
##############################################

molecule,movies,geoms=input_check()
print("Molecule: ",molecule,"\n Geoms: ",geoms)
print("#######################\n")

analyze_geoms  = process_movies(movies,geoms)   # np.array returning time, channel over all geoms
print(analyze_geoms)
statistic      = channel_statistics(analyze_geoms)



#distance_matrix(movies)
