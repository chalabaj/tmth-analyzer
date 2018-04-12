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
##############################################
##############################################
##############################################

def input_check():
    global molecule,natoms   # same number and molecule for all movies and geoms
    movies  = []
    lines   = []
    geoms   = []
    molecule = sys.argv[1]   # TMTH or THMS
    movs = sys.argv[2:]      # list of movies to process
    # IF NEW MOLECULE ADDER, EXTEND HERE NEW PARAMETERS
    if molecule == "tm":
       lines_per_mol = 17
       natoms = 15
    elif molecule == "th":
       lines_per_mol = 13
       natoms = 11
    else:
       print('Wrong system (tm/th), received:  ', molecule)
       sys.exit(1)
       
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

def geoms_check(mov,lines_per_mol):   # checks the integrity of movie files, fast number_of_lines_ reader, exploits limited buffer size                 
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
       print('Nlines is not divisible by l_p_m: ',lines, lines_per_mol)
       sys.exit(1)
    return lines,geoms

# MAIN ROUTINE TO GO THROUGH EACH MOVIE - READ XYZ, CALCULATE DISTANCE, ANALYZE GEOMETRY
def process_movies(movies,geoms):

     for m,mov in enumerate(movies):   # iterate over movies
         #global natoms
         print("Processing ",m+1,"movie: ",mov)
         with open(mov,'r') as f:
          #geoms[m] = 5 
          
                  # comment for real run
          for g in range(1,int(geoms[m])+1):  # iterate over geoms in each mov file, first index is inclusive, last exclusive!
              #natoms = int(f.readline())
              f.readline()                    # atoms
              time   = f.readline().split()[6] # comment + time
              if g == 1:
                 xyz = np.zeros(shape=(natoms,3))
                 if os.path.isfile(os.path.join(os.getcwd(),'dist_mat.dat')): os.remove('dist_mat.dat')     
              print('geometry: ',g)
              for at in range(0,natoms):        # iterate over atoms in each geometry
                  line = f.readline().split()
                  xyz[at]=[float(line[1]),float(line[2]),float(line[3])]
              #print(time,"\n",xyz)
              
              dist_mat = distance_matrix(xyz)
              # Other molecules can be 
              if molecule == "tm"   :  channel  = analyze_tm(dist_mat)
              elif molecule == "th" :  channel  = analyze_th(dist_mat)
              
              print("channel: ",channel)
          f.close()
    
     
# DISTANCE MATRIX
def distance_matrix(xyz):    
# all combination of pairs: list(itertools.combinations(range(natoms),2)) - yet still need to loop over two-indices to call dist func
# create empty dist matrix - matrix is not stored for future
# brute force number of combinations len(list((itertools.combinations(range(natoms),2))))
    dist_mat = np.zeros(shape=(natoms-1,natoms))
    #print(dist_mat)
    for k in range(0,natoms):
          for l in range(k+1,natoms):
                v1, v2 = np.array(xyz[k]), np.array(xyz[l])
                dist = [(a - b)**2 for a, b in zip(v1, v2)]
                dist_mat[k][l] = math.sqrt(sum(dist))
                # print(v1,v2,l,k) # combination check
    with open('dist_mat.dat','a') as file_dist_save:  # save dist_mat in file for check if needed
      np.savetxt(file_dist_save, dist_mat, newline='\n', fmt='%.8e',footer =" ")
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
  0 Sn      
  1    C       Sn-C  = 1,2
  2    C       Sn-C  = 1,3
  3    C       Sn-C  = 1,4
  4    O       Sn-O  = 1,5
  5-14 H
  """
#1) WHERE ARE HYDROGEN ATOMS
  h_diss  = 0                      # number of dissciated hydrogen atom
  h_diss_index = []                # which H atoms are dissociated
  h_bonds = []                     # list of X - H bonds to test for shortest distance 
  h_ats_on_heavies = [0,0,0,0,0]   # how many H atoms are on each heavy atom
  
  for hydrogen_atom in range(5,natoms):
     for heavy_atom in range(0,5):                      #  last index excluded, upper diagonal l matrix, first index < second one
         h_bonds.append(dist_mat[heavy_atom][hydrogen_atom])
         print(hydrogen_atom,heavy_atom," : ",dist_mat[heavy_atom][hydrogen_atom])
    
     shortest_bond = min((j,i) for i,j in enumerate(h_bonds))         # find the smallest bod and print heavy atom related to it, enumerate over heavy atoms 0 - 5
     h_bonds.clear()  # dont need anymore now
     
     #1a) how many hydrogens are on each heavy atom
     if shortest_bond[0] < H_diss_dist: 
       h_ats_on_heavies[shortest_bond[1]] = h_ats_on_heavies[shortest_bond[1]] + 1
     else : 
       h_diss = h_diss + 1 
       h_diss_index.append(shortest_bond[1])
       if h_diss >= 2: 
         print("2 diss H, check movie: ",movies)
  print(h_ats_on_heavies) 
   
#2) Where are the heavy atoms 
  for heavy_atom in range(0,5):  
     if dist_mat[0][heavy_atom] > SnX_bond_dist
 
                                 
  return channel,h_diss         
           
           
            
##############################################
     ##########  MAIN   ##########
##############################################

molecule,movies,geoms=input_check()
print("Movie::\n".join(movies),"\n","Molecule: ",molecule,"\n Geoms: ",geoms)
print("#######################\n")

# create channel, statistics!!!
channel = process_movies(movies,geoms)   # np.array returned as upper diagonal matrix





#distance_matrix(movies)
