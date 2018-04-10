#!/usr/bin/env python

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
#Constant:


def input_check():
    movies  = []
    lines   = []
    geoms   = []

    molecule = sys.argv[1]   # TMTH or THMS
    movs = sys.argv[2:]      # list of movies to process
        
    if molecule == "tm":
       lines_per_mol = 17 
    elif molecule == "th":
       lines_per_mol = 13
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


def distance_matrix(movies,geoms):
     
     for m,mov in enumerate(movies):   # iterate over movies
         print("Processing ",m,"movie: ",mov)
         with open(mov,'r') as f:
          geoms[m] = 1                        # comment for real run
          for g in range(1,int(geoms[m])+1):  # iterate over geoms in each mov file, first index is inclusive, last exclusive!
              natoms = int(f.readline())
              time  = f.readline().split()[6]
              if g == 1:
                 xyz = np.zeros(shape=(natoms,3))
              print('g: ',g)
              for at in range(0,natoms):        # iterate over atoms in each geometry
                  line = f.readline().split()
                  xyz[at]=[float(line[1]),float(line[2]),float(line[3])]
              print(time,"\n",xyz)
              # all combination of pairs, easily by all_pairs = list(itertools.combinations(range(natoms),2)) - yet still need to loop over two-indices to call eucdist func
              # CREATE empty dist matrix  - should be upper triangular (jagged)  but since the memory is stored only for each geometry - no nned
              dist_mat = np.zeros(shape=(natoms-1,natoms))
              print(dist_mat)
              for k in range(0,natoms):
                for l in range(k+1,natoms):
                    v1, v2 = np.array(xyz[k]), np.array(xyz[l])
                    dist = [(a - b)**2 for a, b in zip(v1, v2)]
                    dist_mat[k][l] = math.sqrt(sum(dist))
                    # combination check: print(v1,v2,l,k)
              #print(list(itertools.combinations(range(natoms),2)),len(list((itertools.combinations(range(natoms),2)))))
              np.savetxt('dist_mat.dat', dist_mat, newline='\n', fmt='%.8e')
          f.close()
     dist_mat = []
     return dist_mat
    
def geoms_check(mov,lines_per_mol):   # checks the integry of movie files, fast number_of_lines_ reader, exploits limited buffer size                 
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
    if not (lines % lines_per_mol):  
       geoms = lines / lines_per_mol
    else:
       print('Nlines is not divisible by lpm: ',lines, lines_per_mol)
       sys.exit(1)
    
    return lines,geoms

##############################################
##############################################
##############################################

molecule,movies,geoms=input_check()
print("Movie::\n".join(movies),"\n","Molecule: ",molecule,"\n Geoms: ",geoms)
print("#######################\n")

distance_matrix(movies,geoms)



#distance_matrix(movies)
