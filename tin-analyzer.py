#!/usr/bin/env python

import math
import sys
import numpy as np
import random
import time
import os
import subprocess
import string

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
         print(mov ,'ERROR NOT EXISTS.')    
    print('File check finished',"\n",'##########FILES:############################')
    return molecule,movies,geoms
#end input check
def eudis5(v1, v2):
    dist = [(a - b)**2 for a, b in zip(v1, v2)]
    dist = math.sqrt(sum(dist))
    return dist

def distance_matrix(movies,geoms):
     
     for m in range(len(movies)):
       print(m,movies[m])
       with open(movies[m],'r') as f:
       
        geoms[m] = 3
        for g in range(1,int(geoms[m])):
          natoms = int(f.readline())
          title  = f.readline()
          if g == 1:
            xyz = np.zeros(shape=(natoms,3))
          #print('g: ',g)
          for at in range(natoms):
             line = f.readline().split()
             xyz[at]=[float(line[1]),float(line[2]),float(line[2])]
          print(xyz)                        
        f.close()
     dist_mat = []
     return dist_mat
    
def geoms_check(mov,lines_per_mol):   # checks the integry of movie files                  
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
print("\n".join(movies),"\n",molecule,geoms)
print("#######################\n")

distance_matrix(movies,geoms)



#distance_matrix(movies)