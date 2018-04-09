#!/usr/bin/env python

import math
import sys
import numpy
import random
import time
import os
import subprocess
import string

##############################################
##############################################
##############################################

def input_check():
    #global movies
    molecule = sys.argv[1]   # TMTH or THMS
    if not (molecule == "tm" or molecule == "th"):
       print('Wrong system (tm/th), received:  ', molecule)
       sys.exit(1)
       
    movs = sys.argv[2:]      # list of movies to process
    print('Number of movie files requested:', len(sys.argv)-2)

    movies  = []
    lines   = []
    geoms   = []
    for mov in movs:
       cwd = os.getcwd()
       movie_path = os.path.join(cwd,mov)
       #print(movie_path)
       if os.path.isfile(movie_path):
         movies.append(movie_path)
         gc = geoms_check(movies[-1],molecule) #no need to specify order
         print(gc)
         lines.append(gc[0])
         geoms.append(gc[1])
         print(mov,' OK, Nlines: ', lines[-1], 'Ngeoms: ', geoms[-1])
       else:
         print(mov ,'ERROR NOT EXISTS.')    
    print('File check finished',"\n",'##########FILES:############################')
    return molecule,movies
#end input check

def distance_matrix(movies):

    dist_mat = []
    return dist_mat
    
def geoms_check(mov,molecule):   # checks the integry of movie files                  
    lines = 0
    geoms = 1
    buf_size = 1024 * 1024
    with open(mov,'r') as f:
     read_f = f.read
     buf = read_f(buf_size)
     while buf:
        lines += buf.count('\n')  
        #\n) is left at the end of the string, and is only omitted on the last line of the file
        buf = read_f(buf_size)
    f.close()
    return lines,geoms

##############################################
##############################################
##############################################
molecule,movies=input_check()
print("\n".join(movies))


#distance_matrix(movies)
