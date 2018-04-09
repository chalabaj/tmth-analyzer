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
    n_movies = len(sys.argv) # 0 is script name
    molecule = sys.argv[1]   # TMTH or THMS
    if not (molecule == "tm" or molecule == "th"):
       print('Wrong system (tm/th), received:  ', molecule)
       sys.exit(1)
    
    movs = sys.argv[2:]      # list of movies to process
    print('Number of movie files requested:', n_movies)
    i = 0
    movies = []
    for mov in movs:
       base_dir = os.getcwd()
       movie_path = os.path.join(base_dir, mov)
       #print(movie_path)
       if os.path.isfile(movie_path):
         movies.append(movie_path)
         print(mov,' OK')
       else:
         print(mov ,'ERROR NOT EXISTS.')    
    print('File check finished',"\n",'##########FILES:############################')
    return molecule,movies
#end input check

def distance_matrix(movies):

    dist_mat = []
    return dist_mat
    
def n_lines(mov):
    f = open(mov)                  
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization
    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)

    return lines

##############################################
##############################################
##############################################
molecule,movies=input_check()
print("\n".join(movies))

for mov in movies: 
     nlines=n_lines(mov)
     print(nlines)
#distance_matrix(movies)
