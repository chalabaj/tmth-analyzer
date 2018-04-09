#!/usr/bin/env python

import math
import sys
import numpy
import random
import time
import os

def input_check():
    n_movies = len(sys.argv)  #sys.argv[0] is script name
    movies = sys.argv  # list of movies to process
    print('Number of movie files requested:', n_movies)
    print('Movie: ',movies)
    i = 0
    movie = []
    for mov in movies:
       print(mov)
       movie_path = os.getcwd()+'/'+mov
       print(movie_path)
       if os.path.isfile(movie_path):
         movie.append(movie_path)
         print('File:',i,movie[i])
         i = i + 1
       else:
         print('ERROR: File', movie_path ,'does not exists.')    
   
     
input_check()
print('File check finished!')
