# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 19:06:54 2021

@author: Kush Singhal

This programme can be used to find the group G_{L, q', q} if the group G_{L, Q', Q} has already been computed, where Q is a larger power of the prime p than q. For instance, if the group G_{L, 16, 4} has been computed, then there is no need to run the C++ programme to evaluate G_{L, 16, 2} again, simply run this programme instead.

Input is self-explanatory.
"""
#import numpy as np

def reduceMatrices(a, p, N, new_N): #N and new_N are assumed to be powers of p
  filenameformat = "%d-adic matrices a=(%d,%d,%d) N=%d.txt"
  oldfilename = filenameformat%(p, a[0], a[1],a[2], N)
  newfilename = filenameformat%(p, a[0], a[1],a[2], new_N)
  reduced_matrices = set()
  with open(oldfilename, 'r') as oldf:
    matrix = [[],[],[]]
    eof_reached = False
    while not eof_reached:
      for i in xrange(3):
        row = oldf.readline()
        if row == '': 
          eof_reached = True
          break
        matrix[i] = map(int, row.strip().split())
        for j in xrange(3):
          matrix[i][j] %= new_N
        matrix[i] = tuple(matrix[i])
      reduced_matrices.add(tuple(matrix))
      if not eof_reached: oldf.readline()
    with open(newfilename, 'w') as newf:
      for matrices in reduced_matrices:  
        for i in xrange(3):
          newf.write("%d %d %d\n"%matrices[i])
        newf.write('\n')
        
reduceMatrices((3,4,12), 2, 4, 2)