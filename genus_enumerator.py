"""
@author: Kush Singhal

Code to evaluate the genus of all lattice cosets of the form NL+h, for a given conductor N, and L a ternary quadratic form assumed to be of class number 1. This code can also display the p-adic orbits of a specific lattice coset under O^+(NL_p+h).
The programme first looks for a (suitably formatted and suitably titled) .txt file in the current folder that contains all of the p-adic matrices in G_{L, q', q} (c.f. paper) for all p|N. If such a .txt file is not found, it creates the .txt file on its own, which may be very ineffecient for large dL or N. It is recommended to run the matrix finder C++ file first, which quickly generates the input file. After creating the input file, run this programme.

Based on the algorithm described in Section 4 of my paper.

Input is self-explanatory at the bottom of the file, though some helpful comments are also available.
Uses the 'advanced_math_functions.py' module for some functionality.

Any questions please contact me via email.
"""

import advanced_math_functions as amf
import numpy as np, itertools as it, os

def makeSmaller(a, N):
  if a <= N//2: return a
  return a-N

def getModulus(form, N, p):
  max_ord = -1
  #k = 0
  for i in xrange(3):
    o = amf.order(form[i,i], p)
    if o >= max_ord:
      max_ord = o
      #k = i
  #print max_ord
  M = p ** (max_ord + amf.order(N,p))
  if p == 2: 
    M *= 2
    if N % 4 == 2: 
      M *= 2
  return M 

def getFileNamePath(form, N, p):
  N = p**(amf.order(N,p))
  filename = "%d-adic matrices a=(%d,%d,%d) N=%d.txt"%(p, form[0,0], form[1,1], form[2,2], N)
  return os.path.join(os.getcwd(), filename)

def writeIterTextFile(form, N, p, setX):
  filenamepath = getFileNamePath(form, N, p)
  with open(filenamepath, 'w') as f:
    for X in setX:
      if isinstance(X, tuple) and len(X) == 9:
        for i in xrange(3):
          f.write("%d %d %d\n"%X[3*i:3*i+3])
      elif isinstance(X, np.ndarray) and X.shape == (3, 3):
        for i in xrange(3):
          f.write("%d %d %d\n"%tuple(X[i]))
      else:
        f.write("Bad Input: %s\n"%repr(X))
      f.write('\n')

def existsFile(form, N, p):
  filename = getFileNamePath(form, N, p)
  #print form[0,0], form[1,1], form[2,2], N, p, filename
  return os.path.exists(filename)

def readIterTextFile(form, N, p, debug = False):
  filenamepath = getFileNamePath(form, N, p)
  c = 0
  with open(filenamepath, 'r') as matrixFile:
    X = np.array(((0,0,0), (0,0,0), (0,0,0)))
    for line in matrixFile:
      if c % 4 != 3:
        X[c%4] = map(int, line.strip().split())
      else:
        yield X
        #if debug: print 1+(c//4), ')', X[0], X[1], X[2]
      c += 1

def findAllpAdicMatricesMod(form, N = 8, p = 2, rotations_only = True, debug = False):
  N = p**(amf.order(N, p))
  M = getModulus(form, N, p)
  num_matrices = 0
  if existsFile(form, N, p):
    if debug: print "Reading from file"
    p1 = N
    for X in readIterTextFile(form, N, p, debug):
      if (N == 2 and p == 2) or (rotations_only and amf.det(X) % p1 == 1):
        num_matrices += 1
        yield X
    if debug: print "Number of Matrices =", num_matrices
    return
  cols = [[], [], []]
  for i in xrange(3): 
    for j in xrange(i): #optimize -- prevent repeat calculations
      if (form[i,i] - form[j,j]) % M == 0:
        cols[i] = list(cols[j])
        break
      #print cols[2]
    if cols[i] == []:
      for a in xrange(M):
        for b in xrange(M):
          for c in xrange(M):
            if (form[0,0]*a*a + form[1,1]*b*b + form[2,2]*c*c - form[i,i]) % M == 0:
              cols[i].append((a,b,c))
    if debug: 
      print "Number of possible candidates for column %d: %d"%(i,len(cols[i]))
  pairs_of_col12 = []
  col3 = list(cols[2])
  for col1 in cols[0]:
    for col2 in cols[1]:
      if sum(form[i,i]*col1[i]*col2[i] for i in xrange(3)) % M == 0:
        pairs_of_col12.append((col1, col2))
  del cols
  if debug: print "Number of Iterations needed: %d x %d = %d"%(len(pairs_of_col12), len(col3), len(pairs_of_col12)*len(col3))
  p1 = M#p if p != 2 else 2*p
  #return
  G = set()
  for c1, c2 in pairs_of_col12:
    for c3 in col3:
      X = np.array((c1, c2, c3)).transpose()
      if (sum(form[i,i]*c3[i]*c2[i] for i in xrange(3)) % M == 0) and              (sum(form[i,i]*c3[i]*c1[i] for i in xrange(3)) % M == 0) and amf.det(X) % M == 1:
        G.add(tuple((X % N).flatten()))  
  if debug: 
    print ("Number of %d-adic matrices ="%p), len(G)     
  writeIterTextFile(form, N, p, G)             
  for X in G:
    A = np.array(tuple(tuple(makeSmaller(X[i+3*j], N) for i in xrange(3)) for j in xrange(3)))
    if (not rotations_only):
      num_matrices += 1
      yield A
  if debug: print num_matrices


sign_changes = [np.array(((-1, 0, 0), (0, 1, 0), (0, 0, 1))), np.array(((1, 0, 0), (0, -1, 0), (0, 0, 1))), np.array(((1, 0, 0), (0, 1, 0), (0, 0, -1))), np.array(((-1, 0, 0), (0, -1, 0), (0, 0, -1)))]
for i in xrange(3):
  sign_changes.append(sign_changes[i].dot(sign_changes[3]))


def iterAllLatticeRotations(form, rotations_only = True):
  #print sign_changes
  for X in it.permutations(((1, 0, 0), (0, 1, 0), (0, 0, 1))):
    preserves = True
    #print X
    for i in xrange(3):
      j = 0
      while X[i][j] == 0: j += 1
      if form[j,j] != form[i,i]:
        preserves = False
        break
    #print preserves
    if preserves:
      nX = np.array(X)
      if not rotations_only or amf.det(nX) == 1:
        yield nX
      for S in sign_changes:
        Y = nX.dot(S)
        if not rotations_only or amf.det(Y) == 1:
          yield Y
              
          
def enumerate_pAdic_class(form, h, N, p, opfile = None, rotations_only = True, debug = False):      
  """Prints out all residues (mod N.Z_p) of rotations of the form, along with the image mod N.Z_p of h. It then numerates all lattice cosets of the type N*lattice + h' which are p-adically equivalent to N*lattice + h. Finally, it prints out the orbit of these N*lattice + h' cosets under the isometries of the (global) form
  
  form must be a diagonal 3x3 matrix in a numpy array datatype""" 
  N = p**amf.order(N,p)
  if opfile == None:
    adjective = "proper " if rotations_only else ''
    opfile = 'output %d-adic %scls - form=(%d,%d,%d) h=(%d,%d,%d) N=%d.txt'%(p, adjective, form[0,0], form[1,1], form[2,2], h[0,0], h[1,0], h[2,0], N)
  c = 0
  c_true = 0
  images = set()
  #The matrices X over here are actually [\sigma] mod N, for some \sigma in SO(NL_p) (here [\sigma] is the matrix of \sigma in the basis in which the lattice has the diagonal splitting. It is usually simple enough to compute [\sigma] by solving a couple of simultaneous diophantine equations
  file_mode = 'w'
  for X in findAllpAdicMatricesMod(form, N, p, rotations_only, debug):
    c += 1
    h1 = X.dot(h) % N
    images.add(tuple(h1.flatten()))
    with open(opfile, file_mode) as f:
      f.write(str(c) + ') ' + str(h1.flatten()) + '\n' + str(X) + "\n\n")
#      for i in xrange(3):
#        f.write(str(X[i]) + '\n')
#      f.write('\n')
    if ((h1 - h) % N == 0).all():
      c_true += 1
    file_mode = 'a'
  
  if debug:
    print c, '/', c_true
    print
  with open(opfile, 'a') as f:
    f.write("Group Index = %d / %d = %d\n\n"%(c, c_true, c/c_true)) #c / c_true is the group index [SO(NL_p) : SO(NL_p + h)]
  O_L = [X for X in iterAllLatticeRotations(form)]
  images = list(images)
  orbits = {}
  while len(images) > 0:
    h1 = np.array([[images[0][0]], [images[0][1]], [images[0][2]]])
    key = images[0]
    orbits[key] = set()
    for X in O_L:
      Y = tuple((X.dot(h1) % N).flatten())
      if Y in images: 
        images.remove(Y)
        orbits[key].add(Y)
    #print 
  with open(opfile, 'a') as f:
    for key in orbits.keys():
      f.write(str(key) + ' :')
      for t in orbits[key]:
        f.write(' ' + str(t))
      f.write('\n')

def getAllCosets(N):
  all_cosets = set()
  for h0 in xrange(N):
    for h1 in xrange(N):
      for h2 in xrange(N):
        all_cosets.add((h0,h1,h2))
  return all_cosets

def dotVector(X, h):
  if isinstance(h, tuple):
    h1 = np.array(((h[0],), (h[1],), (h[2],)))
  ans = X.dot(h1)
  if isinstance(h, tuple):
    return (ans[0,0], ans[1,0], ans[2,0])
  return ans

def classifyAllGenera(form, N, proper_classes = True, proper_genera = True, debug = False, opfile = None):
  if opfile == None:
    opfile = "output  (%d,%d,%d) N=%d %s class %s gen.txt"%(form[0,0], form[1,1], form[2,2], N, ('prop' if proper_classes else 'all'), ('prop' if proper_genera else 'all'))
  set_of_primes = amf.primeDivisors(N)
  padic_classes = {p:None for p in set_of_primes}
  stab_count_padic_classes = {p:None for p in set_of_primes}
  group_indices_of_padic_classes = {p:None for p in set_of_primes}
  if debug:
    print "Set of primes:", set_of_primes
    
  #Run through and calculate all p-adic classes
  for p in set_of_primes: 
    MOD = p ** (amf.order(N, p))
    padic_classes[p] = amf.UnionFindDS(getAllCosets(MOD))
    stab_count_padic_classes[p] = {h:0 for h in getAllCosets(MOD)}
    tot_X = 0 
    if debug: 
      print "Prime under consideration = %d\nNumber of %d L_%d + h type cosets = %d"%(p, N, p, padic_classes[p].n)
    #Run through and find orbit of each coset
    for X in findAllpAdicMatricesMod(form, MOD, p, proper_genera, debug):
      tot_X += 1
      for h in padic_classes[p].iterRepresentatives():
        h1 = dotVector(X, h)
        h2 = tuple(h1[i] % MOD for i in xrange(3))
        if h == h2:
          stab_count_padic_classes[p][h] += 1
        else:
          padic_classes[p].union(h, h2)
          
    #Find the group index [O^+(NL_p) : O^+(NL_p + h)] -- used for mass calculations
    group_indices_of_padic_classes[p] = {h: (tot_X / stab_count_padic_classes[p][h]) for h in padic_classes[p].iterRepresentatives()}
    
    if debug: #program tracker
      print "All Classes of %d L_%d + h, for all h:"%(N,p)
      padic_classes[p].display()
      print "Group Indices are:"
      for h in padic_classes[p].iterRepresentatives():
        print h, ':', group_indices_of_padic_classes[p][h]
  
  #Run through all cosets and classify the genera --> two cosets are in the same genus if they are in the same p-adic class for all p.      
  all_cosets = list(getAllCosets(N))
  genera = amf.UnionFindDS(all_cosets)
  for i in xrange(len(all_cosets)):
    for j in xrange(i+1, len(all_cosets)):
      h1 = all_cosets[i]
      h2 = all_cosets[j]
      same_genus = True
      for p in set_of_primes:
        #MOD_M = getModulus(form, N, p)
        MOD_N = p ** (amf.order(N, p))
        h1_modp = tuple(h1[k] % MOD_N for k in xrange(3))
        h2_modp = tuple(h2[k] % MOD_N for k in xrange(3))
        if not padic_classes[p].sameParent(h1_modp, h2_modp):
          same_genus = False
          break
      if same_genus:
        genera.union(h1, h2)
        
  tot_num_units_L = len(tuple(0 for _ in iterAllLatticeRotations(form, proper_classes))) #used for mass calculations
  
  cls_stringmodifier = 'Proper ' if proper_classes else ''
  gen_stringmodifier = 'Proper ' if proper_genera else ''
  if debug: print '\nGenera Classification Done' #clear space for actual output
  with open(opfile, 'w') as f:
    gen_representatives = list(genera.iterRepresentatives())
    gen_representatives.sort() #make things look nicer
    if debug: print gen_representatives
    for H in gen_representatives: #iterate through each genus
      entire_genus = sorted(list(genera.iterSameSet(H)))
      if debug: print H, ':', entire_genus
      f.write(("%sGenus of "%gen_stringmodifier) +  str(entire_genus[0]) + ':\n')
      
      #Classify all the classes in the same genus
      classes = amf.UnionFindDS(entire_genus)
      num_units = {x : 0 for x in entire_genus}
      for h in entire_genus: 
        for X in iterAllLatticeRotations(form, proper_classes):
          h_np = np.array(tuple((h[i],) for i in xrange(3)))
          h1_np = X.dot(h_np) % N
          h1 = tuple(h1_np[i,0] for i in xrange(3))
          if h1 == h:
            num_units[h] += 1
          classes.union(h1, h)
      
      #Mass Calculations for verification
      predicted_mass = amf.RationalNumber(0)
      num_classes = 0
      for h in classes.iterRepresentatives():
        num_classes += 1
        predicted_mass += amf.RationalNumber(1, num_units[h])
        f.write(("%sClass of "%cls_stringmodifier)+str(h)+ ': ')
        for h1 in classes.iterSameSet(h):
          f.write(str(h1) + ' ')
        f.write('\n')
      product_of_indices = 1
      for p in set_of_primes:
        MOD = p ** (amf.order(N, p))
        h_mod = tuple(h[i]%MOD for i in xrange(3))
        h_mod_par = padic_classes[p].getParent(h_mod)
        product_of_indices *= group_indices_of_padic_classes[p][h_mod_par]
      actual_mass = amf.RationalNumber(product_of_indices, tot_num_units_L)
      f.write("Number of %sclasses in the %sGenus = %d\n"%(cls_stringmodifier.lower(), gen_stringmodifier, num_classes))
      f.write("Mass of the %sGenus = %s; calculated mass of cosets listed above = %s\n"%(gen_stringmodifier, str(actual_mass), str(predicted_mass)))
      f.write('\n')
  if debug:
    print "Output text file ready"


if __name__ == '__main__':
  a = (1,8,8)
  form = np.array(((a[0], 0, 0), (0, a[1], 0), (0, 0, a[2])))
  
  enumerate_pAdic_class(form, h = np.array([[4], [1], [1]]), N = 8, p = 2, debug = True) 
  """
  classifyAllGenera(form, N = 8, debug = True) 
  #if debug == True, then auxiliary debugging messages will be printed to console -- useful for keeping track of what the program is doing. It will also print out all the p-adic classes (for p|N) if debug == True"""
  
  
  
  
  
  
  