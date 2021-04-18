import numpy as np, itertools as it

def primeDivisors(N):
  l = []
  p = 2
  while p <= N//2:
    if N % p == 0:
      l.append(p)
      while N % p == 0:
        N /= p
    if p % 2 == 0:
      p += 1
    else:
      p += 2
  if N > 1:
    l.append(N)
  return l

class RationalNumber():
    def __init__(self, p = 0, q = 1):
        if q < 0:
            q = -q
            p = -p
        self.num = p
        self.den = q
        self.lowestTerms()
        
    @staticmethod
    def convert(obj):
        if isinstance(obj, int):
            return RationalNumber(obj, 1)
        if isinstance(obj, float):
            return RationalNumber(int(100000*obj), 100000)
        return obj
    
    def lowestTerms(self):
        g = np.gcd(self.num, self.den)
        if g > 1:
            self.num /= g
            self.den /= g
    def __getcomparitor(self, other):
        return (self.num * other.den - self.den * other.num)
    def __eq__(self, other):
        other = RationalNumber.convert(other)
        return (self.__getcomparitor(other) == 0)
    def __lt__(self, other):
        other = RationalNumber.convert(other)
        return (self.__getcomparitor(other) < 0)
    def __le__(self, other):
        other = RationalNumber.convert(other)
        return (self.__getcomparitor(other) <= 0)
    def __gt__(self, other):
        other = RationalNumber.convert(other)
        return (self.__getcomparitor(other) > 0)
    def __ge__(self, other):
        other = RationalNumber.convert(other)
        return (self.__getcomparitor(other) >= 0)
    def __ne__(self, other):
        other = RationalNumber.convert(other)
        return (self.__getcomparitor(other) != 0)
    def __repr__(self): return str(self)
    def __str__(self):
        self.lowestTerms()
        if self.den == 1:
            return str(self.num)
        else:
            return "%d/%d"%(self.num, self.den)
    def __mul__(self, other):
        other = RationalNumber.convert(other)
        return RationalNumber(self.num * other.num, self.den * other.den)
    def __div__(self, other):
        other = RationalNumber.convert(other)
        return RationalNumber(self.num*other.den, self.den * other.num)
    def __neg__(self):
        return RationalNumber(-self.num, self.den)
    def __add__(self, other):
        other = RationalNumber.convert(other)
        return RationalNumber(self.num * other.den + self.den * other.num, self.den * other.den)
    def __sub__(self, other):
        return self + (-other)
    def approx(self):
        return float(self.num) / float(self.den)
    def __floordiv__(self, other):
        return int((self/other).approx())
    def __pow__(self, other):
        if isinstance(other, int):
            return RationalNumber(self.num**other, self.den**other)
        if isinstance(other, RationalNumber):
            return self.approx() ** other.approx()
    def __float__(self):
        return self


def det(X):
    sign = {(0, 1, 2):1, (0, 2, 1):-1, (1, 0, 2):-1, (1, 2, 0):1, (2, 0, 1):1, (2, 1, 0):-1}
    d = 0
    for p1, p2, p3 in it.permutations(xrange(3)):
        d += sign[(p1, p2, p3)]*X[0][p1]*X[1][p2]*X[2][p3]
    return d

def order(n, p):
    c = 0
    while n % p == 0:
        c += 1
        n /= p
    return c

class UnionFindDS:
  #a Union Find Data Structure -- store and compute orbits extremely efficiently (on average constant time)
  def __init__(self, objects):
    self.indices = {}
    self.objects = []
    self.n = len(objects)
    c = 0
    for x in objects:
      self.indices[x] = c
      self.objects.append(x)
      c += 1
    self.parent = [-1 for _ in xrange(self.n)]
  
  def getParent(self, x):
    i = self.find(x)
    if i == -1:
      raise ValueError("Object %s not in structure"%repr(x))
    return self.objects[i]
  
  def find(self, x):
    if x not in self.indices:
      return -1
    k = self.indices[x]
    stack = []
    i = k
    while self.parent[i] != -1:
      stack.append(i)
      i = self.parent[i]
    ans = i
    for j in stack:
      self.parent[j] = ans
    return ans
  
  def union(self, x, y):
    if x not in self.indices:
      raise ValueError("Object %s not in structure"%repr(x))
    if y not in self.indices:
      raise ValueError("Object %s not in structure"%repr(y))
    i = self.find(x)
    j = self.find(y)
    if i != j:
      try:
        if self.objects[i] < self.objects[j]:
          i, j = j, i
      except:
        pass
      self.parent[i] = j
  
  def sameParent(self, x, y):
    return (self.find(x) == self.find(y))
  
  def iterRepresentatives(self):
    for i in xrange(self.n):
      if self.parent[i] == -1:
        yield self.objects[i]
        
  def iterSameSet(self, rep):
    if self.find(rep) == -1:
      raise ValueError("Representative %s not in structure"%repr(rep))
    for x in self.objects:
      if self.sameParent(x, rep):
        yield x
  
  def getDict(self):
    d = {}
    for x in self.objects:
      p_index = self.find(x)
      parent = self.objects[p_index]
      if parent not in d.keys():
        d[parent] = [x]
      else:
        d[parent].append(x)
    return d
  
  def display(self):
    d = self.getDict()
    for x in d.keys():
      print repr(x), ':', 
      for y in d[x]:
        print repr(y),
      print
      
