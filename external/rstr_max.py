#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''The MIT License (MIT)

Copyright (c) 2011 https://code.google.com/p/py-rstr-max/

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.'''

## Code obtained from: https://code.google.com/p/py-rstr-max/

from tools_karkkainen_sanders import direct_kark_sort
from array import array

class Rstr_max :

  def __init__(self) :
    self.array_str = []

  def add_str(self, str_unicode) :
    self.array_str.append(str_unicode)
    
  def step1_sort_suffix(self) :
    char_frontier = chr(2)

    self.global_suffix = char_frontier.join(self.array_str)
    
    nbChars = len(self.global_suffix)
    init = [-1]*nbChars
    self.idxString = array('i', init)
    self.idxPos = array('i', init)
    self.endAt = array('i', init)
    
    k = idx = 0
    for mot in self.array_str :
      last = k + len(mot)
      for p in xrange(len(mot)) :
        self.idxString[k] = idx
        self.idxPos[k] = p
        self.endAt[k] = last
        k += 1
      idx += 1
      k += 1

    self.res = direct_kark_sort(self.global_suffix)

  def step2_lcp(self) :
    n = len(self.res)
    init = [0]*n
    rank = array('i', init)
    LCP = array('i', init)
    
    s = self.global_suffix
    suffix_array = self.res
    endAt = self.endAt

    for i in xrange(len(self.array_str),n):
      v = self.res[i]
      rank[v] = i

    repeatLength = 0
    for j in xrange(n):
      if(repeatLength > 0) :
        repeatLength -= 1
      i = rank[j]
      j2 = suffix_array[i-1]
      if i:
        while repeatLength + j < endAt[j] and repeatLength + j2 < endAt[j2] and s[j+repeatLength] == s[j2+repeatLength]:
          repeatLength += 1
        LCP[i-1] = repeatLength 
      else:
        repeatLength = 0
    self.lcp = LCP

  def step3_rstr(self) :
    prev_len = 0
    idx = 0
    results = {}
    len_lcp = len(self.lcp) -1

#    lcp = self.lcp
#    res = self.res

    class Stack:
      pass
    stack = Stack()
    stack._top = 0
    stack.lst_max = []

    if len(self.res) == 0 :
      return {}

    pos1 = self.res[0]
    #offset1 = self.idxPos[self.res[0]]
    #idStr1 = self.idxString[self.res[0]]
    for idx in xrange(len_lcp):
      current_len = self.lcp[idx]
      pos2 = self.res[idx+1]
      #offset2 = self.idxPos[pos2]
      #idStr2 = self.idxString[pos2]
      #offset2, idStr2  = self.array_suffix[idx+1]
      end_ = max(pos1, pos2) + current_len# max(pos1, pos2) + current_len 
#      e = max((idStr1, offset1), (idStr2, offset2))
#      end_ = (e[0],e[1]+current_len)
      n = prev_len - current_len
      if n < 0 :
        #pushMany
        stack.lst_max.append([-n, idx, end_])
        stack._top += -n
      elif n > 0:
        self.removeMany(stack, results, n, idx)
      elif stack._top > 0 and end_ > stack.lst_max[-1][-1] :
        #setMax
        stack.lst_max[-1][-1] = end_

      prev_len = current_len
      pos1 = pos2
      #offset1 = offset2
      #idStr1 = idStr2
      
    if(stack._top > 0) :
      self.removeMany(stack, results, stack._top, idx+1)

    return results

  def removeMany(self, stack, results, m, idxEnd):
      prevStart = -1
      while m > 0:
        n, idxStart, maxEnd = stack.lst_max.pop()
        if prevStart != idxStart:
          #idStr = self.idxString[maxEnd-1]
          #pos = self.idxPos[maxEnd-1]
          id_ = (maxEnd, idxEnd-idxStart+1)
          if id_ not in results or results[id_][0] < stack._top:
              results[id_] = (stack._top,idxStart)
          prevStart = idxStart
        m -= n
        stack._top -= n
      if m < 0:
        stack.lst_max.append([-m, idxStart, maxEnd-n-m])
        stack._top -= m    
    
  def go(self) :
#    import time
#    t_start = t0 = time.time()
#    t0 = time.time()
    self.step1_sort_suffix()
#    print time.time() - t0
#    t0 = time.time()
    self.step2_lcp()
#    print time.time() - t0
#    t0 = time.time()
    r = self.step3_rstr()
#    print time.time() - t0
#    print time.time() - t_start
    return r

if (__name__ == '__main__') :
  str1 = 'toto'
  str1_unicode = unicode(str1,'utf-8','replace')
  rstr = Rstr_max()
  rstr.add_str(str1_unicode)
  rstr.add_str(str1_unicode)
  r = rstr.go()
  for ((id_str, end), nb), (repeatLength, startPlaceInSuffixArray) in r.iteritems():
    ss = rstr.array_str[id_str][end-repeatLength:end]
    print '[%s] %d'%(ss.encode('utf-8'), nb)
