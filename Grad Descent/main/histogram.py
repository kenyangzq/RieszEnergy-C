#!/usr/bin/python

import pylab as P
from numpy import *

histo = open('hist.txt','r')
lister = histo.readlines()
li = [int(elem) for elem in lister]

P.figure(1)
P.title('Histogram')
P.xlabel('Bins')
P.ylabel('Counts')
P.hist(li)

#plot=open('plot.txt','r')
#pList=plot.readlines()
#n=len(pList)
#print n
#nList=arange(0,n,1)
#yLi=arange(0,n,1)
#xLi=zeros(n)

#for i in nList:
#	yLi[i]=int(pList[i].split()[1])
#	xLi[i]=float(pList[i].split()[0])

#print xLi
#print yLi

#P.figure(2)
#P.title('Plot of Histogram')
#P.xlabel('Distance Apart')
#P.ylabel('Counts')
#P.plot(xLi,yLi)

P.show()
