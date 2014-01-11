# ScoringXlinks.py
# This script was implemented to run within IMP software package 
# Documentation nd information on how to run it can be found where:
# https://github.com/integrativemodeling/hybrid_ms_method
# This script scores the violation of cross-links in the
# candidate model structures

import re
#import os, numarray
import operator
from operator import itemgetter
from math import *
from Numeric import * # imports numerical python
import csv, sys, os

# get List of theoretical X-links (vxl type of file)
if len(sys.argv) !=3:
   sys.exit(" Must provide the input text file")

# Convert the arguments into strings and number
InFname = str(sys.argv[1])
RefD = float(sys.argv[2])
#filenameA = str(sys.argv[1])


#this function reads the input txt file
def ReadList_Xlinks(InFname):
        #list =[]
        InFile = open(InFname, 'r')
        InFLines = InFile.read().splitlines() # string of lines from pdb file
        InFile.close()
        return InFLines
# Two lists are generated: one with the first crosslink and and another with the second linker
# For inter-protein crosslinks, the links from two different chains (subunits) of the same proteins
# or two different proteins are used

def getList_Xlinks0(lines):
       outList0=[]
       #outList2=[]
       for i in lines[:]:
              l = re.findall ("(\S+)", i) 
              L0 = int(l[4])
	      #L1 = int(l[8])
              outList0.append(L0)
	      #outList2.append(L2)
       return outList0#, outList2


def getList_Xlinks1(lines):
       #outList1=[]
       outList1=[]
       for i in lines[:]:
              l = re.findall ("(\S+)", i) 
              #L0 = int(l[4])
	      L1 = int(l[8])
              #outList1.append(L1)
	      outList1.append(L1)
       return outList1

ls0= getList_Xlinks0(ReadList_Xlinks(InFname))
ls1= getList_Xlinks1(ReadList_Xlinks(InFname))

print ls0, ls1

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The user is prompted to enter the first and the last number of structures
NAstrs=raw_input("enter the first number of structures: ")
NBstrs=raw_input("enter the last number of structures: ")
NA=int(NAstrs)
NB=int(NBstrs)

list0 = []
list1 = []
list2 = []
list3 = []
list4 = []
list5 = []
list6 = []
listNames=[]

##Read input file *.mfj type od structures
for i in range(NA, NB):
       #InputFileName="module1_20k_15Dec."+str(i)+".mfj"
       InputFileName="models."+str(i)+".mfj"
       def ReadFile(InputFile):
              InputFile = open(InputFileName,'r') # Open file to read
              InputFileLines = InputFile.read().splitlines() # string of lines from pdb file
              InputFile.close()
              return InputFileLines

# Store only the coordinates from the files
       def coordinates(lines):
              outX=[]
              outY=[]
              outZ=[]
              outR=[]
              for i in lines[6:]:
                     l = re.findall ("(\S+)", i) 
                     X0 = float(l[0])
                     Y0 = float(l[1])
                     Z0 = float(l[2])
                     R0 = float(l[3])
                     outX.append(X0)
                     outY.append(Y0)
                     outZ.append(Z0)
                     outR.append(R0)
                     #print Scr2
              return outX, outY, outZ, outR
       b= coordinates(ReadFile(InputFileName))

# get the euclidean distances in the stored coordinates that are specified in the input *txt file
       listDlink = []
       def getDistances(lines):
	      for i,j in zip(ls0, ls1):
			Dlink = sqrt(pow(b[0][512+i-6]-b[0][j-1], 2) + pow(b[1][512+i-6]-b[1][j-1], 2) + pow(b[2][512+i-6]-b[2][j-1], 2))
			listDlink.append(Dlink)
	      return listDlink
	
       cA = getDistances(ReadFile(InputFileName))

# Count the violation of the cross-link restraints for a given cut-off distance       
       def getScoreXlinks(lines, cLink):
              if cLink<RefD:
                     score1=0.0
              else:
                     score1=1.0
              Score_Xlinks= score1 
              return Score_Xlinks
              
       
# append the file names into a lit
       listNames.append(InputFileName)
       listS = []

#calculate the sum of violations
       def getScores(lines):
	       for i in range(len(ls0)):
		       S = getScoreXlinks(ReadFile(InputFileName), cA[i])
		       listS.append(S)
	       return listS
# print overall violation score       
       Sc = getScores(ReadFile(InputFileName))
       print Sc
       TotalScore = sum(Sc)
       print TotalScore
       
#%%%%%%%%%%%%%%%%%%%%%%%       
       eSF=str(TotalScore)
       list3.append(eSF)

# write scores into a list (*.csv type of file)
       SummaryFileName = 'outList-Xlinks_Scores.csv'
       SummaryFile = open(SummaryFileName, 'w')
       SummaryFile.write(str(listNames) +'\n')
       SummaryFile.write(str(list3) +'\n')
       

SummaryFile.flush()
SummaryFile.close()

