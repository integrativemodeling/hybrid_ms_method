mobcal.py

# This script reads the mobcal exe file and calculate the PA for all structures in the folder
# Input type of files: *.mfj
# For use in coarse grained type of files

import sys, os, subprocess, shutil,time,glob
import re
#NumTraj = 10000 #Number of trajectories for each simulation.
#NumRepeats = 16 #Number of replication simulations for each structure.

################################################################################
#   Retrieves an array of strings from a text file.
def GetLines(InputFileName):
    # Read input file
    # Creates an array containing elements for each lines of the input file.
    InputFile = open(InputFileName,'r') # Open file to read
   # Copy the specific input file to the general file.
    shutil.copyfile(InputFileName,'t20.mfj')
    Lines = InputFile.read().splitlines() # string of lines from pdb file
    InputFile.close()
    return Lines
####################################################################

####################################################################
#   Return mass for an element.
#   Average atomic weights from http://en.wikipedia.org/wiki/List_of_elements_by_atomic_weight
MassList = {
    'H':1.00794,
    'C':12.0107,
    'N':14.0067,
    'O':15.9994,
    'Na':22.98976928,
    'Mg':24.3050,
    'P':30.973762,
    'S':32.065,
    'Cl':35.453,
    'K':39.0983,
    'Ca':40.078}
def mass(element):
    if element in MassList:
        return MassList[element]
    else:
        print 'element not in list'
        return 40.0 #   To prevent failure for random metal ions.
####################################################################

####################################################################
# NEEDS ATTENTION
# This section calls the compiled fortran code.  Unfortunately,
# this often results in errors, thus the error handling sections.
# Not sure how to address these problems better, but this
# solution does not affect the output.
#
def CalcCCS(text):
    try:
        p = subprocess.Popen('./mobcal_cg',bufsize=-1)
        p.wait()
    except KeyboardInterrupt:
        print
        print "User Exit"
        os._exit(0)
    except:
        print
        print "Error calling PA/EHSS Code.  Will try again."
        print
        time.sleep(0.5)
        CalcCCS(text)
################################################################################

# Loads PDB filenames from command line.  Wild cards are processed.
MFJs = []
for item in sys.argv[1:]:
    MFJs += glob.glob(item)

print
print MFJs
print

# Creates a summary file.
SummaryFileName = 'outPA-summary.csv'
SummaryFile = open(SummaryFileName, 'w')

# Processes MFJ files, one at a time.
for MFJ in MFJs:

    print
    print "Performing calculations for: " +MFJ 
    print

    # Creates an input file with a name based on the PDB file.
    OutPutFileName = MFJ + ".CCS.in"
    #OutPutFileName = "t20.mfj"

    # Copy the specific input file to the general file.
    #shutil.copyfile(InputFileName,'t20.mfj')


    Spheres = []
    for line in GetLines(MFJ):
        if line[0:4] == '    ':
	    l = re.findall ("(\S+)", line)
	    Spheres.append([float(l[0]),float(l[1]),float(l[2]),float(l[3])])
#            Spheres.append([float(line[30:38]),float(line[38:46]),float(line[46:54]),mass(line[77:78])])
    OutPutFile = open(OutPutFileName, 'w')
#    OutPutFile.write(str(NumTraj) + '\n')
#    OutPutFile.write(str(NumRepeats) + '\n')
    OutPutFile.write(str(len(Spheres))+'\n')
    for Sphere in Spheres:
        OutPutFile.write("%e\t%e\t%e\t%e\n" % (Sphere[0],Sphere[1],Sphere[2],Sphere[3]))
    OutPutFile.flush()
    OutPutFile.close()

    # Copy the specific input file to the general file.
    shutil.copyfile(OutPutFileName,'CCSin.txt')

    # Run the compiled fortran code.
    CalcCCS('./mobcal_cg')

    # Copy the generic output file to the specific file.
    shutil.copyfile('t20.out',MFJ + ".CCS.out")
    
    # Average results from replicate simulations.
    PA = []
    EHSS = []
    for line in GetLines(MFJ + ".CCS.out"):
#	if line[0:12] == 'average PA c': 
	if line.find('average PA cross section')>=0:
#        dx1 = lines.index('    set     PA CS       PA MOB^-1      EHSS CS      EHSS MOB^-1    ASYMP')
		PA.append(float(line[27:41])) 
#	PA.append(float(line[0:24]))
#        EHSS.append(float(line[39:]))

    # Report results on screen, and in the summary file.
    #print
    #os.remove('./MFJ + ".CCS.out"')
    SummaryFile.write(MFJ +' ')
    average =PA
    print "PA:   " +str(average) + '\n'
    SummaryFile.write(str(average) + ',')
#    average =sum(EHSS)/len(EHSS)
#    print "EHSS: " +str(average) + '\n'
    SummaryFile.write(str(average) + ',\n')

SummaryFile.flush()
SummaryFile.close()

