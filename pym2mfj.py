pym2mfj.py

# This script converts pym files into mfj for use by mobcal.py and/ or scoringXlinks.py

#import wx
import re

# enetr the number of input pym structures that you want to covert into mfj
Nst=raw_input('enter number of structures: ')
Nstrs=int(Nst)

# Enter the overall umber of spheres in each structure (that can be the number of eresidues, atoms. subunits, domains)
Ns=raw_input("enter number of spheres: ")
Nsphs=str(Ns)

# read input files 
for i in range(Nstrs):
        #InputFileName='cC.pym'
        InputFileName="config."+str(i)+".pym"
	InputFile = open(InputFileName,'r') # Open file to read
	InputFileLines = InputFile.read().splitlines() # string of lines from pdb file

# The output filenames (*mfj) to save
	Fname = 'outConfig.'+str(i)+'.mfj'

        f1 = open(Fname, 'w')
        #f1 = open(Fname, 'w')
        f1.write(Fname+'\n')
        f1.write('1'+'\n')
        f1.write(Nsphs+'\n')
        f1.write('ang'+'\n')
        f1.write('none'+'\n')
        f1.write('1.0000'+'\n')

# do the actual job: to convert the pym files into mfj
	for line in InputFileLines:
		if line.find('SPHERE')>=0:
                        seq = line[:]
                        Newseq = seq.replace(",", " ")

                        l = re.findall ("(\S+)", Newseq) 
                        x = float(l[1])
                        y = float(l[2])
                        z = float(l[3])
                        r = float(l[4])
                        #f1.write('REMARK'+'\n')
                        f1.write("%s%8.3f%s%8.3f%s%8.3f%7.1f%s" % ('    ',x,'      ',y,'      ',z,r,'      .00000')+'\n')
        
        f1.close()

