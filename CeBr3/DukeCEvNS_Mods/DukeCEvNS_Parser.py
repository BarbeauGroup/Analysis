#                                                                         #
# ======================================================================= #
# ========================================================================#
# ================== Python Based Duke CEvNS Parser ===================== #
# ======================================================================= #
# ================== Written by C. Awe - April 2021 ===================== #
# ======================================================================= #
# ======================================================================= #
#                                                                         #
#              Code to compute CeBr3 quenching factors.                   #
#                                                                         #

import sys
import ROOT
import random
from ROOT import TGraph, TH1D, TH1F, TH2F, TCanvas, TFile, TTree
import numpy
import scipy as sp
import matplotlib
from array import array
from tqdm import tqdm
import emcee
import corner
import csv as csvlib
import os
from matplotlib import pyplot as plt

#gROOT.SetStyle("Plain")

#Check that the user supplied a filename.
if len(sys.argv) != 3:
	print( "\nERROR: Incorrect call." )
	print( "Correct usage is python DukeCEvNS_Parser.py <filename> <output file>\n" )
	
#Open the file and retrieve the data.
filename = sys.argv[1]
file = open( filename, 'r' )
text = file.read()
file.close()

#Open the output file.
outputname = sys.argv[2]
outFile = ROOT.TFile( outputname, "RECREATE" )

#Create a TTree.
cevnsTree = ROOT.TTree("cevnsTree","Duke CEvNS Output")
cevnsHist = ROOT.TH1F("cevsnHist","CEvNS Spectrum", 100, 0, 100)

#Set up branches for the TTree.
energy = array('d',[0])
counts = array('d',[0])
cevnsTree.Branch('energy',energy,'energy/D')
cevnsTree.Branch('counts',counts,'counts/D')


#Split text into lines.
lines = text.split('\n')

print( "Parsing..." )

#Step through filling the TTree and Histogram.
for index, line in enumerate(lines[0:len(lines)-1]):

	args = line.split() #Split on whitespace.
	energy[0] = float(args[0]) * 1000 #Convert to keVee.
	counts[0] = float(args[1])
	cevnsTree.Fill()
	cevnsHist.Fill(energy[0],counts[0])
	
print( "Done!" )
print( outputname + " created." )

#Write the output.
cevnsTree.Write()
cevnsHist.Write()
 
#Close the output file.
outFile.Close()















