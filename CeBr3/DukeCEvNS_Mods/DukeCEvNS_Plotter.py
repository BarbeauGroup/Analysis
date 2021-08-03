#                                                                         #
# ======================================================================= #
# ======================================================================= #
# ================= Python Based Duke CEvNS Plotter ===================== #
# ======================================================================= #
# ================= Written by C. Awe - August 2021 ===================== #
# ======================================================================= #
# ======================================================================= #
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

#Check that the user supplied a filename.
if len(sys.argv) != 3:
	print( "\nERROR: Incorrect call." )
	print( "Correct usage is python DukeCEvNS_Plotter.py <Path to parsed Duke CEvNS Files>  <Base name for output files>" )
	
#Open the folder containing all three (parsed) Duke CEvNS Sims
path = sys.argv[1]
avgFilename = path + "SNS_AvgQF_Quenched.root"
ceFilename = path + "SNS_CeQF_Quenched.root"
brFilename = path + "SNS_BrQF_Quenched.root"
avgFile = ROOT.TFile( avgFilename )
ceFile = ROOT.TFile( ceFilename )
brFile = ROOT.TFile( brFilename )
avgHist = avgFile.Get("cevnsHist")
ceHist = ceFile.Get("cevnsHist")
brHist = brFile.Get("cevnsHist")

#Make a canvas to plot on.
c1=ROOT.TCanvas("c1","c1")
avgHist.GetXaxis().SetTitle("Recoil Energy [keVee]")
avgHist.GetXaxis().CenterTitle()
avgHist.GetYaxis().SetTitle("\\mbox{Counts / 100 eV } $\cdot$ \\mbox{kg} $\cdot$ \\mbox{year}")
avgHist.GetYaxis().CenterTitle()
avgHist.SetLineWidth(3)
avgHist.SetLineColor(46)
avgHist.SetTitle("SNS CEvNS Spectra (CeBr_{3})")
avgHist.SetStats(0)
ceHist.GetXaxis().SetTitle("Recoil Energy [keVee]")
ceHist.GetXaxis().CenterTitle()
ceHist.GetYaxis().SetTitle("\\mbox{Counts / 100 eV } $\cdot$ \\mbox{kg} $\cdot$ \\mbox{year}")
ceHist.GetYaxis().CenterTitle()
ceHist.SetLineWidth(3)
ceHist.SetLineColor(30)
ceHist.SetTitle("SNS CEvNS Spectra (CeBr_{3})")
ceHist.SetStats(0)
brHist.GetXaxis().SetTitle("Recoil Energy [keVee]")
brHist.GetXaxis().CenterTitle()
brHist.GetYaxis().SetTitle("\\mbox{Counts / 100 eV } $\cdot$ \\mbox{kg} $\cdot$ \\mbox{year}")
brHist.GetYaxis().CenterTitle()
brHist.SetLineWidth(3)
brHist.SetLineColor(4)
brHist.SetTitle("SNS CEvNS Spectra (CeBr_{3})")
brHist.SetStats(0)

#Build a line to show threshold.
thresh = ROOT.TLine(0.55,0,0.55,40)
thresh.SetLineStyle(9)
thresh.SetLineWidth(3)
thresh.SetLineColor(12)

#Build a legend.
leg = ROOT.TLegend(0.5,0.55,0.87,0.85); 
leg.AddEntry(avgHist,"Quenching Factors derived from both nuclei.","l")
leg.AddEntry(ceHist,"Quenching Factors derived from cerium recoils.","l")
leg.AddEntry(brHist,"Quenching Factors derived from bromine recoils.","l")
leg.AddEntry(thresh,"10 PE Threshold","l")
#gStyle.SetLegendTextSize(3.);

#Draw.
c1.cd()
brHist.Draw("hist")
ceHist.Draw("hist,same")
avgHist.Draw("hist,same")
leg.Draw()
thresh.Draw()
c1.Modified()
c1.Update()

#Save the output. 
outputFilename = path + sys.argv[2]
c1.SaveAs(outputFilename+".png")
rootFile = ROOT.TFile( outputFilename + ".root", "recreate" )
rootFile.Write()

#Close stuff.
avgFile.Close()
ceFile.Close()
brFile.Close()
rootFile.Close()







