#                                                                                            #
# ========================================================================================== #
# ========================================================================================== #
# ================== Python Based Backing Detector Calibration - CeBr3 ===================== #
# ========================================================================================== #
# ============================ Written by C. Awe - March 2021 ============================== #
# ========================================================================================== #
# ========================================================================================== #
#                                                                                            #
#                        Code to calibrate CeBr3 Backing Detectors.                          #
#                                                                                            #
#                          Based on conversations with S. Hedges                             #

import ROOT
import random
from ROOT import TGraph, TH1D, TH1F, TH2F, TCanvas, TFile, TTree
import numpy as np
import scipy as sp
import matplotlib as mpl
from array import array
from tqdm import tqdm
import emcee

ls_channelList = [8,9,10,11,12,13,14,15] #Converts BD number to channel number.
bd_cellList = [306,314,322,330,326,318,310,302] #Converts BD number to MCNP cell number.

def main():

	####################
	## Get User Input ##
	####################

	calibrationFolder = sys.argv[1]
	if not calibrationFolder.endswith('/'):
		calibrationFolder += "/"
	#Specify where to write plots.
	plotPath = calibrationFolder

	calibrationChannel = sys.argv[2]

	calibrationSources = ["na22","co60"]

	#Specify data tree name.
	dataTreeName = "analysisTree"
	dataTreeChannelBranchName = "LS_channel"
	dataTreeChannelBranchType = 'i'
	dataTreeIntegralBranchName = "LS_integral"
	dataTreeIntegralBranchType = 'd'

	simTreeName = "totalEnergyTree"
	simTreeEnergyBranchName = "energy"
	simTreeEnergyBranchType = 'd'
	simTreeChannelBranchName = "cellNum"
	simTreeChannelBranchType = 'd'


	#################
	## Emcee Setup ##
	#################

	ndim = 5 #Alpha, Beta, Gamma, Slope, Offset
	for source in calibrationSources:
		ndim += 1 #Add an amplitude for each source.

	#Parameter names for plots.
	labels = ['alpha','beta','gamma','slope','offset']
	for source in calibrationSources:
		labels.append("amp_"+source)
	
	#Minimum and maximum values for the parameters.
	mins = [5.5,0.3,0.0,20,-10]
	maxes = [7.5,0.7,0.5,30,10]
	for source in calibrationSources:
		mins.append(0.90)
		maxes.append(1.01)

	nwalkers = 500 #Based on Sam's code.
	nBurnInSteps = 100
	nSteps = 1000
	
	#########################
	## Simulation Settings ##
	#########################
	
	valsToGen = 10 #How many smeared data points to generate for each real one.
	bgndPdfType = "keys" #Either keys or binned, keys ensures no zero bins.
	nBins = 230 #specify binning if bgndPdfType is "binned".

	###############################
	## Specify Import/Fit Ranges ##
	###############################
	fitRangeMin = 2000
	fitRangeMax = 25000
	lowerIntegralBound = 1
	upperIntegralBound = 30000
	lowerEnergyBound = 1
	UpperEnergyBound = 200000
	
	#####################
	## RooFit Settings ##
	#####################
	ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
	#Set the number of steps used for integration.
	ROOT.RooAbsReal.defaultIntegratorConfig().getConfigSection("RooIntegrator1D").setRealValue("maxSteps",30);
	#Set integration step sizes.
	ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-7)
	ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-7)

	###############################
	## Specify Data/Format Names ##
	###############################
	bdNum = int(calibrationChannel)
	dataTreeChannelNum = ls_channelList[ bdNum - 1 ]
	bgndFilename = calibrationFolder + "bgnd.root"
	sourceFilenames = []
	simFilenames = []
	for source in calibrationSources:
		sourceFilenames.append( calibrationFolder + "/" + source + ".root" )
		print( "Source file is: " + sourceFilenames[-1] )
		simFilenames.append( calibrationFolder + "/" + source + "_sim.root" )
		print( "Simulation file is: " + simFilenames[-1] )
		
	dataTreeIntegral = array.array(dataTreeIntegralBranchType,[0])
	dataTreeChannel = array.array(dataTreeChannelBranchType,[0])
	
	########################
	## Set up observables ##
	########################
	integralVar = ROOT.RooRealVar( "integralVar","integralVar", lowerIntegralBound, upperIntegralBound )
	binning = ROOT.RooBinning( nBins, lowerIntegralBound, upperIntegralBound, "binning" )
	integralVar.setBinning( "binning" )
	integralVar.setRange( "fitRange", fitRangeMin, fitRangeMax )
	argSet = ROOT.RooArgSet( integralVar )
	argList = ROOT.RooArgList( integralVar )
	#Cut used to determine amplitudes.
	cut = "integralVar >= " + str( fitRangeMin ) + " && integralVar <= " + str( fitRangeMax )
	
	##########################
	## Load Background Data ##
	##########################
	bgndFile = ROOT.TFile( bgndFilename, "READ" )
	bgndTree = bgndFile.Get( dataTreeName )
	bgndTree.SetBranchAddress( dataTreeIntegralBranchName, dataTreeIntegral )
	bgndTree.SetBranchAddress( dataTreeChannelBranchName, dataTreeChannel )
	bgndDataSet = ROOT.RooDataSet( "bgndDataSet", "bgndDataSet", argSet )
	nBgndEntries = bgndTree.GetEntries()
	print( "\nFound "+str(nBgndEntries)+" in bgnd tree" )
	for entry in range( 0, nBgndEntries ):
			bgndTree.GetEntry(entry)
			if dataTreeChannel[0] == dataTreeChannelNum:
				if ( dataTreeIntegral[0] >= lowerIntegralBound ) and ( dataTreeIntegral[0] <= upperIntegralBound ):
					integralVar.setVal( dataTreeIntegral[0] )
					bgndDataSet.add( argSet )
					
	#Get counts in fit range.				
	reducedBgndDataSet = bgndDataSet.reduce(ROOT.RooFit.Cut(cut))
	bgndCountsInFitRange = reducedBgndDataSet.numEntries()
	print( "Found "+str(bgndCountsInFitRange)+" bgnd entries in the fit range" )
	
	#Make a background PDF.
	if bgndPdfType == "binned":
		print( "Making binned bgnd pdf.\n" )
		(bgndDataSet.get().find("integralVar")).setBins(nBins)
		bgndDataHist = bgndDataSet.binnedClone()
		bgndDataPdf = ROOT.RooHistPdf("bgndDataPdf","bgndDataPdf",argSet,bgndDataHist,1) #1 specifies interpolation order.
	else:
		print( "Making RooKeys bgnd pdf.\n" )
		bgndDataPdf = ROOT.RooKeysPdf("bgndDataPdf","bgndDataPdf",integralVar,bgndDataSet,ROOT.RooKeysPdf.NoMirror,1.2)
	
	######################
	## Load Source Data ##
	######################
	sourceDataSets = []
	sourceDataHists = []
	sourceCountsInFitRanges = []
	srcDataHistPdfs = []
	for sourceNum in range(0,len(calibrationSources)):
		source = calibraitonSources[ sourceNum ]
		sourceFile = ROOT.TFile( sourceFilenames[ sourceNum ], "READ" )
		sourceTree = sourceFile.Get( dataTreeName )
		sourceTree.SetBranchAddress( dataTreeIntegralBranchName, dataTreeIntegral )
		sourceTree.SetBranchAddress( dataTreeChannelBranchName, dataTreeChannel )
		sourceDataSets.append(ROOT.RooDataSet( "sourceDataSet_"+source, "sourceDataSet_"+source, argSet )
		nSourceEntries = sourceTree.GetEntries()
		print( "Found " + str( nSourceEntries ) + " entries in " + source + " data set." )
		for entry in range( 0, nSourceEntries ):
			sourceTree.GetEntry(entry)
			if dataTreeChannel[0] == dataTreeChannelNum:
				if ( dataTreeIntegral[0] >= lowerIntegralBound ) and ( dataTreeIntegral[0] <= upperIntegralBound ):
					integralVar.setVal( dataTreeIntegral[0] )
					sourceDataSets[ sourceNum ].add( argSet )
		sourceCountsInFitRanges.append( sourceDataSets[ sourceNum ].numEntries() )
		print( "Found "+str(sourceCountsInFitRanges[-1])+" entries for " + source + " in the fit range" )
		( sourceDataSets[ sourceNum ].get().find("integralVar") ).setBins( nBins )
		sourceDataHists.append( sourceDataSets[ sourceNum ].binnedClone() )
		
	##########################
	## Load Simulation Data ##
	##########################
	simTreeChannelNum = bd_cellList[ bdNum - 1 ]
	simTreeEnergy = array.array( simTreeEnergyBranchType, [0] )
	simTreeChannel = array.array( simTreeChannelBranchType, [0] )
	simDataSets = []
	simArrays = []
	for sourceNum in range(0,len(calibrationSources)):
		source = calibraitonSources[ sourceNum ]
		simFile = ROOT.TFile( simFilenames[ sourceNum ], "READ" )
		simTree = simFile.Get( simTreeName )
		simTree.SetBranchAddress( simTreeEnergyBranchName, simTreeEnergy )
		simTree.SetBranchAddress( simTreeChannelBranchName, simTreeChannel )
		simDataSets.append(ROOT.RooDataSet( "simDataSet_"+source, "simDataSet_"+source, argSet )
		nSimEntries = simTree.GetEntries()
		print( "Found " + str( nSimEntries ) + " entries in " + source + " sim data set." )
		#Numpy way
		simList = []
		for entry in range( 0, nSimEntries ):
			simTree.GetEntry(entry)
			if simTreeChannel[0] == simTreeChannelNum:
				if ( simTreeEnergy[0] >= lowerEnergyBound ) and ( simTreeEnergy[0] <= upperEnergyBound ):
					simList.append( simTreeEnergy[0] )
		simArrays.append(numpy.array( simList )
		
	#################
	## Emcee Setup ##
	#################
	PLOT = 0
	#Convert limits into numpy arrays.
	pos_min = numpy.array( mins )
	pos_max = numpy.array( maxes )
	
	#Make ranges for corner plots.
	ranges = []
	for i in range( 0, len( mins ) ):
		entry = []
		entry.append( mins[i] )
		entry.append( maxes[i] )
		ranges.append( entry )
		
	#Size of each parameter space.
	psize = mos_max - pos_min
	
	#Generate random values within those spaces for each walker. 
	pos = [ pos_min + psize * numpy.random.rand( ndim ) for i in range( nwalkers ) ]
	
	#Define Emcee functions.
	f = lambda x,alpha,beta,gamma: numpy.sqrt( numpy.power(alpha*x,2) + numpy.power(beta,2)*x + numpy.power(gamma,2) )
	
	#Returns 0 if all parameters are in their allowed range, otherwise returns -infinity.
	def lnprior( theta ):
		allParsInRange = 1
		for i in range( 0, len( theta ) ):
			if not mins[i] < theta[i] < maxes[i]:
				allParsInRange = 0
		
		if allParsInRange == 1:
			return 0.0
		else:
			return -numpy.inf
	
	#RooFit function that does positive log likelihood.
	
	

#Execute main function 	
if __name__== "__main__":
  main()		
