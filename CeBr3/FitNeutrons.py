#                                                                                            #
# ========================================================================================== #
# ========================================================================================== #
# ================== Python Based Backing Detector Calibration - CeBr3 ===================== #
# ========================================================================================== #
# ============================ Written by C. Awe - April 2021 ============================== #
# ========================================================================================== #
# ========================================================================================== #
#                                                                                            #
#                         Code to compute CeBr3 quenching factors.                           #
#                                                                                            #
#                                Based on code by S. Hedges                                  #

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

ls_channelList = [8,9,10,11,12,13,14,15] #Converts BD number to channel number.
bd_slopeList = [18.15,18.40,16.51,17.53,16.91,18.20,17.56,18.15] #Calibration constants for each BD.
bd_offsetList = [-1.42,10.17,-1.82,8.14,-5.81,11.95,-3.96,4.13] #Calibration offsets for each BD.

lowerEnergyBound = 0 #keVee
upperEnergyBound = 1000 #keVee
energyVar = ROOT.RooRealVar( "energyVar","energyVar", lowerEnergyBound, upperEnergyBound )

def fitNeutrons( channel, folder ):

	#Specify where to write plots.
	plotPath = folder

	recoilChannel = channel

	calibrationSources = ["neutrons"]

	#Specify data tree name and branches.
	dataTreeName = "analysisTree"
	dataTreeChannelBranchName = "LS_channel"
	dataTreeChannelBranchType = 'i'
	dataTreeIntegralBranchName = "LS_integral"
	dataTreeIntegralBranchType = 'd'
  	dataTreePSDBranchName = "LS_psd"
	dataTreePSDBranchType = 'd'
  	dataTreeTimingBranchName = "LS_timeToBPM"
	dataTreeTiminglBranchType = 'd'
  	dataTreeSignalBranchName = "scatterer_integral"
	dataTreeSignalBranchType = 'd'
  	dataTreeBgndBranchName = "scatterer_noise"
	dataTreeBgndBranchType = 'd'
  
  	#Cuts to be applied to data tree.
  	time_Low = 306
  	time_High = 336
  	psd_Low = 0.26
  	psd_High = 0.6
  	integral_Low = 10000
  	integral_High = 35000

  	#Specify sim tree name and branches.
	simTreeName = "totalEnergyTree"
	simTreeEnergyBranchName = "energy"
	simTreeEnergyBranchType = 'd'
	simTreeChannelBranchName = "cellNum"
	simTreeChannelBranchType = 'd'
  	simTreeChannelBranchName = "tof"
	simTreeChannelBranchType = 'd'


	#################
	## Emcee Setup ##
	#################

	ndim = 5 #Alpha, Beta, Gamma, QF, Offset
	for source in calibrationSources:
		ndim += 1 #Add an amplitude for each source.

	#Parameter names for plots.
	labels = ['alpha','beta','gamma','qf','offset']
	for source in calibrationSources:
		labels.append("amp_"+source)
	
	#Minimum and maximum values for the parameters.
	mins = [0,0,0,0.01,-200]
	maxes = [7.5,3,1,20]
	for source in calibrationSources:
		mins.append(0.01)
		maxes.append(1.01)

	nwalkers = 100 #Based on Sam's code.
	nBurnInSteps = 500
	nSteps = 3500
	
	#########################
	## Simulation Settings ##
	#########################
	
	valsToGen = 10 #How many smeared data points to generate for each real one.
	bgndPdfType = "keys" #Either keys or binned, keys ensures no zero bins.
	nBins = 99 #specify binning if bgndPdfType is "binned".

	###############################
	## Specify Import/Fit Ranges ##
	###############################
	fitRangeMin = 1
	fitRangeMax = 20

	
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
	bdNum = int(channel)
	dataTreeChannelNum = ls_channelList[ bdNum - 1 ]
	sourceFilenames = []
	simFilenames = []
	for source in calibrationSources:
		sourceFilenames.append( folder + source + ".txt" )
		print( "Source file is: " + sourceFilenames[-1] )
		simFilenames.append( folder + source + "-" + str(recoilChannel) + "-sim.root" )
		print( "Simulation file is: " + simFilenames[-1] )
		
	dataTreeIntegral = array(dataTreeIntegralBranchType,[0])
	dataTreeChannel = array(dataTreeChannelBranchType,[0])
	dataTreePSD = array(dataTreePSDBranchType,[0])
	dataTreeSignal = array(dataTreeSignalBranchType,[0])
	DataTreeTiming = array(dataTreeTimingBranchType,[0])
	
	########################
	## Set up observables ##
	########################
	global energyVar
	binning = ROOT.RooBinning( nBins, lowerEnergyBound, upperEnergyBound, "binning" )
	energyVar.setBinning( binning )
	energyVar.setRange( "fitRange", fitRangeMin, fitRangeMax )
	argSet = ROOT.RooArgSet( energyVar )
	argList = ROOT.RooArgList( energyVar )
	#Cut used to determine amplitudes.
	cut = "energyVar >= " + str( fitRangeMin ) + " && energyVar <= " + str( fitRangeMax )
	
	######################
	## Load Source Data ##
	######################
	adc_to_keV = 9.65e-4
	expectedSourceCounts=[]
	scaledBgndEntriesVars=[]
	sourceDataSets = []
	sourceDataHists = []
	sourceCountsInFitRanges = []
	srcDataHistPdfs = []
	bgndDataSet = ROOT.RooDataSet( "bgndDataSet_"+source, "bgndDataSet_"+source, argSet ) )
	for sourceNum in range(0,len(calibrationSources)):
		source = calibrationSources[ sourceNum ]
		sourceFile = ROOT.TFile( sourceFilenames[ sourceNum ], "READ" )
		sourceTree = sourceFile.Get( dataTreeName )
		sourceTree.SetBranchAddress( dataTreeIntegralBranchName, dataTreeIntegral )
		sourceTree.SetBranchAddress( dataTreeChannelBranchName, dataTreeChannel )
		sourceTree.SetBranchAddress( dataTreePSDBranchName, dataTreePSD )
		sourceTree.SetBranchAddress( dataTreeSignalBranchName, dataTreeSignal )
		sourceTree.SetBranchAddress( dataTreeTimingBranchName, dataTreeTiming )
		sourceDataSets.append(ROOT.RooDataSet( "sourceDataSet_"+source, "sourceDataSet_"+source, argSet ) )
		nSourceEntries = sourceTree.GetEntries()
		print( "Found " + str( nSourceEntries ) + " entries in " + source + " data set." )
		for entry in range( 0, nSourceEntries ):
			sourceTree.GetEntry(entry)
			if dataTreeChannel[0] == dataTreeChannelNum:
				if ( dataTreePSD[0] >= psd_Low ) and ( dataTreePSD[0] <= psd_High ):
					if ( dataTreeTiming[0] >= time_Low ) and ( dataTreeTiming[0] <= time_High ):
						if ( dataTreeIntegral[0] >= integral_Low ) and ( dataTreeIntegral[0] <= integral_High ):
							if ( dataTreeSignal[0] * adc_to_keV >= lowerEnergyBound ) and ( dataTreeSignal[0] * adc_to_keV <= upperEnergyBound ):
								energyVar.setVal( dataTreeSignal[0] * adc_to_keV )
								sourceDataSets[ sourceNum ].add( argSet )
								energyVar.setVal( dataTreeNoise[0] * adc_to_keV )
								bgndDataSet.add( argSet )
		reducedBgndDataSet = bgndDataSet.reduce(ROOT.RooFit.Cut(cut))
		bgndCountsInFitRange = reducedBgndDataSet.numEntries() )
		print( "Found "+str(bgndCountsInFitRange)+" bgnd entries in the fit range" )
		sourceCountsInFitRanges.append( sourceDataSets[ sourceNum ].numEntries() )
		print( "Found "+str(sourceCountsInFitRanges[-1])+" entries for " + source + " in the fit range" )
		( sourceDataSets[ sourceNum ].get().find("integralVar") ).setBins( nBins )
		sourceDataHists.append( sourceDataSets[ sourceNum ].binnedClone() )
	
					
	#Make a background PDF.
	if bgndPdfType == "binned":
		print( "Making binned bgnd pdf.\n" )
		(bgndDataSet.get().find("integralVar")).setBins(nBins)
		bgndDataHist = bgndDataSet.binnedClone()
		bgndDataPdf = ROOT.RooHistPdf("bgndDataPdf","bgndDataPdf",argSet,bgndDataHist,1) #1 specifies interpolation order.
	else:
		print( "Making RooKeys bgnd pdf.\n" )
		bgndDataPdf = ROOT.RooKeysPdf("bgndDataPdf","bgndDataPdf",integralVar,bgndDataSet,ROOT.RooKeysPdf.NoMirror,1.2)
	
		
	##########################
	## Load Simulation Data ##
	##########################
	simTreeChannelNum = bd_cellList[ bdNum - 1 ]
	simTreeEnergy = array( simTreeEnergyBranchType, [0] )
	simTreeChannel = array( simTreeChannelBranchType, [0] )
	simDataSets = []
	simArrays = []
	for sourceNum in range(0,len(calibrationSources)):
		source = calibrationSources[ sourceNum ]
		simFile = ROOT.TFile( simFilenames[ sourceNum ], "READ" )
		simTree = simFile.Get( simTreeName )
		simTree.SetBranchAddress( simTreeEnergyBranchName, simTreeEnergy )
		simTree.SetBranchAddress( simTreeChannelBranchName, simTreeChannel )
		simDataSets.append(ROOT.RooDataSet( "simDataSet_"+source, "simDataSet_"+source, argSet ) )
		nSimEntries = simTree.GetEntries()
		print( "Found " + str( nSimEntries ) + " entries in " + source + " sim data set." )
		#Numpy way
		simList = []
		for entry in range( 0, nSimEntries ):
			simTree.GetEntry(entry)
			if simTreeChannel[0] == simTreeChannelNum:
				if ( simTreeEnergy[0] >= lowerEnergyBound ) and ( simTreeEnergy[0] <= upperEnergyBound ):
					simList.append( simTreeEnergy[0] )
		simArrays.append(numpy.array( simList ))
		
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
	psize = pos_max - pos_min
	
	#Generate random values within those spaces for each walker. 
	pos = [ pos_min + psize * numpy.random.rand( ndim ) for i in range( nwalkers ) ]
	
	#Define Emcee functions.
	f = lambda x,alpha,beta: numpy.sqrt( numpy.power(alpha*x,2) + numpy.power(beta,2)*x )
	
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
	#This is our RooFit function that returns the POSITIVE log likelihood
	def lnlike(theta):
 	 	#energyVar gets modified in loop. So far hasn't caused any issues with multiprocessing
		global integralVar

  		# Reset nllVal to zero because we'll add to this for each source we generate
  		# an nll value for
		nllVal=0
  
  		#Load theta values into parameters
		alpha=theta[0]
		beta=theta[1]
		slope=theta[2]
		offset=theta[3]
  
		#Load remaining theta values for each source, generate smeared pdf, calculate 
 		#nll value
		for sourceNum in range(0,len(calibrationSources)):
    
			amplitude=theta[4+sourceNum]

    			#Make ampltitude var 
			sourceCountsVar = ROOT.RooRealVar("sourceCountsVar","sourceCountsVar",amplitude*expectedSourceCounts[sourceNum])
			sourceCountsVar.setConstant(1)
    
			#Make array of sigmas from res params and sim values
			sigmas=f(simArrays[sourceNum],alpha,beta)
    
			#Generate valsToGen random values for every entry in sim
			smearedArrayList=[]
			for i in range(0,valsToGen):
				smearedArrayList.append(slope*(numpy.random.normal(simArrays[sourceNum],sigmas)-offset))
    
			smearedArrayArray=numpy.array(smearedArrayList)
			flatArray=smearedArrayArray.flatten()
    
			#Make smeared data set
			smearedSimDataSet=ROOT.RooDataSet("smearedSimDataSet","smearedSimDataSet",argSet)
    
			#Numpy array->TH1->RooDataHist->RooHistPdf
			#~0.03 seconds, much faster than iterating through array to fill
			w=numpy.full(flatArray.size,1.)
			h=ROOT.TH1D("h","",nBins,lowerEnergyBound,upperEnergyBound)
			h.FillN(flatArray.size,flatArray,w)
			smearedSimDataHist=ROOT.RooDataHist("smearedSimDataHist","smearedSimDataHist",argList,h)
			simPdf = ROOT.RooHistPdf("simPdf","simPdf",argSet,smearedSimDataHist,0) #1 specifies interpolation order
			h.Delete()
			del h
    
			##Make Model
			pdfList = ROOT.RooArgList(bgndDataPdf,simPdf)
			ampList = ROOT.RooArgList(scaledBgndEntriesVars[sourceNum],sourceCountsVar)
			model = ROOT.RooAddPdf("model","model",pdfList,ampList)
			model.fixCoefRange("fitRange")
    
			#Compute nll
			nll = model.createNLL(sourceDataHists[sourceNum],
			 ROOT.RooFit.Extended(1),
			 ROOT.RooFit.Verbose(0),
			 ROOT.RooFit.Range("fitRange"),
			 ROOT.RooFit.SumCoefRange("fitRange"),
			 ROOT.RooFit.NumCPU(1)
			)

			#Make NLL positive
			nllVal += (-1*nll.getVal())
    
			if (PLOT==1):
				try:
					c1
				except NameError:
					c1=ROOT.TCanvas("c1","c1")
      
				#Reduce integrator for plotting, massively speeds plotting
				ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-5)
				ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-5)
				frame = integralVar.frame(lowerIntegralBound,upperIntegralBound,nBins)
				frame.SetTitle(calibrationSources[sourceNum]+": alpha="+str(alpha)+", beta="+str(beta)+",  slope="+str(slope)+", offset="+str(offset))
      
				#Plot source data
				sourceDataSets[sourceNum].plotOn(
				 frame,
				 ROOT.RooFit.Name("Source"),
				 ROOT.RooFit.MarkerColor(1),
				 ROOT.RooFit.FillColor(0)
				)
				#Plot components
				model.plotOn(
				 frame,
				 ROOT.RooFit.Name("Bgnd"),
				 ROOT.RooFit.Components("bgndDataPdf"),
				 ROOT.RooFit.LineColor(ROOT.kSolid),
				 ROOT.RooFit.FillColor(0),
				 ROOT.RooFit.ProjWData(sourceDataSets[sourceNum])
				)
				model.plotOn(
				 frame,ROOT.RooFit.Name("Sim"),
				 ROOT.RooFit.Components("simPdf"),
				 ROOT.RooFit.LineColor(ROOT.kRed),
				 ROOT.RooFit.FillColor(0),
				 ROOT.RooFit.ProjWData(sourceDataSets[sourceNum]),
				 ROOT.RooFit.AddTo("Bgnd")
				)
      
				#Draw
				frame.Draw()
      
				#Add legend
				leg = ROOT.TLegend(0.65,0.65,0.95,0.92); 
				sourceObj = frame.findObject("Source");
				bgndObj = frame.findObject("Bgnd");
				simObj = frame.findObject("Sim");
				leg.AddEntry(sourceObj,"Source","P")
				leg.AddEntry(bgndObj,"Background","L")
				leg.AddEntry(simObj,"Sim","L")
				leg.Draw("same")
      
				#Draw
				c1.SetLogy()
				c1.Modified()
				c1.Update()
				c1.SaveAs(plotPath+"bestFit_simultaneous_"+calibrationSources[sourceNum]+"_ch"+str(calibrationChannel)+".pdf")
      
				#Reset integrator for step size
				ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-7)
				ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-7)
      
				#Memory management for plotting
				frame.Delete()
				#del frame()
      
			#Memory management
			smearedSimDataSet.reset()
			#smearedSimDataSet.Delete()
			del smearedSimDataSet
			#pdfList.Delete()
			del pdfList
			#ampList.Delete()
			del ampList
			#sourceCountsVar.Delete()
			del sourceCountsVar
			#smearedSimDataHist.Delete()
			del smearedSimDataHist
			#simPdf.Delete()
			del simPdf
			del model
			nll.Delete()
			del nll
			#gc.collect()
    
    
		#Return total nllval from all sources
		return nllVal
  
	#Calls out roofit function, makes sure output is not infinite and parameters in allowed range 
	def lnprob(theta):
		lp = lnprior(theta)
		if not numpy.isfinite(lp):
			return -numpy.inf
		return lp + lnlike(theta)



##########################MC RUN STARTS HERE####################################

	#Single-threaded
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
	print("Starting burn in...")
	pos, prob, state  = sampler.run_mcmc(pos, nBurnInSteps, progress=True)
	sampler.reset()
	print("Burn-in complete!")
	pos, prob, state  = sampler.run_mcmc(pos, nSteps, progress=True)

	#Parallel processing
	#with Pool() as pool:
 	#sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,pool=pool)
  	#Burn in
  	#print("Starting burn in...")
  	#pos, prob, state  = sampler.run_mcmc(pos, nBurnInSteps, progress=True)
  	#print("Burn-in complete! Mean acceptance fraction: {0:.3f}".format(numpy.mean(sampler.acceptance_fraction)))
 	#sampler.reset()
 	# pos, prob, state  = sampler.run_mcmc(pos,nSteps,progress=True)

###########################MC PLOTTING HERE#####################################
#My computer doesn't have python-tk, so I can't view plots and had to save them
#as a PDF to look at the results 
	matplotlib.use('PDF')

	#GET THE LL VALUES--From grayson's code--do this first in case plotting fails.
	samples=sampler.flatchain
	lnprobs = sampler.lnprobability[:,:]
	flatLnprobs = lnprobs.reshape(-1)
	with open(plotPath+"sampler_simultaneous_ch"+str(calibrationChannel)+".csv", 'w') as sampleOutFile:
		theWriter = csvlib.writer(sampleOutFile, delimiter=',')
		for sampleLine, llval in zip(samples, flatLnprobs):
			theWriter.writerow(numpy.append(sampleLine,llval))

	#MAKE TRACE PLOTS OF EACH PARAMATER
	fig = plt.figure(figsize=(10,ndim*2))
	gs = fig.add_gridspec(ndim,1)
	plt.subplots_adjust(hspace=0.4)
	for i in range(0,ndim):
		axes = fig.add_subplot(gs[i,:])
		axes.plot(sampler.chain[:,:,i].T, '-', color='k', alpha=0.3)
		axes.set_title(labels[i])
	plt.savefig(plotPath+"traceplots_simultaneous_ch"+str(calibrationChannel)+".pdf")

	#CORNER PLOT HERE 
	samples=sampler.flatchain
	fig = corner.corner(samples, labels=labels, ranges=ranges, quantiles=[0.16,0.5,0.84],show_titles=True,title_kwargs={'fontsize':12})
	fig.savefig(plotPath+"corner_simultaneous_ch"+str(calibrationChannel)+".pdf")

	#CALCULATE QUANTILES HERE
	bestFitValues=[]
	for i in range(ndim):
		mcmc = numpy.percentile(sampler.chain[:,:,i],[16, 50, 84])
		q = numpy.diff(mcmc)
		print(labels[i]+": "+str(mcmc[1])+"+"+str(q[0])+" -"+str(q[1])+"\n")
		bestFitValues.append(mcmc[1])

	#Plot best fit 
	PLOT=1
	lnlike(bestFitValues) #Unpack list to arguments for our call to lnlike

	#Print out stats
	print("Mean acceptance fraction: {0:.3f}".format(numpy.mean(sampler.acceptance_fraction))+"\n")
	print("Mean autocorrelation time: {0:.3f} steps".format(numpy.mean(sampler.get_autocorr_time(c=1,quiet=True))))

#Step through all BDs, invoking fitBDs at each step.
def main():

	#Get the location of the data and sims.
	calibrationFolder = sys.argv[1]
	if not calibrationFolder.endswith('/'):
		calibrationFolder += "/"
	#Step through backing detectors.
	for i in range(1,8):
		calibrationChannel = i
		fitBD( calibrationChannel, calibrationFolder )

#Execute main function 	
if __name__== "__main__":
  main()		
