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
bd_cellList = [306,314,322,330,326,318,310,302] #Converts BD number to MCNP cell number.
bd_slopeList = [18.15,18.40,16.51,17.53,16.91,18.20,17.56,18.15] #Calibration constants for each BD.
bd_offsetList = [-1.42,10.17,-1.82,8.14,-5.81,11.95,-3.96,4.13] #Calibration offsets for each BD.
time_lowList = [321,314,308,306,308,311,317,325] #Timing cut lower bounds for each BD in ns.
time_HighList = [351,344,338,336,338,341,347,355] #Timing cut high bounds for each BD in ns.
bd_timingAdjustments = [315.18,305.68,299.7,328,299.34,303.18,309.05,317.61] #Adjusts from simulated tof to "real" tof. Temporary. 

lowerEnergyBound = 0 #keVee
upperEnergyBound = 200 #keVee
energyVar = ROOT.RooRealVar( "energyVar","energyVar", lowerEnergyBound, upperEnergyBound )

def main():
	
	#Get the location of the data and sims.
	folder = sys.argv[1]
	if not folder.endswith('/'):
		folder += "/"

	#Get the backing detector to calibrate.
	bdNum = sys.argv[2]
	
	#Specify where to write plots.
	plotPath = folder


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
  
  	#Cuts to be applied to the data.
	psd_Low = 0.26
	psd_High = 0.6
	integral_Low = 10000
	integral_High = 35000

  	#Specify sim tree name and branches.
	simTreeName = "totalEnergyTree"
	simTreeEnergyBranchName = "energy"
	simTreeEnergyBranchType = 'd'
	simTreeLS_EnergyBranchName = "ls_energy"
	simTreeLS_EnergyBranchType = 'd'
	simTreeTOFBranchName = "tof"
	simTreeTOFBranchType = 'd'


	#################
	## Emcee Setup ##
	#################

	ndim = 6 #Alpha, Beta, Gamma, Offset, QF, Amplitude.
	#for ch in channels: 
		#ndim += 2 #Add a qeunching factor and amplitude for each source.

	#Parameter names for plots.
	labels = ['alpha','beta','gamma','offset','qf','amplitude']
	#for ch in channels:
		#labels.append("qf_"+str(ch)) #QF for each channel.
		#labels.append("amp_"+str(ch))
	
	#Minimum and maximum values for the parameters.
	mins = [0,0,0,-200,0.001,0.01]
	maxes = [7.5,3,5,20,0.05,1.5]
	#for ch in channels:
		#mins.append(0.01) #QF mins
		#maxes.append(0.3) #QF maxes
		#mins.append(0.01) #Amplitude mins
		#maxes.append(1.5) #Amplitude maxes

	nwalkers = 55 #Based on Sam's code.
	nBurnInSteps = 100
	nSteps = 350
	
	#########################
	## Simulation Settings ##
	#########################
	
	valsToGen = 100 #How many smeared data points to generate for each real one.
	bgndPdfType = "keys" #Either keys or binned, keys ensures no zero bins.
	nBins = 2000 #specify binning if bgndPdfType is "binned".

	###############################
	## Specify Import/Fit Ranges ##
	###############################
	fitRangeMin = 0 #keVee
	fitRangeMax = 12 #keVee

	
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
	sourceFilename = folder + "neutrons.txt"
	simFilenames = folder + str(bdNum) + "-sim.root"
	#for ch in channels:
		#simFilenames.append( folder + str(ch) + "-sim.root" )
		#print( "Simulation file is: " + simFilenames[-1] )
	
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
	#sourceDataSets = []
	#sourceDataHists = []
	#expectedSourceCounts=[]
	#scaledBgndEntriesVars=[]
	#chDataSets = []
	#chDataHists = []
	#bgndCountsInFitRanges = []
	#chDataHistPdfs = []
	#bgndDataSets = []
	#bgndDataPdfs = []
	sourceChain = ROOT.TChain( dataTreeName )
	file = open( sourceFilename, 'r' )
	for line in file:
		print( "Adding file " + line )
		sourceChain.AddFile( line[:-1] )
	print( "TChain Built." )
	#for chNum in range(0,len(channels)):
		
	#print("On BD " +str(chNum)+" of "+str(len(channels)))
	time_Low = time_lowList[int(bdNum) - 1]
	time_High = time_HighList[int(bdNum) - 1]
	#sourceDataSets.append(ROOT.RooDataSet("sourceDataSet_"+str(chNum),"sourceDataSet_"+str(chNum),argSet))
	sourceDataSet = ROOT.RooDataSet("sourceDataSet","sourceDataSet",argSet)
	#bgndDataSets.append(ROOT.RooDataSet( "bgndDataSet_"+str(chNum), "bgndDataSet_"+str(chNum), argSet ))
	bgndDataSet = ROOT.RooDataSet( "bgndDataSet", "bgndDataSet", argSet )
	numEntries = sourceChain.GetEntries()
	count = 0;
	for entry in sourceChain:
		if count % 10000 == 0:
			print("On entry " + str(count) + " of " + str(numEntries))
		if entry.LS_channel == ls_channelList[ int(bdNum) - 1 ]:
			if ( entry.LS_psd >= psd_Low ) and ( entry.LS_psd <= psd_High ):
				if ( entry.LS_timeToBPM >= time_Low ) and ( entry.LS_timeToBPM <= time_High ):
					if ( entry.LS_integral >= integral_Low ) and ( entry.LS_integral <= integral_High ):
						if ( entry.scatterer_integral * adc_to_keV >= lowerEnergyBound ) and ( entry.scatterer_integral * adc_to_keV <= upperEnergyBound ):
							energyVar.setVal( entry.scatterer_integral * adc_to_keV )
							sourceDataSet.add( argSet )
							energyVar.setVal( entry.scatterer_noise * adc_to_keV )
							bgndDataSet.add( argSet )
		count += 1
	#reducedBgndDataSet = bgndDataSets[ chNum ].reduce(ROOT.RooFit.Cut(cut))
	reducedBgndDataSet = bgndDataSet.reduce(ROOT.RooFit.Cut(cut))
	#bgndCountsInFitRanges.append(reducedBgndDataSet.numEntries())
	bgndCountsInFitRange = reducedBgndDataSet.numEntries()
	print( "Found "+str(bgndCountsInFitRange)+" bgnd entries in the fit range" )
	scaledBgndEntries = bgndCountsInFitRange
	scaledBgndEntriesVar = ROOT.RooRealVar("scaledBgndEntriesVar","scaledBgndEntriesVar",scaledBgndEntries)
	scaledBgndEntriesVar.setConstant()
	expectedSourceCounts = sourceDataSet.numEntries()
	print( "Found "+str(expectedSourceCounts)+" entries for BD " + str(bdNum) + " in the fit range" )
	( sourceDataSet.get().find("energyVar") ).setBins( nBins )
	sourceDataHist = sourceDataSet.binnedClone()
	
					
	#Make a background PDF.
	if bgndPdfType == "binned":
		print( "Making binned bgnd pdf.\n" )
		(bgndDataSet.get().find("energyVar")).setBins(nBins)
		bgndDataHist = bgndDataSet.binnedClone()
		bgndDataPdfs = ROOT.RooHistPdf("bgndDataPdf","bgndDataPdf",argSet,bgndDataHist,1) #1 specifies interpolation order.
	else:
		print( "Making RooKeys bgnd pdf.\n" )
		bgndDataPdf = ROOT.RooKeysPdf("bgndDataPdf","bgndDataPdf",energyVar,bgndDataSet,ROOT.RooKeysPdf.NoMirror,1.2)
		
	##########################
	## Load Simulation Data ##
	##########################
	simTreeEnergy = array( simTreeEnergyBranchType, [0] ) #Energy in the CeBr3 crystal.
	simTreeLS_Energy = array( simTreeLS_EnergyBranchType, [0] ) #Energy in the LS Backing Detectors.
	simTreeTOF = array( simTreeTOFBranchType, [0] ) #TOF from creation to BD.
	#simDataSets = []
	simArray = []
	#for chNum in range(0,len(channels)): 
	timingOffset = bd_timingAdjustments[int(bdNum)-1]
	timing_Low = time_lowList[int(bdNum)-1]
	timing_High = time_HighList[int(bdNum)-1]
	simFile = ROOT.TFile( simFilenames, "READ" )
	simTree = simFile.Get( simTreeName )
	simTree.SetBranchAddress( simTreeEnergyBranchName, simTreeEnergy )
	simTree.SetBranchAddress( simTreeTOFBranchName, simTreeTOF )
	simDataSet = ROOT.RooDataSet( "simDataSet", "simDataSet", argSet )
	nSimEntries = simTree.GetEntries()
	print( "Found " + str( nSimEntries ) + " entries in BD " + str(bdNum) + " sim data set." )
	#Numpy way
	simList = []
	for entry in range( 0, nSimEntries, 2 ): #Events are time-ordered and will go scatterer - bd - scatterer - bd etc.
		simTree.GetEntry(entry) #Get the scatterer energy.
		scatterer_Enrg = simTreeEnergy[0]
		simTOF = simTreeTOF[0] + timingOffset
		#print( "TOF is: " + str(simTOF) )
		#if ( simTOF > timing_Low ) and ( simTOF < timing_High ):
			#print("Passed TOF cut.")
			#simTree.GetEntry(entry+1) #Get BD energy.
			#dummyIntegral = (simTreeEnergy[0] * bd_slopeList[chNum-1])+bd_offsetList[chNum-1] #Use BD calibrations to mimic our cuts.
			#if ( dummyIntegral >= integral_Low ) and ( dummyIntegral <= integral_High ):
				#print("Passed integral cut.")
		simList.append( scatterer_Enrg )
		#print("Found a simulated event with energy "+str(scatterer_enrg))
	simArray = numpy.array( simList )
		
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
	#This is our RooFit function that returns the POSITIVE log likelihood
	def lnlike(theta):
 	 	#energyVar gets modified in loop. So far hasn't caused any issues with multiprocessing
		global energyVar

  		# Reset nllVal to zero because we'll add to this for each source we generate
  		# an nll value for
		nllVal=0
  
  		#Load theta values into parameters
		alpha=theta[0]
		beta=theta[1]
		gamma=theta[2]
		offset=theta[3]
  
		#Load remaining theta values for each source, generate smeared pdf, calculate 
 		#nll value
		#for chNum in range(0,len(channels)):
    
		qf=theta[4]
		amplitude=theta[5]

    		#Make ampltitude var 
		sourceCountsVar = ROOT.RooRealVar("sourceCountsVar","sourceCountsVar",amplitude*expectedSourceCounts)
		#sourceCountsVar.setConstant(1)
    
		#Make array of sigmas from res params and sim values
		sigmas=f(simArray,alpha,beta,gamma)
    
		#Generate valsToGen random values for every entry in sim
		smearedArrayList=[]
		for i in range(0,valsToGen):
			#print(str(qf*(numpy.random.normal(simArrays[chNum],sigmas)-offset)))
			smearedArrayList.append(qf*(numpy.random.normal(simArray,sigmas)-offset))
    
		#print("Smeared array is: " + smearedArrayList  ) 
		smearedArrayArray=numpy.array(smearedArrayList)
		flatArray=smearedArrayArray.flatten()
		#print( "Flat array is: " + flatArray )
    
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
    
		#Make Model
		pdfList = ROOT.RooArgList(bgndDataPdf,simPdf)
		ampList = ROOT.RooArgList(scaledBgndEntriesVar,sourceCountsVar)
		model = ROOT.RooAddPdf("model","model",pdfList,ampList)
		model.fixCoefRange("fitRange")
    
		#Compute nll
		nll = model.createNLL(sourceDataHist,
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
			frame = energyVar.frame(0,20,20)
			frame.SetTitle("Channel "+str(bdNum)+" alpha="+str(alpha)+", beta="+str(beta)+", gamma="+str(gamma)+", offset="+str(offset))
			#Plot source data
			sourceDataSet.plotOn(
				frame,
				ROOT.RooFit.Name("Source"),
				ROOT.RooFit.MarkerColor(1),
				ROOT.RooFit.FillColor(0),
				ROOT.RooFit.Range("fitRange"),
				ROOT.RooFit.Binning(binning)
			)
			#Plot components
			bgndName = "bgndDataPdf"
			model.plotOn(
				frame,
				ROOT.RooFit.Name("Bgnd"),
				ROOT.RooFit.Components(bgndName),
				ROOT.RooFit.LineColor(ROOT.kSolid),
				ROOT.RooFit.FillColor(0),
				ROOT.RooFit.ProjWData(sourceDataSet),
				ROOT.RooFit.Range("fitRange")
			)
			model.plotOn(
				frame,ROOT.RooFit.Name("Sim"),
				ROOT.RooFit.Components("simPdf"),
				ROOT.RooFit.LineColor(ROOT.kRed),
				ROOT.RooFit.FillColor(0),
				ROOT.RooFit.ProjWData(sourceDataSet),
				ROOT.RooFit.AddTo("Bgnd"),
				ROOT.RooFit.Range("fitRange")
			)
      
			#Draw
			#frame.GetYaxis.SetRangeUser(0.1,1e5)
			frame.Draw()
			frame.SetMinimum(0.1)
			frame.SetMaximum(1e4)
      
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
			c1.SaveAs(plotPath+"bestFit_simultaneous_ch"+str(bdNum)+".pdf")
      
			#Reset integrator for step size
			ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-7)
			ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-7)
      
			#Memory management for plotting
			frame.Delete("a")
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
	with open(plotPath+"sampler_simultaneous_ch"+str(bdNum)+".csv", 'w') as sampleOutFile:
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
	plt.savefig(plotPath+"traceplots_simultaneous"+".pdf")

	#CORNER PLOT HERE 
	samples=sampler.flatchain
	fig = corner.corner(samples, labels=labels, ranges=ranges, quantiles=[0.16,0.5,0.84],show_titles=True,title_kwargs={'fontsize':12})
	fig.savefig(plotPath+"corner_simultaneous"+".pdf")

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


#Execute main function 	
if __name__== "__main__":
  main()		
