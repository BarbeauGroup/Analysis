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

import sys
import ROOT
import random
from ROOT import TGraph, TH1D, TH1F, TH2F, TCanvas, TFile, TTree
import numpy
import scipy as sp
import matplotlib as mpl
from array import array
from tqdm import tqdm
import emcee
import corner

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
	upperEnergyBound = 200000
	
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
		sourceFilenames.append( calibrationFolder + source + ".root" )
		print( "Source file is: " + sourceFilenames[-1] )
		simFilenames.append( calibrationFolder + source + "-sim.root" )
		print( "Simulation file is: " + simFilenames[-1] )
		
	dataTreeIntegral = array(dataTreeIntegralBranchType,[0])
	dataTreeChannel = array(dataTreeChannelBranchType,[0])
	
	########################
	## Set up observables ##
	########################
	integralVar = ROOT.RooRealVar( "integralVar","integralVar", lowerIntegralBound, upperIntegralBound )
	binning = ROOT.RooBinning( nBins, lowerIntegralBound, upperIntegralBound, "binning" )
	integralVar.setBinning( binning )
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
	expectedSourceCounts=[]
	scaledBgndEntriesVars=[]
	sourceDataSets = []
	sourceDataHists = []
	sourceCountsInFitRanges = []
	srcDataHistPdfs = []
	for sourceNum in range(0,len(calibrationSources)):
		source = calibrationSources[ sourceNum ]
		scaledBgndEntries = bgndCountsInFitRange
		scaledBgndEntriesVars.append(ROOT.RooRealVar("scaledBgndEntriesVar_"+source,"scaledBgndEntriesVar_"+source,scaledBgndEntries))
		scaledBgndEntriesVars[-1].setConstant()
		sourceFile = ROOT.TFile( sourceFilenames[ sourceNum ], "READ" )
		sourceTree = sourceFile.Get( dataTreeName )
		sourceTree.SetBranchAddress( dataTreeIntegralBranchName, dataTreeIntegral )
		sourceTree.SetBranchAddress( dataTreeChannelBranchName, dataTreeChannel )
		sourceDataSets.append(ROOT.RooDataSet( "sourceDataSet_"+source, "sourceDataSet_"+source, argSet ) )
		nSourceEntries = sourceTree.GetEntries()
		print( "Found " + str( nSourceEntries ) + " entries in " + source + " data set." )
		for entry in range( 0, nSourceEntries ):
			sourceTree.GetEntry(entry)
			if dataTreeChannel[0] == dataTreeChannelNum:
				if ( dataTreeIntegral[0] >= lowerIntegralBound ) and ( dataTreeIntegral[0] <= upperIntegralBound ):
					integralVar.setVal( dataTreeIntegral[0] )
					sourceDataSets[ sourceNum ].add( argSet )
		sourceCountsInFitRanges.append( sourceDataSets[ sourceNum ].numEntries() )
		expectedSourceCounts.append(sourceCountsInFitRanges[sourceNum]-bgndCountsInFitRange)
		print( "Found "+str(sourceCountsInFitRanges[-1])+" entries for " + source + " in the fit range" )
		( sourceDataSets[ sourceNum ].get().find("integralVar") ).setBins( nBins )
		sourceDataHists.append( sourceDataSets[ sourceNum ].binnedClone() )
		
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
	PLOT = 1
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
		global integralVar

  		# Reset nllVal to zero because we'll add to this for each source we generate
  		# an nll value for
		nllVal=0
  
  		#Load theta values into parameters
		alpha=theta[0]
		beta=theta[1]
		gamma=theta[2]
		slope=theta[3]
		offset=theta[4]
  
		#Load remaining theta values for each source, generate smeared pdf, calculate 
 		#nll value
		for sourceNum in range(0,len(calibrationSources)):
    
			amplitude=theta[5+sourceNum]

    			#Make ampltitude var 
			sourceCountsVar = ROOT.RooRealVar("sourceCountsVar","sourceCountsVar",amplitude*expectedSourceCounts[sourceNum])
			sourceCountsVar.setConstant(1)
    
			#Make array of sigmas from res params and sim values
			sigmas=f(simArrays[sourceNum],alpha,beta,gamma)
    
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
				frame.SetTitle(calibrationSources[sourceNum]+": alpha="+str(alpha)+", beta="+str(beta)+", gamma="+str(gamma)+", slope="+str(slope)+", offset="+str(offset))
      
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
				c1.SaveAs(plotPath+"bestFit_simultaneous_"+calibrationSources[sourceNum]+"_ch"+calibrationChannel+".pdf")
      
				#Reset integrator for step size
				ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-7)
				ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-7)
      
				#Memory management for plotting
				frame.Delete()
				#del Frame()
      
			#Memory management
			smearedSimDataSet.reset()
			smearedSimDataSet.Delete()
			#del smearedSimDataSet
			pdfList.Delete()
			#del pdfList
			ampList.Delete()
			#del ampList
			sourceCountsVar.Delete()
			#del sourceCountsVar
			smearedSimDataHist.Delete()
			#del smearedSimDataHist
			simPdf.Delete()
			#del simPdf
			#del model
			nll.Delete()
			#del nll
			gc.collect()
    
    
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
	with open(plotPath+"sampler_simultaneous_ch"+calibrationChannel+".csv", 'w') as sampleOutFile:
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
	plt.savefig(plotPath+"traceplots_simultaneous_ch"+calibrationChannel+".pdf")

	#CORNER PLOT HERE 
	samples=sampler.flatchain
	fig = corner.corner(samples, labels=labels, ranges=ranges, quantiles=[0.16,0.5,0.84],show_titles=True,title_kwargs={'fontsize':12})
	fig.savefig(plotPath+"corner_simultaneous_ch"+calibrationChannel+".pdf")

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
