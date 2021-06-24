#                                                                                     #
# =================================================================================== #
# =================================================================================== #
# ================== Python Based Nuclear Recoil Fitter - CeBr3 ===================== #
# =================================================================================== #
# ============================ Written by C. Awe - April 2021 ======================= #
# =================================================================================== #
# =================================================================================== #
#                                                                                     #
#                         Code to compute CeBr3 quenching factors.                    #
#                                                                                     #

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
time_lowList = [321,314,308,306,308,311,317,325] #Timing cut lower bounds for each BD in ns.
time_HighList = [351,344,338,336,338,341,347,355] #Timing cut high bounds for each BD in ns.
bd_timingAdjustments = [315.18,305.68,299.7,328,299.34,303.18,309.05,317.61] #Adjusts from simulated tof to "real" tof. Temporary. 
enrgList = [71.6,43.1,18.0,2.1,8.9,25.9,57.7,86.7] #Recoil energy for each backing detector.

lowerEnergyBound = 0 #keVee
upperEnergyBound = 100 #keVee
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

	#################
	## Emcee Setup ##
	#################

	ndim = 4 #Gamma, Beta, Signal Amplitude, Background Amplitude.
	nBins = 500

	#Parameter names for plots.
	labels = ['Gamma','Beta','Signal','Bgnd']

	#Minimum and maximum values for the parameters.
	mins = [0,0,400,10000]
	maxes = [20,1,20000,100000]

	#MC Stats
	nwalkers = 55 #Based on Sam's code.
	nBurnInSteps = 300
	nSteps = 1800

	###############################
	## Specify Import/Fit Ranges ##
	###############################
	fitRangeMin = 0 #keVee
	fitRangeMax = 8 #keVee

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
	sourceChain = ROOT.TChain( dataTreeName )
	file = open( sourceFilename, 'r' )
	for line in file:
		print( "Adding file " + line )
		sourceChain.AddFile( line[:-1] )
	print( "TChain Built." )
	time_Low = time_lowList[int(bdNum) - 1]
	time_High = time_HighList[int(bdNum) - 1]
	sourceDataSet = ROOT.RooDataSet("sourceDataSet","sourceDataSet",argSet)
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
						energyVar.setVal( entry.scatterer_integral * adc_to_keV )
						sourceDataSet.add( argSet )
						energyVar.setVal( entry.scatterer_noise * adc_to_keV )
						bgndDataSet.add( argSet )
		count += 1
	reducedBgndDataSet = bgndDataSet.reduce(ROOT.RooFit.Cut(cut))
	bgndCountsInFitRange = reducedBgndDataSet.numEntries()
	print( "Found "+str(bgndCountsInFitRange)+" bgnd entries in the fit range" )
	scaledBgndEntries = bgndCountsInFitRange
	#scaledBgndEntriesVar = ROOT.RooRealVar("scaledBgndEntriesVar","scaledBgndEntriesVar",scaledBgndEntries)
	#scaledBgndEntriesVar.setConstant()
	expectedSourceCounts = sourceDataSet.numEntries()
	print( "Found "+str(expectedSourceCounts)+" entries for BD " + str(bdNum) + " in the fit range" )
	( sourceDataSet.get().find("energyVar") ).setBins( nBins )
	sourceDataHist = sourceDataSet.binnedClone()

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
		Gamma=theta[0]
		Beta=theta[1]
		Signal=theta[2]
		Background=theta[3]

    		#Make vars and pdfs.
		sourceCountsVar = ROOT.RooRealVar("sourceCountsVar","Signal",Signal)
		sourceCountsVar.setConstant(1)
		scaledBgndEntriesVar = ROOT.RooRealVar("scaledBgndEntriesVar","Background",Background)
		scaledBgndEntriesVar.setConstant(1)
		( bgndDataSet.get().find("energyVar") ).setBins( nBins )
		bgndDataHist = bgndDataSet.binnedClone()
		bgndDataPdf = ROOT.RooHistPdf("bgndDataPdf","Background",energyVar,bgndDataHist,0) #1 specifies interpolation order
		gammaVar = ROOT.RooRealVar("gammaVar","Gamma",Gamma,1,20)
		gammaVar.setConstant(1)
		betaVar = ROOT.RooRealVar("betaVar","Beta",Beta,0.03,1)
		betaVar.setConstant(1)
		muVar = ROOT.RooRealVar("muVar","Mu",0,0,1)
		muVar.setConstant(1)
		gammaPdf = ROOT.RooGamma("gammaPdf","Gamma PDF",energyVar, gammaVar, betaVar, muVar)
		#Make Model
		pdfList = ROOT.RooArgList(bgndDataPdf,gammaPdf)
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

		#Make a binned version of the model.
		gammaData = gammaPdf.generate(energyVar,100000)
		gammaHist = gammaData.binnedClone()
		binnedGamma = ROOT.RooHistPdf("binnedGamma","Binned Gamma", energyVar, gammaHist, 0)
		binnedPdfList = ROOT.RooArgList(bgndDataPdf,binnedGamma)
		binnedModel = ROOT.RooAddPdf("binnedModel","Binned Model", binnedPdfList,ampList)
    
		if (PLOT==1):
			try:
				c1
			except NameError:
				c1=ROOT.TCanvas("c1","c1")
      
			#Reduce integrator for plotting, massively speeds plotting
			ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-5)
			ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-5)
			frame = energyVar.frame(0,8,8)
			frame.SetTitle("Backing Detector "+str(bdNum)+", Recoil Energy = "+str(enrgList[bdNum-1])+" keVnr")
			frame.GetXaxis().SetTitle("Energy [keVee]")
			frame.GetYaxis().SetTitle("Counts / 0.2 keV")
			frame.GetXaxis().CenterTitle();
			frame.GetYaxis().CenterTitle();
			#Plot source data
			sourceDataSet.plotOn(
				frame,
				ROOT.RooFit.Name("Data"),
				ROOT.RooFit.MarkerColor(1),
				ROOT.RooFit.FillColor(0),
				ROOT.RooFit.Range("fitRange"),
				ROOT.RooFit.Binning(binning)
			)
			#Plot components
			bgndName = "bgndDataPdf"
			binnedModel.plotOn(
				frame,
				ROOT.RooFit.Name("Background"),
				ROOT.RooFit.Components("bgndDataPdf"),
				ROOT.RooFit.LineColor(33),
				ROOT.RooFit.FillColor(0),
				ROOT.RooFit.ProjWData(sourceDataSet),
			)
			binnedModel.plotOn(
				frame,ROOT.RooFit.Name("Nuclear Recoils"),
				ROOT.RooFit.Components("binnedGamma"),
				ROOT.RooFit.LineColor(46),
				ROOT.RooFit.FillColor(0),
				ROOT.RooFit.ProjWData(sourceDataSet),
			)
			binnedModel.plotOn(
				frame,
				ROOT.RooFit.LineColor(30),
				ROOT.RooFit.FillColor(0),
				ROOT.RooFit.ProjWData(sourceDataSet)
			)
			resHist = frame.residHist()
      
			#Draw
			#frame.GetYaxis.SetRangeUser(0.1,1e5)
			frame.Draw()
			frame.SetMinimum(0.001)
			frame.SetMaximum(10e4)
      
			#Add legend
			leg = ROOT.TLegend(0.62,0.67,0.87,0.85); 
			sourceObj = frame.findObject("Data");
			bgndObj = frame.findObject("Background");
			recObj = frame.findObject("Nuclear Recoils");
			leg.AddEntry(sourceObj,"Source","P")
			leg.AddEntry(bgndObj,"Background","L")
			leg.AddEntry(recObj,"Nuclear Recoils","L")
			leg.Draw("same")
      
			#Draw
			c1.SetLogy()
			c1.Modified()
			c1.Update()
			c1.SaveAs(plotPath+"bestFit_simultaneous_ch"+str(bdNum)+".pdf")
				
			#Make Residuals Plot.
			#resFrame = energyVar.frame(0,8,8)
			#res.SetTitle("Channel "+str(bdNum)+" Residuals")
			#resFrame.GetXaxis().CenterTitle();
			#resFrame.GetXaxis().SetRangeUser(fitMin,fitMax);
			#resFrame.GetXaxis().SetTitle("keVee");
			#resFrame.GetYaxis().CenterTitle();
			#resFrame.GetYaxis().SetTitle("Counts / 0.2 keV");
			#resFrame.addPlotable(resHist,"p");
			#resFrame.Draw(); 
      
			#Reset integrator for step size
			ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-7)
			ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-7)
      
			#Memory management for plotting
			frame.Delete("a")
			#del frame()
      
			#Memory management
			del pdfList
			del ampList
			del sourceCountsVar
			del scaledBgndEntriesVar
			del bgndDataHist
			del bgndDataPdf
			del model
			del gammaData
			del gammaHist
			del binnedGamma
			del binnedPdfList
			del binnedModel
			nll.Delete()
			del nll
			del gammaVar
			del betaVar
			del muVar
			del gammaPdf
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
	plt.savefig(plotPath+"traceplots_simultaneous_bd"+str(bdNum)+".pdf")

	#CORNER PLOT HERE 
	samples=sampler.flatchain
	fig = corner.corner(samples, labels=labels, ranges=ranges, quantiles=[0.16,0.5,0.84],show_titles=True,title_kwargs={'fontsize':12})
	fig.savefig(plotPath+"corner_simultaneous_bd"+str(bdNum)+".pdf")

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










