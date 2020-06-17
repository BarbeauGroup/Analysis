#Python port of grayson's code for CsI[Na]. Some modifications.
#
#Usage: python genToyCeBr3Data.py <file containing traces to harvest noise from> <output file name>
#
# Notes:
#   - You can potentially change the way fake pulses are generated if fast and slow decay time parameters
#     are not known for your scintillator. This is done in the getPETimes function.
#   - getPETimes has been modified to create pulses with the passed in integral rather than number of SPE
#   - Be sure the ranges of the RooFit variables in getPEtimes are large enough to cover the input parameters
#   - The code assumes that the majority of the scintillator waveforms we're harvesting noise from are blank.
#     It calculates the mode of the waveform and considers that the baseline.
#   - In the function gathering noise pulses, you can set a limit for the max pulse height allowed (after
#     subtracting baseline) for a trace to be considered empty.
#
import ROOT
import sys
import array
from scipy import stats
import random
import gc

######################
##Crystal parameters##
######################
#From https://www.berkeleynucleonics.com/cerium-bromide, using
#uncertainty as sigma
fastDecayTime=18. #ns
fastDecayTimeSigma=1. #ns
slowDecayTime=1040 #ns
slowDecayTimeSigma=10 #ns
slowFraction=0.0 #Some energy dependence to this
#From https://arxiv.org/pdf/1307.1398.pdf
riseTime = 0.7 #ns
riseTimeSigma = 0.1 #Unc. used as sigma

##################
##DAQ PARAMETERS##
##################
waveformSamples=500
nsPerSample=4

#######################
##Analysis parameters##
#######################
scatterer_preOnsetIntegralSamples=10
scatterer_postOnsetIntegralSamples=100

##################
##Get input file##
##################
inpFilename=sys.argv[1]
outputFilename=sys.argv[2]

#Make RooFit observable global so we can use it inside and outside of cutnions
sampleTime=ROOT.RooRealVar("sampleTime","sampleTime",0,waveformSamples*nsPerSample)

onsetLoc=800 #Time in ns where we generate sample pulses

#############
##Make PDFs##
#############
print("Making PDFs")
#Don't do inside a function so we can use the PDFs globally##
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

#Rise time
#Mean corresponds to onset of pulse, not sure of limits @grayson
riseGauss_mean = ROOT.RooRealVar("riseGauss_mean","riseGauss_mean",onsetLoc,-1000,2500)
#Not sure how to properly take into account uncertainty on the rise time, but this
#corresponds to the rise time @grayson
riseGauss_sigma = ROOT.RooRealVar("riseGauss_sigma","riseGauss_sigma",riseTime,0.,1000.)

#Decay time
fastDecayTime_mean = ROOT.RooRealVar("fastDecayTime_mean","Mean of fast decay time",fastDecayTime)
fastDecayTime_sigma = ROOT.RooRealVar("fastDecayTime_sigma","Sigma of fast decay time",fastDecayTimeSigma)
slowDecayTime_mean = ROOT.RooRealVar("slowDecayTime_mean","Mean of slow decay time",slowDecayTime)
slowDecayTime_sigma = ROOT.RooRealVar("slowDecayTime_sigma","Sigma of slow decay time",slowDecayTimeSigma)

#GaussModels to form the basis for the fast and slow components
riseGaussianFast = ROOT.RooGaussModel("riseGaussianFast","riseGaussianFast",
sampleTime,riseGauss_mean,riseGauss_sigma)
riseGaussianSlow = ROOT.RooGaussModel("riseGaussianSlow","riseGaussianSlow",
sampleTime,riseGauss_mean,riseGauss_sigma)

#RooDecays for the actual fast and slow components. @grayson why these limits for the RooRealVars?
decayTime_fast = ROOT.RooRealVar("decayTime_fast","decayTime_fast",fastDecayTime,0,10000)
pdf_fastComponent = ROOT.RooDecay("pdf_fastComponent","pdf_fastComponent",
sampleTime,decayTime_fast,riseGaussianFast,ROOT.RooDecay.SingleSided)
decayTime_slow = ROOT.RooRealVar("decayTime_slow","decayTime_slow",slowDecayTime,0,20000)
pdf_slowComponent = ROOT.RooDecay("pdf_slowComponent","pdf_slowComponent",
  sampleTime,decayTime_slow,riseGaussianSlow,ROOT.RooDecay.SingleSided)
  
#Fraction of fast component
fastComponentFraction=ROOT.RooRealVar("fastComponentFraction","fastComponentFraction",1-slowFraction,0,1)
#fastComponentFraction.SetConstant();

#Make model without systematics by combining fast and slow fractions
pdf_signal_preSystematic=ROOT.RooAddPdf("pdf_signal_preSystematic","pdf_signal_preSystematic",
ROOT.RooArgList(pdf_fastComponent, pdf_slowComponent),ROOT.RooArgList(fastComponentFraction),ROOT.kTRUE)

#Make uncertainty PDFs
fastDecayTime_uncPdf = ROOT.RooGaussian("fastDecayTime_uncPdf","Uncertainty PDF on fast decay time",
decayTime_fast,fastDecayTime_mean,fastDecayTime_sigma)
slowDecayTime_uncPdf = ROOT.RooGaussian("slowDecayTime_uncPdf","Uncertainty PDF on slow decay time",
decayTime_slow,slowDecayTime_mean,slowDecayTime_sigma)
decayTimeSystematics = ROOT.RooProdPdf("decayTimeSystematics","Systematic uncertainty PDF of decay times",
slowDecayTime_uncPdf, fastDecayTime_uncPdf)

#Make arg sets because otherwise this is a memory leak when passed into a RooFit function
fastArgSet=ROOT.RooArgSet(decayTime_fast)
slowArgSet=ROOT.RooArgSet(decayTime_slow)
argSet=ROOT.RooArgSet(sampleTime)

#Make signal by multiplying systematic-free signal with uncertainty PDF
pdf_signal = ROOT.RooProdPdf("pdf_signal", "pdf_signal",pdf_signal_preSystematic,decayTimeSystematics)


#Making these global, trying to fix memory leaks  
dataSamples=ROOT.RooDataSet("dataSamples","dataSamples",argSet)
protoData_fast=ROOT.RooDataSet("protoData_fast","protoData_fast",argSet)
protoData_decayTimes=ROOT.RooDataSet("protoData_decayTimes","protoData_decayTimes",argSet)


#########
##pause##
#########
#Helper function to avoid repeated code. Pauses until user types a key on the terminal.
#TCanvases are interactive while this is going on
def pause():
  try:
    input("Press enter to continue")
  except SyntaxError:
    pass


############
##plotList##
############
#Plots a list of samples as a waveform
def plotList(samples):
  nSamples=len(samples)
  hist=ROOT.TH1D("hist","hist",nSamples,0,nSamples)
  c1=ROOT.TCanvas("c1","c1")
  for i in range(0,nSamples):
    hist.SetBinContent(i+1,samples[i])
    
  hist.Draw()
  c1.Modified()
  c1.Update()
  pause()
  hist.Delete()
  
  
##############
##GetPEtimes##
##############
#Output is the same as Grayson's, but takes in the PDF and fast and slow decay time RooRealVars
def getPEtimes(nPEs):
  global protoData_fast
  global protoData_decayTimes
  global dataSamples
  
  #Sample from fast data
  protoData_fast = fastDecayTime_uncPdf.generate(fastArgSet,nPEs)
  #Use that data set as a prototype to generate data with fast and slow decay times
  protoDataTemp = ROOT.RooFit.ProtoData(protoData_fast)
  protoData_decayTimes = decayTimeSystematics.generate(slowArgSet,protoDataTemp)
  
  #Use that data to generate samples
  protoDataTemp2=ROOT.RooFit.ProtoData(protoData_decayTimes)
  dataSamples = pdf_signal.generate(argSet,protoDataTemp2)
  
  protoData_fast.reset()
  protoData_decayTimes.reset()
  protoDataTemp.Delete()
  del protoDataTemp
  protoDataTemp2.Delete()
  del protoDataTemp2
  gc.collect
  
  
  
  #Return RooDataSet containing PE times
  return dataSamples


################
##genToyPulses##
################
#Generates toy pulses
def genToyPulses(nPulsesToGenerate):
  global dataSamples
  
  #Open root file containing histogram of fixedIntegral from CeBr3 channel. We'll
  #use this to randomly sample from our integral distribution.
  histFile=ROOT.TFile("Processed_SIS3316Raw_20200207153627_16.root","READ")
  cebr3IntegralHist=ROOT.TH1F("cebr3IntegralHist","Scatterer Spectrum",100000,0,10000)
  tree=histFile.Get("analysisTree")
  for entry in tree:
    if( entry.scatterer_foundPulse == 1):
  		sample=entry.scatterer_integral
  		#print( "Sample is: " + str( sample ) )
  		cebr3IntegralHist.Fill(sample)
  
  #How many toy pulses to generate
  pulses=[]
  for i in range(0,nPulsesToGenerate):
    #Get random integral from distribution
    integral=cebr3IntegralHist.GetRandom()
    #print("Integral is: " + str(integral))
    #Generate a RooDataSet with a fake pulse based on that integral
    dataSamples.reset()
    dataSamples=getPEtimes(integral)

    #Make a histogram version of the data set
    hist=ROOT.TH1D("hist","hist",waveformSamples,0,waveformSamples*nsPerSample)
    nEntries=dataSamples.numEntries()
    for entry in range(0,nEntries):
      hist.Fill(dataSamples.get(entry).getRealValue("sampleTime"))
      #print( str(dataSamples.get(entry).getRealValue("sampleTime")))
    
    #Make the histogram into a list, append to list of pulses
    pulse=[]
    for bin in range(1,hist.GetNbinsX()+1):
      pulse.append(int(hist.GetBinContent(bin)))
    pulses.append(pulse)
    
    #Memory management
    hist.Delete()
    
    if i%10==0:
      print("Generated pulse "+str(i)+" of "+str(nPulsesToGenerate))
    
  return pulses
  

####################################
##getFarmedScattererBaselinePulses##
####################################
#Modified version of Grayson's code. This returns a list
#of empty pulses we can draw randomly from to get noise
#pulses to add to our fake pulses. It operates on the processed
#pulses from the reconstruction code, and requires the CeBr3 waveforms
#to be saved. The criteria for a trace being empty are:
#  - fir_foundPulse==0 - We didn't find a pulse with our FIR filter
#  - maxHeight<50
def getFarmedScattererBaselinePulses(inpFilename):
  #Get input file, ttree, and number of entries
  inpFile = ROOT.TFile(inpFilename)
  tree = inpFile.Get("analysisTree")
  nEntries = tree.GetEntries()
  
  #Will hold our empty traces
  emptyTraces=[]
  
  #Step through tree. If waveforms meet our conditions for being considered empty, add them
  #to the list
  k = 0
  for entry in tree:
    if entry.scatterer_foundPulse==0:
      #print "Found and empty pulse."
      waveform=[]
      for i in range(0,len(entry.scatterer_waveform)):
        waveform.append(entry.scatterer_waveform[i])
        #print(str(entry.scatterer_waveform[i]))
      baseline=stats.mode(waveform)[0]
      maxHeight=max(waveform)-baseline[0]
      
      #To check if pulse is empty
      if maxHeight<75:
        emptyTraces.append(waveform)
        
    if (k%1000==0):
      print("On entry "+str(k)+" of "+str(nEntries))
    k += 1
  
  #Close the input file
  inpFile.Close()
  
  #Return the list
  return emptyTraces

#######################
##generateToyDataTree##
#######################
#Takes a list of toy pulses, noise traces, and combines
#to make fake pulses. Writes these to a tree to mimic
#the sis3316 output format
def generateToyDataTree(pulses,noiseTraces,outputFilename):

  numSamples=waveformSamples

  #Output file
  outFile=ROOT.TFile(outputFilename,"RECREATE")
  
  #Output tree
  sis3316tree=ROOT.TTree("sis3316tree","Unsorted events")
  
  #Output tree branches that match sis3316tree
  channelID=array.array('H',[0])
  timestamp=array.array('L',[0])
  peakHighIndex=array.array('H',[0])
  peakHighValue=array.array('H',[0])
  pileupFlag=array.array('H',[0])
  nSamples=array.array('i',[0])
  accumulatorSum=array.array('d',8*[0])
  waveform=array.array('H',500*[0])

  sis3316tree.Branch('channelID',channelID,'channelID/s')
  sis3316tree.Branch('timestamp',timestamp,'timestamp/l')
  sis3316tree.Branch('peakHighIndex',peakHighIndex,'peakHighIndex/s')
  sis3316tree.Branch('peakHighValue',peakHighValue,'peakHighValue/s')
  sis3316tree.Branch('pileupFlag',pileupFlag,'pileupFlag/O')
  sis3316tree.Branch('nSamples',nSamples,'nSamples/i')
  sis3316tree.Branch('waveform',waveform,'waveform[500]/s')
  sis3316tree.Branch('accumulatorSum',accumulatorSum,'accumulatorSum[8]/i')
  
  #Fake data branch we added
  trueIntegral=array.array('d',[0])
  sis3316tree.Branch('trueIntegral',trueIntegral,'trueIntegral/D')

  channelID[0]=20
  timestamp[0]=0
  peakHighIndex[0]=0
  peakHighValue[0]=0
  pileupFlag[0]=0
  nSamples[0]=numSamples

  for pulseNum in range(0,len(pulses)):
    pulse=pulses[pulseNum]
    noiseTrace=noiseTraces[random.randrange(len(noiseTraces))]
    realPulse=[(a + b) for a, b in zip(pulse, noiseTrace)]
    
    #plotList(noiseTrace)
    #plotList(pulse)
    #plotList(realPulse)
    
    #Check whether the event is saturated
    saturatedEvents=[i for i in realPulse if i >= 16384]
    
    #Only process non-saturating events
    if len(saturatedEvents)==0:
    
      #Fill branch
      for sample in range(0,len(realPulse)):
        waveform[sample]=realPulse[sample]
        #print realPulse[sample]
      
      #Get "true integral"
      integratedSection=[pulse[i] for i in range(onsetLoc/nsPerSample-scatterer_preOnsetIntegralSamples,onsetLoc/nsPerSample+scatterer_postOnsetIntegralSamples)]
      trueIntegral[0]=sum(integratedSection)
      sis3316tree.Fill()
      
  sis3316tree.Write()
  outFile.Close()
  

#############
##MAIN CODE##
#############
print("Generating toy pulses...")
nPulsesToGenerate=60000
pulses=genToyPulses(nPulsesToGenerate)
print("Done!")

print("Generating noise pulses...")
noisePulses=getFarmedScattererBaselinePulses(inpFilename)
print("Done! Generated "+str(len(noisePulses))+" noise pulses")

print("Making fake pulse tree...")
generateToyDataTree(pulses,noisePulses,outputFilename)
print("Done!")
