#Python code to optimize CMA filter parameters on toy data.
#
#Usage: python optimizeCMA.py <file containing toy data> <output file name>
#
# Notes:
#
import ROOT
import sys
from array import array
from scipy import stats
import random
import gc
import math

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
  
#Gets a guess at the baseline value of a waveform and its rms.
def getBaselineGuess( waveform ):

  baselineInfo = [0]*2
  variationGuess = 10
  baselineHist = ROOT.TH1F("baselineHist","baselineHist",16000,0,16000)
    
  #Fill the histogram.
  for i in range(0, len( waveform )):
    baselineHist.Fill( waveform[i] )
        
  #Find the most frequent value.
  binMax = baselineHist.GetMaximumBin()
  #print binMax
  
  baseline = 0.
  weight = 0.
  stdDev = 0.
  n = 0.
  
  if (binMax > variationGuess and binMax + variationGuess < 16384):
      for i in range(binMax - 1 - variationGuess, binMax + variationGuess - 1):
          baseline += baselineHist.GetBinCenter(i) * baselineHist.GetBinContent(i)
          weight += baselineHist.GetBinContent(i)
          stdDev += baselineHist.GetBinContent(i)*pow(baselineHist.GetBinCenter(i)-binMax,2)
          n += baselineHist.GetBinContent(i)
          
  else:
    print("Error computing the baseline, baseline too close to the edges of ADC.")
    return baselineInfo
    
  #Calculate baseline.
  if weight > 0:
  	baseline = baseline / weight
    
  #Calculate std. dev.
  if n > 0:
      stdDev = math.sqrt(1./n * stdDev)
  
  baselineInfo[0] = baseline
  baselineInfo[1] = stdDev
  return baselineInfo
  
#Gets the CMA filter of a waveform. Implementation of C++ code by G.C. Rich.
def getCMAFilter( waveform, halfWidth, preloadValue, rejectThreshold ):

  #Initialize our filter.
  cmaFilter = []

  if halfWidth > len(waveform):
    print("Halfwidth greater than waveform length! Setting filter to zero.")
    return cmaFilter
    
  #Preload the filter.
  movingBaselineSum = 0.
  movingBaselineFilter = []
  for i in range(0,int(halfWidth)):
    movingBaselineFilter.append(preloadValue)
    movingBaselineSum += preloadValue
  
  #Now step through the first halfWidth of values, adding to the average if within reject threshold.
  for j in range(0, int(halfWidth)):
    #Calculate the current moving average.
    movingAverage = movingBaselineSum / len( movingBaselineFilter )
    #Check that we're within the reject threshold.
    if abs( waveform[j] - movingAverage ) < rejectThreshold:
        movingBaselineFilter.append( waveform[j] )
        movingBaselineSum += waveform[j]

  #print len( movingBaselineFilter )
  #Now moving average full. Start at element 0 and compute moving average.
  for k in range(0, len( waveform )):
	#print("k is " + str(k))
	#print("halfwidth is " + str( halfWidth ))
	#print("waveform legth is " + str( len( waveform ) ))
	#Check if there is an element k + halfWidth ahead. If out of scope, start popping off elements.
	#Push back the moving average.
	#print(str(len( movingBaselineFilter )))
	movingAverage = movingBaselineSum / len( movingBaselineFilter )
	cmaFilter.append( movingAverage )
	if ((k + int(halfWidth)) < len( waveform )):
		#print("k + halfwidth < len( waveform )")
		movingAverage = movingBaselineSum / len( movingBaselineFilter )
		#Check that the value at k + halfWidth is within reject threshold.
		if abs( waveform[k + int(halfWidth)] - movingAverage ) < rejectThreshold:
			movingBaselineFilter.append( waveform[k + int(halfWidth)] )
			movingBaselineSum += waveform[k + int(halfWidth)]
			#Check if the filter is larger than 2*halfWidth + 1. If so, pop off the front element.
		if len( movingBaselineFilter ) > 2 * int(halfWidth) + 1:
			#print("2*halfwidth+1 < len( movingbaselinefilter ), popping off elements.")
			movingBaselineSum -= movingBaselineFilter[0]
			del movingBaselineFilter[0]
	else:
		#print("k + halfwidth > len( waveform ), popping off elements.")
		movingBaselineSum -= movingBaselineFilter[0]
		del movingBaselineFilter[0]
    
    
  return cmaFilter


#Code to integrate waveform pulses. Adapted from C++ code written by S. Hedges.
def getIntegral( waveform, baseline, startSample, numSamples ):
    
  if startSample + numSamples > len( waveform ):
      numSamples = len( waveform ) - startSample - 1
  
  integral = 0.
  for i in range(startSample, startSample + numSamples):
      integral += waveform[i] - baseline
      
  return integral
  
#Code to find pulses within a waveform. Supply this with filtered waves.
def getPulsesFromFIR( waveform, baseline, windowSize, gapSize, firThresh, fraction, holdOffSamples, FIRSamples):
    
    #Holds the locations of found pulses.
    pulseTimes = []
    
    #Holds the moving average of two windows.
    leadingWindow = []
    trailingWindow = []
    sumLeadingWindow = 0.
    sumTrailingWindow = 0.
    
    #Holds the difference between the moving averages.
    firVector = []
    triggerSample = 0
    
    #Calculate the FIR vector.
    for i in range(0, len( waveform )):
        #For the first window+gapSize samples, load baseline values.
        if i < windowSize + gapSize:
            leadingWindow.append( baseline )
            sumLeadingWindow += baseline
            trailingWindow.append( baseline )
            sumTrailingWindow += baseline
            
        else:
            #Push new values and sum.
            leadingWindow.append( waveform[i] )
            sumLeadingWindow += waveform[i]
            trailingWindow.append( waveform[i - windowSize - gapSize] )
            sumTrailingWindow += waveform[i - windowSize - gapSize]
            
        #If vector sizes are > windowSize, remove oldest values.
        if len( leadingWindow ) > windowSize:
            sumLeadingWindow -= leadingWindow[0]
            del leadingWindow[0]
        if len( trailingWindow ) > windowSize:
            sumTrailingWindow -= trailingWindow[0]
            del trailingWindow[0]
            
        #Calculate the difference.
        firVector.append( sumLeadingWindow - sumTrailingWindow )
        
    #Now look for stuff above threshold.
    j = 0
    while j < len( waveform ):
        if firVector[j] >= firThresh:
            maxValue = 0
            maxSample = 0
            triggerSample = 0
            #Search for max value up to FIRSamples.
            for k in range(j, j + FIRSamples):
                #Make sure we don't extend  beyond the range of the waveform.
                if k > len( waveform ):
                    break
                if firVector[k] > maxValue:
                    maxValue = firVector[k]
                    maxSample = k
                    
            #Check if we have a maxValue at position > 0
            if maxSample > 0:
                for k in range(maxSample, j + FIRSamples):
                    if k > len( waveform ):
                        break
                    if firVector[k] <= fraction*maxValue:
                        triggerSample = k
                        break
                #Make sure we can interpolate and then do so.
                if triggerSample > 0:
                    valueAtTrigger = firVector[triggerSample]
                    valueBeforeTrigger = firVector[triggerSample-1]
                    
                    #Interpolation.
                    diffMaxDiv2AndMaxAtTrig = maxValue / 2. - valueAtTrigger
                    diffBeforeAndAtTrig = valueBeforeTrigger - valueAtTrigger
                    corrector = diffMaxDiv2AndMaxAtTrig / diffBeforeAndAtTrig
                    
                    pulseTimes.append( triggerSample - corrector )
            j += holdOffSamples
            
        else:
            j += 1
            
    return pulseTimes

#Main function. This will handle running CMA filters with different parameters on the data and comparing
#the resulting integral distributions to the "trueIntegral" distribution of the toy data.  
def main():

  #Analysis values.
  scatterer_preOnsetIntegralSamples=10
  scatterer_postOnsetIntegralSamples=100
  scatterer_riseTime = 4
  scatterer_gapTime = 4
  scatterer_firThresh = 400
  scatterer_holdOffSamples = 100
  scatterer_FIRSamples = 4
    
  #Initializations.
  integral = array( 'd', [0] )
  integral_delta = array( 'd', [0] )
  deltaMean = array( 'd', [0] )
  deltaRMS = array( 'd', [0] )
  pulseTime = array( 'd', [0] )
  halfWidth = array( 'd', [0] )
  rejectThreshold = array( 'd', [0] )
  halfWidthMin = 38. #Units of baseline samples.
  halfWidthMax = 48. #Units of samples.
  halfWidthSteps = 10
  halfWidthStepSize = float(( halfWidthMax - halfWidthMin ) / halfWidthSteps)
  print "Half width step size will be " + str( halfWidthStepSize )
  rejectThresholdMin = 1.85 #Units of baseline std. dev.
  rejectThresholdMax = 2.15 #Units of baseline std. dev.
  rejectThresholdSteps = 6
  rejectThresholdStepSize = float(( rejectThresholdMax - rejectThresholdMin ) / rejectThresholdSteps)
  print "Reject Threshold step size will be " + str( rejectThresholdStepSize )

  #Get run parameters from the user.
  #filename = input( "Enter the path to the toy data file: " )
  filename = "CeBr3_FakePulses.root"
  #output = input( "Enter the path and name of the output file you wish to produce: " )
  output = "CeBr3_CMAOptimization.root"
  inputFile=ROOT.TFile(filename,"READ")
  inputTree=inputFile.Get("sis3316tree")
    
  #Set up root stuff.
  outputFile = ROOT.TFile( output, "recreate" )
  cmaTree = ROOT.TTree( "cmaTree", "Output Tree" )
  cmaTree.Branch( "deltaMean", deltaMean, "deltaMean/D" )
  cmaTree.Branch( "deltaRMS", deltaRMS, "deltaRMS/D")
  cmaTree.Branch( "halfWidth", halfWidth, "halfWidth/D" )
  cmaTree.Branch( "rejectThreshold", rejectThreshold, "rejectThreshold/D" )  
  optimizationHist = ROOT.TH2F("optimizationHist","Half Width and Reject Threshold",halfWidthSteps,halfWidthMin,halfWidthMax,rejectThresholdSteps,rejectThresholdMin,rejectThresholdMax)
  deltaHist = ROOT.TH1F("deltaHist","deltaHist",16000,0,16000) #Keeps track of how we do relative to the true integral.
  #Step through our grid and process.
  #minDelta = 1000 #Start at a large value.
  for i in range( 0, halfWidthSteps):
  	print "On " + str( i + 1 ) + " of " + str( halfWidthSteps ) + " half width steps."
	halfWidth[0] = float(halfWidthMin + i * halfWidthStepSize)
	print "halfWidth is " + str(halfWidth[0])
	for j in range( 0, rejectThresholdSteps ):
		print "On " + str( j + 1 ) + " of " + str( rejectThresholdSteps ) + " reject threshold steps."
		rejectThreshold[0] = rejectThresholdMin + j * rejectThresholdStepSize
		#deltaHist = ROOT.TH1F("deltaHist","deltaHist",16000,0,16000) #Keeps track of how we do relative to the true integral.
		for entry in inputTree:
			baselineInfo = getBaselineGuess( entry.waveform )
			baseline = baselineInfo[0]
			#print baseline
			#print baselineInfo[0]
			stdDev = baselineInfo[1]
			cmaFilter = getCMAFilter( entry.waveform, halfWidth[0], baseline, rejectThreshold[0] * stdDev )
			waveList = []
			filteredWave = []
			for k in range( 0, len( entry.waveform ) ):
				filteredWave.append(entry.waveform[k] - cmaFilter[k])
				waveList.append(entry.waveform[k])
			#plotList( waveList )
			#plotList( cmaFilter )
			#plotList( filteredWave )
			pulses = getPulsesFromFIR( filteredWave, 0., scatterer_riseTime, scatterer_gapTime, scatterer_firThresh, 0.5, scatterer_holdOffSamples, scatterer_FIRSamples)
			if len( pulses ) > 0:
				startSample = int( pulses[0] - scatterer_preOnsetIntegralSamples )
				endSample = int( pulses[0] + scatterer_postOnsetIntegralSamples )
				numSamples = int( scatterer_preOnsetIntegralSamples + scatterer_postOnsetIntegralSamples )
				if (( startSample > 0 ) and ( endSample < len( filteredWave ))):
					integral[0] = getIntegral( filteredWave, 0., startSample, numSamples )
					integral_delta[0] = abs( integral[0] - entry.trueIntegral )
					deltaHist.Fill( integral_delta[0] )
		deltaMean[0] = deltaHist.GetMean()
		deltaRMS[0] = deltaHist.GetRMS()
		optimizationHist.Fill(halfWidth[0],rejectThreshold[0],deltaMean[0])
		cmaTree.Fill()
		deltaHist.Reset()
        
  print("Routine compelte, writing the output file.")
  outputFile.Write()
  outputFile.Close()
  
#Execute main function
if __name__== "__main__":
  main()
    