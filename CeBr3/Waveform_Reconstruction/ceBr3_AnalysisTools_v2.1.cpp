/* Analysis tools package for CeBr3 
Written by C. Awe based on existing code by S. Hedges et. al.
https://github.com/schedges/cebr3-qf-analysis
4/28/2020

This toolset is loaded directly into root rather than compiling.
While slower, this style is perhaps easier to use across differing
systems.

Usage:
*/

//C++ include files
#include <iostream>
#include <string>
#include <cstdlib>
#include <list>
#include <time.h>
#include <algorithm>
#include <vector>
#include "stdio.h"
#include <dirent.h>
#include <sys/stat.h>
#include <deque>
#include <stdio.h> 
#include <string.h> 

//ROOT include files
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1D.h"
#include "TSpline.h"
#include "TF1.h"
#include "TApplication.h"
#include "TCanvas.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////HELPER FUNCTIONS///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//Plot all values in waveform, get weighted avg. in region around mode

pair<Double_t,Double_t> getBaselineGuess(vector<Double_t> waveform,Int_t baselineVariationGuess=10) {

  //Make a histogram for storing baseline values
	TH1D* baselineHist = new TH1D("baselineHist","baselineHist",16000,0,16000);
	
	//Fill histogram
	for (Int_t i=0; i < waveform.size(); i++) {
		baselineHist->Fill(waveform[i]);
	}
	
	//Find the most frequent value
	Double_t binMax = static_cast<Double_t>(baselineHist->GetMaximumBin());
	
	//Deterministic approach for determining baseline--Grayson's algorithm
	//Around baselineVariationGuess bins of maximum value, get bin centers and
	//multiply by bin height for a weighted average of the baseline. This is 
	//replacing the previous approach of using a TF1 fit for getting the baseline 
	Double_t baseline=0;
	Double_t weight=0;
	Double_t stdDev=0;
	Double_t n=0;
	//Make sure the max isn't too close to the edges of the histogram
	if ((binMax > baselineVariationGuess) && (binMax+baselineVariationGuess<16384)) {
    for(Int_t i = binMax-1-baselineVariationGuess; i<=binMax+baselineVariationGuess-1; i++) {
      baseline += baselineHist->GetBinCenter(i)*baselineHist->GetBinContent(i);
      weight += baselineHist->GetBinContent(i);
      
      stdDev += baselineHist->GetBinContent(i)*pow(baselineHist->GetBinCenter(i)-binMax,2);
      n += baselineHist->GetBinContent(i);
    }
  }
  else {
    cout<<"Error computing baseline, baseline too close to edges of ADC"<<endl;
    delete baselineHist;
    return make_pair(0,0);
  }
  
  delete baselineHist;
   
  //Calculate baseline 
  if (weight>0) {
    baseline = baseline/weight;
  }
  else {
    cout<<"Error calculating baseline, weight of 0"<<endl;
    return make_pair(0,0);
  }

  //Calculate std. dev
  if (n>0) {
    stdDev = sqrt(1./n * stdDev);
  }
  else {
    cout<<"Error calculating stdev, zero entries"<<endl;
    return make_pair(0,0);
  }

  return make_pair(baseline,stdDev);
}


/*
//Fitting version
//Plot all values in waveform, fit region within 10 samples of the mode
pair<Double_t,Double_t> getBaselineGuess(vector<Double_t> waveform,Int_t baselineVariationGuess=10) {

  //Make a histogram for storing baseline values
	TH1D* baselineHist = new TH1D("baselineHist","baselineHist",16384,0,16384);
	
	//Fill histogram
	for (Int_t i=0; i < waveform.size(); i++) {
		baselineHist->Fill(waveform[i]);
	}
	
	//Find the most frequent value
	Double_t binMax = static_cast<Double_t>(baselineHist->GetMaximumBin());
	
	
	if ((binMax > baselineVariationGuess) && (binMax+baselineVariationGuess < 16384)) {
	
    TF1* fit = new TF1("fit","gaus",binMax-baselineVariationGuess,binMax+baselineVariationGuess);
    fit->SetParameter(1,binMax);
    baselineHist->Fit("fit","QR0","",binMax-baselineVariationGuess,binMax+baselineVariationGuess);
    
    Double_t baseline=fit->GetParameter(1);
    Double_t baselineRMS=fit->GetParameter(2);
   
    baselineHist->Delete();
    fit->Delete();
    
    return make_pair(baseline,baselineRMS);
	}
	else {
	  delete baselineHist;
	  return make_pair(0,0);
	}
}
*/

	
	
	
//Plots waveform and a Tline at a time passed in
void plotWaveformAndPulse(vector<Double_t> waveVec,Double_t nsPerSample,Double_t pulseSample) {

	TH1D* waveformPlot = new TH1D("waveformPlot","waveformPlot; time (ns); ADC",waveVec.size(),0,nsPerSample*waveVec.size());
	waveformPlot->SetStats(0);
	for (Int_t i=0; i < waveVec.size(); i++) {
		waveformPlot->SetBinContent(i+1,waveVec.at(i));
	}
	TCanvas* c1 = new TCanvas("c1","c1");
	waveformPlot->Draw();

  TLine pulseLine(pulseSample*nsPerSample,0,pulseSample*nsPerSample,170.);
  pulseLine.SetLineWidth(3);
	pulseLine.Draw();
	
	c1->Modified();
	c1->Update();
	c1->WaitPrimitive("A");
	
	delete c1;
	delete waveformPlot;
}


//Software implementation of the Struck 3316 trigger. Returns a vector of pulses
//that satisfy this trigger condition 
vector<Double_t> getPulsesFromFIR(vector<Double_t> waveVec,Double_t waveformBaseline=0,Int_t windowSize=5,Int_t gapSize=5,Int_t firThresh=100,Double_t fraction=0.5, Int_t holdOffSamples=100, Int_t FIRSamples=100) {

  //Holds vectors of times of found pulses in the waveform
  vector<Double_t> pulseTimes;
  pulseTimes.clear();

  //Deque's will hold moving averages in two windows
  deque<Double_t> leadingWindow;
  deque<Double_t> trailingWindow;
  Double_t sumInLeadingWindow;
  Double_t sumInTrailingWindow;
  
  //FIR vector will hold difference between two moving averages
  vector<Double_t> firVector;
  firVector.clear();
  
  Int_t triggerSample;
  
  //Calculate FIR vector
  for (Int_t i=0; i < waveVec.size(); i++) {
  
    //For first windowSize+gapSize samples, load baseline values
    if (i < windowSize+gapSize) {
      leadingWindow.push_back(waveformBaseline);
      sumInLeadingWindow += waveformBaseline;
    
      trailingWindow.push_back(waveformBaseline);
      sumInTrailingWindow += waveformBaseline;
    }
    else {
      //Push new values to deques & sum
      leadingWindow.push_back(static_cast<Double_t>(waveVec.at(i)));
      sumInLeadingWindow += static_cast<Double_t>(waveVec.at(i) );
    
      trailingWindow.push_back(static_cast<Double_t>(waveVec.at(i-(windowSize+gapSize))));
      sumInTrailingWindow += static_cast<Double_t>(waveVec.at(i-(windowSize+gapSize))); 
    }
    
    //If the vector sizes are > the window size, remove oldest value & from sum 
    if (leadingWindow.size() > windowSize) {
      //Remove older values from deques & sum 
      sumInLeadingWindow -= leadingWindow.front();
      leadingWindow.pop_front();
    }
    if (trailingWindow.size() > windowSize) {
      sumInTrailingWindow -= trailingWindow.front();
      trailingWindow.pop_front();
    }
    
    //Calculate difference 
    firVector.push_back((sumInLeadingWindow-sumInTrailingWindow));
  }
    
  //Now look for instances above threshold
  Int_t i=0;
  while (i < waveVec.size()) {
  
    //We found an event above threshold. 
    if (firVector.at(i) >= firThresh) {

      Double_t maxValue=0;
      Int_t maxSample=0;
      triggerSample=0;
      
      //Search for max value up to FIRSamples after this 
      for (Int_t j=i; j < i+FIRSamples; j++) {
      
        //Make sure we don't extend beyond edge of waveform looking for FIR max
        if (j >= waveVec.size()) {
          break;
        }
        if (firVector.at(j) > maxValue) {
          maxValue=firVector.at(j);
          maxSample=j;
        }
      }
      
      //Check if we actually found a maxValue at position > 0. If not, ignore
      if (maxSample > 0) {
      
        //Now go from maxSample up to i+holdOffSample searching for where fir drops below fraction of max
        for (Int_t j=maxSample; j < i+FIRSamples; j++) {
        
          //If we go beyond edge of waveform looking for fraction of max
          if (j >= waveVec.size()) {
            break;
          }
          if (firVector.at(j) <= fraction*maxValue) {
            triggerSample=j;
            break;
          }
        }
        
        //Make sure we can do interpolation (also checks we actually set triggerSample),
        //and do interpolation to get onset time to sub-ns precision
        if (triggerSample>0) {
          //Get value at and before FIR trigger
          Double_t valueAtTrigger = firVector.at(triggerSample);
          Double_t valueBeforeTrigger = firVector.at(triggerSample-1);
          
          //Interpolation
          Double_t differenceMaxDiv2AndMaxAtTrig = maxValue/2.-valueAtTrigger;
          Double_t differenceBeforeTrigAndAtTrig = valueBeforeTrigger-valueAtTrigger;
          Double_t corrector = differenceMaxDiv2AndMaxAtTrig/differenceBeforeTrigAndAtTrig;

          pulseTimes.push_back(static_cast<Double_t>(triggerSample)-corrector);
        }
      }
      
      i+=holdOffSamples;
    }
    else {
      i++;
    }
  }
  
  /*
  //Plot FIR for debugging
  cout<<"Trigger sample"<<triggerSample<<endl;
  TCanvas* c1 = new TCanvas("c1","c1");
  TH1D* firHist = new TH1D("firHist","",firVector.size(),0,firVector.size());
  for (Int_t i=0; i < firVector.size(); i++) {
    firHist->SetBinContent(i+1,firVector.at(i));
  }
  firHist->Draw();
  c1->Modified();
  c1->Update();
  c1->WaitPrimitive("ABA");
  delete firHist;
  delete c1;
  */
  
  return pulseTimes;
}




//Takes in a vector representing a waveform, a baseline guess, and a rejectThreshold.
//Looks through next baselineSamples following startSample, histograms values within
//rejectThreshold of the guess, and fits them. Returns zero if there's something 
//weird with this baseline
Double_t getLocalBaseline(vector<Double_t> waveVec,Double_t baselineGuess,Int_t rejectThreshold,Int_t startSample,Int_t baselineSamples) {

  if (startSample+baselineSamples >= waveVec.size()) {
    cout<<"Not enough samples left in waveform to compute baseline!"<<endl;	  
    cout<<"Start sample: "<<startSample<<", baselineSamples: "<<baselineSamples<<", waveform size: "<<waveVec.size()<<endl;
	  return 0.;
  }
  
  //We'll histogram waveform values into this
  TH1D* baselineHist = new TH1D("baselineHist","",rejectThreshold,baselineGuess-rejectThreshold,baselineGuess+rejectThreshold);
  
  //Fill values in hist located within rejectThreshold of baseline
  Int_t nEntries=0;
  for (Int_t i=startSample; i < startSample+baselineSamples; i++) {
    if ( fabs(static_cast<Double_t>(waveVec.at(i))-baselineGuess) < rejectThreshold) {
      baselineHist->Fill(static_cast<Double_t>(waveVec.at(i)));
      nEntries++;
    }
  }
  
  
  //Draw hist, for debugging
  TCanvas* c1 = new TCanvas("c1","c1");
  baselineHist->Draw();
  c1->Modified();
  c1->Update();
  c1->WaitPrimitive("ABA");
  delete c1;
  
  
  //Require a minimum of 10 entries to be found in the expected waveform region.
  //If there aren't that many entries, we're not going to use this event. Return
  //baseline of zero
  if (nEntries<20) {
    delete baselineHist;
    return 0;
  }
  
  TF1* fit = new TF1("fit","gaus",baselineGuess-rejectThreshold,baselineGuess+rejectThreshold);
  baselineHist->Fit("fit","QR0","",baselineGuess-rejectThreshold,baselineGuess+rejectThreshold);
  
  Double_t baseline=fit->GetParameter(1);
  
  delete baselineHist;
  delete fit;
  
  return baseline;
}


//Pass in a vector holding the waveform, a baseline value, a start sample, a 
//number of samples to look through, and a fraction. Finds max value in waveform,
//and then starts at the max value sample and traverses backwards until we get 
//below cfdFraction*maxValue. Returns that sample (no interpolation done)
Double_t getPHVCFDSample(vector<Double_t> wave,Double_t baseline,Int_t startSample,Int_t numSamples,Double_t cfdFraction) {
	
	if (startSample+numSamples > wave.size()) {
	  cout<<"Region passed in to PHV CFD algorithm extends beyond waveform!"<<endl;
	  cout<<"Start sample: "<<startSample<<", numSamples: "<<numSamples<<", waveform size: "<<wave.size()<<endl;
	  return 0.;
	}
   
  Double_t maxValue=0; 
  Int_t maxSample=0;
  Int_t cfdSample=0;
  
	for (Int_t i=startSample; i < startSample+numSamples; i++) {
		if (wave.at(i)-baseline > maxValue) {
			maxSample=i;
			maxValue=wave[i]-baseline;
		}
	}
	
  //Now starting at the max sample, go back until we find the first sample less than fraction of the max values
  /*
  for (Int_t i=maxSample; i >= startSample; i--) {
    if (wave.at(i)-baseline <= cfdFraction*maxValue) {
      cfdSample=i;
      break;
    }
  }
  */
	
  //Now starting at first sample, go forward until we find the first sample >= the fraction the max value
  for (Int_t i=startSample; i < startSample+numSamples; i++) {
    if (wave.at(i)-baseline >= cfdFraction*maxValue) {
      cfdSample=i;
      break;
    }
  }
	
	return cfdSample;
}


//Pass in a vector representing the waveform, a baseline, a startSample, and a 
//number of samples to look through. Finds the peak high value in the region
//passed in, baseline-subtracted.
Double_t getPHV(vector<Double_t> waveVec,Double_t baseline,Int_t startSample,Int_t numSamples,Int_t searchBackwards) {

	Double_t maxVal=0;
		
	if (searchBackwards==0) {
		if (startSample+numSamples >= waveVec.size()) {
			numSamples = waveVec.size()-startSample-1;
		}
		
		for (Int_t i=startSample; i < startSample+numSamples; i++) {
			if (waveVec.at(i) > maxVal) {
					maxVal=waveVec.at(i);
			}
		}
		return maxVal-baseline;
	}
	else {
		if (startSample-numSamples < 0) {
			numSamples = startSample;
		}
		for (Int_t i=startSample; i > startSample-numSamples; i--) {
			if (waveVec.at(i) > maxVal) {
				maxVal=waveVec.at(i);
			}
		}
		return maxVal-baseline;
	}
}


//Passing in a vector representing the waveform, a baseline value, a start sample,
//and a number of samples to integrate. Computes the integral (sum of ADC counts
//in region, NOT taking bin size into account)
Double_t getIntegral(vector<Double_t> waveVec,Double_t baseline,Int_t startSample,Int_t numSamples) {

	if (startSample+numSamples >= waveVec.size()) {
	  numSamples = waveVec.size()-startSample-1;
	}
	
	Double_t integral=0.;
	for (Int_t i=startSample; i < startSample+numSamples; i++) {
		integral += (waveVec.at(i)-baseline);
	}
	
	return integral;
}



//Finds pulses by edge counting. Pass in a vector representing the waveform, 
//a value of the waveform baseline, a threshold, and a number of hold-off samples.
//Returns a vector of times when that threshold is first exceeded, not searching
//for pulses again until after holdOffSamples has passed.
vector<Double_t> getPulsesFromEdgeCounting(vector<Double_t> waveVec,Double_t waveformBaseline,Int_t threshold,Int_t holdOffSamples) {

  vector<Double_t> pulses;
  pulses.clear();
  
  Int_t sampleNum=0;
  while (sampleNum < waveVec.size()) {
    if (static_cast<Double_t>(waveVec.at(sampleNum)) - waveformBaseline >= threshold) {
      pulses.push_back(static_cast<Double_t>(sampleNum));
      sampleNum+=holdOffSamples;
    }
    else {
      sampleNum++;
    }
  }
  
  return pulses;
}


//getTsincVector + helper functions:
//sinc definition
double sinc(double x) {
  if (x == 0.0) {
	return 1.0;
  }
  else {
	return sin(x)/x;
  }
}

//tstinc definition
double tsinc(double x,int n, int t) {
  return sinc(x*M_PI/static_cast<Double_t>(n))*exp(-1*pow(x/static_cast<Double_t>(t),2));
}


//Returns an interpolated vector
vector<Double_t> getTsincVector(vector<Double_t> wave,Double_t baseline,Int_t startSample,Int_t numSamples, Int_t n=8, Int_t t=30, Int_t l=6) {

  vector<Double_t> sincVector;
  sincVector.clear();
  if (startSample < l) {
	cout<<"Not enough samples before startSample to do interpolation"<<endl;
	cout<<"startSample: "<<startSample<<", samples required for interpolation: "<<l<<endl;
	return sincVector;
  }

  if (startSample + numSamples +l >= wave.size()) {
	cout<<"Not enough samples remaining in waveform to do interpolation"<<endl;
	cout<<"startSample: "<<startSample<<", nSamples: "<<numSamples<<", samples required for interpolation: "<<l<<", waveformSize: "<<wave.size()<<endl;
	return sincVector;
  }

  //TH1D* originalVector = new TH1D("originalVector","",numSamples,0,numSamples);
  for (Int_t j=startSample; j < startSample+numSamples; j++) {

  //originalVector->SetBinContent(j+1,wave.
      
	//Push non-interpolated point to vector
	sincVector.push_back(wave.at(j));
			
	for (Int_t k=1; k < n; k++) {
		Double_t interpVal = 0.;

		for (Int_t i=0; i < l-1; i++) {
			interpVal += wave.at(j-i)*tsinc(i*n+k,n,t)+wave.at(j+1+i)*tsinc((i+1)*n-k,n,t);
		}
		//cout<<j+k<<","<<interpVal<<endl;
		sincVector.push_back(interpVal);
	}
  }
  return sincVector;
}



//Vector version of Grayson Rich's CMA filter 
vector<Double_t> getCMAFilter( vector<Double_t> wave,Int_t halfWidth,Double_t preloadValue, Double_t rejectThreshold) {
  
    if (halfWidth > wave.size()) {
      cout<<"Half-width greater than waveform length!"<<endl;
      vector<Double_t> a;
      a.clear();
      return a;
    }
  
    deque<Double_t> movingBaselineFilter;
    movingBaselineFilter.clear();
    Double_t movingBaselineSum=0;

    vector<Double_t>CMAWave;
    CMAWave.clear();

    //Preload filter with baselines guess
    for (Int_t i=0; i < halfWidth; i++) {
      movingBaselineFilter.push_back(preloadValue);
      movingBaselineSum += preloadValue;
    }
  
    //Now step through first halfWidth values of waveform, adding to moving average 
    //if within rejectThreshold of the current moving average
    for (Int_t i=0; i < halfWidth; i++) {
      //Calculate current moving average
      Double_t movingAverage = movingBaselineSum/movingBaselineFilter.size();
      
      //Check if current value within rejectThreshold of current moving average 
      if ( fabs(wave.at(i)-movingAverage) < rejectThreshold) {
        movingBaselineFilter.push_back(wave.at(i));
        movingBaselineSum += wave.at(i);
      }
    }
    
    //Now moving average full. Start at element 0, and compute moving average 
    //(which is the average of the half-width elements on either side). 
    for (Int_t i=0; i < wave.size(); i++) {
    
      //Check if there is an element i+halfWidth ahead. If it's beyond the scope
      //of the waveform, start popping off elements as we advance update CMA
      if (i+halfWidth < wave.size()) {
      
        //Calculate current moving average
        Double_t movingAverage = movingBaselineSum/movingBaselineFilter.size();
        
        //Check the value at i+halfWidth is within the threshold
        if ( fabs(wave.at(i+halfWidth)-movingAverage) < rejectThreshold) {

          movingBaselineFilter.push_back(wave.at(i+halfWidth));
          movingBaselineSum += wave.at(i+halfWidth);
        }
          
        //Check if size of filter is greater than 2*halfWidth+1. If so, remove oldest element 
        if (movingBaselineFilter.size() >= 2*halfWidth+1) {
          movingBaselineSum -= movingBaselineFilter.front();
          movingBaselineFilter.pop_front();
        }
      }
      else {
        movingBaselineSum -= movingBaselineFilter.front();
        movingBaselineFilter.pop_front();
      }
      
      //push back moving average
      Double_t movingAverage = movingBaselineSum/movingBaselineFilter.size();
      CMAWave.push_back(movingAverage);
    }
    
    //Debugging, plot 
    /*
    TH1D* waveHist = new TH1D("wave","",wave.size(),0,wave.size());
    TH1D* cmaHist = new TH1D("cmaHist","",CMAWave.size(),0,CMAWave.size());
    cmaHist->SetLineColor(2);
    for (Int_t i=0; i < wave.size(); i++) {
      waveHist->SetBinContent(i+1,wave.at(i));
      cmaHist->SetBinContent(i+1,CMAWave.at(i));
    }
    TCanvas* c1 = new TCanvas("c1","c1");
    waveHist->Draw();
    cmaHist->Draw("same");
    c1->Modified();
    c1->Update();
    c1->WaitPrimitive("ABA");
    delete waveHist;
    delete cmaHist;
    delete c1;
    */
    return CMAWave;
} 


Double_t getMeanTime(vector<Double_t> waveVec,Double_t baseline,Int_t startSample,Int_t nSamples,Double_t maxValue) {

  Int_t meanTimeSamples = 15;
  Double_t meanTimeNumerator=0;
  Double_t meanTimeDenominator=0;
  
  //Find 50% of maxValue
  Int_t halfMax=0;
  for (Int_t i=startSample; i < startSample+nSamples; i++) {
    if (waveVec.at(i) >= 0.5*maxValue) {
      halfMax=i;
      break;
    }
  }
  
  //Get mean time 
  if ((halfMax == 0)||(halfMax+meanTimeSamples > waveVec.size())) {
    return 0.;
  }
  else {
    for (Int_t i=0; i < meanTimeSamples; i++) {
      meanTimeNumerator += i*waveVec.at(halfMax+i);
      meanTimeDenominator += waveVec.at(halfMax+i);
    }
  }
  
  if (meanTimeDenominator != 0) {
    return meanTimeNumerator/meanTimeDenominator;
  }
  else {
    return 0;
  }
}


Int_t checkIfWaveformSaturates(vector<Double_t> waveVec) {
	for (Int_t i=0; i < waveVec.size(); i++) {
		if (waveVec.at(i) >= pow(2,14)-1) {
			return 1;
		}
		else {
			return 0;
		}
	}
}


//Find the max, then find the first sample below the baseline after that max
//Doesn't do interpolation, designed to be run on interpolated waveform 
Double_t getBPMZeroCrossing(vector<Double_t> waveVec,Double_t baseline,Int_t startSample,Int_t numSamples) {

	Double_t maxVal=0;
	Int_t maxIndex=0;
	for (Int_t i=0; i < numSamples; i++) {
		if (waveVec.at(i) > maxVal) {
			maxVal=waveVec.at(i);
			maxIndex=i;
		}
	}
	
	Double_t zeroCrossingIndex=0;
	for (Int_t i=maxIndex; i < numSamples; i++) {
		if (waveVec.at(i) <= baseline) {
			zeroCrossingIndex=static_cast<Double_t>(i);
			return zeroCrossingIndex;		
		}
	}
}

Double_t getSplineZeroCrossing(vector<Double_t> waveVec, Double_t baseline, Int_t startSample, Int_t numSamples) {
	
	Double_t maxVal=0;
	Int_t maxIndex=0;
	for (Int_t i=startSample; i < startSample+numSamples; i++) {
		if (waveVec.at(i) > maxVal) {
			maxVal=waveVec.at(i);
			maxIndex=i;
		}
	}
	Double_t minVal=16000;
	Int_t minIndex=0;
	for (Int_t i=maxIndex; i < maxIndex+numSamples; i++) {
		if (waveVec.at(i) < minVal) {
			minVal=waveVec.at(i);
			minIndex=i;
		}
	}
	TH1D* splHist = new TH1D("splHist","splHist",(minIndex-maxIndex),maxIndex,minIndex);
	for (Int_t i=maxIndex; i < minIndex; i++) {
		splHist->SetBinContent(splHist->FindBin(i),waveVec.at(i));
	}
	
	TSpline3 spl(splHist);
	for (Double_t i=maxIndex; i < minIndex; i+=0.1) {
		if (spl.Eval(i) <=baseline) {
			delete splHist;
			return i;
		}
	}
}


void processData(TString inputFilename, TString outputFilename){

  //-------------FLAGS----------------//
  
  //Save waveforms of the scatterer
  Int_t saveScattererWaveform=1;
  
  //Info about the run so we know what kind of channels to process 
  Int_t hasBPM=1;
  Int_t hasBDs=1;
  Int_t hasScatterer=1;
	
  //0=NIM, 1=BIPOLAR
  Int_t BPMType=0;
  
  //Set to 1 for TOF run, used different pre-trigger delays
  Int_t isTOFRun=0;
  
  
//-------------DAQ SETTINGS----------------//
  Double_t nsPerSample=4.;
  
  //BD trigger settings, used for production run
  Int_t BD_riseTime=4.;
  Int_t BD_gapTime=4.;
  Int_t BD_firThresh=100.;
  Int_t BD_preTriggerDelay;
  if (isTOFRun==1) {
	BD_preTriggerDelay=155; //155 used for TOF data, 
  }
  else {
	BD_preTriggerDelay=195; //195 for standard run & calibrations
  }
  
  //Scatterer trigger settings, used for source/bgnd run
  //Also used for identifying pulses
  Int_t scattererTrigger_riseTime=4;
  Int_t scattererTrigger_gapTime=4;
  Int_t scattererTrigger_firThresh=400;
  //Used both in scatterer calibrations and production run
  Int_t scatterer_preTriggerDelay=195;

  
  //BPM settings
  Int_t BPM_preTriggerDelay;
  if (isTOFRun==1) {
	BPM_preTriggerDelay=155; //155 used for TOF data, 
  }
  else {
	BPM_preTriggerDelay=195; //195 for standard run & calibrations
  }
  
  //Beam settings - Samples
  Int_t beam_samplesBetweenBPMPulses=100;
  
  //Channel settings
  vector<Int_t> LSChannels;
  if (hasBDs==1) {
  	LSChannels = {8,9,10,11,12,13,14,15};
  }
  Int_t zeroDegreeChannel=2;
  Int_t scattererChannel = 0;
  Int_t BPMChannel = 1;


  //-------------ANALYSIS SETTINGS----------------//
  
  
  //BD parameters
  Int_t BD_baselineVariationGuess = 10; //Guess of how much the baseline varies by, in ADC
  Int_t BD_preTraceIntegralSamples = 25; //Samples to use for calculating the preTraceIntegral
  Int_t BD_preOnsetIntegralSamples = 10; //Samples to integrate before found onset
  Int_t BD_preOnsetSamplesRequired = BD_preTraceIntegralSamples+BD_preOnsetIntegralSamples; //Samples required before onset in waveform for pulse to be valid
  Int_t BD_preOnsetInterpolationSamples = 10;
  Int_t BD_postOnsetIntegralSamples = 90; //Samples to integrate
  Int_t BD_tailIntegralSamples = 85; //Tail samples to integrate
  Int_t BD_postOnsetIntegralSamplesRequired = BD_postOnsetIntegralSamples; //Samples required after onset in waveform for pulse to be valid
  Int_t BD_holdOffSamples = 250+BD_preOnsetSamplesRequired; //Hold-off samples when looking for subsequent pulses
  Double_t BD_cfdFraction=0.2; //CFD fraction for onset definition 
  Int_t BD_FIRSamples=20; //Samples after onset to look for zero-crossing, only for Bipolar
  
  /////////////////////////////////
  //Scatterer Analysis Parameters//
  /////////////////////////////////
  Int_t scatterer_preTraceIntegralSamples=10;
  Int_t scatterer_baselineVariationGuess=10;
  Int_t scatterer_preOnsetIntegralSamples=10;
  Int_t scatterer_preOnsetSamplesRequired=scatterer_preTraceIntegralSamples+scatterer_preOnsetIntegralSamples;
  Int_t scatterer_preOnsetInterpolationSamples=10;
  Int_t scatterer_postOnsetIntegralSamples=100;
  Int_t scatterer_postOnsetIntegralSamplesRequired=scatterer_postOnsetIntegralSamples;
  Int_t scatterer_holdOffSamples=100;
  Double_t scatterer_cfdFraction=0.2;
  Int_t scatterer_FIRsamples=4;
  
  //BPM
  Int_t BPM_baselineVariationGuess=50;
  Int_t BPM_threshold;
  if (BPMType==0) {
    BPM_threshold=4000;
  }
  else {
    BPM_threshold=1500;
  }
  
  Int_t BPM_preOnsetSamplesRequired = 4;
  Int_t BPM_postOnsetSamplesRequired = 9;
  Double_t BPM_cfdFraction = 0.4;
  Int_t BPM_holdOffSamples=80; //Nominally 100 samples apart

  //-------------CMA SETTINGS----------------//
  
  Int_t BD_halfWidth=50;
  Int_t BD_rejectThresholdSigma=2.;
  Int_t scatterer_halfWidth=45;
  Double_t scatterer_rejectThresholdSigma=2.5;
  

  //-------------TSINC SETTINGS----------------//
  Int_t l=6;
  Int_t fTaper=30;
  Int_t nInterPoints=16; //# points between samples, 16=0.25ns precision

  

  //-------------FILE I/O----------------//


  /////////////////////////////////////////////////////////////
  //Read in file name, make output file/folder (if necessary)//
  /////////////////////////////////////////////////////////////
  //Get filenames passed in
  TFile* input = new TFile(inputFilename,"READ");
  cout<<"Output file is "<<outputFilename<<endl;
  TFile* output = new TFile(outputFilename,"RECREATE");

  //////////////////////////////////
  //Create TTree, set up for read//
  //////////////////////////////////
  TTree* sis3316tree = (TTree*)input->Get("sis3316tree");
	
  //From parsed output 
  UShort_t channelID;
  Bool_t pileupFlag;
  ULong64_t timestamp;
  UShort_t peakHighValue;
  UShort_t peakHighIndex;
  UInt_t nSamples;
  UInt_t accumulatorSum[8];
  UShort_t shortWaveform[65536];
    
  //Signal info
  sis3316tree->SetBranchAddress("channelID",&channelID);
  sis3316tree->SetBranchAddress("pileupFlag",&pileupFlag);
  sis3316tree->SetBranchAddress("timestamp",&timestamp);
  sis3316tree->SetBranchAddress("peakHighValue",&peakHighValue);
  sis3316tree->SetBranchAddress("peakHighIndex",&peakHighIndex);
  sis3316tree->SetBranchAddress("nSamples",&nSamples);
  sis3316tree->SetBranchAddress("accumulatorSum",accumulatorSum);
  sis3316tree->SetBranchAddress("waveform",shortWaveform);

  Long64_t numSignals=int(sis3316tree->GetEntries()/1);
  cout<<"Found "<<numSignals<<" signals"<<endl;



///////////////////////////PROCESS BPM SIGNALS FIRST////////////////////////////



  ////////////////////////////
  //Prepare BPM output trees//
  ////////////////////////////
  //BPM Tree//
  Double_t BPM_waveformStartTime;
  Double_t BPM_onsetTime;
  Double_t BPM_peakHeight=0;
  Double_t BPM_timeToPrevBPM=0;
  Double_t BPM_prevTime=0;
  TTree* BPMTree = new TTree("BPMTree","BPMTree tree");
  BPMTree->SetMaxTreeSize(2000000000000LL);
  BPMTree->Branch("BPM_waveformStartTime", &BPM_waveformStartTime, "BPM_waveformStartTime/D");
  BPMTree->Branch("BPM_onsetTime", &BPM_onsetTime, "BPM_onsetTime/D");
  BPMTree->Branch("BPM_peakHeight", &BPM_peakHeight, "BPM_peakHeight/D");
  BPMTree->Branch("BPM_timeToPrevBPM", &BPM_timeToPrevBPM, "BPM_timeToPrevBPM/D");

  /////////////////////////////////////////////////////////////
  //Optimizations for stepping through and finding BPM pulses//
  /////////////////////////////////////////////////////////////
  TBranch* channelIDBranch = sis3316tree->GetBranch("channelID");
  TBranch* timestampBranch = sis3316tree->GetBranch("timestamp");
  TBranch* nSamplesBranch = sis3316tree->GetBranch("nSamples");
  TBranch* waveformBranch = sis3316tree->GetBranch("waveform");
  //Tree read optimizations
  Int_t cachesize = 500000000;
  sis3316tree->SetCacheSize(cachesize);
  sis3316tree->AddBranchToCache(channelIDBranch,1);
  sis3316tree->AddBranchToCache(timestampBranch,1);
  sis3316tree->AddBranchToCache(nSamplesBranch,1);
  sis3316tree->AddBranchToCache(waveformBranch,1);
  sis3316tree->StopCacheLearningPhase();

  //////////////////////////////////////////////////
  //Prepare timing vector to store BPM pulse times//
  /////////////////////////////////////////////////
  vector<Double_t> BPM_timingVector;
  //Step through BPM pulses
  if (hasBPM==1) {
    for (Int_t entry=0; entry < numSignals; entry++) {
    
      if (entry%10000==0) {
        cout<<"Looking for BPM signals, on event "<<entry<<" of "<<numSignals<<endl;
      }
    
      sis3316tree->LoadTree(entry);
      channelIDBranch->GetEntry(entry);
      
      if (channelID==BPMChannel) {
        
          nSamplesBranch->GetEntry(entry);
          timestampBranch->GetEntry(entry);
          waveformBranch->GetEntry(entry);
          
          //Corresponds to trigger time 
          BPM_waveformStartTime = static_cast<Double_t>(timestamp)*nsPerSample;
          //Subtract off the preTriggerDelay to get the start time of the waveform
          BPM_waveformStartTime -= static_cast<Double_t>(BPM_preTriggerDelay)*nsPerSample;
    
          //Convert array of shorts into a vector
          vector<Double_t> waveformVector(shortWaveform,shortWaveform+nSamples);
    
          //Get baseline (mode of entire waveform)
          //pair<Double_t,Double_t> baselineContainer = getBaselineGuess(waveformVector,BPM_baselineVariationGuess);
          //Double_t waveformBaseline = baselineContainer.first;
          //Double_t waveformBaselineRMS = baselineContainer.second;
          
          //Set baseline to zero, not enough non-NIM pulse samples to get accurate calculation
          Double_t waveformBaseline=0;
				
          //Find instances where BPM went above threshold
          vector<Double_t> pulses = getPulsesFromEdgeCounting(waveformVector,waveformBaseline,BPM_threshold,BPM_holdOffSamples);
  
	  
        
          //Interpolate to get better timing resolution
          for (Int_t i=0; i < pulses.size(); i++) {
            
            //plotWaveformAndPulse(waveformVector,nsPerSample,pulses.at(i)); //Plot BPM waveform w/pulse
                
            //Require BPM_preOnsetSamplesRequired+l samples before onset, BPM_postOnsetSamplesRequired+l samples after 
            if ((pulses.at(i) >= BPM_preOnsetSamplesRequired+l) && (pulses.at(i)+BPM_postOnsetSamplesRequired+l < waveformVector.size())) {
              
              //interpolate region where we expect to find onset 
              vector<Double_t> interpWaveVec;
              interpWaveVec.clear();
              interpWaveVec = getTsincVector(
                waveformVector,
                0,
                pulses.at(i)-BPM_preOnsetSamplesRequired,
                BPM_preOnsetSamplesRequired+BPM_postOnsetSamplesRequired,
                nInterPoints,
                fTaper,
                l
               );
            
            Double_t interpOnsetSample=0;
            if (BPMType==0) {
              //Get subsample where interpolated waveform crosses 20% max value 
              interpOnsetSample = getPHVCFDSample(interpWaveVec,waveformBaseline,0,interpWaveVec.size(),BPM_cfdFraction);
            }
            
            if (BPMType==1) {
              interpOnsetSample = getBPMZeroCrossing(interpWaveVec,waveformBaseline,0.,interpWaveVec.size());
            }
            
            //plotWaveformAndPulse(interpWaveVec,nsPerSample/static_cast<Double_t>(nInterPoints),interpOnsetSample);	 //Plot interpolated waveform + pulse
          
            //Convert from the interpolated samples to uninterpolated samples
            interpOnsetSample /= static_cast<Double_t>(nInterPoints);
            //Add interpolated onset to where we started out interpolation 
            interpOnsetSample += (pulses.at(i)-BPM_preOnsetSamplesRequired);
            //Convert from samples to time
            BPM_onsetTime = interpOnsetSample*nsPerSample;
               
            //Get onset sample of refined onset time to make sure we won't go out of bounds when looking for the peak height
            Int_t onsetSample = round(BPM_onsetTime/nsPerSample);
            
            //Make sure we have enough samples before and after onset sample
            if ((onsetSample > BPM_preOnsetSamplesRequired) && (onsetSample + BPM_postOnsetSamplesRequired < nSamples)) {
              if (BPMType==1) {
                BPM_peakHeight=getPHV(waveformVector,0,onsetSample,BPM_postOnsetSamplesRequired,1);
              }
              else {
                BPM_peakHeight=getPHV(waveformVector,0,onsetSample,BPM_postOnsetSamplesRequired,0);
              }
              //plotWaveformAndPulse(waveformVector,nsPerSample,onsetSample);
              
              BPM_timingVector.push_back(BPM_waveformStartTime+BPM_onsetTime);
              if (BPM_prevTime!=0) {
                BPMTree->Fill();      
              }
            }
            BPM_timeToPrevBPM=BPM_waveformStartTime+BPM_onsetTime-BPM_prevTime;
            BPM_prevTime=BPM_waveformStartTime+BPM_onsetTime;
          }
        }
        
        /*
        //Non-interpolation version
        for (Int_t i=0; i < pulses.size(); i++) {
          if ((pulses.at(i) > BPM_preOnsetSamplesRequired) && (pulses.at(i) + BPM_postOnsetSamplesRequired < nSamples)) {

            BPM_onsetTime = pulses.at(i)*nsPerSample;
            BPM_timingVector.push_back(BPM_waveformStartTime+BPM_onsetTime);
            BPM_peakHeight=getPHV(waveformVector,0,pulses.at(i),BPM_postOnsetSamplesRequired);
            BPMTree->Fill();    
          }
        }
        */
       
      }
    }
  }

  
////////////////////////////PROCESS LS SIGNALS NEXT/////////////////////////////



  ////////////////////////
  //Set up analysis tree//
  ////////////////////////
  Int_t LS_channel;
  Int_t LS_pileUpInScatterIntegrationWindow;
  Int_t LS_saturated;
  Double_t LS_waveformStartTime;
  Double_t LS_onsetTime;
  Double_t LS_timeToBPM;
  Double_t LS_preTraceIntegral;
  Double_t LS_preTraceIntegralDev;
  Double_t LS_integral;
  Double_t LS_peakHeight;
  Double_t LS_psd;
  
  Int_t scatterer_foundPulse;
  Double_t scatterer_waveformStartTime;
  Double_t scatterer_onsetTime;
  Double_t scatterer_timeToBPM;
  Double_t scatterer_timeToBD;
  Double_t scatterer_waveformStartTimeToNextBD;
  Double_t scatterer_preTraceIntegral;
  Double_t scatterer_integral;
  Double_t scatterer_peakHeight;
  Double_t scatterer_meanTime;
  Double_t scatterer_baseline;
  Double_t scatterer_noise;
  vector<UShort_t> scatterer_waveform;

  //*NEW* Parameters for a fixed window integration.
  Int_t scatterer_startSample = 165; //Calibration runs.
  Double_t bpm_Offset = 120.0; //Production runs. Distance of "good" pulses from the next BPM.

  TTree* analysisTree = new TTree("analysisTree","analysisTree");
  analysisTree->SetMaxTreeSize(2000000000000LL);
  
  //If we have BDs, we'll need these branches
  if (hasBDs==1) {
    analysisTree->Branch("LS_channel", &LS_channel, "LS_channel/I");
    analysisTree->Branch("LS_waveformStartTime", &LS_waveformStartTime, "LS_waveformStartTime/D");
    analysisTree->Branch("LS_onsetTime", &LS_onsetTime, "LS_onsetTime/D");
    analysisTree->Branch("LS_timeToBPM", &LS_timeToBPM, "LS_timeToBPM/D");
	  analysisTree->Branch("LS_preTraceIntegral", &LS_preTraceIntegral, "LS_preTraceIntegral/D");
	  analysisTree->Branch("LS_preTraceIntegralDev", &LS_preTraceIntegralDev, "LS_preTraceIntegralDev/D");
  	analysisTree->Branch("LS_integral", &LS_integral, "LS_integral/D");
  	analysisTree->Branch("LS_peakHeight", &LS_peakHeight, "LS_peakHeight/D");
  	analysisTree->Branch("LS_psd", &LS_psd, "LS_psd/D");
    analysisTree->Branch("LS_saturated", &LS_saturated, "LS_saturated/I");
  }
  //If we have a scatterer, we'll need these branches
  if (hasScatterer==1) {
    if (hasBDs==1) {
      analysisTree->Branch("LS_pileUpInScatterIntegrationWindow", &LS_pileUpInScatterIntegrationWindow, "LS_pileUpInScatterIntegrationWindow/I");
    }
		analysisTree->Branch("scatterer_foundPulse", &scatterer_foundPulse, "scatterer_foundPulse/I");
		analysisTree->Branch("scatterer_waveformStartTime", &scatterer_waveformStartTime, "scatterer_waveformStartTime/D");
		analysisTree->Branch("scatterer_onsetTime", &scatterer_onsetTime, "scatterer_onsetTime/D");
		analysisTree->Branch("scatterer_timeToBPM", &scatterer_timeToBPM, "scatterer_timeToBPM/D");
		analysisTree->Branch("scatterer_timeToBD", &scatterer_timeToBD, "scatterer_timeToBD/D");
		analysisTree->Branch("scatterer_waveformStartTimeToNextBD", &scatterer_waveformStartTimeToNextBD, "scatterer_waveformStartTimeToNextBD/D");
		analysisTree->Branch("scatterer_preTraceIntegral", &scatterer_preTraceIntegral, "scatterer_preTraceIntegral/D");
		analysisTree->Branch("scatterer_integral", &scatterer_integral, "scatterer_integral/D");
		analysisTree->Branch("scatterer_noise", &scatterer_noise, "scatterer_noise/D");
		analysisTree->Branch("scatterer_peakHeight", &scatterer_peakHeight, "scatterer_peakHeight/D");
		analysisTree->Branch("scatterer_meanTime", &scatterer_meanTime, "scatterer_meanTime/D");
		analysisTree->Branch("scatterer_baseline", &scatterer_baseline, "scatterer_baseline/D");
		if (saveScattererWaveform == 1) {
			analysisTree->Branch("scatterer_waveform", &scatterer_waveform);
		}
	}

	
  //Prepare vectors for storing LS info. We'll preload the found pulses here and later
  //Fill in the tree with the appropriate scatterer signal if found 
  vector<Int_t> LS_channelVector;
  vector<Int_t> LS_saturatedVector;
  vector<Double_t> LS_waveformStartTimeVector;
  vector<Double_t> LS_onsetTimeVector;
  vector<Double_t> LS_timeToBPMVector;
  vector<Double_t> LS_preTraceIntegralVector;
  vector<Double_t> LS_preTraceIntegralDevVector;
  vector<Double_t> LS_integralVector;
  vector<Double_t> LS_peakHeightVector;
  vector<Double_t> LS_psdVector;

	Int_t bpmIndex=0;
	
  if (hasBDs==1) {
    
    //Step through LS pulses
    for (Int_t entry=0; entry < numSignals; entry++) {
    
      if (entry%10000==0) {
        cout<<"Looking for LS signals, on event "<<entry<<" of "<<numSignals<<endl;
      }
    
      sis3316tree->LoadTree(entry);
      channelIDBranch->GetEntry(entry);
      
        //cout<<channelID<<endl;
      //Check if this is an LS channel we're processing 
			if ( (find(LSChannels.begin(), LSChannels.end(), channelID) != LSChannels.end()) || (channelID==zeroDegreeChannel)) {
        //If so load other branches
        nSamplesBranch->GetEntry(entry);
        timestampBranch->GetEntry(entry);
        waveformBranch->GetEntry(entry);
      
        //Convert waveform array to vector for easier processing 
        vector<Double_t> waveformVector(shortWaveform,shortWaveform+nSamples);
       
        //Corresponds to trigger time 
        LS_waveformStartTime = static_cast<Double_t>(timestamp)*nsPerSample;
        //Subtract off the preTriggerDelay to get the start time of the waveform
        LS_waveformStartTime -= static_cast<Double_t>(BD_preTriggerDelay)*nsPerSample;
				
				//Get baseline guess
				pair<Double_t,Double_t> baselineContainer = getBaselineGuess(waveformVector,BD_baselineVariationGuess);
				Double_t waveformBaseline = baselineContainer.first;
				Double_t waveformBaselineRMS = baselineContainer.second;
				
				//baseline of zero indicates we had problems, most likely due to pile-up. Ignore these pulses.
				if (waveformBaseline != 0) {
			
					//Check raw waveform to see if we have saturation
					LS_saturated=checkIfWaveformSaturates(waveformVector);
			
					//Get vector of CMA for baseline subtraction
					vector<Double_t> cmaFilter = getCMAFilter(waveformVector,BD_halfWidth,waveformBaseline,BD_rejectThresholdSigma*waveformBaselineRMS);
					//Subtract CMA to subtract baseline 
					for (Int_t i=0; i < waveformVector.size(); i++) {
						waveformVector.at(i) = waveformVector.at(i)-cmaFilter.at(i);
					}
					
				 
					//Get vector of pulses found in waveform. For now, we copy the DAQ trigger settings
					vector<Double_t> pulses = getPulsesFromFIR(waveformVector,0,BD_riseTime,BD_gapTime,BD_firThresh,0.5,BD_holdOffSamples,BD_FIRSamples); 
					//Now step through found pulses to do timing refinement
					for (Int_t i=0; i < pulses.size(); i++) {       
						
						//cout<<channelID<<endl;
						//plotWaveformAndPulse(waveformVector,nsPerSample,pulses.at(i)); //Plot BPM waveform w/pulse
						
						//Check we can traverse BD_preOnsetIntegralSamples+l samples back in time and l samples forward in time to do interpolation
						if ((pulses.at(i)-l > BD_preOnsetInterpolationSamples)&&(pulses.at(i)+l < waveformVector.size())) {
						
							//interpolate region where we expect to find onset 
							vector<Double_t> interpWaveVec;
							interpWaveVec.clear();
							interpWaveVec = getTsincVector(waveformVector,0,pulses.at(i)-BD_preOnsetInterpolationSamples,BD_preOnsetInterpolationSamples,nInterPoints,fTaper,l);
						
							//Get subsample where waveform crosses 20% max value 
							Double_t interpOnsetSample = getPHVCFDSample(interpWaveVec,0,0,interpWaveVec.size(),BD_cfdFraction);
							//plotWaveformAndPulse(interpWaveVec,0.25,interpOnsetSample);
							interpOnsetSample /= static_cast<Double_t>(nInterPoints);
							interpOnsetSample += (pulses.at(i)-BD_preOnsetInterpolationSamples);
							LS_onsetTime = interpOnsetSample*nsPerSample; //only used for timing
								 
							//Get onset sample, only use sub-sample timing not for integrating, etc.
							Int_t onsetSample = round(LS_onsetTime/nsPerSample);
							
							//Make sure we have enough samples before and after onset sample
							if ((onsetSample > BD_preOnsetSamplesRequired) && (onsetSample + BD_postOnsetIntegralSamplesRequired < nSamples)) {
							
						    //plotWaveformAndPulse(waveformVector,nsPerSample,onsetSample);
						
								//Find PHV within same region we're going to integrate for 
								LS_peakHeight = getPHV(waveformVector,0,onsetSample,BD_postOnsetIntegralSamples,0);
									
								LS_preTraceIntegral = getIntegral(waveformVector,0,onsetSample-BD_preOnsetSamplesRequired,BD_preTraceIntegralSamples);
									
								//Get integral 
								LS_integral = getIntegral(waveformVector,0,onsetSample-BD_preOnsetIntegralSamples,BD_preOnsetIntegralSamples+BD_postOnsetIntegralSamples);
									
								//Get tail integral 
								Double_t LS_tailIntegral = getIntegral(waveformVector,0,onsetSample+BD_postOnsetIntegralSamples-BD_tailIntegralSamples,BD_tailIntegralSamples);
								LS_psd = LS_tailIntegral/LS_integral;
								
								// Check BPM pulses exist
								if (BPM_timingVector.size() != 0) {
									//Get time to previous BPM pulse
									while (bpmIndex < BPM_timingVector.size()) {
										if (BPM_timingVector[bpmIndex] > LS_waveformStartTime+LS_onsetTime) {
											if (bpmIndex > 0) {
												//LS_timeToBPM = LS_waveformStartTime+LS_onsetTime-BPM_timingVector[bpmIndex-1];
												LS_timeToBPM = BPM_timingVector[bpmIndex]-(LS_waveformStartTime+LS_onsetTime);
												
												//Fill vectors 
												LS_channelVector.push_back(channelID);
												LS_saturatedVector.push_back(LS_saturated);
												LS_waveformStartTimeVector.push_back(LS_waveformStartTime);
												LS_onsetTimeVector.push_back(LS_onsetTime);
												LS_timeToBPMVector.push_back(LS_timeToBPM);
												LS_preTraceIntegralVector.push_back(LS_preTraceIntegral);
												LS_preTraceIntegralDevVector.push_back(waveformBaselineRMS);
												LS_integralVector.push_back(LS_integral);
												LS_peakHeightVector.push_back(LS_peakHeight);
												LS_psdVector.push_back(LS_psd);
												
												//If we don't have a scatterer, then fill tree here.
												if (hasScatterer == 0) {
													analysisTree->Fill();
												}
												
												break;
											}
										}
										bpmIndex++;
									} //end while
								} // end check BPM pulses exist
								else {
								  //If there are no BPM pulses, we still may want to add this pulse to tree
                  if (hasScatterer == 0) {
										analysisTree->Fill();
									}
								}//end check for BPM
							}//end check have enough samples to integrate
						}//end check have enough samples to interpolate 
					}//end loop through found pulses
				}//end check for valid baseline
      }//end check is a LS channel
    }//end loop through entries
  }//end has backing detectors


////////////////////////FINALLY DEAL WITH SCATTERER PULSES////////////////////////////


	//Reset index counters 
	bpmIndex=0;
	Int_t lsIndex=0;
	if (hasScatterer==1) {
	
		//Step through scatterer pulses
		for (Int_t entry=0; entry < numSignals; entry++) {
		
		  if (entry%10000==0) {
		    cout<<"Looking for scatterer signals, on event "<<entry<<" of "<<numSignals<<endl;
		  }
		
			sis3316tree->LoadTree(entry);
			channelIDBranch->GetEntry(entry);
			
			//Check if this is an NaI channel we're processing 
			if (channelID==scattererChannel) {
		
		    nSamplesBranch->GetEntry(entry);
		    timestampBranch->GetEntry(entry);
		    waveformBranch->GetEntry(entry);
		    
				//Convert array of shorts into a vector 
				scatterer_waveform.assign(shortWaveform, shortWaveform + nSamples);
				//Convert vector of shorts into a vector of doubles
				vector<Double_t> waveformVector(scatterer_waveform.begin(),scatterer_waveform.end());
				
				//Corresponds to trigger time 
				scatterer_waveformStartTime = static_cast<Double_t>(timestamp)*nsPerSample;
				//Subtract off the preTriggerDelay to get the start time of the waveform
				scatterer_waveformStartTime -= static_cast<Double_t>(scatterer_preTriggerDelay)*nsPerSample;
			
				//Get baseline guess--fit gaussian near mode of waveform, get mean and sigma
				pair<Double_t,Double_t> baselineContainer = getBaselineGuess(waveformVector,scatterer_baselineVariationGuess);
				Double_t waveformBaseline = baselineContainer.first;
				scatterer_baseline = waveformBaseline;
				Double_t waveformBaselineRMS = baselineContainer.second;
				
				//plotWaveformAndPulse(waveformVector,nsPerSample,0);
				//If we found a valid NaI pulse, set this to one. If not, still record the waveform, but set all other quantities to zero
				LS_pileUpInScatterIntegrationWindow=0;
				scatterer_foundPulse=0;
				scatterer_onsetTime=0;
				scatterer_timeToBPM=0;
				scatterer_timeToBD=0;
				scatterer_preTraceIntegral=0;
				scatterer_integral=0;
				scatterer_peakHeight=0;
				scatterer_meanTime=0;

				//baseline of zero indicates we had problems, most likely due to pile-up. Don't search for pulses in these 
				if (waveformBaseline != 0) {
				
					//Get vector of CMA for baseline subtraction
					vector<Double_t> cmaFilter = getCMAFilter(waveformVector,scatterer_halfWidth,waveformBaseline,scatterer_rejectThresholdSigma*waveformBaselineRMS);
					//Subtract CMA to subtract baseline 
					for (Int_t i=0; i < waveformVector.size(); i++) {
						waveformVector.at(i) = waveformVector.at(i)-cmaFilter.at(i);
					}

								//Locate our fixed integration region.
							        Int_t onsetSample;
								if (BPM_timingVector.size()>0) {
								  while (bpmIndex < BPM_timingVector.size()) {
									  if (BPM_timingVector[bpmIndex] > scatterer_waveformStartTime+bpm_Offset) {
										  if (bpmIndex > 0) {
											  scatterer_onsetTime = BPM_timingVector[bpmIndex]-bpm_Offset;
										    	  scatterer_onsetTime -= scatterer_waveformStartTime;
					    						  onsetSample = round( scatterer_onsetTime / nsPerSample );
											  break;
										  }
									  }
									  bpmIndex++;
								  }
								}
		        	    
					    //Make sure we have enough samples before and after onset sample
					    if ((onsetSample > scatterer_preOnsetSamplesRequired) && (onsetSample + scatterer_postOnsetIntegralSamplesRequired < nSamples)) {
						
								scatterer_foundPulse=1;
						
							  scatterer_preTraceIntegral = getIntegral(waveformVector,0,onsetSample-scatterer_preOnsetSamplesRequired,scatterer_preTraceIntegralSamples);
							    
					      //Get integral and noise.
							  scatterer_integral = getIntegral(waveformVector,0,onsetSample-scatterer_preOnsetIntegralSamples,scatterer_preOnsetIntegralSamples+scatterer_postOnsetIntegralSamples);

							  scatterer_noise = getIntegral(waveformVector,0,onsetSample+scatterer_postOnsetIntegralSamples,scatterer_preOnsetIntegralSamples+scatterer_postOnsetIntegralSamples);
							    
				        //Find PHV within same region we're going to integrate for 
				        scatterer_peakHeight = getPHV(waveformVector,0,onsetSample,scatterer_postOnsetIntegralSamples,0);

//interpolate region where we expect to find onset 
					vector<Double_t> interpWaveVec;
					interpWaveVec.clear();
					interpWaveVec = getTsincVector(waveformVector,0,onsetSample-BD_preOnsetInterpolationSamples,BD_preOnsetInterpolationSamples,nInterPoints,fTaper,l);

					//Find the interpolated onset sample for computing timing.
					Double_t interpOnsetSample = getPHVCFDSample(interpWaveVec,0,0,interpWaveVec.size(),BD_cfdFraction);
					interpOnsetSample /= nInterPoints;
					interpOnsetSample += onsetSample-scatterer_preOnsetInterpolationSamples;
					scatterer_onsetTime = interpOnsetSample*nsPerSample;
					scatterer_timeToBPM = BPM_timingVector[bpmIndex]-(scatterer_waveformStartTime+scatterer_onsetTime);        
							  //Get mean time 
							  scatterer_meanTime = getMeanTime(waveformVector,0,onsetSample-scatterer_preOnsetIntegralSamples,scatterer_preOnsetIntegralSamples+scatterer_postOnsetIntegralSamples,scatterer_peakHeight);
							  
					  	} //end check had enough samples to integrate 
						//} //end check had enough samples to interpolate  
					//} //end check we found at least one pulse 
				} //end check we found a valid baseline
				
				//If LS detectors were active (not a NaI calibration run), check if we have LS-pileup present
				if (LS_onsetTimeVector.size() != 0) {
					while (lsIndex < LS_onsetTimeVector.size()) {
					
						Double_t lsTime=LS_waveformStartTimeVector[lsIndex]+LS_onsetTimeVector[lsIndex];
						
						//Find first LS event after the start of the NaI waveform 
						if (lsTime >= scatterer_waveformStartTime) {
						
							scatterer_waveformStartTimeToNextBD=lsTime-scatterer_waveformStartTime;
						
							//Let's check if the following LS pulse occurs within scatterer_preOnsetIntegralSamples+scatterer_postOnsetIntegralSamples
							//of the current one--i.e. there is pileup. First make sure we won't go out of bounds.
							if (lsIndex+1 < LS_onsetTimeVector.size()) {
								Double_t next_lsTime = LS_waveformStartTimeVector[lsIndex+1]+LS_onsetTimeVector[lsIndex+1];
								if (next_lsTime - lsTime < scatterer_preTriggerDelay+scatterer_preOnsetIntegralSamples+scatterer_postOnsetIntegralSamples) {
									LS_pileUpInScatterIntegrationWindow=1;
								}
							}
							
							//If we found a pulse in scatterer, calculate this time 
							if (scatterer_foundPulse==1) {
								scatterer_timeToBD = (scatterer_waveformStartTime+scatterer_onsetTime)-lsTime;
							}
							
							//Fill TTree
							LS_channel=LS_channelVector[lsIndex];
							LS_waveformStartTime=LS_waveformStartTimeVector[lsIndex];
							LS_onsetTime=LS_onsetTimeVector[lsIndex];
							LS_timeToBPM=LS_timeToBPMVector[lsIndex];
							LS_preTraceIntegral=LS_preTraceIntegralVector[lsIndex];
							LS_preTraceIntegralDev=LS_preTraceIntegralDevVector[lsIndex];
							LS_integral=LS_integralVector[lsIndex];
							LS_peakHeight=LS_peakHeightVector[lsIndex];
							LS_psd=LS_psdVector[lsIndex];
							LS_saturated=LS_saturatedVector[lsIndex];
							analysisTree->Fill();
							
							break;
						}
						
						//Increment counter if we haven't gotten to an LS after the NaI waveform start time 
						lsIndex++;
					}//End while
				} // End check LS pulses exist

				else {
				  //Still fill waveform branch even if we don't find any pulses
					analysisTree->Fill();
				}
			} //end check this was a NaI channel 
		} //end for loop over all entries
	} //end check we had scatterer channel
  
  //Write 
  output->cd();
  if (BPMTree->GetEntries()>0) {
	BPMTree->Write("BPMTree", TObject::kOverwrite);
  }
  if (analysisTree->GetEntries()>0) {
	analysisTree->Write("analysisTree",TObject::kOverwrite);
  }

  //Close files.
  output->Close();
  input->Close();

  //Delete stuff.
  LS_channelVector.clear();
  LS_saturatedVector.clear();
  LS_waveformStartTimeVector.clear();
  LS_onsetTimeVector.clear();
  LS_timeToBPMVector.clear();
  LS_preTraceIntegralVector.clear();
  LS_preTraceIntegralDevVector.clear();
  LS_integralVector.clear();
  LS_peakHeightVector.clear();
  LS_psdVector.clear();

}

//Supply this with a path to a .txt file containing the names of every root file you want to process
//and the path to the folder you want output files written to.
//Only works for files of the same type i.e. a TOF Run or a Production Run.
void processMultipleFiles( TString inputFilename, TString pathToOutputFiles ){

  //Step through the .txt file, processing as we go.
  ifstream file( inputFilename );
  string line;
  TString lineTString;
  Int_t lineNum = 0;
  while( file.good() ){
	if (getline ( file, line, '\n' )){
		lineTString = line;
		cout << "Processing file " << line << endl;
		TString delimiter1 = "SIS3316Raw";
		TString delimiter2 = ".root";
		TString outputName = "Processed_";
		outputName += line.substr(line.find(delimiter1),line.find(delimiter2));
		TString outputFilename = pathToOutputFiles + outputName;
		//cout <<"The processed file will be " << outputFilename << endl;
		processData( lineTString, outputFilename);
	}
  }
  cout << "Routine complete!" << endl;
}








