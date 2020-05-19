/* Plotting and Fitting tools package for CeBr3 
Written by C. Awe
5/19/2020

This toolset is loaded directly into root rather than compiling.
While slower, this style is perhaps easier to use across differing
systems. The code here runs on data that has already gone through initial 
processing with ceBr3_AnalysisTools.cpp.

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

//Supply this with the path to a .txt file containing the full paths
//to all your processed files.
void makePlots( TString filename ){

  //Variables to assign to ttree branches.
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

  //Create a TChain to hold all of our data.
  TChain* chain = new TChain("analysisTree");
  chain->SetBranchAddress("scatterer_foundPulse",&scatterer_foundPulse);
  chain->SetBranchAddress("LS_psd",&LS_psd);
  chain->SetBranchAddress("LS_integral",&LS_integral);
  chain->SetBranchAddress("scatterer_integral",&scatterer_integral);
  chain->SetBranchAddress("scatterer_timeToBPM",&scatterer_timeToBPM);
  chain->SetBranchAddress("scatterer_timeToBD",&scatterer_timeToBD);
  chain->SetBranchAddress("scatterer_onsetTime",&scatterer_onsetTime);
  ifstream file( filename );
  string line;
  TString lineTString;
  Int_t lineNum = 0;
  while( file.good() ){
	if (getline ( file, line, '\n' )){
		lineTString = line;
		chain->Add( lineTString );
	}
  }

  //Prepare histograms in order to plot our data.
  TH1D* LS_psd_Hist = new TH1D("LS_psd_Hist","Backing Detector PSD",100,0,1);
  TH1D* LS_integral_Hist = new TH1D("LS_integral_Hist","Backing Detector Integral",20000,0,200000);
  TH2D* LS_2D_Hist = new TH2D("LS_2D_Hist","Backing Detector Integral and PSD",20000,0,200000,100,0,1);
  TH1D* scatterer_integral_Hist = new TH1D("scatterer_integral_Hist","Scatterer Integral",1000,0,10000);
  TH1D* scatterer_timeToBPM_Hist = new TH1D("scatterer_timeToBPM_Hist","Scatterer Time to BPM",500,0,500);
  TH1D* scatterer_timeToBD_Hist = new TH1D("scatterer_timeToBD_Hist","Scatterer Time to Backing Detector",1100,0,1100);

  //Fill the histograms.
  Int_t numEntries = chain->GetEntries();
  for( Int_t i = 0; i < numEntries; i++ ){
	chain->GetEntry(i);
	if( scatterer_foundPulse == 1 ){
		if( ( LS_psd > 0.258 ) && ( LS_psd < 0.6 ) ){
			if( ( LS_integral > 6000 ) && ( LS_integral < 40000 ) ){
				LS_psd_Hist->Fill(LS_psd);
				LS_integral_Hist->Fill(LS_integral);
				LS_2D_Hist->Fill(LS_integral,LS_psd);
				scatterer_integral_Hist->Fill(scatterer_integral);
				scatterer_timeToBPM_Hist->Fill(scatterer_timeToBPM);
				scatterer_timeToBD_Hist->Fill(scatterer_timeToBD);
			}
		}
	}
  }

  //Prepare TCanvases and make plots.
  TCanvas* LS_psd_Canvas = new TCanvas("LS_psd_Canvas","Liquid Scintillator PSD");
  TCanvas* LS_integral_Canvas = new TCanvas("LS_integral_Canvas","Backing Detector Integral");
  TCanvas* LS_2D_Canvas = new TCanvas("LS_2D_Canvas","Backing Detector Integral and PSD");
  TCanvas* scatterer_integral_Canvas = new TCanvas("scatterer_integral_Canvas","Scatterer Pulse Integral");
  TCanvas* scatterer_timeToBPM_Canvas = new TCanvas("scatterer_timeToBPM_Canvas","Scatterer Time to Prev. BPM");
  TCanvas* scatterer_timeToBD_Canvas = new TCanvas("scatterer_timeToBD_Canvas","Scatterer Time to Prev. BD");
  LS_psd_Canvas->cd();
  LS_psd_Hist->Draw();
  LS_integral_Canvas->cd();
  LS_integral_Hist->Draw();
  LS_2D_Canvas->cd();
  //Make a box showing our PSD cut.
  TLine* line1 = new TLine(6000,0.258,40000,0.258);
  TLine* line2 = new TLine(40000,0.258,40000,0.6);
  TLine* line3 = new TLine(6000,0.6,40000,0.6);
  TLine* line4 = new TLine(6000,0.258,6000,0.6);
  LS_2D_Hist->Draw("COLZ");
  line1->Draw("SAME");
  line2->Draw("SAME");
  line3->Draw("SAME");
  line4->Draw("SAME");
  scatterer_integral_Canvas->cd();
  scatterer_integral_Hist->Draw();
  scatterer_timeToBPM_Canvas->cd();
  scatterer_timeToBPM_Hist->Draw();
  scatterer_timeToBD_Canvas->cd();
  scatterer_timeToBD_Hist->Draw();

}




