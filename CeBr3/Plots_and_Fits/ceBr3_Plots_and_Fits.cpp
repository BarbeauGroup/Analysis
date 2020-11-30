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
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooGamma.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooAddition.h"
#include "RooMinimizer.h"
#include "RooNumConvPdf.h"
#include "RooFitResult.h"

using namespace std;
using namespace RooFit;

//Supply this with the path to a .txt file containing the full paths
//to all your processed files and a flag if the scatterer is present.
void fitPEs( TString filename ){

  //Silence function evaluation errors - weird quirk with RooGamma.
  gErrorIgnoreLevel = kFatal;
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  //Set the fit range.
  Long64_t fitMin = 0;
  Long64_t fitMax = 10e3;

  //Declare our observable.
  TCanvas* fitCanvas = new TCanvas("fitCanvas","fitCanvas");
  RooRealVar integral("integral","integral",0,10e3);
  RooArgSet observables("observables");
  observables.add(integral);
  RooDataSet noise("noise","Noise Pulses",integral);

  // Set #bins to be used for FFT sampling to 10000
  integral.setBins(10000,"cache"); 

  //Variables to assign to ttree branches.
  Double_t scatterer_noise;
  //Create a TChain to hold all of our data.
  TChain* chain = new TChain("analysisTree");
  chain->SetBranchAddress("scatterer_noise",&scatterer_noise);
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
  //Fill our data set.
  Long64_t numEntries = chain->GetEntries();
  for( Int_t i = 0; i < numEntries; i++ ){
	chain->GetEntry( i );
	integral = scatterer_noise;
	noise.add( integral );

  }

  //Declare our RooFit Variables for the Gaussian Pedestal.
  RooRealVar pedMean("pedMean","pedMean",0,-1,20);
  RooRealVar pedSigma("pedSigma","pedSigma",10,0,15);
  RooFormulaVar pedMean2("pedMean2","2*pedMean",RooArgList(pedMean));
  RooFormulaVar pedMean3("pedMean3","3*pedMean",RooArgList(pedMean));


  //Declare our RooFit Variables for the Polya (RooGamma).
  //mu is set to be equal to the pedestal mean, and beta2 and gamma 2 ( beta3 and gamma 3, etc.)
  //are set such that the mean of the 2nd Polya is an appropriate multiple of that of the first Polya.
  RooRealVar gamma1("gamma1","gamma1",9,2,15);
  RooRealVar beta1("beta1","beta1",8,2,25);
  RooFormulaVar gamma2("gamma2","2*gamma1",RooArgList(gamma1));
  RooFormulaVar gamma3("gamma3","3*gamma1",RooArgList(gamma1));

  //Declare our RooFit Variables for Exponential Noise.
  RooRealVar noiseDecay("noiseDecay","noiseDecay",200,100,300);

  //Declare PDFs.
  RooGaussian pedestal("pedestal","Gaussian Pedestal PDF",integral,pedMean,pedSigma);
  RooGaussModel pedestalResponse("pedestalResponse","pedestalResponse",integral,pedMean,pedSigma);
  RooDecay decay("decay","Exponential Noise Decay PDF",integral,noiseDecay,pedestalResponse,RooDecay::SingleSided);
  RooGamma spe("spe","Polya PDF for SPEs",integral,gamma1,beta1,pedMean);
  RooGamma dpe("dpe","Polya PDF for DPEs",integral,gamma2,beta1,pedMean2);
  RooGamma tpe("tpe","Polya PDF for TPEs",integral,gamma3,beta1,pedMean3);
  //RooGamma qpe("qpe","Polya PDF for QPEs",integral,gamma4,beta4,pedMean);
  //RooGamma fpe("fpe","Polya PDF for FPEs",integral,gamma5,beta5,pedMean);
  //("spe","Polya PDF for SPEs",integral,gamma1,beta1,pedMean);
  //Add PDFs.
  RooRealVar pedEvents("pedEvents","pedEvents",2200,0,1000000);
  RooRealVar decayEvents("decayEvents","decayEvents",300,0,100000);
  RooRealVar peEvents("peEvents","peEvents",25000,0,1000000);
  RooRealVar peMu("peMu","peMu",10,1e-3,1000);
  RooRealVar speEvents("speEvents","speEvents",10000,0,2000000);
  RooRealVar dpeEvents("dpeEvents","dpeEvents",1000,0,1000000);
  RooRealVar tpeEvents("tpeEvents","tpeEvents",25000,0,500000);
  //RooRealVar qpeEvents("qpeEvents","qpeEvents",1000,0,100000);
  //RooRealVar fpeEvents("fpeEvents","fpeEvents",100,0,10000);
  RooAddPdf speModel("speModel","Polya + Pedestal",RooArgList(pedestal,decay,spe),RooArgList(pedEvents,decayEvents,speEvents));

  //Fit the data.
  speModel.fitTo(noise);
  //Plot.
  fitCanvas->cd();
  RooPlot* frame = integral.frame();
  frame->SetTitle("SPE Calibration");
  noise.plotOn(frame,MarkerStyle(20));
  speModel.plotOn(frame,Precision(1e-6));
  frame->Draw();

}

//Supply this with the path to a .txt file containing the full paths
//to all your processed files and a flag if the scatterer is present.
void makePlots( TString filename, Int_t hasScatterer ){

  //Flag to plot average waveforms.
  Int_t plotAvgWaves = 0;

  //Calibration constant.
  Double_t adc_to_keV = 9.65e-4;

  //Cuts
  Double_t time_Low = 306;
  Double_t time_High = 336;
  Double_t psd_Low = 0.26;
  Double_t psd_High = 0.6;
  Double_t integral_Low = 10000;
  Double_t integral_High = 35000;
  Int_t BD_Choice =11; //LS Channel you want to look at.

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
  Double_t scatterer_noise;
  Double_t scatterer_peakHeight;
  Double_t scatterer_baseline; 
  Double_t scatterer_meanTime;
  vector<UShort_t>* scatterer_waveform = new vector<UShort_t>;

  //Create a TChain to hold all of our data.
  TChain* chain = new TChain("analysisTree");
  chain->SetBranchAddress("LS_channel", &LS_channel);
  chain->SetBranchAddress("LS_psd",&LS_psd);
  chain->SetBranchAddress("LS_integral",&LS_integral);
  chain->SetBranchAddress("LS_timeToBPM",&LS_timeToBPM);
  if( hasScatterer == 1 ){
	chain->SetBranchAddress("scatterer_foundPulse",&scatterer_foundPulse);
	chain->SetBranchAddress("scatterer_integral",&scatterer_integral);
	chain->SetBranchAddress("scatterer_noise",&scatterer_noise);
	chain->SetBranchAddress("scatterer_timeToBPM",&scatterer_timeToBPM);
	chain->SetBranchAddress("scatterer_timeToBD",&scatterer_timeToBD);
	chain->SetBranchAddress("scatterer_onsetTime",&scatterer_onsetTime);
	chain->SetBranchAddress("scatterer_waveform",&scatterer_waveform);
	chain->SetBranchAddress("scatterer_peakHeight",&scatterer_peakHeight);
	chain->SetBranchAddress("scatterer_baseline",&scatterer_baseline);
  }
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
  TH1D* LS_timing_Hist = new TH1D("LS_timing_Hist","Backing Detector Time to BPM",400,0,400);
  TH2D* LS_2D_Hist = new TH2D("LS_2D_Hist","Backing Detector Integral and PSD",20000,0,200000,100,0,1);
  TH2D* LS_psd_timing_Hist = new TH2D("LS_psd_timing_Hist","Backing Detector PSD and Time to BPM",400,0,400,100,0,1);
  TH2D* LS_integral_timing_Hist = new TH2D("LS_integral_timing_Hist","Backing Detector integral and Time to BPM",400,0,400,20000,0,200000);
  TH1D* scatterer_energy_Hist = new TH1D("scatterer_energy_Hist","Scatterer Energy",50,0,5000*adc_to_keV);
  TH1D* scatterer_timeToBPM_Hist = new TH1D("scatterer_timeToBPM_Hist","Scatterer Time to BPM",400,0,400);
  TH1D* scatterer_timeToBD_Hist = new TH1D("scatterer_timeToBD_Hist","Scatterer Time to Backing Detector",1100,0,1100);
  UShort_t neutronArr[500];
  TH1D* avgNeutronHist = new TH1D("avgNeutronHist","Average Neutron Waveform",500,0,500);
  TH2D* scatterer_neutronHist = new TH2D("scatterer_neutronHist","scatterer_neutronHist",5000,0,500,1000,0,1);
  TH1D* avgGammaHist = new TH1D("avgGammaHist","Average Gamma Waveform",500,0,500);
  TH2D* scatterer_gammaHist = new TH2D("scatterer_gammaHist","scatterer_wfHist",5000,0,500,1000,0,1);
  TH1D* scatterer_noise_Hist = new TH1D("scatterer_noise_Hist","Scatterer Noise",50,0,5000*adc_to_keV);

  //Fill the histograms.
  Int_t numEntries = chain->GetEntries();
  Int_t numNeutrons = 0;
  Int_t numGammas = 0;
  UShort_t wfSample;
  Double_t peak;
  for( Int_t i = 0; i < numEntries; i++ ){
	chain->GetEntry(i);
	if ( plotAvgWaves == 1 && LS_channel == BD_Choice ){
		if ( hasScatterer == 1 && ( LS_psd > psd_Low ) && ( LS_psd < psd_High ) ){
			if( scatterer_integral > integral_Low && scatterer_integral < integral_High ){
				peak = scatterer_peakHeight;
				for( Int_t j = 0; j < 500; j++ ){
					wfSample = scatterer_waveform->at(j);
					//wfSample -= scatterer_baseline;
					avgNeutronHist->Fill(j,wfSample / ( peak + scatterer_baseline ));
					scatterer_neutronHist->Fill(j,wfSample / ( peak + scatterer_baseline ));
					numNeutrons += 1;
				}
			}
		}
		else if ( hasScatterer == 1 && ( LS_psd < psd_Low ) ){
			if( scatterer_integral > integral_Low && scatterer_integral < integral_High ){
				peak = scatterer_peakHeight;
				for( Int_t j = 0; j < 500; j++ ){
					wfSample = scatterer_waveform->at(j);
					//wfSample -= scatterer_baseline;
					avgGammaHist->Fill(j,wfSample / ( peak + scatterer_baseline ));
					scatterer_gammaHist->Fill(j,wfSample / ( peak + scatterer_baseline ));
					numGammas += 1;
				}
			}
		}
	}
	if(hasScatterer == 1 && LS_channel == BD_Choice){
		LS_psd_timing_Hist->Fill(LS_timeToBPM,LS_psd);
		LS_integral_timing_Hist->Fill(LS_timeToBPM,LS_integral);
		LS_2D_Hist->Fill(LS_integral,LS_psd);
		if( ( LS_psd > psd_Low ) && ( LS_psd < psd_High ) ){
			if( ( LS_integral > integral_Low ) && ( LS_integral < integral_High ) ){
				LS_timing_Hist->Fill(LS_timeToBPM);
				if( ( LS_timeToBPM > time_Low ) && ( LS_timeToBPM < time_High ) ){
					LS_psd_Hist->Fill(LS_psd);
					LS_integral_Hist->Fill(LS_integral);
					scatterer_energy_Hist->Fill(scatterer_integral*adc_to_keV);
					scatterer_noise_Hist->Fill(scatterer_noise*adc_to_keV);
					scatterer_timeToBPM_Hist->Fill(scatterer_timeToBPM);
					scatterer_timeToBD_Hist->Fill(scatterer_timeToBD);
				}
			}
		}
	}
	else if( hasScatterer == 0 ){
		LS_psd_timing_Hist->Fill(LS_timeToBPM,LS_psd);
		LS_integral_timing_Hist->Fill(LS_timeToBPM,LS_integral);
		LS_2D_Hist->Fill(LS_integral,LS_psd);
		if( ( LS_psd > psd_Low ) && ( LS_psd < psd_High ) ){
			if( ( LS_integral > integral_Low ) && ( LS_integral < integral_High ) ){
				LS_timing_Hist->Fill(LS_timeToBPM);
				if( ( LS_timeToBPM > time_Low ) && ( LS_timeToBPM < time_High ) ){
					LS_psd_Hist->Fill(LS_psd);
					LS_integral_Hist->Fill(LS_integral);
				}
			}
		}
	}
  }

  //Normalize our avergae waveforms.
  Double_t norm = 1.;
  Double_t neutronScale = norm/avgNeutronHist->Integral();
  avgNeutronHist->Scale(neutronScale);
  Double_t gammaScale = norm/avgGammaHist->Integral();
  avgGammaHist->Scale(gammaScale);	

  //Prepare TCanvases and make plots.
  //TCanvas* LS_psd_Canvas = new TCanvas("LS_psd_Canvas","Liquid Scintillator PSD");
  //TCanvas* LS_integral_Canvas = new TCanvas("LS_integral_Canvas","Backing Detector Integral");
  //TCanvas* LS_timing_Canvas = new TCanvas("LS_timing_Canvas","Backing Detector Time to BPM");
  //TCanvas* LS_2D_Canvas = new TCanvas("LS_2D_Canvas","Backing Detector Integral and PSD");
  //TCanvas* LS_psd_bpm_Canvas = new TCanvas("LS_psd_bpm_Canvas","Backing Detector PSD and Time to BPM");
  //TCanvas* LS_integral_bpm_Canvas = new TCanvas("LS_integral_bpm_Canvas","Backing Detector Integral and Time to BPM");
  if( hasScatterer == 1 ){
	TCanvas* scatterer_energy_Canvas = new TCanvas("scatterer_energy_Canvas","Scatterer Pulse Energy");
	TCanvas* scatterer_timeToBPM_Canvas = new TCanvas("scatterer_timeToBPM_Canvas","Scatterer Time to Prev. BPM");
	TCanvas* scatterer_timeToBD_Canvas = new TCanvas("scatterer_timeToBD_Canvas","Scatterer Time to Prev. BD");
	if ( plotAvgWaves == 1 ){
		TCanvas* avgNeutronCanvas = new TCanvas("avgNeutronCanvas","Neutron Waveforms");
		avgNeutronCanvas->Divide(2,1);
		avgNeutronCanvas->cd(1);
		avgNeutronHist->SetMarkerColor(2);
		avgNeutronHist->Draw("C");
		avgNeutronCanvas->cd(2);
		scatterer_neutronHist->Draw("COLZ");
		TCanvas* avgGammaCanvas = new TCanvas("avgGammaCanvas","Gamma Waveforms");
		avgGammaCanvas->Divide(2,1);
		avgGammaCanvas->cd(1);
		avgGammaHist->Draw("C");
		//avgNeutronHist->DrawCopy("SAME");
		avgGammaCanvas->cd(2);
		scatterer_gammaHist->Draw("COLZ");
	}

	scatterer_energy_Canvas->cd();
	//scatterer_energy_Hist->Draw();
	scatterer_noise_Hist->SetLineColor(2);
	//scatterer_noise_Hist->Draw("SAME");
	auto rp = new TRatioPlot(scatterer_energy_Hist, scatterer_noise_Hist,"diff");
	//rp->GetLowerRefYaxis()->SetTitle("Residuals");
	//rp->GetUpperRefYaxis()->SetTitle("Counts");
	rp->Draw();
	scatterer_timeToBPM_Canvas->cd();
	scatterer_timeToBPM_Hist->Draw();
	scatterer_timeToBD_Canvas->cd();
	scatterer_timeToBD_Hist->Draw();
  }
  //LS_psd_Canvas->cd();
  //LS_psd_Hist->Draw();
  //LS_integral_Canvas->cd();
  //LS_integral_Hist->Draw();
  //LS_timing_Canvas->cd();
  //LS_timing_Hist->Draw();
  //LS_2D_Canvas->cd();
  //LS_2D_Hist->Draw("COLZ");
  //Make a box showing our PSD cut if we're looking at a production run.
  //if( hasScatterer == 1 ){
	//TLine* line1 = new TLine(integral_Low,psd_Low,integral_High,psd_Low);
	//TLine* line2 = new TLine(integral_High,psd_Low,integral_High,psd_High);
	//TLine* line3 = new TLine(integral_Low,psd_High,integral_High,psd_High);
	//TLine* line4 = new TLine(integral_Low,0.26,integral_Low,psd_High);
	//line1->SetLineWidth(3);
	//line2->SetLineWidth(3);
	//line3->SetLineWidth(3);
	//line4->SetLineWidth(3);
	//line1->SetLineColor(kRed);
	//line2->SetLineColor(kRed);
	//line3->SetLineColor(kRed);
	//line4->SetLineColor(kRed);
	//line1->Draw("SAME");
	//line2->Draw("SAME");
	//line3->Draw("SAME");
	//line4->Draw("SAME");
  //}
  //LS_psd_bpm_Canvas->cd();
  /*LS_psd_timing_Hist->Draw("COLZ");
  if( hasScatterer == 1 ){
	TLine* line5 = new TLine(time_Low,psd_Low,time_High,psd_Low);
	TLine* line6 = new TLine(time_High,psd_Low,time_High,psd_High);
	TLine* line7 = new TLine(time_Low,psd_High,time_High,psd_High);
	TLine* line8 = new TLine(time_Low,psd_Low,time_Low,psd_High);
	line5->SetLineWidth(3);
	line6->SetLineWidth(3);
	line7->SetLineWidth(3);
	line8->SetLineWidth(3);
	line5->SetLineColor(kRed);
	line6->SetLineColor(kRed);
	line7->SetLineColor(kRed);
	line8->SetLineColor(kRed);
	line5->Draw("SAME");
	line6->Draw("SAME");
	line7->Draw("SAME");
	line8->Draw("SAME");
  }
  LS_integral_bpm_Canvas->cd();
  LS_integral_timing_Hist->Draw("COLZ");
  if( hasScatterer == 1 ){
	TLine* line9 = new TLine(time_Low,integral_Low,time_High,integral_Low);
	TLine* line10 = new TLine(time_High,integral_Low,time_High,integral_High);
	TLine* line11 = new TLine(time_Low,integral_High,time_High,integral_High);
	TLine* line12 = new TLine(time_Low,integral_Low,time_Low,integral_High);
	line9->SetLineWidth(3);
	line10->SetLineWidth(3);
	line11->SetLineWidth(3);
	line12->SetLineWidth(3);
	line9->SetLineColor(kRed);
	line10->SetLineColor(kRed);
	line11->SetLineColor(kRed);
	line12->SetLineColor(kRed);
	line9->Draw("SAME");
	line10->Draw("SAME");
	line11->Draw("SAME");
	line12->Draw("SAME");
  }
*/
}

//Fits the peaks in the calibration spectra.
void calibrate( ){

  //Silence function evaluation errors - weird quirk with RooGamma.
  gErrorIgnoreLevel = kFatal;
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  //Set the fit range.
  Long64_t fitMin = 140e3;
  Long64_t fitMax = 155e3;

  //Declare our observables.
  RooRealVar integral("integral","integral",0,185e3);
  integral.setBins(fitMax - fitMin);
  RooArgSet observables("observables");
  observables.add(integral);

  //Set up data containers.
  RooDataSet data("data","Ba-133 Calibration Data",observables);
  RooDataSet bgnd("bgnd","Room Backgrounds",observables);

  //Read in background data for a PDF.
  Double_t scatterer_integral;
  TString bgndFilename = "/var/phy/project/phil/cma46/CeBr3/processedRuns_FixedIntegral/CalibrationRun_1/Processed_SIS3316Raw_20200207142051_1.root";
  TFile bgndFile( bgndFilename );
  TTree* bgndTree = (TTree*)bgndFile.Get("analysisTree");
  bgndTree->SetBranchAddress("scatterer_integral",&scatterer_integral);
  Long64_t numBGNDEntries = bgndTree->GetEntries();
  for( Int_t i = 0; i < numBGNDEntries; i++ ){
	bgndTree->GetEntry( i );
	if( scatterer_integral > fitMin && scatterer_integral < fitMax ){
		integral = scatterer_integral;
		bgnd.add( integral );
	}
  }

  //Read in calibration data.
  TString filename = "/var/phy/project/phil/cma46/CeBr3/processedRuns_FixedIntegral/CalibrationRun_1/Processed_SIS3316Raw_20200207141804_1.root";
  TFile dataFile( filename );
  TTree* dataTree = (TTree*)dataFile.Get("analysisTree");
  dataTree->SetBranchAddress("scatterer_integral",&scatterer_integral);
  Long64_t numDataEntries = dataTree->GetEntries();
  for( Int_t i = 0; i < numDataEntries; i++ ){
	dataTree->GetEntry( i );
	if( scatterer_integral > fitMin && scatterer_integral < fitMax ){
		integral = scatterer_integral;
		data.add( integral );
	}
  }

  //Data health check.
  data.Print();
  bgnd.Print();

  //Set up binning.
  ((RooRealVar*)data.get()->find("integral"))->setBins(fitMax - fitMin);
  ((RooRealVar*)bgnd.get()->find("integral"))->setBins(fitMax - fitMin);

  //Create PDFs.
  RooRealVar mean1("mean1","mean1",15185,14500,15400);
  RooRealVar sigma1("sigma1","sigma1",1e3,0,5e3);
  RooRealVar mean2("mean2","mean2",30681,30500,30900);
  RooRealVar sigma2("sigma2","sigma2",1e3,0,5e3);
  RooRealVar mean3("mean3","mean3",41012,40000,41200);
  RooRealVar sigma3("sigma3","sigma3",1e3,0,5e3);
  RooRealVar mean4("mean4","mean4",56508,56300,56700);
  RooRealVar sigma4("sigma4","sigma4",1e3,0,5e3);
  RooRealVar mean5("mean5","mean5",82334,82100,82500);
  RooRealVar sigma5("sigma5","sigma5",1e3,0,5e3);
  RooRealVar mean6("mean6","mean6",148400,148200,148600);
  RooRealVar sigma6("sigma6","sigma6",1e3,0,5e3);
  RooRealVar mean7("mean7","mean7",161364,161100,161500);
  RooRealVar sigma7("sigma7","sigma7",1e3,0,5e3);
  RooRealVar mean8("mean8","mean8",179442,179200,179600);
  RooRealVar sigma8("sigma8","sigma8",1e3,0,5e3);
  RooRealVar mean9("mean9","mean9",190289,190100,190500);
  RooRealVar sigma9("sigma9","sigma9",1e3,0,5e3);
  RooRealVar mean10("mean10","mean10",203202,203000,203500);
  RooRealVar sigma10("sigma10","sigma10",1e3,0,5e3);
  RooGaussian peak1("peak1","First Peak",integral,mean1,sigma1);
  RooGaussian peak2("peak2","Second Peak",integral,mean2,sigma2);
  RooGaussian peak3("peak3","Third Peak",integral,mean3,sigma3);
  RooGaussian peak4("peak4","Fourth Peak",integral,mean4,sigma4);
  RooGaussian peak5("peak5","Fifth Peak",integral,mean5,sigma5);
  RooGaussian peak6("peak6","Sixth Peak",integral,mean6,sigma5);
  RooGaussian peak7("peak7","Seventh Peak",integral,mean7,sigma7);
  RooGaussian peak8("peak8","Eigth Peak",integral,mean8,sigma8);
  RooGaussian peak9("peak9","Ninth Peak",integral,mean9,sigma9);
  RooGaussian peak10("peak10","Tenth Peak",integral,mean10,sigma10);
  RooDataHist* bgndHist = bgnd.binnedClone();
  RooHistPdf bgndPDF("bgndPDF","Room Backgrounds",observables,*bgndHist,0);
  RooRealVar amp1("amp1","Amplitude of the first peak",5e7,0,100e6);
  RooRealVar amp2("amp2","Amplitude of the second peak",1e7,0,100e6);
  RooRealVar amp3("amp3","Amplitude of the third peak",3e7,0,100e6);
  RooRealVar amp4("amp4","Amplitude of the fourth peak",0,100e6);
  RooRealVar amp5("amp5","Amplitude of the fifth peak",0,100e6);
  RooRealVar amp6("amp6","Amplitude of the sixth peak",0,100e6);
  RooRealVar amp7("amp7","Amplitude of the seventh peak",0,100e6);
  RooRealVar amp8("amp8","Amplitude of the eigth peak",0,100e6);
  RooRealVar amp9("amp9","Amplitude of the ninth peak",0,100e6);
  RooRealVar amp10("amp10","Amplitude of the tenth peak",0,100e6);
  RooRealVar bgndAmp("bgndAmp","Amplitude of the background pdf",5e7,0,100e8);
  RooRealVar a("a","a",221600.57,1e3,1e5);
  RooRealVar b("b","b",-0.00014,-10,10);
  RooPolynomial simpleBGND("simpleBGND","Simple linear background model",integral,RooArgSet(a,b), 0);
  RooAddPdf model("model","Peaks + BGND",RooArgList(peak7,peak8,simpleBGND),RooArgList(amp7,amp8,bgndAmp));

  //Construct a binned version of our data.
  RooAbsData* binnedData = data.binnedClone();

  //Do the fit.
  model.fitTo(data,Save(1),Range(fitMin,fitMax));

  //Create a canvas.
  TCanvas* c1 = new TCanvas("c1","c1");
    
  //Produce Plots and Residuals.
  c1->cd(1);
  c1->SetLogy();
  RooPlot* frame = integral.frame(Title("Ba-133 Calibration"));
  binnedData->plotOn(frame,Name("Data"),MarkerColor(1),FillColor(0),Binning(185));
  model.plotOn(frame,Name("Total Model"),LineColor(9),FillColor(0),Binning(fitMax - fitMin));
  //RooHist* resHist = frame->residHist();
  //resHist->SetMarkerColorAlpha(30,1);
  model.plotOn(frame,Name("Backgrounds"),Components(bgndPDF),LineColor(2),LineStyle(kDashed),FillColor(0));
  //recModel.plotOn(frame,Name("Poisson Residuals"),Components(recPoisson),LineColor(3),LineStyle(kDashed),FillColor(0),Binning(fitMax - fitMin));
  frame->GetXaxis()->CenterTitle();
  frame->GetXaxis()->SetRangeUser(fitMin,fitMax);
  frame->GetXaxis()->SetTitle("Integral [Arb.]");
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitle("Counts / Bin");
  frame->GetYaxis()->SetRangeUser(0.1,10000000);
  //frame->addPlotable(resHist,"p");
  frame->Draw();

  c1->Modified();
  c1->Update();

  c1->WaitPrimitive();

  //Clean up.
  dataFile.Close();
  bgndFile.Close();

}

