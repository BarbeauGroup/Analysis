//Simplified code to fit CeBr3 Calibration Data
//Written by C. Awe - 06/07/2021

#include "TString.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include <iostream>
#include <fstream> 
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
#include "TLegend.h"
#include <iostream>
#include <limits>

using namespace std;
using namespace RooFit;

void fitNeutrons( Int_t bdNum ){

  //////////////////////////
  // Setting Up Variables //
  //////////////////////////

  //Path to our data.
  TString filename = "/var/phy/project/phil/cma46/CeBr3/processedRuns_FixedIntegral/ProductionRun_1/processedFilelist.txt";

  //Don't plot canvas on screen.
  gStyle->SetCanvasPreferGL(kTRUE);
  gROOT->SetBatch(1);

  //Fitter options.
  ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(20000); 
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(2000);

  //Energy range to fit.
  Long64_t fitMin=0;
  Long64_t fitMax=10;
  Long64_t keVMin = 0;
  Long64_t keVMax = 100;
  Double_t binning = 500;

  //Conversion from BD Number to LS Channel.
  Int_t ls_channelList[] = { 8, 9, 10, 11, 12, 13, 14, 15 };

  //Declare our observable.
  RooRealVar keV("keV","keV",keVMin,keVMax);
  RooArgSet observables("observables");
  observables.add(keV);

  //Unit conversion.
  Double_t adc_to_keV = 2.00114445451e-3;

  //Data Cuts.
  Double_t time_Low = 325;
  Double_t time_High = 355;
  Double_t psd_Low = 0.26;
  Double_t psd_High = 0.6;
  Double_t integral_Low = 10000;
  Double_t integral_High = 35000;

  /////////////////////////
  // Loading Source Data //
  /////////////////////////

  //Create empty data sets.
  RooDataSet sourceData("sourceData","sourceData",keV);
  RooDataSet noise("noise","noise",keV);	

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

  //Create a TChain to hold all of our data.
  TChain* chain = new TChain("analysisTree");
  chain->SetBranchAddress("LS_channel", &LS_channel);
  chain->SetBranchAddress("LS_psd",&LS_psd);
  chain->SetBranchAddress("LS_integral",&LS_integral);
  chain->SetBranchAddress("LS_timeToBPM",&LS_timeToBPM);
  chain->SetBranchAddress("scatterer_foundPulse",&scatterer_foundPulse);
  chain->SetBranchAddress("scatterer_integral",&scatterer_integral);
  chain->SetBranchAddress("scatterer_noise",&scatterer_noise);
  chain->SetBranchAddress("scatterer_timeToBPM",&scatterer_timeToBPM);
  chain->SetBranchAddress("scatterer_timeToBD",&scatterer_timeToBD);
  chain->SetBranchAddress("scatterer_onsetTime",&scatterer_onsetTime);
  chain->SetBranchAddress("scatterer_peakHeight",&scatterer_peakHeight);
  chain->SetBranchAddress("scatterer_baseline",&scatterer_baseline);
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

  //Now step through our TChain and pick out the right data.
  Int_t numSourceEntries = chain->GetEntries();
  Int_t bd_channel = ls_channelList[bdNum-1]; //Converts BD_Num to a channel number.
  for (Long64_t entry=0; entry < numSourceEntries; entry++) {
	chain->GetEntry(entry);
	if( ( LS_psd > psd_Low ) && ( LS_psd < psd_High ) ){
		if( ( LS_integral > integral_Low ) && ( LS_integral < integral_High ) ){
			if( ( LS_timeToBPM > time_Low ) && ( LS_timeToBPM < time_High ) ){
					Double_t enrg = scatterer_noise * adc_to_keV;
					keV = enrg;
					noise.add( keV );
					if(LS_channel == bd_channel){
						enrg = scatterer_integral * adc_to_keV;
						keV = enrg;
						sourceData.add( keV );
					
					}
			}
		}
	}
		
	//Give user information on screen
	if (entry % 10000 == 0) {
		printf("Processing entry %lu of %lu entries\n", entry, numSourceEntries);
	}
  }

  cout << "Filled!" << endl;

  cout << "Printing source data info." << endl;
  sourceData.Print("V");

  cout << "Generating RooPDFs..." << endl;

  /////////////////
  // Making PDFs //
  /////////////////

  //Build our model for the fit - noise pedestal + gaussian.
  //RooRealVar tau("tau","tau",3,1.5,5);
  //RooTruthModel idealRes("idealRes","Ideal Resolution Model",keV);
  //RooDecay decay("decay","decay",keV,tau,idealRes,RooDecay::SingleSided);
  //Make a noise hist pdf.
  ((RooRealVar*)noise.get()->find("keV"))->setBins(binning);
  RooDataHist* noiseHist = noise.binnedClone();

  RooHistPdf noiseHistPdf("noiseHistPdf","Noise Histogram",keV,*noiseHist,0);
  RooRealVar nNoise("nNoise","nNoise",10000,5000,100000);
  //RooRealVar mean("mean","mean",3.5,0.1,6);
  //RooRealVar sigma("sigma","sigma",0.2,0.01,1); //sigma probably won't ever be greater than fitMax.
  //RooGaussian gauss("gauss","gauss",keV,mean,sigma);
  RooRealVar mu("mu","mu",0); //Offset built into gamma distribution.
  mu.setConstant();
  RooRealVar beta("beta","beta",0.5,0.03,1); //Decay constant built into gamma pdf, roughly corresponds to the spread.
  RooRealVar gamma("gamma","gamma",5,1.0,20); //Another gamma pdf constant, harder to describe. gamma * beta gives the mean of the distribution.
  RooGamma signal("signal","Polya PDF for Nuclear Recoils",keV,gamma,beta,mu);
  RooRealVar nSignal("nSignal","nSignal",3000,100,20000); //Pdf amplitude.
  RooAddPdf model("model","model",RooArgSet(noiseHistPdf,signal),RooArgSet(nNoise,nSignal));

  //////////////////////////
  // Fitting and Plotting //
  //////////////////////////

  //Create canvas
  TCanvas* c1 = new TCanvas("c1","c1");
  //c1->Divide(1,2);
  c1->cd();
  //c1->SetLogx();
  TString bd1_title = "Neutron Recoil Fits - BD1";
  TString bd2_title = "Neutron Recoil Fits - BD2";
  TString bd3_title = "Neutron Recoil Fits - BD3";
  TString bd4_title = "Neutron Recoil Fits - BD4";
  TString bd5_title = "Neutron Recoil Fits - BD5";
  TString bd6_title = "Neutron Recoil Fits - BD6";
  TString bd7_title = "Neutron Recoil Fits - BD7";
  TString bd8_title = "Neutron Recoil Fits - BD8";
  TString plotTitle_List[] = { bd1_title, bd2_title, bd3_title, bd4_title, bd5_title, bd6_title, bd7_title, bd8_title };
  TString plotTitle = plotTitle_List[bdNum-1];
  RooPlot* keVFrame = keV.frame(Title(plotTitle));

  //Fit.
  model.fitTo(sourceData,Save(1),Verbose(3),Range(fitMin,fitMax),PrintEvalErrors(10));

  //Plot Model
  keV.setBins(binning);
  nNoise.setConstant();
  nSignal.setConstant();
  RooDataSet* gammaData = signal.generate(keV,100000);
  RooDataHist* gammaHist = gammaData->binnedClone();
  RooHistPdf binnedGamma("binnedGamma","Binned Gamma", keV, *gammaHist, 0);
  RooAddPdf binnedModel("binnedModel","Binned Model", RooArgSet(noiseHistPdf, binnedGamma),RooArgSet(nNoise,nSignal));
  sourceData.plotOn(keVFrame,Name("sourceData"),MarkerColor(1),FillColor(0),Binning(binning),Range(fitMin,fitMax));
  binnedModel.plotOn(keVFrame,Name("noiseHistPdf"),Components(noiseHistPdf),LineColor(4),FillColor(0),LineWidth(3),Range(fitMin,fitMax));
  binnedModel.plotOn(keVFrame,Name("binnedGamma"),Components(binnedGamma),LineColor(2),FillColor(46),LineWidth(3),Range(fitMin,fitMax));
  binnedModel.plotOn(keVFrame,Name("binnedModel"),LineColor(8),FillColor(0),LineWidth(3),Range(fitMin,fitMax));
  RooHist* resHist = keVFrame->residHist();
  keVFrame->GetXaxis()->SetTitle("Energy [keVee]");
  keVFrame->GetYaxis()->SetTitle("Counts / 0.2 keV");
  keVFrame->GetXaxis()->CenterTitle();
  keVFrame->GetYaxis()->CenterTitle();
  keVFrame->SetTitle("");
  keVFrame->SetMinimum(0.001);
  keVFrame->SetMaximum(10e4);
  keVFrame->GetXaxis()->SetRangeUser(fitMin,fitMax);

  //Draw
  keVFrame->Draw();
  c1->SetLogy();
  c1->Modified();
  c1->Update();
  TString imagePath = "/var/phy/project/phil/cma46/CeBr3/cebr3-qf-analysis/Plots/BestFit-Neutrons-BD8.pdf";
  c1->Print(imagePath);

  //Draw residuals.
  TCanvas* c2 = new TCanvas("c2","c2");
  c2->cd();
  TString bd1_res = "Neutron Recoil Residuals - BD1";
  TString bd2_res = "Neutron Recoil Residuals - BD2";
  TString bd3_res = "Neutron Recoil Residuals - BD3";
  TString bd4_res = "Neutron Recoil Residuals - BD4";
  TString bd5_res = "Neutron Recoil Residuals - BD5";
  TString bd6_res = "Neutron Recoil Residuals - BD6";
  TString bd7_res = "Neutron Recoil Residuals - BD7";
  TString bd8_res = "Neutron Recoil Residuals - BD8";
  TString resTitle_List[] = { bd1_res, bd2_res, bd3_res, bd4_res, bd5_res, bd6_res, bd7_res, bd8_res };
  TString resTitle = resTitle_List[bdNum-1];
  RooPlot* resFrame = keV.frame(Title(resTitle));
  resFrame->GetXaxis()->CenterTitle();
  resFrame->GetXaxis()->SetRangeUser(fitMin,fitMax);
  resFrame->GetXaxis()->SetTitle("keVee");
  resFrame->GetYaxis()->CenterTitle();
  resFrame->GetYaxis()->SetTitle("Counts / 0.2 keV");
  resFrame->addPlotable(resHist,"p");
  resFrame->Draw(); 
  c2->Modified();
  c2->Update();

  TString resPath = "/var/phy/project/phil/cma46/CeBr3/cebr3-qf-analysis/Plots/Residuals-Neutrons-BD8.pdf";
  c2->Print(resPath);

}
















//Plot all background PDFs together to compare.

void compareBackgrounds(){

  //////////////////////////
  // Setting Up Variables //
  //////////////////////////

  //Path to our data.
  TString filename = "/var/phy/project/phil/cma46/CeBr3/processedRuns_FixedIntegral/ProductionRun_1/processedFilelist.txt";

  //Don't plot canvas on screen.
  gStyle->SetCanvasPreferGL(kTRUE);
  gROOT->SetBatch(1);

  //Energy range to fit.
  Long64_t fitMin=0;
  Long64_t fitMax=5;
  Long64_t keVMin = 0;
  Long64_t keVMax = 100;
  Double_t binning = 1000;

  //Conversion from BD Number to LS Channel.
  Int_t ls_channelList[] = { 8, 9, 10, 11, 12, 13, 14, 15 };

  //Declare our observable.
  RooRealVar keV("keV","keV",keVMin,keVMax);
  RooArgSet observables("observables");
  observables.add(keV);

  //Unit conversion.
  Double_t adc_to_keV = 2.00114445451e-3;

  //Data Cuts.
  Double_t time_Low = 325;
  Double_t time_High = 355;
  Double_t psd_Low = 0.26;
  Double_t psd_High = 0.6;
  Double_t integral_Low = 10000;
  Double_t integral_High = 35000;

  /////////////////////////
  // Loading Source Data //
  /////////////////////////

  //Create empty data sets.
  RooDataSet noise1("noise1","noise1",keV);	
  RooDataSet noise2("noise2","noise2",keV);
  RooDataSet noise3("noise3","noise3",keV);
  RooDataSet noise4("noise4","noise4",keV);
  RooDataSet noise5("noise5","noise5",keV);
  RooDataSet noise6("noise6","noise6",keV);
  RooDataSet noise7("noise7","noise7",keV);
  RooDataSet noise8("noise8","noise8",keV);
  RooDataSet noise[] = {noise1,noise2,noise3,noise4,noise5,noise6,noise7,noise8};

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

  //Create a TChain to hold all of our data.
  TChain* chain = new TChain("analysisTree");
  chain->SetBranchAddress("LS_channel", &LS_channel);
  chain->SetBranchAddress("LS_psd",&LS_psd);
  chain->SetBranchAddress("LS_integral",&LS_integral);
  chain->SetBranchAddress("LS_timeToBPM",&LS_timeToBPM);
  chain->SetBranchAddress("scatterer_foundPulse",&scatterer_foundPulse);
  chain->SetBranchAddress("scatterer_integral",&scatterer_integral);
  chain->SetBranchAddress("scatterer_noise",&scatterer_noise);
  chain->SetBranchAddress("scatterer_timeToBPM",&scatterer_timeToBPM);
  chain->SetBranchAddress("scatterer_timeToBD",&scatterer_timeToBD);
  chain->SetBranchAddress("scatterer_onsetTime",&scatterer_onsetTime);
  chain->SetBranchAddress("scatterer_peakHeight",&scatterer_peakHeight);
  chain->SetBranchAddress("scatterer_baseline",&scatterer_baseline);
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

  //Now step through our TChain and pick out the right data.
  Int_t numSourceEntries = chain->GetEntries();
  for (Long64_t entry=0; entry < numSourceEntries; entry++) {
	chain->GetEntry(entry);
	for( Int_t i = 0; i < 8; i++ ){
		Int_t bd_channel = ls_channelList[i];
		if(LS_channel == bd_channel){
			if( ( LS_psd > psd_Low ) && ( LS_psd < psd_High ) ){
				if( ( LS_integral > integral_Low ) && ( LS_integral < integral_High ) ){
					if( ( LS_timeToBPM > time_Low ) && ( LS_timeToBPM < time_High ) ){
						enrg = scatterer_noise * adc_to_keV;
						keV = enrg;
						noise[i].add( keV );

					}
				}
			}
		}
	}
		
	//Give user information on screen
	if (entry % 10000 == 0) {
		printf("Processing entry %lu of %lu entries\n", entry, numSourceEntries);
	}
  }

  cout << "Filled!" << endl;
  cout << "Generating RooPDFs..." << endl;

  /////////////////
  // Making PDFs //
  /////////////////
  for( Int_t i = 0; i < 8; i++ ){
  	((RooRealVar*)noise[i].get()->find("keV"))->setBins(binning);
  }
  RooDataHist* noiseHist1 = noise[0].binnedClone();
  RooDataHist* noiseHist2 = noise[1].binnedClone();
  RooDataHist* noiseHist3 = noise[2].binnedClone();
  RooDataHist* noiseHist4 = noise[3].binnedClone();
  RooDataHist* noiseHist5 = noise[4].binnedClone();
  RooDataHist* noiseHist6 = noise[5].binnedClone();
  RooDataHist* noiseHist7 = noise[6].binnedClone();
  RooDataHist* noiseHist8 = noise[7].binnedClone();
  RooHistPdf noiseHistPdf1("noiseHistPdf1","Noise1 Histogram",keV,*noiseHist1,0);
  RooHistPdf noiseHistPdf2("noiseHistPdf2","Noise2 Histogram",keV,*noiseHist2,0);
  RooHistPdf noiseHistPdf3("noiseHistPdf3","Noise3 Histogram",keV,*noiseHist3,0);
  RooHistPdf noiseHistPdf4("noiseHistPdf4","Noise4 Histogram",keV,*noiseHist4,0);
  RooHistPdf noiseHistPdf5("noiseHistPdf5","Noise5 Histogram",keV,*noiseHist5,0);
  RooHistPdf noiseHistPdf6("noiseHistPdf6","Noise6 Histogram",keV,*noiseHist6,0);
  RooHistPdf noiseHistPdf7("noiseHistPdf7","Noise7 Histogram",keV,*noiseHist7,0);
  RooHistPdf noiseHistPdf8("noiseHistPdf8","Noise8 Histogram",keV,*noiseHist8,0);

  //////////////
  // Plotting //
  //////////////

  //Create canvas
  TCanvas* c1 = new TCanvas("c1","c1");
  //c1->Divide(1,2);
  c1->cd();
  //c1->SetLogx();
  RooPlot* keVFrame = keV.frame(Title("Background Spectra - All Detectors"));

  //Plot Model
  noiseHistPdf1.plotOn(keVFrame,Name("noiseHistPdf1"),LineColor(1),FillColor(0),LineWidth(1),Range(fitMin,fitMax));
  noiseHistPdf2.plotOn(keVFrame,Name("noiseHistPdf2"),LineColor(2),FillColor(0),LineWidth(1),Range(fitMin,fitMax));
  noiseHistPdf3.plotOn(keVFrame,Name("noiseHistPdf3"),LineColor(3),FillColor(0),LineWidth(1),Range(fitMin,fitMax));
  noiseHistPdf4.plotOn(keVFrame,Name("noiseHistPdf4"),LineColor(4),FillColor(0),LineWidth(1),Range(fitMin,fitMax));
  noiseHistPdf5.plotOn(keVFrame,Name("noiseHistPdf5"),LineColor(5),FillColor(0),LineWidth(1),Range(fitMin,fitMax));
  noiseHistPdf6.plotOn(keVFrame,Name("noiseHistPdf6"),LineColor(6),FillColor(0),LineWidth(1),Range(fitMin,fitMax));
  noiseHistPdf7.plotOn(keVFrame,Name("noiseHistPdf7"),LineColor(7),FillColor(0),LineWidth(1),Range(fitMin,fitMax));
  noiseHistPdf8.plotOn(keVFrame,Name("noiseHistPdf8"),LineColor(8),FillColor(0),LineWidth(1),Range(fitMin,fitMax));
  keVFrame->GetXaxis()->SetTitle("Energy [keVee]");
  keVFrame->GetYaxis()->SetTitle("Counts / 0.1 keV");
  keVFrame->GetXaxis()->CenterTitle();
  keVFrame->GetYaxis()->CenterTitle();
  keVFrame->SetTitle("");
  keVFrame->SetMinimum(0.00001);
  keVFrame->SetMaximum(10);
  keVFrame->GetXaxis()->SetRangeUser(fitMin,fitMax);

  //Draw
  keVFrame->Draw();
  c1->SetLogy();
  c1->Modified();
  c1->Update();
  TString imagePath = "/var/phy/project/phil/cma46/CeBr3/cebr3-qf-analysis/Plots/AllBGNDs.pdf";  
  c1->Print(imagePath);
  
}


