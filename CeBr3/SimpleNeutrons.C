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
  Long64_t fitMax=6;
  Long64_t keVMin = 0;
  Long64_t keVMax = 100;
  Double_t binning = 1000;

  //Conversion from BD Number to LS Channel.
  Int_t ls_channelList[] = { 8, 9, 10, 11, 12, 13, 14, 15 };

  //Bin adjustments - this is just an excersion exercise to help set error bars on the QF values.
  Int_t bin_flag = 0; //Set to 0 for standard fits, 1 to adjust for afterglow in the noise pdf.
  Int_t bin1_list[] = {1150,1175,1400,1325,1650,1350,1200,1450};
  Int_t bin2_list[] = {-800,-900,-1175,-1725,-1650,-1020,-750,-1000};
  Int_t bin3_list[] = {-175,-200,-200,400,0,-250,-250,200};
  Int_t bin4_list[] = {-50,-100,-10,400,100,-100,-100,-100};

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
	if(LS_channel == bd_channel){
		if( ( LS_psd > psd_Low ) && ( LS_psd < psd_High ) ){
			if( ( LS_integral > integral_Low ) && ( LS_integral < integral_High ) ){
				if( ( LS_timeToBPM > time_Low ) && ( LS_timeToBPM < time_High ) ){
					Double_t enrg = scatterer_integral * adc_to_keV;
					//if( enrg < keVMax && enrg > keVMin ){
					keV = enrg;
					sourceData.add( keV );
					//}
					enrg = scatterer_noise * adc_to_keV;
					//if( enrg < keVMax && enrg > keVMin ){
					keV = enrg;
					noise.add( keV );
					//}
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
  if( bin_flag == 1 ){
	cout <<"Adjusting noise bins! Make sure this is what you want." << endl;
	noiseHist.Fill(0.0, bin1_list[bdNum - 1]);
	noiseHist.Fill(0.1, bin2_list[bdNum - 1]);
	noiseHist.Fill(0.2, bin3_list[bdNum - 1]);
	noiseHist.Fill(0.3, bin4_list[bdNum - 1]);
  }

  RooHistPdf noiseHistPdf("noiseHistPdf","Noise Histogram",keV,*noiseHist,0);
  RooRealVar nNoise("nNoise","nNoise",10000,5000,100000);
  //RooRealVar mean("mean","mean",3.5,0.1,6);
  //RooRealVar sigma("sigma","sigma",0.2,0.01,1); //sigma probably won't ever be greater than fitMax.
  //RooGaussian gauss("gauss","gauss",keV,mean,sigma);
  RooRealVar mu("mu","mu",0); //Offset built into gamma distribution.
  mu.setConstant();
  RooRealVar beta("beta","beta",0.5,0,1); //Decay constant built into gamma pdf, roughly corresponds to the spread.
  RooRealVar gamma("gamma","gamma",5,0,20); //Another gamma pdf constant, harder to describe. gamma * beta gives the mean of the distribution.
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
  sourceData.plotOn(keVFrame,Name("sourceData"),MarkerColor(1),FillColor(0),Binning(binning),Range(fitMin,fitMax));
  model.plotOn(keVFrame,Name("noiseHistPdf"),Components(noiseHistPdf),LineColor(4),FillColor(0),LineWidth(3),Range(fitMin,fitMax));
  model.plotOn(keVFrame,Name("signal"),Components(signal),LineColor(2),FillColor(46),LineWidth(3),Binning(binning),Range(fitMin,fitMax));
  model.plotOn(keVFrame,Name("model"),LineColor(8),FillColor(0),LineWidth(3),Range(fitMin,fitMax));
  RooHist* resHist = keVFrame->residHist();
  keVFrame->GetXaxis()->SetTitle("Energy [keVee]");
  keVFrame->GetYaxis()->SetTitle("Counts / 0.1 keV");
  keVFrame->GetXaxis()->CenterTitle();
  keVFrame->GetYaxis()->CenterTitle();
  keVFrame->SetTitle("");
  keVFrame->SetMinimum(0.001);
  keVFrame->SetMaximum(10e4);
  keVFrame->GetXaxis()->SetRangeUser(0,10);

  //Draw
  keVFrame->Draw();
  c1->SetLogy();
  c1->Modified();
  c1->Update();
  if( bin_flag == 1 ){
  	TString imagePath = "/var/phy/project/phil/cma46/CeBr3/cebr3-qf-analysis/Plots/BestFit-Afterglow-BD8.pdf"; //I apologize for hard-coding this. Fuck C.
  }
  else{
	TString imagePath = "/var/phy/project/phil/cma46/CeBr3/cebr3-qf-analysis/Plots/BestFit-Neutrons-BD1.pdf";
  }
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
  resFrame->GetXaxis()->SetRangeUser(0,10);
  resFrame->GetXaxis()->SetTitle("keVee");
  resFrame->GetYaxis()->CenterTitle();
  resFrame->GetYaxis()->SetTitle("Counts / 0.1 keV");
  resFrame->addPlotable(resHist,"p");
  resFrame->Draw(); 
  c2->Modified();
  c2->Update();

  if( bin_flag == 1 ){
  	TString resPath = "/var/phy/project/phil/cma46/CeBr3/cebr3-qf-analysis/Plots/Residuals-Afterglow-BD8.pdf"; 
  }
  else{
	TString resPath = "/var/phy/project/phil/cma46/CeBr3/cebr3-qf-analysis/Plots/Residuals-Neutrons-BD1.pdf";
  } //Ditto the above.Just tryng to finish my damn thesis.
  c2->Print(resPath);

}





