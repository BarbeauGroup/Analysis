//Code to fit MCNP Simulations to CeBr3 Backing Detector Calibration Data
//Written by C. Awe (based heavily on work by S. Hedges)
//03/02/2021

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

//Create Gaussian random numbers using a Box Muller distribution.
//Code from https://stackoverflow.com/questions/19944111/creating-a-gaussian-random-generator-with-a-mean-and-standard-deviation
double rand_normal(double mean, double stddev )
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

double genRandomData( double alpha, double beta, double enrg ){

	Double_t sdev = alpha + beta * sqrt( enrg );
	Double_t sample = rand_normal(enrg, sdev );
	sample = sample;
	return sample;
}

//Pass in alpha coefficient for resolution, return RooFitResult from fit
RooFitResult getNll(Double_t alpha=4.7, Double_t beta=0.015, Int_t fitNum=0, Int_t bdNum=8 ) {

//Don't plot canvas on screen.
gStyle->SetCanvasPreferGL(kTRUE);
gROOT->SetBatch(1);

//Fitter options.
ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(20000); 
ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(2000);

//Integral range to fit.
Long64_t fitMin=0;
Long64_t fitMax=30e3;
Double_t binning = 300;

//Conversion from BD Number to LS Channel.
Int_t ls_channelList[] = { 8, 9, 10, 11, 12, 13, 14, 15 };

//Converts BD Number to MCNP Cell.
Int_t bd_cellList[] = { 306, 314, 322, 330, 326, 318, 310, 302 };

//Simulation file paths.
TString na22_sim = "~/Desktop/cebr3-qf-analysis/analysis/BD_Sims/dumn1-bdcalibration-na22.root";
TString co60_sim = "~/Desktop/cebr3-qf-analysis/analysis/BD_Sims/dumn1-bdcalibration-co60.root";

//Declare our observable.
TCanvas* fitCanvas = new TCanvas("fitCanvas","fitCanvas");
RooRealVar integral("integral","integral",fitMin,fitMax);
RooArgSet observables("observables");
observables.add(integral);

//Load the simulations.
TFile* na22_simFile = new TFile( na22_sim );
TFile* co60_simFile = new TFile( co60_sim );

//////////////////////////////////
//////Load simulation ttrees//////
//////////////////////////////////

TTree *na22_simTree = (TTree*)na22_simFile->Get("totalEnergyTree");
TTree *co60_simTree = (TTree*)co60_simFile->Get("totalEnergyTree");

//Load energy branches.
Double_t energy;
Double_t cellNum;
na22_simTree->SetBranchAddress("energy",&energy);
na22_simTree->SetBranchAddress("cellNum",&cellNum);
co60_simTree->SetBranchAddress("energy",&energy);
co60_simTree->SetBranchAddress("cellNum",&cellNum);

//Get number of entries in both trees.
Long64_t na22_numSimEntries=na22_simTree->GetEntries();
Long64_t co60_numSimEntries=co60_simTree->GetEntries();

//Create gaussian to generate random data
RooRealVar mean("mean","mean",fitMin,fitMax);
RooRealVar sigma("sigma","sigma",fitMin,fitMax); //sigma probably won't ever be greater than fitMax.
RooGaussian gauss("gauss","gauss",integral,mean,sigma);

//Create data sets to hold smeared simulation data.
RooDataSet na22_smearedData("na22_smearedData","na22_smearedData",integral);
RooDataSet co60_smearedData("co60_smearedData","co60_smearedData",integral);
Int_t valsToGen=10; //create 10 gaussian-distributed data points based on alpha for each point in simulation
default_random_engine generator;

/////////////////////////////////////
/////////MAKING SMEARED DATA/////////
/////////////////////////////////////

Double_t adc_to_keV = 20; //An initial guess for keV->ADC.
Int_t bd_cell = bd_cellList[bdNum-1]; //Converts BD_Num to a channel number.

//Start with Na22.
cout<<"Generating smeared Na22 data..."<<endl;
for (Long64_t entry=0; entry < na22_numSimEntries; entry++) {

	//Get entry
	na22_simTree->GetEntry(entry);
	if( cellNum == bd_cell ){
		Double_t enrg = energy;
		for (Int_t val=0; val < valsToGen; val++) {
			Double_t sample = genRandomData( alpha, beta, enrg );
			sample = sample * adc_to_keV; //Convert it.
			integral = sample;
			na22_smearedData.add( integral );
		}
	}
	//Give user information on screen
	if (entry % 10000 == 0) {
		printf("Processing entry %lu of %lu entry\n", entry, na22_numSimEntries);
	}
}

cout << "Generated!" << endl;
//Now do Cobalt 60.
cout<<"Generating smeared Co60 data..."<<endl;
for (Long64_t entry=0; entry < co60_numSimEntries; entry++) {

	//Get entry
	co60_simTree->GetEntry(entry);
	if( cellNum == bd_cell ){
		Double_t enrg = energy;
		for (Int_t val=0; val < valsToGen; val++) {
			Double_t sample = genRandomData( alpha, beta, enrg );
			sample = sample * adc_to_keV; //Convert it.
			integral = sample;
			co60_smearedData.add( integral );
		}
	}
	//Give user information on screen
	if (entry % 10000 == 0) {
		printf("Processing entry %lu of %lu entry\n", entry, co60_numSimEntries);
	}
}

cout << "Generated!" << endl;

na22_smearedData.Print("V");
co60_smearedData.Print("V");

/////////////////////////////////////
/////////LOADING SOURCE DATA/////////
/////////////////////////////////////

//Create empty data sets
RooDataSet na22_sourceData("na22_sourceData","na22_sourceData",integral);
RooDataSet co60_sourceData("co60_sourceData","co60_sourceData",integral);
RooDataSet na22_bgnd("na22_bgnd","na22_bgnd",integral);	
RooDataSet co60_bgnd("co60_bgnd","co60_bgnd",integral);	
//Long64_t numSrcDataInRange=0;

//Path to our data.
TString na22_filename = "~/Desktop/cebr3-qf-analysis/analysis/BD_Sims/Processed_SIS3316Raw_20200207140100_1.root";
TString co60_filename = "~/Desktop/cebr3-qf-analysis/analysis/BD_Sims/Processed_SIS3316Raw_20200207140645_1.root";
TString bgnd_filename = "~/Desktop/cebr3-qf-analysis/analysis/BD_Sims/Processed_SIS3316Raw_20200207141225_1.root";

//Variables to assign to ttree branches.
Double_t LS_integral;
Int_t LS_channel;

//Create a TChain to hold all of our data.
TFile* na22_File = new TFile(na22_filename);
TFile* co60_File = new TFile(co60_filename);
TFile* bgnd_File = new TFile(bgnd_filename);
TTree* na22_Tree = (TTree*) na22_File->Get("analysisTree");
TTree* co60_Tree = (TTree*) co60_File->Get("analysisTree");
TTree* bgnd_Tree = (TTree*) bgnd_File->Get("analysisTree");
na22_Tree->SetBranchAddress("LS_integral",&LS_integral);
co60_Tree->SetBranchAddress("LS_integral",&LS_integral);
bgnd_Tree->SetBranchAddress("LS_integral",&LS_integral);
na22_Tree->SetBranchAddress("LS_channel",&LS_channel);
co60_Tree->SetBranchAddress("LS_channel",&LS_channel);
bgnd_Tree->SetBranchAddress("LS_channel",&LS_channel);

//Now step through our TFiles and pick out the right data.
Int_t bd_channel = ls_channelList[bdNum-1]; //Converts BD_Num to a channel number.

//Na22 first.
Int_t na22_numSourceEntries = na22_Tree->GetEntries();
for (Long64_t entry=0; entry < na22_numSourceEntries; entry++) {
	na22_Tree->GetEntry(entry);
	if( LS_channel == bd_channel ){
		if( ( LS_integral > fitMin ) && ( LS_integral < fitMax ) ){
			integral = LS_integral;
			na22_sourceData.add( integral );
		}
	}	
	//Give user information on screen
	if (entry % 10000 == 0) {
		printf("Processing entry %lu of %lu entries\n", entry, na22_numSourceEntries);
	}
}
cout << "Filled!" << endl;

//Now do Co60.
Int_t co60_numSourceEntries = co60_Tree->GetEntries();
for (Long64_t entry=0; entry < co60_numSourceEntries; entry++) {
	co60_Tree->GetEntry(entry);
	//cout << LS_channel << endl;
	//cout << bd_channel << endl;
	if( LS_channel == bd_channel ){
		if( ( LS_integral > fitMin ) && ( LS_integral < fitMax ) ){
			//cout << "Should be adding data" << endl;
			integral = LS_integral;
			co60_sourceData.add( integral );
		}
	}	
	//Give user information on screen
	if (entry % 10000 == 0) {
		printf("Processing entry %lu of %lu entries\n", entry, co60_numSourceEntries);
	}
}
cout << "Filled!" << endl;

//And finally the background.
Int_t bgnd_numSourceEntries = bgnd_Tree->GetEntries();
for (Long64_t entry=0; entry < bgnd_numSourceEntries; entry++) {
	bgnd_Tree->GetEntry(entry);
	if( LS_channel == bd_channel ){
		if( ( LS_integral > fitMin ) && ( LS_integral < fitMax ) ){
			integral = LS_integral;
			na22_bgnd.add( integral );
			co60_bgnd.add( integral );
		}
	}	
	//Give user information on screen
	if (entry % 10000 == 0) {
		printf("Processing entry %lu of %lu entries\n", entry, bgnd_numSourceEntries);
	}
}
cout << "Filled!" << endl;

cout << "Printing source data info." << endl;
na22_sourceData.Print("V");
co60_sourceData.Print("V");
na22_bgnd.Print("V");
co60_bgnd.Print("V");

cout << "Generating RooPDFs..." << endl;

////////////////////////////////
/////////Make Hist PDFs/////////
////////////////////////////////

//Set binnings
((RooRealVar*)na22_smearedData.get()->find("integral"))->setBins(binning);
((RooRealVar*)co60_smearedData.get()->find("integral"))->setBins(binning);

//Represent simulations data as PDFs.
RooDataHist* na22_simHist = na22_smearedData.binnedClone();
RooDataHist* co60_simHist = co60_smearedData.binnedClone();

//Make a noise hist pdf.
((RooRealVar*)na22_bgnd.get()->find("integral"))->setBins(binning);
((RooRealVar*)co60_bgnd.get()->find("integral"))->setBins(binning);
RooDataHist* na22_bgndHist = na22_bgnd.binnedClone();
RooHistPdf na22_bgndHistPdf("na22_bgndHistPdf","Background Histogram for Na22",integral,*na22_bgndHist,0);
RooDataHist* co60_bgndHist = co60_bgnd.binnedClone();
RooHistPdf co60_bgndHistPdf("co60_bgndHistPdf","Background Histogram for Co60",integral,*co60_bgndHist,0);

//Create scaling factors.
RooRealVar tempSlope("tempSlope","tempSlope",1,0.1,10); //Converts keV to ADC.
RooFormulaVar tempScaling("tempScaling","integral*tempSlope",RooArgSet(integral,tempSlope));
RooHistPdf na22_tempSimHist("na22_tempSimHist","na22_tempSimHist",tempScaling,integral,*na22_simHist,1);
RooRealVar na22_bgndEvents("na22_bgndEvents","Background Counts for Na22",2e6,0,1e7);
RooRealVar na22_gammaEvents("na22_gammaEvents","Gamma Scatter Counts for Na22",30e5,1,1e7);
RooAddPdf na22_tempPdf("na22_tempPdf","Simulation + Background",RooArgList(na22_bgndHistPdf,na22_tempSimHist),RooArgList(na22_bgndEvents,na22_gammaEvents));
RooHistPdf co60_tempSimHist("co60_tempSimHist","co60_tempSimHist",tempScaling,integral,*co60_simHist,1);
RooRealVar co60_bgndEvents("co60_bgndEvents","Background Counts for co60",2e6,0,1e7);
RooRealVar co60_gammaEvents("co60_gammaEvents","Gamma Scatter Counts for co60",30e5,1,1e7);
RooAddPdf co60_tempPdf("co60_tempPdf","Simulation + Background",RooArgList(co60_bgndHistPdf,co60_tempSimHist),RooArgList(co60_bgndEvents,co60_gammaEvents));

//Categroize data sets to perform simultaneous fits.
RooCategory source("source","source");
source.defineType("na22");
source.defineType("co60");

//Build combined data set.
RooDataSet combData("combData","Combined Data Sets",RooArgSet(integral),Index(source),Import("na22",na22_sourceData),Import("co60",co60_sourceData));

//Build a simultaneous PDF.
RooSimultaneous simPdf("simPdf","Simultaneous PDF",source);
simPdf.addPdf(na22_tempPdf,"na22");
simPdf.addPdf(co60_tempPdf,"co60");

//Fit.
RooFitResult* tempRes = simPdf.fitTo(combData,Save(1),Verbose(3),PrintEvalErrors(10));
auto pars = tempRes->floatParsFinal();

///////////////////////////////////
/////////Plotting/Fits/////////////
///////////////////////////////////

//Create canvas
TCanvas* c1 = new TCanvas("c1","c1");
c1->Divide(1,2);
c1->cd(1);
c1->SetLogy();
TString na22_plotTitle = "Backing Detector Calibration - Na22";
RooPlot* na22_Frame = integral.frame(Title(na22_plotTitle));
//Plot Model
na22_sourceData.plotOn(na22_Frame,Name("na22_sourceData"),MarkerColor(1),FillColor(0),Binning(binning));
simPdf.plotOn(na22_Frame,Slice(source, "na22"),ProjWData(source, combData),Name("Total Model"),LineColor(38),FillColor(0),LineWidth(3),Range(fitMin,fitMax));
simPdf.plotOn(na22_Frame,Slice(source, "na22"),Components("na22_bgndHistPdf"),ProjWData(source, combData),Name("Background"),LineColor(46),FillColor(0),LineWidth(3),Range(fitMin,fitMax));
simPdf.plotOn(na22_Frame,Slice(source, "na22"),Components("na22_tempSimHist"),ProjWData(source, combData),Name("Simulation"),LineColor(29),FillColor(0),LineWidth(3),Range(fitMin,fitMax));
na22_Frame->GetXaxis()->SetTitle("Integral [Arb.]");
na22_Frame->GetYaxis()->SetTitle("Counts / 100 ADC");
na22_Frame->GetXaxis()->CenterTitle();
na22_Frame->GetYaxis()->CenterTitle();
na22_Frame->SetTitle("");
na22_Frame->GetYaxis()->SetRangeUser(0.01,5e2);
na22_Frame->GetXaxis()->SetRangeUser(0,30e3);
//Draw
na22_Frame->Draw();
c1->cd(2);
c1->SetLogy();
TString co60_plotTitle = "Backing Detector Calibration - Co60";
RooPlot* co60_Frame = integral.frame(Title(co60_plotTitle));
//Plot Model
co60_sourceData.plotOn(co60_Frame,Name("co60_sourceData"),MarkerColor(1),FillColor(0),Binning(binning));
simPdf.plotOn(co60_Frame,Slice(source, "co60"),ProjWData(source, combData),Name("Total Model"),LineColor(38),FillColor(0),LineWidth(3),Range(fitMin,fitMax));
simPdf.plotOn(co60_Frame,Slice(source, "co60"),Components("co60_bgndHistPdf"),ProjWData(source, combData),Name("Background"),LineColor(46),FillColor(0),LineWidth(3),Range(fitMin,fitMax));
simPdf.plotOn(co60_Frame,Slice(source, "co60"),Components("co60_tempSimHist"),ProjWData(source, combData),Name("Simulation"),LineColor(29),FillColor(0),LineWidth(3),Range(fitMin,fitMax));
co60_Frame->GetXaxis()->SetTitle("Integral [Arb.]");
co60_Frame->GetYaxis()->SetTitle("Counts / 100 ADC");
co60_Frame->GetXaxis()->CenterTitle();
co60_Frame->GetYaxis()->CenterTitle();
co60_Frame->SetTitle("");
co60_Frame->GetYaxis()->SetRangeUser(0.01,5e2);
co60_Frame->GetXaxis()->SetRangeUser(0,30e3);
//Draw
na22_Frame->Draw();
c1->Modified();
c1->Update();

TString imagePath = "~/Desktop/cebr3-qf-analysis/Plots/";
imagePath.Append("BestFit_BD");
imagePath.Append(".pdf");
c1->Print(imagePath);
//c1->WaitPrimitive();

//return the fit.
return *tempRes;

}

void scanResolution( Int_t bdNum ) {

	//Select appropriate file name to save output
	TString fileName;
	TString fname1 = "~/Desktop/cebr3-qf-analysis/analysis/BD_Sims/calibraitonScan_bd1.root";
	TString fname2 = "~/Desktop/cebr3-qf-analysis/analysis/BD_Sims/calibrationScan_bd2.root";
	TString fname3 = "~/Desktop/cebr3-qf-analysis/analysis/BD_Sims/calibrationScan_bd3.root";
	TString fname4 = "~/Desktop/cebr3-qf-analysis/analysis/BD_Sims/calibrationScan_bd4.root";
	TString fname5 = "~/Desktop/cebr3-qf-analysis/analysis/BD_Sims/calibrationScan_bd5.root";
	TString fname6 = "~/Desktop/cebr3-qf-analysis/analysis/BD_Sims/calibrationScan_bd6.root";
	TString fname7 = "~/Desktop/cebr3-qf-analysis/analysis/BD_Sims/calibrationScan_bd7.root";
	TString fname8 = "~/Desktop/cebr3-qf-analysis/analysis/BD_Sims/calibrationScan_bd8.root";
	TString fname_List[] = { fname1, fname2, fname3, fname4, fname5, fname6, fname7, fname8 };
	fileName = fname_List[bdNum-1];

	//Create file
	TFile* f = new TFile(fileName,"RECREATE");

	//Set number of alphas to run fits for, minimum alpha value, alpha step value 
	const Int_t numAlphas=10;
	const Int_t numBetas=10;
	const Int_t numFits=numAlphas*numBetas;
	Double_t alphaMin=4.0;
	Double_t alphaStep=0.1;
	Double_t alpha; //stores present alpha value
	Double_t betaMin=0.01;
	Double_t betaStep=0.001;
	Double_t beta; //stores present alpha value

	//Dummies for generating TTree.
	Double_t nll;
	Double_t slope;
	Double_t status;

	//Arrays for holding data until we normalize the NLL.
	Double_t alphas[numFits]={};
	Double_t betas[numFits]={}; 
	Double_t nlls[numFits]={};
	Double_t slopes[numFits]={};

	//Other dummy variables.
	Int_t numGoodFits=0; //Stores number of successful fits
	Int_t isFirstFit=1; //Flag to determine if this is the first fit, for normalizing NLL value
	Double_t minNll; //Stores min NLL value, for normalizing

	cout <<"Starting For loop..." << endl;
	for (Int_t alphaNum=0; alphaNum < numAlphas; alphaNum++) {
		//Update alpha
		alpha=alphaMin+alphaNum*alphaStep;
		for( Int_t betaNum=0; betaNum < numBetas; betaNum++ ){
			cout << "Beta is: " << beta << endl;
			beta=betaMin+betaNum*betaStep;
			
			//Do the fits.
			RooFitResult res = getNll(alpha,beta,alphaNum + betaNum, bdNum);	

			//Get fit status.
			status = res.status();	

			if (status == 0) { //status=0 means fit successful
				if (isFirstFit==1) {
					isFirstFit=0; //Set first fit to 0
					minNll=res.minNll(); //Set initial value for the minNll
				}
				else {
					if (res.minNll() < minNll) { //Check to see if current NLL < minNll, if so update minNll
						minNll=res.minNll();
					}
				}
		
				//Get fit parameters
				auto pars = res.floatParsFinal();
				alphas[ numGoodFits ] = alpha;
				betas[ numGoodFits ] = beta;
				//Store slopes from successful fits.
				slopes[ numGoodFits ] = ((RooRealVar*) pars.find("tempSlope"))->getVal();
				//Store Nlls from successful fits
				nlls[ numGoodFits ] = res.minNll();
				//Increment number of good fits found
				numGoodFits++;
			}
		}
	}
	//Outside of fitting loop
	cout << "The minNll was: " << minNll << endl;

	//Declare a new TTree and set the branches.
	cout << "Generating new TTree..." << endl;
    	TTree* nllTree = new TTree("nllTree","Tree Containing NLL Fit Values");
    	nllTree->Branch("alpha",&alpha,"alpha/D");
	nllTree->Branch("beta",&beta,"beta/D");
	nllTree->Branch("nll",&nll,"nll/D");
	nllTree->Branch("slope",&slope,"slope/D");
    	cout << "TTree generated!" << endl;
    	nllTree->SetDirectory(0);
	TH2D* nllHist = new TH2D("nllHist","nllHist",numAlphas,alphaMin,alphaMin + numAlphas*alphaStep, numBetas, betaMin, betaMin + numBetas*betaStep);

	//Subtract minNll value from all nlls so lowest one has a value of 0 and populate the Tree.
	for (Int_t i=0; i < numGoodFits; i++) {
		nlls[i] -= minNll;
		nll = nlls[i];
		alpha = alphas[i];
		beta = betas[i];
		slope = slopes[i];
		nllHist->SetBinContent( nllHist->FindBin( alpha, beta ), nll );
		nllTree->Fill();
		if( nll == 0 ){
			cout << "The minimizing alpha is: " << alpha << endl;
			cout << "The minimizing beta is: " << beta << endl;
		}
		
	}

	//Write TTree and Histogram.
	f->cd();
	nllTree->Write();
	nllHist->SetMaximum(5);
	nllHist->Write();

	//Close file and TTree.
	//nllTree->Delete();
	f->Close();
}






























    

