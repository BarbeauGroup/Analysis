//Code to fit MCNP Simulations to CeBr3 Calibration Data
//Written by C. Awe (based heavily on work by S. Hedges)
//12/14/2020

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
double rand_normal(double mean, double stddev)
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

//Pass in alpha coefficient for resolution, return RooFitResult from fit
RooFitResult getNll(Double_t alpha=4.8, Double_t beta=0.021, Int_t fitNum=0) {

//Don't plot canvas on screen
gStyle->SetCanvasPreferGL(kTRUE);
gROOT->SetBatch(1);

//Fitter options
ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(20000); 
ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(2000);

//Integral range to fit
Long64_t fitMin=36e3;
Long64_t fitMax=46e3;

//File paths.
TString simFilename = "~/Desktop/Analysis/CeBr3/Calibrations/dumn1-cebr3-gammacalib.root";
TString dataFilename = "~/Desktop/Analysis/CeBr3/Calibrations/Processed_SIS3316Raw_20200207141804_1.root";

//Load files.
TFile simFile( simFilename );
TFile dataFile( dataFilename );

//Create canvas
TCanvas* c1 = new TCanvas("c1","c1");
c1->cd();
c1->SetLogy();

/////////////////////////////////
//////Load simulation ttree//////
/////////////////////////////////

TTree *simTree = (TTree*)simFile.Get("totalEnergyTree");

//Load energy branch
Double_t energy;
simTree->SetBranchAddress("energy",&energy);

//Get number of entries
Long64_t numSimEntries=simTree->GetEntries();

//Load energy values into nubeIntegral from 0 to fitMax
Int_t minIntegral=6e3;
Int_t maxIntegral=250e3;
RooRealVar integral("integral","integral",minIntegral,maxIntegral);

//Create gaussian to generate random data
RooRealVar mean("mean","mean",minIntegral,maxIntegral);
RooRealVar sigma("sigma","sigma",minIntegral,maxIntegral); //sigma probably won't ever be greater than maxIntegral
RooGaussian gauss("gauss","gauss",integral,mean,sigma);

//Create data set to hold smeared data
RooDataSet smearedData("smearedData","smearedData",integral);

Int_t valsToGen=10; //create 10 gaussian-distributed data points based on alpha for each point in simulation
default_random_engine generator;

/////////////////////////////////////
/////////MAKING SMEARED DATA/////////
/////////////////////////////////////

cout<<"Generating smeared data..."<<endl;
Double_t calibGuess = 550; //Initial guess for the scaling to avoid the need for ultra fine binning. Normally 494

for (Long64_t entry=0; entry < numSimEntries; entry++) {

	//Get entry
	simTree->GetEntry(entry);
	Double_t enrg = energy;

	//Generate sigma, mean for random point generation
	sigma = alpha + beta * sqrt( enrg ); //Normally alpa + beta * sqrt (energy )
	Double_t sdev = alpha + beta * sqrt( enrg );

	//Generate valsToGen random gaussian values based on this distribution, add to smeared data set
	for (Int_t val=0; val < valsToGen; val++) {
		Double_t sample = rand_normal(enrg, sdev );
		sample = sample * calibGuess;
		integral = sample;
		smearedData.add(integral);

	}

	//Give user information on screen
	if (entry % 10000 == 0) {
		printf("Processing entry %lu of %lu entry\n", entry, numSimEntries);
	}
}
cout << "Generated!" << endl;

smearedData.Print("V");

/////////////////////////////////
////////Load source ttree////////
/////////////////////////////////

cout << "Loading source data..." << endl;

TTree *sourceTree = (TTree*)dataFile.Get("analysisTree");

//Get integral data.
Double_t scatterer_integral;
sourceTree->SetBranchAddress("scatterer_integral",&scatterer_integral);

//Get number of entries
Long64_t numSourceEntries=sourceTree->GetEntries();

cout << "Source data loaded!" << endl;

///////////////////////////////////////
//////////Fill source data set/////////
///////////////////////////////////////

//Create empty data set

RooDataSet sourceData("sourceData","sourceData",integral);	

Long64_t numSrcDataInRange=0;

//Fill source data set
cout << "Filling RooDataSet with source data..." << endl;

for (Long64_t entry=0; entry < numSourceEntries; entry++) {
	sourceTree->GetEntry(entry);
	integral=scatterer_integral;
    //cout << scatterer_integral << endl;
	if ( (scatterer_integral>minIntegral) && (scatterer_integral <maxIntegral)){
		numSrcDataInRange++;
		sourceData.add(integral);
		//cout << "Adding data point" << endl;
	}
		
	//Give user information on screen
	if (entry % 10000 == 0) {
		printf("Processing entry %lu of %lu entries\n", entry, numSourceEntries);
	}
}

cout << "Filled!" << endl;

cout << "Generating RooPDFs..." << endl;

////////////////////////////////
/////////Make Hist PDFs/////////
////////////////////////////////

//Set binnings
//((RooRealVar*)sourceData.get()->find("nubeIntegral"))->setBins(100); //Only necessary it want to fit to binned source data
((RooRealVar*)smearedData.get()->find("integral"))->setBins(244);

//Create simulated smeared data pdf with scaling
//Represent simulation data as a PDF
RooDataHist* simHist = smearedData.binnedClone();

//Add a Gaussian noise pedestal.
RooRealVar pedMean("pedMean","pedMean",0,-10,200);
RooRealVar pedSigma("pedSigma","pedSigma",500,0,1500);
RooGaussian pedestal("pedestal","Gaussian Pedestal PDF",integral,pedMean,pedSigma);

//Declare our RooFit Variables for Exponential Noise.
RooRealVar noiseDecay("noiseDecay","noiseDecay",20e3,100,300e5);
RooGaussModel pedestalResponse("pedestalResponse","pedestalResponse",integral,pedMean,pedSigma);
RooDecay decay("decay","Exponential Noise Decay PDF",integral,noiseDecay,pedestalResponse,RooDecay::SingleSided);

cout << "Printing Sim Hist" << endl;
simHist->Print( "V" );
//Create scaling factors
RooRealVar tempSlope("tempSlope","tempSlope",1,0.9,1.1);
//RooRealVar offset("offset","offset",50,0,2000);
RooFormulaVar tempScaling("tempScaling","integral/tempSlope",RooArgSet(integral,tempSlope));
RooHistPdf tempHistPdf("tempSimHistPdf","tempSimHistPdf",tempScaling,integral,*simHist,1);
RooRealVar pedEvents("pedEvents","pedEvents",2e6,1e4,1e7);
RooRealVar decayEvents("decayEvents","decayEvents",30e5,1e4,1e7);
RooRealVar simEvents("simEvents","simEvents",1e7,1e2,1e10);
RooAddPdf tempSimHistPdf("tempSimHistPdf","Simulation + Pedestal",RooArgList(decay,tempHistPdf),RooArgList(decayEvents,simEvents));

///////////////////////////////////
/////////Plotting/Fits/////////////
///////////////////////////////////

RooPlot* integralFrame = integral.frame(); //Use 100 bins for displaying source data

RooFitResult* tempRes = tempHistPdf.fitTo(sourceData,Save(1),Range(fitMin,fitMax));
auto pars = tempRes->floatParsFinal();

//Generate a new model based on the first pass.
Double_t slopeGuess=((RooRealVar*) pars.find("tempSlope"))->getVal();
RooRealVar slope("slope","slope",slopeGuess,0.999*slopeGuess,1.001*slopeGuess);
RooFormulaVar scaling("scaling","integral/slope",RooArgSet(integral,slope));
RooHistPdf histPdf("histPdf","histPdf",scaling,integral,*simHist,1);
RooAddPdf simHistPdf("simHistPdf","Simulation + Pedestal",RooArgList(decay,histPdf),RooArgList(decayEvents,simEvents));

//Fit our improved model. This helps avoid failed fits.
//simHistPdf.fitTo(sourceData,Save(1),Minimizer("Minuit2","simplex"),Range(fitMin,fitMax));
//simHistPdf.fitTo(sourceData,Save(1),Minimizer("Minuit2","migrad"),Range(fitMin,fitMax));
RooFitResult* res = histPdf.fitTo(sourceData,Save(1),Verbose(3),Range(fitMin,fitMax),PrintEvalErrors(10));

//Plot Model
sourceData.plotOn(integralFrame,Name("sourceData"),MarkerColor(1),FillColor(0),Binning(244));
//smearedData.plotOn(integralFrame,Name("smearedData"),MarkerColor(2),FillColor(0),Binning(244));
histPdf.plotOn(integralFrame,Name("sim"),Components(histPdf),LineColor(38),FillColor(33),LineWidth(3));
integralFrame->GetXaxis()->SetTitle("Integral [Arb.]");
integralFrame->GetYaxis()->SetTitle("Counts / 10 Bins");
integralFrame->GetXaxis()->CenterTitle();
integralFrame->GetYaxis()->CenterTitle();
integralFrame->SetTitle("");


//Draw
integralFrame->Draw();
c1->Modified();
c1->Update();
TString imagePath = "~/Desktop/Analysis/CeBr3/Plots/";
imagePath.Append("BestFit");
imagePath.Append(".pdf");
c1->Print(imagePath);
//c1->WaitPrimitive();

//return the fit.
return *res;

}

void scanResolution() {

	//Select appropriate file name to save output
	TString fileName;
	fileName = "~/Desktop/Analysis/CeBr3/Calibrations/resolutionScan.root";

	//Create file
	TFile* f = new TFile(fileName,"RECREATE");

	//Set number of alphas to run fits for, minimum alpha value, alpha step value 
	const Int_t numAlphas=30;
	const Int_t numBetas=30;
	const Int_t numFits=numAlphas*numBetas;
	Double_t alphaMin=3.0;
	Double_t alphaStep=0.1;
	Double_t alpha; //stores present alpha value
	Double_t betaMin=0.00;
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
			RooFitResult res = getNll(alpha,beta,alphaNum + betaNum);	

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
				slopes[ numGoodFits ] = ((RooRealVar*) pars.find("slope"))->getVal();
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
    //nllTree->SetDirectory(f);
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






























    
