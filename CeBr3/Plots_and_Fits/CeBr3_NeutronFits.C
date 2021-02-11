//Code to fit MCNP Simulations to CeBr3 Calibration Data
//Written by C. Awe (based heavily on work by S. Hedges)
//02/09/2021

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
RooFitResult getNll(Double_t alpha=0.0000084, Double_t beta=0.00000004, Int_t fitNum=0, Int_t bdNum=8 ) {

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
Long64_t keVMax = 10;

//Simulation file paths.
TString bd1_sim = "~/Desktop/cebr3-qf-analysis/analysis/Neutron_Sims/dumn1-cebr3-neutrons-bd1.root";
TString bd2_sim = "~/Desktop/cebr3-qf-analysis/analysis/Neutron_Sims/dumn1-cebr3-neutrons-bd2.root";
TString bd3_sim = "~/Desktop/cebr3-qf-analysis/analysis/Neutron_Sims/dumn1-cebr3-neutrons-bd3.root";
TString bd4_sim = "~/Desktop/cebr3-qf-analysis/analysis/Neutron_Sims/dumn1-cebr3-neutrons-bd4.root";
TString bd5_sim = "~/Desktop/cebr3-qf-analysis/analysis/Neutron_Sims/dumn1-cebr3-neutrons-bd5.root";
TString bd6_sim = "~/Desktop/cebr3-qf-analysis/analysis/Neutron_Sims/dumn1-cebr3-neutrons-bd6.root";
TString bd7_sim = "~/Desktop/cebr3-qf-analysis/analysis/Neutron_Sims/dumn1-cebr3-neutrons-bd7.root";
TString bd8_sim = "~/Desktop/cebr3-qf-analysis/analysis/Neutron_Sims/dumn1-cebr3-neutrons-bd8.root";
TString sim_Filelist[] = { bd1_sim, bd2_sim, bd3_sim, bd4_sim, bd5_sim, bd6_sim, bd7_sim, bd8_sim };

//Declare our observable.
TCanvas* fitCanvas = new TCanvas("fitCanvas","fitCanvas");
RooRealVar keV("keV","keV",keVMin,keVMax);
RooArgSet observables("observables");
observables.add(keV);

// Set #bins to be used for FFT sampling to 10000
//keV.setBins(100,"cache"); 

//Load the proper simulation.
TString sim_Filename = sim_Filelist[ bdNum - 1 ];
TFile simFile( sim_Filename );

/////////////////////////////////
//////Load simulation ttree//////
/////////////////////////////////

TTree *simTree = (TTree*)simFile.Get("totalEnergyTree");

//Load energy branch
Double_t energy;
simTree->SetBranchAddress("energy",&energy);

//Get number of entries
Long64_t numSimEntries=simTree->GetEntries();

//Create gaussian to generate random data
RooRealVar mean("mean","mean",keVMin,keVMax);
RooRealVar sigma("sigma","sigma",keVMin,keVMax); //sigma probably won't ever be greater than fitMax.
RooGaussian gauss("gauss","gauss",keV,mean,sigma);

//Create data set to hold smeared simulation data.
RooDataSet smearedData("smearedData","smearedData",keV);

Int_t valsToGen=10; //create 10 gaussian-distributed data points based on alpha for each point in simulation
default_random_engine generator;

/////////////////////////////////////
/////////MAKING SMEARED DATA/////////
/////////////////////////////////////

cout<<"Generating smeared data..."<<endl;
Double_t quenchingGuess = 0.02; //Inital guess for the quenching factor.
for (Long64_t entry=0; entry < numSimEntries; entry++) {

	//Get entry
	simTree->GetEntry(entry);
	Double_t enrg = energy;

	//Generate valsToGen random gaussian values based on this distribution, add to smeared data set
	for (Int_t val=0; val < valsToGen; val++) {
		Double_t sample = genRandomData( alpha, beta, enrg );
		sample = sample * quenchingGuess; //Quench it.
		keV = sample;
		smearedData.add( keV );

	}

	//Give user information on screen
	if (entry % 10000 == 0) {
		printf("Processing entry %lu of %lu entry\n", entry, numSimEntries);
	}
}

cout << "Generated!" << endl;

smearedData.Print("V");

/////////////////////////////////////
/////////LOADING SOURCE DATA/////////
/////////////////////////////////////

//Create empty data set
RooDataSet sourceData("sourceData","sourceData",keV);
RooDataSet noise("noise","noise",keV);	
Long64_t numSrcDataInRange=0;
Double_t adc_to_keV = 9.65e-4;

//Cuts
Double_t time_Low = 306;
Double_t time_High = 336;
Double_t psd_Low = 0.26;
Double_t psd_High = 0.6;
Double_t integral_Low = 10000;
Double_t integral_High = 35000;

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
	if(LS_channel == bdNum){
		if( ( LS_psd > psd_Low ) && ( LS_psd < psd_High ) ){
			if( ( LS_integral > integral_Low ) && ( LS_integral < integral_High ) ){
				if( ( LS_timeToBPM > time_Low ) && ( LS_timeToBPM < time_High ) ){
					Double_t enrg = scatterer_integral * adc_to_keV;
					if( enrg < keVMax && enrg > keVMin ){
						keV = enrg;
						sourceData.add( keV );
					}
					enrg = scatterer_noise * adc_to_keV;
					if( enrg < keVMax && enrg > keVMin ){
						keV = enrg;
						noise.add( keV );
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

cout << "Printing source data info." << endl;
sourceData.Print("V");

noise.Print("V");

cout << "Generating RooPDFs..." << endl;

////////////////////////////////
/////////Make Hist PDFs/////////
////////////////////////////////

//Set binnings
((RooRealVar*)smearedData.get()->find("keV"))->setBins(keVMax-keVMin);

//Represent simulation data as a PDF.
RooDataHist* simHist = smearedData.binnedClone();

//Make a noise hist pdf.
((RooRealVar*)noise.get()->find("keV"))->setBins(keVMax-keVMin);
RooDataHist* noiseHist = noise.binnedClone();
RooHistPdf noiseHistPdf("noiseHistPdf","Noise Histogram",keV,*noiseHist,0);

//Create scaling factors.
RooRealVar tempQF("tempQF","tempQF",1,2,0.1); //Converts keVnr to keVee.
RooFormulaVar tempScaling("tempScaling","keV*tempQF",RooArgSet(keV,tempQF));
RooHistPdf tempSimHist("tempSimHist","tempSimHist",tempScaling,keV,*simHist,1);
RooRealVar noiseEvents("noiseEvents","noiseEvents",2e6,1,1e7);
RooRealVar recoilEvents("recoilEvents","recoilEvents",30e5,1,1e7);
RooAddPdf tempPdf("tempPdf","Simulation + Noise",RooArgList(noiseHistPdf,tempSimHist),RooArgList(noiseEvents,recoilEvents));

///////////////////////////////////
/////////Plotting/Fits/////////////
///////////////////////////////////

//Create canvas
TCanvas* c1 = new TCanvas("c1","c1");
//c1->Divide(1,2);
c1->cd();
c1->SetLogy();
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

RooFitResult* tempRes = tempPdf.fitTo(sourceData,Save(1),Range(fitMin,fitMax));
auto pars = tempRes->floatParsFinal();

//Generate a new model based on the first pass.
Double_t qfGuess=((RooRealVar*) pars.find("tempQF"))->getVal();
RooRealVar qf("qf","qf",qfGuess,0.999*qfGuess,1.001*qfGuess);
RooFormulaVar scaling("scaling","keV*qf",RooArgSet(keV,qf));
RooHistPdf simHistPdf("simHistPdf","simHistPdf",scaling,keV,*simHist,1);
RooAddPdf  totalPdf("totalPdf","Simulation + Noise",RooArgList(noiseHistPdf,simHistPdf),RooArgList(noiseEvents,recoilEvents));

//Fit our improved model. This helps avoid failed fits.
//simHistPdf.fitTo(sourceData,Save(1),Minimizer("Minuit2","simplex"),Range(fitMin,fitMax));
//simHistPdf.fitTo(sourceData,Save(1),Minimizer("Minuit2","migrad"),Range(fitMin,fitMax));
RooFitResult* res = totalPdf.fitTo(sourceData,Save(1),Verbose(3),Range(fitMin,fitMax),PrintEvalErrors(10));

//Plot Model
sourceData.plotOn(keVFrame,Name("sourceData"),MarkerColor(1),FillColor(0),Binning(100),Range(fitMin,fitMax));
totalPdf.plotOn(keVFrame,Name("model"),LineColor(8),FillColor(0),LineWidth(3),Range(fitMin,fitMax));
//RooHist* resHist = keVFrame->residHist();
totalPdf.plotOn(keVFrame,Name("sim"),Components(simHistPdf),LineColor(4),FillColor(0),LineWidth(3),Range(fitMin,fitMax));
totalPdf.plotOn(keVFrame,Name("noise"),Components(noiseHistPdf),LineColor(2),FillColor(46),LineWidth(3),Range(fitMin,fitMax));
keVFrame->GetXaxis()->SetTitle("Energy [keVee]");
keVFrame->GetYaxis()->SetTitle("Counts / keV");
keVFrame->GetXaxis()->CenterTitle();
keVFrame->GetYaxis()->CenterTitle();
keVFrame->SetTitle("");
keVFrame->GetYaxis()->SetRangeUser(0.1,1000);

//Draw
keVFrame->Draw();

//Draw residuals.
/*
c1->cd(2);
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
resFrame->GetYaxis()->SetTitle("Counts / keV");
resFrame->addPlotable(resHist,"p");
resFrame->Draw(); 
*/
c1->Modified();
c1->Update();

TString imagePath = "~/Desktop/cebr3-qf-analysis/Plots/";
imagePath.Append("BestFit");
imagePath.Append(".pdf");
c1->Print(imagePath);
//c1->WaitPrimitive();

//return the fit.
return *res;

}

void scanResolution( Int_t bdNum ) {

	//Select appropriate file name to save output
	TString fileName;
	TString fname1 = "~/Desktop/cebr3-qf-analysis/Plots/resolutionScan_bd1.root";
	TString fname2 = "~/Desktop/cebr3-qf-analysis/Plots/resolutionScan_bd2.root";
	TString fname3 = "~/Desktop/cebr3-qf-analysis/Plots/resolutionScan_bd3.root";
	TString fname4 = "~/Desktop/cebr3-qf-analysis/Plots/resolutionScan_bd4.root";
	TString fname5 = "~/Desktop/cebr3-qf-analysis/Plots/resolutionScan_bd5.root";
	TString fname6 = "~/Desktop/cebr3-qf-analysis/Plots/resolutionScan_bd6.root";
	TString fname7 = "~/Desktop/cebr3-qf-analysis/Plots/resolutionScan_bd7.root";
	TString fname8 = "~/Desktop/cebr3-qf-analysis/Plots/resolutionScan_bd8.root";
	TString fname_List[] = { fname1, fname2, fname3, fname4, fname5, fname6, fname7, fname8 };
	fileName = fname_List[bdNum-1];

	//Create file
	TFile f(fileName,"RECREATE");

	//Set number of alphas to run fits for, minimum alpha value, alpha step value 
	const Int_t numAlphas=20;
	const Int_t numBetas=20;
	const Int_t numFits=numAlphas*numBetas;
	Double_t alphaMin=2.0;
	Double_t alphaStep=0.5;
	Double_t alpha; //stores present alpha value
	Double_t betaMin=0.005;
	Double_t betaStep=0.001;
	Double_t beta; //stores present alpha value

	//Dummies for generating TTree.
	Double_t nll;
	Double_t qf;
	Double_t status;

	//Arrays for holding data until we normalize the NLL.
	Double_t alphas[numFits]={};
	Double_t betas[numFits]={};
	Double_t nlls[numFits]={};
	Double_t qfs[numFits]={};

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
				//Store qfs from successful fits.
				qfs[ numGoodFits ] = ((RooRealVar*) pars.find("qf"))->getVal();
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
	nllTree->Branch("qf",&qf,"qf/D");
    	cout << "TTree generated!" << endl;
    	nllTree->SetDirectory(0);
	TH2D* nllHist = new TH2D("nllHist","nllHist",numAlphas,alphaMin,alphaMin + numAlphas*alphaStep, numBetas, betaMin, betaMin + numBetas*betaStep);

	//Subtract minNll value from all nlls so lowest one has a value of 0 and populate the Tree.
	for (Int_t i=0; i < numGoodFits; i++) {
		nlls[i] -= minNll;
		nll = nlls[i];
		alpha = alphas[i];
		beta = betas[i];
		qf = qfs[i];
		nllHist->SetBinContent( nllHist->FindBin( alpha, beta ), nll );
		nllTree->Fill();
		if( nll == 0 ){
			cout << "The minimizing alpha is: " << alpha << endl;
			cout << "The minimizing beta is: " << beta << endl;
		}
		
	}

	//Write TTree and Histogram.
	nllTree->Write();
	nllHist->SetMaximum(5);
	nllHist->Write();

	//Close file and TTree.
	//nllTree->Delete();
	f.Close();
}






























    

