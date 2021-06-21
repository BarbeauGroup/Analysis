/*

Simple code to plot quenching factor data for CeBr3

C. Awe - 06/21/2021
                               
*/

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
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooAddition.h"
#include "RooMinimizer.h"
#include "RooNumConvPdf.h"
#include "RooFitResult.h"
#include "TLegend.h"
#include "TRandom.h"
#include <iostream>
#include <limits>

using namespace RooFit; 

static const numPoints = 8; //Number of points in our data set.

void makePlots(){
	
	//Declare arrays to hold our data.
	double eVals[numPoints]; //Recoil energy.
	double eErrs[numPoints]; 
	double qfVals[numPoints]; //Quenching factor.
	double qfErrs[numPoints]; 
	double gVals[numPoints]; //Gamma values.
	double gErrs[numPoints];
	double bVals[numPoints]; //Beta values.
	double bErrs[numPoints];
	double sigVals[numPoints]; //Number of signal events.
	double sigErrs[numPoints];
	double noiseVals[numPoints]; //Number of noise/pedestal events.
	double noiseErrs[numPoints];
	double mean = 0.; //Mean effective recoil energy, recalculated at each step.
	double dMean = 0.; //Uncertainty.
	
	//Fill our data set.
	eVals[0] = 2.1;
	eErrs[0] = 0.;
	gVals[0] = 9.20;
	gErrs[0] = 0.42;
	bVals[0] = 0.04;
	bErrs[0] = 0.01;
	mean = bVals[0] * gVals[0];
	dMean = sqrt( pow(bVals[0]*gErrs[0],2) + pow(bErrs[0]*gVals[0],2) );
	qfVals[0] = mean / eVals[0] * 100.;
	qfErrs[0] = dMean / eVals[0] * 100.;
	sigVals[0] = 2945;
	sigErrs[0] = 78;
	noiseVals[0] = 62326;
	noiseErrs[0] = 256;
	
	eVals[1] = 8.9;
	eErrs[1] = 0.;
	gVals[1] = 5.25;
	gErrs[1] = 0.15;
	bVals[1] = 0.10;
	bErrs[1] = 0.01;
	mean = bVals[1] * gVals[1];
	dMean = sqrt( pow(bVals[1]*gErrs[1],2) + pow(bErrs[1]*gVals[1],2) );
	qfVals[1] = mean / eVals[1] * 100.;
	qfErrs[1] = dMean / eVals[1] * 100.;
	sigVals[1] = 4592;
	sigErrs[1] = 88;
	noiseVals[1] = 42742;
	noiseErrs[1] = 214;
	
	eVals[2] = 18.0;
	eErrs[2] = 0.;
	gVals[2] = 4.56;
	gErrs[2] = 0.23;
	bVals[2] = 0.18;
	bErrs[2] = 0.01;
	mean = bVals[2] * gVals[2];
	dMean = sqrt( pow(bVals[2]*gErrs[2],2) + pow(bErrs[2]*gVals[2],2) );
	qfVals[2] = mean / eVals[2] * 100.;
	qfErrs[2] = dMean / eVals[2] * 100.;
	sigVals[2] = 2122;
	sigErrs[2] = 59;
	noiseVals[2] = 26515;
	noiseErrs[2] = 167;
	
	eVals[3] = 25.9;
	eErrs[3] = 0.;
	gVals[3] = 4.42;
	gErrs[3] = 0.30;
	bVals[3] = 0.22;
	bErrs[3] = 0.01;
	mean = bVals[3] * gVals[3];
	dMean = sqrt( pow(bVals[3]*gErrs[3],2) + pow(bErrs[3]*gVals[3],2) );
	qfVals[3] = mean / eVals[3] * 100.;
	qfErrs[3] = dMean / eVals[3] * 100.;
	sigVals[3] = 1192;
	sigErrs[3] = 45;
	noiseVals[3] = 23800;
	noiseErrs[3] = 157;
	
	eVals[4] = 43.1;
	eErrs[4] = 0.;
	gVals[4] = 3.19;
	gErrs[4] = 0.41;
	bVals[4] = 0.36;
	bErrs[4] = 0.05;
	mean = bVals[4] * gVals[4];
	dMean = sqrt( pow(bVals[4]*gErrs[4],2) + pow(bErrs[4]*gVals[4],2) );
	qfVals[4] = mean / eVals[4] * 100.;
	qfErrs[4] = dMean / eVals[4] * 100.;
	sigVals[4] = 487;
	sigErrs[4] = 33;
	noiseVals[4] = 18034;
	noiseErrs[4] = 136;
	
	eVals[5] = 57.7;
	eErrs[5] = 0.;
	gVals[5] = 2.95;
	gErrs[5] = 0.29;
	bVals[5] = 0.71;
	bErrs[5] = 0.07;
	mean = bVals[5] * gVals[5];
	dMean = sqrt( pow(bVals[5]*gErrs[5],2) + pow(bErrs[5]*gVals[5],2) );
	qfVals[5] = mean / eVals[5] * 100.;
	qfErrs[5] = dMean / eVals[5] * 100.;
	sigVals[5] = 561;
	sigErrs[5] = 31;
	noiseVals[5] = 17369;
	noiseErrs[5] = 133;
	
	eVals[6] = 71.6;
	eErrs[6] = 0.;
	gVals[6] = 6.65;
	gErrs[6] = 0.55;
	bVals[6] = 0.44;
	bErrs[6] = 0.04;
	mean = bVals[6] * gVals[6];
	dMean = sqrt( pow(bVals[6]*gErrs[6],2) + pow(bErrs[6]*gVals[6],2) );
	qfVals[6] = mean / eVals[6] * 100.;
	qfErrs[6] = dMean / eVals[6] * 100.;
	sigVals[6] = 672;
	sigErrs[6] = 29;
	noiseVals[6] = 17602;
	noiseErrs[6] = 133;
	
	eVals[7] = 86.7;
	eErrs[7] = 0.;
	gVals[7] = 9.21;
	gErrs[7] = 0.66;
	bVals[7] = 0.39;
	bErrs[7] = 0.03;
	mean = bVals[7] * gVals[7];
	dMean = sqrt( pow(bVals[7]*gErrs[7],2) + pow(bErrs[7]*gVals[7],2) );
	qfVals[7] = mean / eVals[7] * 100.;
	qfErrs[7] = dMean / eVals[7] * 100.;
	sigVals[7] = 848;
	sigErrs[7] = 31;
	noiseVals[7] = 21244;
	noiseErrs[7] = 146;
	
	//Make our TGraphs.
	TGraph* qf_graph = new TGraphErrors(numPoints, eVals, qfVals, eErrs, qfErrs);
	TGraph* gamma_graph = new TGraphErrors(numPoints, eVals, gVals, eErrs, gErrs);
	TGraph* beta_graph = new TGraphErrors(numPoints, eVals, bVals, eErrs, bErrs);
	TGraph* signal_graph = new TGraphErrors(numPoints, eVals, sigVals, eErrs, sigErrs);
	TGraph* noise_graph = new TGraphErrors(numPoints, eVals, noiseVals, eErrs, noiseErrs);
	qf_graph->SetMarkerColor(45);
	qf_graph->SetMarkerStyle(21);
	qf_graph->SetTitle("Nuclear Quenching Factors in CeBr3");
	qf_graph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	qf_graph->GetXaxis()->CenterTitle();
	qf_graph->GetYaxis()->SetTitle("Nuclear Quenching Factor (%)");
	qf_graph->GetYaxis()->CenterTitle();
	gamma_graph->SetMarkerColor(9);
	gamma_graph->SetMarkerStyle(21);
	gamma_graph->SetTitle("Gamma vs. Recoil Energy");
	gamma_graph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	gamma_graph->GetXaxis()->CenterTitle();
	gamma_graph->GetYaxis()->SetTitle("Gamma");
	gamma_graph->GetYaxis()->CenterTitle();
	beta_graph->SetMarkerColor(41);
	beta_graph->SetMarkerStyle(21);
	beta_graph->SetTitle("Beta vs. Recoil Energy");
	beta_graph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	beta_graph->GetXaxis()->CenterTitle();
	beta_graph->GetYaxis()->SetTitle("Beta");
	beta_graph->GetYaxis()->CenterTitle();
	signal_graph->SetMarkerStyle(21);
	signal_graph->SetMarkerColor(29);
	signal_graph->SetTitle("Signal Counts vs. Recoil Energy");
	signal_graph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	signal_graph->GetXaxis()->CenterTitle();
	signal_graph->GetYaxis()->SetTitle("Signal Counts");
	signal_graph->GetYaxis()->CenterTitle();
	noise_graph->SetMarkerStyle(13);
	noise_graph->SetMarkerColor(46);
	noise_graph->SetTitle("Background Counts vs. Recoil Energy");
	noise_graph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	noise_graph->GetXaxis()->CenterTitle();
	noise_graph->GetYaxis()->SetTitle("Background Counts");
	noise_graph->GetYaxis()->CenterTitle();
	
	//Make a canvas and plot.
	TCanvas* c1 = new TCanvas("c1","c1");
	c1->Divide(2,2);
	c1->cd(1);
	gamma_graph->Draw("AP");
	c1->Modified();
	c1->Update();
	c1->cd(2);
	beta_graph->Draw("AP");
	c1->Modified();
	c1->Update();
	c1->cd(3);
	signal_graph->Draw("AP");
	c1->Modified();
	c1->Update();
	c1->cd(4);
	noise_graph->Draw("AP");
	c1->Modified();
	c1->Update();
	
	TCanvas* c2 = new TCanvas("c2","Quenching Factor in CeBr3");
	c2->cd();
	qf_graph->Draw("AP");
	c2->Modified();
	c2->Update();
	
	TString imagePath = "/var/phy/project/phil/cma46/CeBr3/cebr3-qf-analysis/Plots/";
	c1->Print(imagePath + "Parameter_Plots.png");
	c2->Print(imagePath + "QF_Plot.png");
	
}
	
	
	
	
	
	
	
