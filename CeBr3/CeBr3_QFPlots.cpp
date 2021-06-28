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

static const int numPoints = 7; //Number of points in our data set. 

void makePlots(){
	
	//Declare arrays to hold our data.
	double eVals[numPoints]; //Recoil energy.
	double eHighErrs[numPoints]; 
	double eLowErrs[numPoints]; 
	double qfVals[numPoints]; //Quenching factor.
	double qfHighErrs[numPoints]; 
	double qfLowErrs[numPoints]; 
	double gVals[numPoints]; //Gamma values.
	double gHighErrs[numPoints];
	double gLowErrs[numPoints];
	double bVals[numPoints]; //Beta values.
	double bHighErrs[numPoints];
	double bLowErrs[numPoints];
	double sigVals[numPoints]; //Number of signal events.
	double sigHighErrs[numPoints];
	double sigLowErrs[numPoints];
	double noiseVals[numPoints]; //Number of noise/pedestal events.
	double noiseHighErrs[numPoints];
	double noiseLowErrs[numPoints];
	double mean = 0.; //Mean effective recoil energy, recalculated at each step.
	double dMeanHigh = 0.; //Uncertainty.
	double dMeanLow = 0.;
	
	//Fill our data set.
	eVals[0] = 86.7;
	eLowErrs[0] = 0.;
	eHighErrs[0] = 0.;
	gVals[0] = 9.12;
	gHighErrs[0] = 0.67;
	gLowErrs[0] = 0.63;
	bVals[0] = 0.40;
	bHighErrs[0] = 0.03;
	bLowErrs[0] = 0.03;
	mean = bVals[0] * gVals[0];
	dMeanHigh = ((bVals[0]+bHighErrs[0]) * (gVals[0]+gHighErrs[0])) - mean;
	dMeanLow = mean - ((bVals[0]-bLowErrs[0]) * (gVals[0]-gLowErrs[0]));
	qfVals[0] = mean / eVals[0] * 100.;
	qfHighErrs[0] = dMeanHigh / eVals[0] * 100.;
	qfLowErrs[0] = dMeanLow / eVals[0] * 100.;
	sigVals[0] = 838.02;
	sigHighErrs[0] = 30.60;
	sigLowErrs[0] = 29.65;
	noiseVals[0] = 21229.37;
	noiseHighErrs[0] = 146.17;
	noiseHighErrs[0] = 145.60;
	
	eVals[1] = 71.6;
	eLowErrs[1] = 0.;
	eHighErrs[1] = 0.;
	gVals[1] = 6.40;
	gHighErrs[1] = 0.53;
	gLowErrs[1] = 0.50;
	bVals[1] = 0.46;
	bHighErrs[1] = 0.04;
	bLowErrs[1] = 0.03;
	mean = bVals[1] * gVals[1];
	dMeanHigh = ((bVals[1]+bHighErrs[1]) * (gVals[1]+gHighErrs[1])) - mean;
	dMeanLow = mean - ((bVals[1]-bLowErrs[1]) * (gVals[1]-gLowErrs[1]));
	qfVals[1] = mean / eVals[1] * 100.;
	qfHighErrs[1] = dMeanHigh / eVals[1] * 100.;
	qfLowErrs[1] = dMeanLow / eVals[1] * 100.;
	sigVals[1] = 716.74;
	sigHighErrs[1] = 30.23;
	sigLowErrs[1] = 29.61;
	noiseVals[1] = 18580.99;
	noiseHighErrs[1] = 136.68;
	noiseHighErrs[1] = 135.13;
	
	eVals[2] = 57.7;
	eLowErrs[2] = 0.;
	eHighErrs[2] = 0.;
	gVals[2] = 2.37;
	gHighErrs[2] = 0.27;
	gLowErrs[2] = 0.20;
	bVals[2] = 0.88;
	bHighErrs[2] = 0.07;
	bLowErrs[2] = 0.08;
	mean = bVals[2] * gVals[2];
	dMeanHigh = ((bVals[2]+bHighErrs[2]) * (gVals[2]+gHighErrs[2])) - mean;
	dMeanLow = mean - ((bVals[2]-bLowErrs[2]) * (gVals[2]-gLowErrs[2]));
	qfVals[2] = mean / eVals[2] * 100.;
	qfHighErrs[2] = dMeanHigh / eVals[2] * 100.;
	qfLowErrs[2] = dMeanLow / eVals[2] * 100.;
	sigVals[2] = 661.20;
	sigHighErrs[2] = 35.62;
	sigLowErrs[2] = 34.53;
	noiseVals[2] = 19631.98;
	noiseHighErrs[2] = 142.25;
	noiseHighErrs[2] = 141.45;
	
	eVals[3] = 43.1;
	eLowErrs[3] = 0.;
	eHighErrs[3] = 0.;
	gVals[3] = 1.51;
	gHighErrs[3] = 0.18;
	gLowErrs[3] = 0.14;
	bVals[3] = 0.65;
	bHighErrs[3] = 0.06;
	bLowErrs[3] = 0.06;
	mean = bVals[3] * gVals[3];
	dMeanHigh = ((bVals[3]+bHighErrs[3]) * (gVals[3]+gHighErrs[3])) - mean;
	dMeanLow = mean - ((bVals[3]-bLowErrs[3]) * (gVals[3]-gLowErrs[3]));
	qfVals[3] = mean / eVals[3] * 100.;
	qfHighErrs[3] = dMeanHigh / eVals[3] * 100.;
	qfLowErrs[3] = dMeanLow / eVals[3] * 100.;
	sigVals[3] = 736.94;
	sigHighErrs[3] = 55.65;
	sigLowErrs[3] = 53.73;
	noiseVals[3] = 21316.77;
	noiseHighErrs[3] = 153.10;
	noiseHighErrs[3] = 147.78;
	
	eVals[4] = 25.9;
	eLowErrs[4] = 0.;
	eHighErrs[4] = 0.;
	gVals[4] = 3.35;
	gHighErrs[4] = 0.38;
	gLowErrs[4] = 0.44;
	bVals[4] = 0.28;
	bHighErrs[4] = 0.03;
	bLowErrs[4] = 0.03;
	mean = bVals[4] * gVals[4];
	dMeanHigh = ((bVals[4]+bHighErrs[4]) * (gVals[4]+gHighErrs[4])) - mean;
	dMeanLow = mean - ((bVals[4]-bLowErrs[4]) * (gVals[4]-gLowErrs[4]));
	qfVals[4] = mean / eVals[4] * 100.;
	qfHighErrs[4] = dMeanHigh / eVals[4] * 100.;
	qfLowErrs[4] = dMeanLow / eVals[4] * 100.;
	sigVals[4] = 1472.12;
	sigHighErrs[4] = 75.08;
	sigLowErrs[4] = 64.24;
	noiseVals[4] = 29143.40;
	noiseHighErrs[4] = 181.75;
	noiseHighErrs[4] = 181.65;
	
	eVals[5] = 18.0;
	eLowErrs[5] = 0.;
	eHighErrs[5] = 0.;
	gVals[5] = 1.82;
	gHighErrs[5] = 0.13;
	gLowErrs[5] = 0.10;
	bVals[5] = 0.36;
	bHighErrs[5] = 0.02;
	bLowErrs[5] = 0.02;
	mean = bVals[5] * gVals[5];
	dMeanHigh = ((bVals[5]+bHighErrs[5]) * (gVals[5]+gHighErrs[5])) - mean;
	dMeanLow = mean - ((bVals[5]-bLowErrs[5]) * (gVals[5]-gLowErrs[5]));
	qfVals[5] = mean / eVals[5] * 100.;
	qfHighErrs[5] = dMeanHigh / eVals[5] * 100.;
	qfLowErrs[5] = dMeanLow / eVals[5] * 100.;
	sigVals[5] = 3140.23;
	sigHighErrs[5] = 118.91;
	sigLowErrs[5] = 122.32;
	noiseVals[5] = 30710.24;
	noiseHighErrs[5] = 197.74;
	noiseHighErrs[5] = 197.64;
	
	eVals[6] = 8.9;
	eLowErrs[6] = 0.;
	eHighErrs[6] = 0.;
	gVals[6] = 1.66;
	gHighErrs[6] = 0.03;
	gLowErrs[6] = 0.03;
	bVals[6] = 0.21;
	bHighErrs[6] = 0.01;
	bLowErrs[6] = 0.01;
	mean = bVals[6] * gVals[6];
	dMeanHigh = ((bVals[6]+bHighErrs[6]) * (gVals[6]+gHighErrs[6])) - mean;
	dMeanLow = mean - ((bVals[6]-bLowErrs[6]) * (gVals[6]-gLowErrs[6]));
	qfVals[6] = mean / eVals[6] * 100.;
	qfHighErrs[6] = dMeanHigh / eVals[6] * 100.;
	qfLowErrs[6] = dMeanLow / eVals[6] * 100.;
	sigVals[6] = 9366.12;
	sigHighErrs[6] = 211.28;
	sigLowErrs[6] = 211.48;
	noiseVals[6] = 47396.07;
	noiseHighErrs[6] = 275.92;
	noiseHighErrs[6] = 282.52;

	//Outlier data point, this will be plotted on its own as hollow.
	double outlier_e[1];
	outlier_e[0] = 2.1;
	double out_e_HighErrs[1];
	out_e_HighErrs[0] = 0.;
	double out_e_LowErrs[1];
	out_e_LowErrs[0] = 0.;
	double outlier_g[1];
	outlier_g[0] = 11.32;
	double out_g_HighErr[1];
	out_g_HighErr[0] = 0.18;
	double out_g_LowErr[1];
	out_g_LowErr[0] = 0.18;
	double outlier_b[1];
	outlier_b[0] = 0.02;
	double out_b_HighErr[1];
	out_b_HighErr[0] = 0.01;
	double out_b_LowErr[1];
	out_b_LowErr[0] = 0.01;
	mean = outlier_b[0] * outlier_g[0];
	dMeanHigh = ((outlier_b[0]+out_b_HighErr[0]) * (outlier_g[0]+out_g_HighErr[0])) - mean;
	dMeanLow = mean - (outlier_b[0]-out_b_LowErr[0]) * (outlier_g[0]-out_g_LowErr[0]);
	double outlier_qfVal[1];
	outlier_qfVal[0] = mean / outlier_e[0] * 100.;
	double out_qfHighErr[1];
	out_qfHighErr[0] = dMeanHigh / outlier_e[0] * 100.;
	double out_qfLowErr[1];
	out_qfLowErr[0] = dMeanLow / outlier_e[0] * 100.;
	double outlier_sigVal[1];
	outlier_sigVal[0] = 4167.65;
	double out_sigHighErr[1];
	out_sigHighErr[0] = 89.48;
	double out_sigLowErr[1];
	out_sigLowErr[0] = 88.73;
	double outlier_noiseVal[1];
	outlier_noiseVal[0] = 119742.74;
	double out_noiseHighErr[1];
	out_noiseHighErr[0] = 13490.56;
	double out_noiseLowErr[1];
	out_noiseLowErr[0] = 13214.49;
	
	//Make our TGraphs.
	TGraph* qf_graph = new TGraphAsymmErrors(numPoints, eVals, qfVals, eLowErrs, eHighErrs, qfLowErrs, qfHighErrs);
	TGraph* outlier = new TGraphAsymmErrors(1, outlier_e, outlier_qfVal, out_e_LowErrs, out_e_HighErrs, out_qfLowErr, out_qfHighErr);
	TMultiGraph* qf_MultiGraph = new TMultiGraph();
	qf_MultiGraph->Add(qf_graph,"AP");
	qf_MultiGraph->Add(outlier,"AP");
	TGraph* gamma_graph = new TGraphAsymmErrors(numPoints, eVals, eLowErrs, eHighErrs, gLowErrs, gHighErrs);
	TGraph* g_outlier = new TGraphAsymmErrors(1, outlier_e, outlier_g, out_e_LowErrs, out_e_HighErrs, out_g_LowErr, out_g_HighErr);
	TMultiGraph* g_MultiGraph = new TMultiGraph();
	g_MultiGraph->Add(gamma_graph,"AP");
	g_MultiGraph->Add(g_outlier,"AP");
	TGraph* beta_graph = new TGraphAsymmErrors(numPoints, eVals, bVals, eLowErrs, eHighErrs, bLowErrs, bHighErrs);
	TGraph* b_outlier = new TGraphAsymmErrors(1, outlier_e, outlier_b, out_e_LowErrs, out_e_HighErrs, out_b_LowErr, out_b_HighErr);
	TMultiGraph* b_MultiGraph = new TMultiGraph();
	b_MultiGraph->Add(beta_graph,"AP");
	b_MultiGraph->Add(b_outlier,"AP");
	TGraph* signal_graph = new TGraphAsymmErrors(numPoints, eVals, sigVals, eLowErrs, eHighErrs, sigLowErrs, sigHighErrs);
	TGraph* sig_outlier = new TGraphAsymmErrors(1, outlier_e, outlier_sigVal, out_e_LowErrs, out_e_HighErrs, out_sigLowErr, out_sigHighErr);
	TMultiGraph* sig_MultiGraph = new TMultiGraph();
	sig_MultiGraph->Add(signal_graph,"AP");
	sig_MultiGraph->Add(sig_outlier,"AP");
	TGraph* noise_graph = new TGraphAsymmErrors(numPoints, eVals, noiseVals, eLowErrs, eHighErrs, noiseLowErrs, noiseHighErrs);
	TGraph* noise_outlier = new TGraphAsymmErrors(1, outlier_e, outlier_noiseVal, out_e_LowErrs, out_e_HighErrs, out_noiseLowErr, out_noiseHighErr);
	TMultiGraph* noise_MultiGraph = new TMultiGraph();
	noise_MultiGraph->Add(noise_graph,"AP");
	noise_MultiGraph->Add(noise_outlier,"AP");
	qf_graph->SetMarkerColor(45);
	qf_graph->SetMarkerStyle(21);
	qf_MultiGraph->SetTitle("Nuclear Quenching Factors in CeBr3");
	qf_MultiGraph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	qf_MultiGraph->GetXaxis()->CenterTitle();
	qf_MultiGraph->GetYaxis()->SetTitle("Nuclear Quenching Factor (%)");
	qf_MultiGraph->GetYaxis()->CenterTitle();
	outlier->SetMarkerColor(45);
	outlier->SetMarkerStyle(24);
	gamma_graph->SetMarkerColor(9);
	gamma_graph->SetMarkerStyle(21);
	g_MultiGraph->SetTitle("Gamma vs. Recoil Energy");
	g_MultiGraph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	g_MultiGraph->GetXaxis()->CenterTitle();
	g_MultiGraph->GetYaxis()->SetTitle("Gamma");
	g_MultiGraph->GetYaxis()->CenterTitle();
	g_outlier->SetMarkerColor(9);
	g_outlier->SetMarkerStyle(24);
	beta_graph->SetMarkerColor(41);
	beta_graph->SetMarkerStyle(21);
	b_MultiGraph->SetTitle("Beta vs. Recoil Energy");
	b_MultiGraph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	b_MultiGraph->GetXaxis()->CenterTitle();
	b_MultiGraph->GetYaxis()->SetTitle("Beta");
	b_MultiGraph->GetYaxis()->CenterTitle();
	b_outlier->SetMarkerColor(41);
	b_outlier->SetMarkerStyle(24);
	signal_graph->SetMarkerStyle(21);
	signal_graph->SetMarkerColor(29);
	sig_MultiGraph->SetTitle("Signal Counts vs. Recoil Energy");
	sig_MultiGraph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	sig_MultiGraph->GetXaxis()->CenterTitle();
	sig_MultiGraph->GetYaxis()->SetTitle("Signal Counts");
	sig_MultiGraph->GetYaxis()->CenterTitle();
	sig_outlier->SetMarkerColor(29);
	sig_outlier->SetMarkerStyle(24);
	noise_graph->SetMarkerStyle(13);
	noise_graph->SetMarkerColor(46);
	noise_MultiGraph->SetTitle("Background Counts vs. Recoil Energy");
	noise_MultiGraph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	noise_MultiGraph->GetXaxis()->CenterTitle();
	noise_MultiGraph->GetYaxis()->SetTitle("Background Counts");
	noise_MultiGraph->GetYaxis()->CenterTitle();
	noise_outlier->SetMarkerColor(46);
	noise_outlier->SetMarkerStyle(24);
	
	//Make a canvas and plot.
	TCanvas* c1 = new TCanvas("c1","c1");
	c1->Divide(2,2);
	c1->cd(1);
	g_MultiGraph->Draw("a");
	//g_outlier->Draw("SAME");
	c1->Modified();
	c1->Update();
	c1->cd(2);
	b_MultiGraph->Draw("a");
	//b_outlier->Draw("SAME");
	c1->Modified();
	c1->Update();
	c1->cd(3);
	sig_MultiGraph->Draw("a");
	//sig_outlier->Draw("SAME");
	c1->Modified();
	c1->Update();
	c1->cd(4);
	noise_MultiGraph->Draw("a");
	//noise_outlier->Draw("SAME");
	c1->Modified();
	c1->Update();
	
	TCanvas* c2 = new TCanvas("c2","Quenching Factor in CeBr3");
	c2->cd();
	qf_MultiGraph->Draw("a");
	//outlier->Draw("SAME");
	c2->Modified();
	c2->Update();
	
	TString imagePath = "/var/phy/project/phil/cma46/CeBr3/cebr3-qf-analysis/Plots/";
	c1->Print(imagePath + "Parameter_Plots.png");
	c2->Print(imagePath + "QF_Plot.png");
	
}
	
	
	
	
	
	
	
