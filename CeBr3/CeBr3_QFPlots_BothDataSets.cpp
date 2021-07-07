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

	///////////////////////////////
	/// Production Run 1 Values ///
	///////////////////////////////
	
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
	eVals[0] = 76.20;
	eLowErrs[0] = 2.70;
	eHighErrs[0] = 2.70;
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
	noiseLowErrs[0] = 145.60;
	//cout << "QF Vals are " << qfVals[0] << endl;
	//cout << " + " << qfHighErrs[0] << endl;
	//cout << " - " << qfLowErrs[0] << endl;
	
	eVals[1] = 60.33;
	eLowErrs[1] = 2.62;
	eHighErrs[1] = 2.62;
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
	noiseLowErrs[1] = 135.13;
	//cout << "QF Vals are " << qfVals[1] << endl;
	//cout << " + " << qfHighErrs[1] << endl;
	//cout << " - " << qfLowErrs[1] << endl;
	
	eVals[2] = 43.59;
	eLowErrs[2] = 1.91;
	eHighErrs[2] = 1.91;
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
	noiseLowErrs[2] = 141.45;
	//cout << "QF Vals are " << qfVals[2] << endl;
	//cout << " + " << qfHighErrs[2] << endl;
	//cout << " - " << qfLowErrs[2] << endl;
	
	eVals[3] = 27.46;
	eLowErrs[3] = 1.67;
	eHighErrs[3] = 1.67;
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
	noiseLowErrs[3] = 147.78;
	//cout << "QF Vals are " << qfVals[3] << endl;
	//cout << " + " << qfHighErrs[3] << endl;
	//cout << " - " << qfLowErrs[3] << endl;
	
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
	noiseLowErrs[4] = 181.65;
	//cout << "QF Vals are " << qfVals[4] << endl;
	//cout << " + " << qfHighErrs[4] << endl;
	//cout << " - " << qfLowErrs[4] << endl;
	
	eVals[5] = 12.26;
	eLowErrs[5] = 0.97;
	eHighErrs[5] = 0.97;
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
	noiseLowErrs[5] = 197.64;
	//cout << "QF Vals are" << qfVals[5] << endl;
	//cout << " + " << qfHighErrs[5] << endl;
	//cout << " - " << qfLowErrs[5] << endl;
	
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
	noiseLowErrs[6] = 282.52;
	//cout << "QF Vals are " << qfVals[6] << endl;
	//cout << " + " << qfHighErrs[6] << endl;
	//cout << " - " << qfLowErrs[6] << endl;

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
	//cout << outlier_qfVal[0] << endl;
	//cout << out_qfHighErr[0] << endl;
	//cout << out_qfLowErr[0] << endl;

	///////////////////////////////
	/// Production Run 2 Values ///
	///////////////////////////////
	
	//Declare arrays to hold our data.
	double eVals_PR2[4];
	double eHighErrs_PR2[4]; 
	double eLowErrs_PR2[4]; 
	double qfVals_PR2[4]; //Quenching factor.
	double qfHighErrs_PR2[4]; 
	double qfLowErrs_PR2[4]; 
	double gVals_PR2[4]; //Gamma values.
	double gHighErrs_PR2[4];
	double gLowErrs_PR2[4];
	double bVals_PR2[4]; //Beta values.
	double bHighErrs_PR2[4];
	double bLowErrs_PR2[4];
	double sigVals_PR2[4]; //Number of signal events.
	double sigHighErrs_PR2[4];
	double sigLowErrs_PR2[4];
	double noiseVals_PR2[4]; //Number of noise/pedestal events.
	double noiseHighErrs_PR2[4];
	double noiseLowErrs_PR2[4];
	
	//Fill our data set.
	eVals_PR2[0] = eVals[0];
	eLowErrs_PR2[0] = eLowErrs[0];
	eHighErrs_PR2[0] = eHighErrs[0];
	gVals_PR2[0] = 7.95;
	gHighErrs_PR2[0] = 1.03;
	gLowErrs_PR2[0] = 1.09;
	bVals_PR2[0] = 0.52;
	bHighErrs_PR2[0] = 0.07;
	bLowErrs_PR2[0] = 0.05;
	mean = bVals_PR2[0] * gVals_PR2[0];
	dMeanHigh = ((bVals_PR2[0]+bHighErrs_PR2[0]) * (gVals_PR2[0]+gHighErrs_PR2[0])) - mean;
	dMeanLow = mean - ((bVals_PR2[0]-bLowErrs_PR2[0]) * (gVals_PR2[0]-gLowErrs_PR2[0]));
	qfVals_PR2[0] = mean / eVals[0] * 100.;
	qfHighErrs_PR2[0] = dMeanHigh / eVals[0] * 100.;
	qfLowErrs_PR2[0] = dMeanLow / eVals[0] * 100.;
	sigVals_PR2[0] = 1046.24;
	sigHighErrs_PR2[0] = 86.81;
	sigLowErrs_PR2[0] = 69.86;
	noiseVals_PR2[0] = 29413.22;
	noiseHighErrs_PR2[0] = 184.70;
	noiseLowErrs_PR2[0] = 196.29;
	//cout << "QF Vals are " << qfVals_PR2[0] << endl;
	//cout << " + " << qfHighErrs_PR2[0] << endl;
	//cout << " - " << qfLowErrs_PR2[0] << endl;
	
	eVals_PR2[1] = eVals[1];
	eLowErrs_PR2[1] = eLowErrs[1];
	eHighErrs_PR2[1] = eHighErrs[1];
	gVals_PR2[1] = 5.93;
	gHighErrs_PR2[1] = 0.83;
	gLowErrs_PR2[1] = 0.75;
	bVals_PR2[1] = 0.59;
	bHighErrs_PR2[1] = 0.07;
	bLowErrs_PR2[1] = 0.06;
	mean = bVals_PR2[1] * gVals_PR2[1];
	dMeanHigh = ((bVals_PR2[1]+bHighErrs_PR2[1]) * (gVals_PR2[1]+gHighErrs_PR2[1])) - mean;
	dMeanLow = mean - ((bVals_PR2[1]-bLowErrs_PR2[1]) * (gVals_PR2[1]-gLowErrs_PR2[1]));
	qfVals_PR2[1] = mean / eVals[1] * 100.;
	qfHighErrs_PR2[1] = dMeanHigh / eVals[1] * 100.;
	qfLowErrs_PR2[1] = dMeanLow / eVals[1] * 100.;
	sigVals_PR2[1] = 1147.22;
	sigHighErrs_PR2[1] = 88.32;
	sigLowErrs_PR2[1] = 82.51;
	noiseVals_PR2[1] = 25576.79;
	noiseHighErrs_PR2[1] = 175.50;
	noiseLowErrs_PR2[1] = 176.58;
	//cout << "QF Vals are " << qfVals_PR2[1] << endl;
	//cout << " + " << qfHighErrs_PR2[1] << endl;
	//cout << " - " << qfLowErrs_PR2[1] << endl;
	
	eVals_PR2[2] = eVals[2];
	eLowErrs_PR2[2] = eLowErrs[2];
	eHighErrs_PR2[2] = eHighErrs[2];
	gVals_PR2[2] = 1.27;
	gHighErrs_PR2[2] = 0.04;
	gLowErrs_PR2[2] = 0.04;
	bVals_PR2[2] = 0.93;
	bHighErrs_PR2[2] = 0.05;
	bLowErrs_PR2[2] = 0.05;
	mean = bVals_PR2[2] * gVals_PR2[2];
	dMeanHigh = ((bVals_PR2[2]+bHighErrs_PR2[2]) * (gVals_PR2[2]+gHighErrs_PR2[2])) - mean;
	dMeanLow = mean - ((bVals_PR2[2]-bLowErrs_PR2[2]) * (gVals_PR2[2]-gLowErrs_PR2[2]));
	qfVals_PR2[2] = mean / eVals[2] * 100.;
	qfHighErrs_PR2[2] = dMeanHigh / eVals[2] * 100.;
	qfLowErrs_PR2[2] = dMeanLow / eVals[2] * 100.;
	sigVals_PR2[2] = 4795.33;
	sigHighErrs_PR2[2] = 668.39;
	sigLowErrs_PR2[2] = 661.50;
	noiseVals_PR2[2] = 23639.24;
	noiseHighErrs_PR2[2] = 676.17;
	noiseLowErrs_PR2[2] = 678.26;
	//cout << "QF Vals are " << qfVals_PR2[2] << endl;
	//cout << " + " << qfHighErrs_PR2[2] << endl;
	//cout << " - " << qfLowErrs_PR2[2] << endl;

	eVals_PR2[3] = eVals[3];
	eLowErrs_PR2[3] = eLowErrs[3];
	eHighErrs_PR2[3] = eHighErrs[3];
	gVals_PR2[3] = 1.44;
	gHighErrs_PR2[3] = 0.03;
	gLowErrs_PR2[3] = 0.03;
	bVals_PR2[3] = 0.61;
	bHighErrs_PR2[3] = 0.02;
	bLowErrs_PR2[3] = 0.02;
	mean = bVals_PR2[3] * gVals_PR2[3];
	dMeanHigh = ((bVals[3]+bHighErrs[3]) * (gVals[3]+gHighErrs[3])) - mean;
	dMeanLow = mean - ((bVals_PR2[3]-bLowErrs_PR2[3]) * (gVals_PR2[3]-gLowErrs_PR2[3]));
	qfVals_PR2[3] = mean / eVals[3] * 100.;
	qfHighErrs_PR2[3] = dMeanHigh / eVals[3] * 100.;
	qfLowErrs_PR2[3] = dMeanLow / eVals[3] * 100.;
	sigVals_PR2[3] = 8828.29;
	sigHighErrs_PR2[3] = 890.40;
	sigLowErrs_PR2[3] = 978.05;
	noiseVals_PR2[3] = 21315.47;
	noiseHighErrs_PR2[3] = 990.12;
	noiseLowErrs_PR2[3] = 901.06;
	//cout << "QF Vals are " << qfVals_PR2[3] << endl;
	//cout << " + " << qfHighErrs_PR2[3] << endl;
	//cout << " - " << qfLowErrs_PR2[3] << endl;
	
	//Make our TGraphs.
	TGraph* qf_graph = new TGraphAsymmErrors(numPoints, eVals, qfVals, eLowErrs, eHighErrs, qfLowErrs, qfHighErrs);
	qf_graph->SetTitle("Production Run 1");
	TGraph* outlier = new TGraphAsymmErrors(1, outlier_e, outlier_qfVal, out_e_LowErrs, out_e_HighErrs, out_qfLowErr, out_qfHighErr);
	outlier->SetTitle("Production Run 1 Outlier");
	TGraph* qf_graph_pr2 = new TGraphAsymmErrors(4, eVals_PR2, qfVals_PR2, eLowErrs_PR2, eHighErrs_PR2, qfLowErrs_PR2, qfHighErrs_PR2);
	qf_graph_pr2->SetTitle("Production Run 2");
	TMultiGraph* qf_MultiGraph = new TMultiGraph();
	qf_MultiGraph->Add(qf_graph,"AP");
	qf_MultiGraph->Add(outlier,"AP");
	qf_MultiGraph->Add(qf_graph_pr2,"AP");
	TGraph* gamma_graph = new TGraphAsymmErrors(numPoints, eVals, gVals, eLowErrs, eHighErrs, gLowErrs, gHighErrs);
	gamma_graph->SetTitle("Production Run 1");
	TGraph* g_outlier = new TGraphAsymmErrors(1, outlier_e, outlier_g, out_e_LowErrs, out_e_HighErrs, out_g_LowErr, out_g_HighErr);
	g_outlier->SetTitle("Production Run 1 Outlier");
	TGraph* gamma_graph_pr2 = new TGraphAsymmErrors(4, eVals_PR2, gVals_PR2, eLowErrs_PR2, eHighErrs_PR2, gLowErrs_PR2, gHighErrs_PR2);
	gamma_graph_pr2->SetTitle("Production Run 2");
	TMultiGraph* g_MultiGraph = new TMultiGraph();
	g_MultiGraph->Add(gamma_graph,"AP");
	g_MultiGraph->Add(g_outlier,"AP");
	g_MultiGraph->Add(gamma_graph_pr2,"AP");
	TGraph* beta_graph = new TGraphAsymmErrors(numPoints, eVals, bVals, eLowErrs, eHighErrs, bLowErrs, bHighErrs);
	beta_graph->SetTitle("Production Run 1");
	TGraph* b_outlier = new TGraphAsymmErrors(1, outlier_e, outlier_b, out_e_LowErrs, out_e_HighErrs, out_b_LowErr, out_b_HighErr);
	b_outlier->SetTitle("Production Run 1 Outlier");
	TGraph* beta_graph_pr2 = new TGraphAsymmErrors(4, eVals_PR2, bVals_PR2, eLowErrs_PR2, eHighErrs_PR2, bLowErrs_PR2, bHighErrs_PR2);
	beta_graph_pr2->SetTitle("Production Run 2");
	TMultiGraph* b_MultiGraph = new TMultiGraph();
	b_MultiGraph->Add(beta_graph,"AP");
	b_MultiGraph->Add(b_outlier,"AP");
	b_MultiGraph->Add(beta_graph_pr2,"AP");
	TGraph* signal_graph = new TGraphAsymmErrors(numPoints, eVals, sigVals, eLowErrs, eHighErrs, sigLowErrs, sigHighErrs);
	signal_graph->SetTitle("Production Run 1");
	TGraph* sig_outlier = new TGraphAsymmErrors(1, outlier_e, outlier_sigVal, out_e_LowErrs, out_e_HighErrs, out_sigLowErr, out_sigHighErr);
	sig_outlier->SetTitle("Production Run 1 Outlier");
	TGraph* signal_graph_pr2 = new TGraphAsymmErrors(4, eVals_PR2, sigVals_PR2, eLowErrs, eHighErrs, sigLowErrs_PR2, sigHighErrs_PR2);
	signal_graph_pr2->SetTitle("Production Run 2");
	TMultiGraph* sig_MultiGraph = new TMultiGraph();
	sig_MultiGraph->Add(signal_graph,"AP");
	sig_MultiGraph->Add(sig_outlier,"AP");
	sig_MultiGraph->Add(signal_graph_pr2,"AP");
	TGraph* noise_graph = new TGraphAsymmErrors(numPoints, eVals, noiseVals, eLowErrs, eHighErrs, noiseLowErrs, noiseHighErrs);
	noise_graph->SetTitle("Production Run 1");
	TGraph* noise_outlier = new TGraphAsymmErrors(1, outlier_e, outlier_noiseVal, out_e_LowErrs, out_e_HighErrs, out_noiseLowErr, out_noiseHighErr);
	noise_outlier->SetTitle("Production Run 1 Outlier");
	TGraph* noise_graph_pr2 = new TGraphAsymmErrors(4, eVals_PR2, noiseVals_PR2, eLowErrs, eHighErrs, noiseLowErrs_PR2, noiseHighErrs_PR2);
	noise_graph_pr2->SetTitle("Production Run 2");
	TMultiGraph* noise_MultiGraph = new TMultiGraph();
	noise_MultiGraph->Add(noise_graph,"AP");
	noise_MultiGraph->Add(noise_outlier,"AP");
	noise_MultiGraph->Add(noise_graph_pr2,"AP");
	qf_graph->SetMarkerColor(45);
	qf_graph->SetMarkerStyle(21);
	qf_graph_pr2->SetMarkerColor(40);
	qf_graph_pr2->SetMarkerStyle(21);
	qf_MultiGraph->SetTitle("Nuclear Quenching Factors in CeBr3");
	qf_MultiGraph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	qf_MultiGraph->GetXaxis()->CenterTitle();
	qf_MultiGraph->GetYaxis()->SetTitle("Nuclear Quenching Factor (%)");
	qf_MultiGraph->GetYaxis()->CenterTitle();
	qf_MultiGraph->GetYaxis()->SetRangeUser(0,20);
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
	gamma_graph_pr2->SetMarkerColor(28);
	gamma_graph_pr2->SetMarkerStyle(21);
	beta_graph->SetMarkerColor(41);
	beta_graph->SetMarkerStyle(21);
	beta_graph_pr2->SetMarkerColor(45);
	beta_graph_pr2->SetMarkerStyle(21);
	b_MultiGraph->SetTitle("Beta vs. Recoil Energy");
	b_MultiGraph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	b_MultiGraph->GetXaxis()->CenterTitle();
	b_MultiGraph->GetYaxis()->SetTitle("Beta");
	b_MultiGraph->GetYaxis()->CenterTitle();
	b_outlier->SetMarkerColor(41);
	b_outlier->SetMarkerStyle(24);
	signal_graph->SetMarkerStyle(21);
	signal_graph->SetMarkerColor(29);
	signal_graph_pr2->SetMarkerStyle(21);
	signal_graph_pr2->SetMarkerColor(49);
	sig_MultiGraph->SetTitle("Signal Counts vs. Recoil Energy");
	sig_MultiGraph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	sig_MultiGraph->GetXaxis()->CenterTitle();
	sig_MultiGraph->GetYaxis()->SetTitle("Signal Counts");
	sig_MultiGraph->GetYaxis()->CenterTitle();
	sig_outlier->SetMarkerColor(29);
	sig_outlier->SetMarkerStyle(24);
	noise_graph->SetMarkerStyle(21);
	noise_graph->SetMarkerColor(13);
	noise_graph_pr2->SetMarkerStyle(21);
	noise_graph_pr2->SetMarkerColor(42);
	noise_MultiGraph->SetTitle("Background Counts vs. Recoil Energy");
	noise_MultiGraph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	noise_MultiGraph->GetXaxis()->CenterTitle();
	noise_MultiGraph->GetYaxis()->SetTitle("Background Counts");
	noise_MultiGraph->GetYaxis()->CenterTitle();
	noise_outlier->SetMarkerColor(13);
	noise_outlier->SetMarkerStyle(24);
	
	//Make a canvas and plot.
	TCanvas* c1 = new TCanvas("c1","c1");
	c1->Divide(2,2);
	c1->cd(1);
	g_MultiGraph->Draw("a");
	c1->BuildLegend(0.62,0.67,0.87,0.85);
	c1->Modified();
	c1->Update();
	c1->cd(2);
	b_MultiGraph->Draw("a");
	c1->BuildLegend(0.62,0.67,0.87,0.85);
	c1->Modified();
	c1->Update();
	c1->cd(3);
	sig_MultiGraph->Draw("a");
	c1->BuildLegend(0.62,0.67,0.87,0.85);
	c1->Modified();
	c1->Update();
	c1->cd(4);
	noise_MultiGraph->Draw("a");
	c1->BuildLegend(0.62,0.67,0.87,0.85);
	c1->Modified();
	c1->Update();
	
	TCanvas* c2 = new TCanvas("c2","Quenching Factor in CeBr3");
	c2->cd();
	qf_MultiGraph->Draw("a");
	c2->BuildLegend(0.62,0.67,0.87,0.85);
	c2->Modified();
	c2->Update();
	
	TString imagePath = "/var/phy/project/phil/cma46/CeBr3/cebr3-qf-analysis/Plots/";
	c1->Print(imagePath + "Parameter_Plots.png");
	c2->Print(imagePath + "QF_Plot.png");
	
}
	
	
	
	
	
	
	
