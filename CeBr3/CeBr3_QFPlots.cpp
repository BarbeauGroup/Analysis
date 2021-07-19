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
	double eVals[numPoints]; //Recoil energy (best guess)
	double eHighErrs[numPoints]; 
	double eLowErrs[numPoints]; 
	double brVals[numPoints]; //Recoil energy assuming only bromine recoils.
	double brHighErrs[numPoints]; 
	double brLowErrs[numPoints]; 
	double ceVals[numPoints]; //Recoil energy assuming only cerium recoils.
	double ceHighErrs[numPoints]; 
	double ceLowErrs[numPoints]; 
	double yield[numPoints];
	double yieldLowErrs[numPoints];
	double yieldHighErrs[numPoints];
	double qfVals[numPoints]; //Quenching factor assuming best guess recoil energies.
	double qfHighErrs[numPoints]; 
	double qfLowErrs[numPoints]; 
	double br_qfVals[numPoints]; //Quenching factor assuming bromine recoil energies.
	double br_qfHighErrs[numPoints]; 
	double br_qfLowErrs[numPoints]; 
	double ce_qfVals[numPoints]; //Quenching factor assuming cerium recoil energies.
	double ce_qfHighErrs[numPoints]; 
	double ce_qfLowErrs[numPoints]; 
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
	brVals[0] = 78.84;
	brLowErrs[0] = 2.71;
	brHighErrs[0] = 2.71;
	ceVals[0] = 45.34;
	ceLowErrs[0] = 2.63;
	ceHighErrs[0] = 2.63;
	gVals[0] = 9.12;
	gHighErrs[0] = 0.67;
	gLowErrs[0] = 0.63;
	bVals[0] = 0.40;
	bHighErrs[0] = 0.03;
	bLowErrs[0] = 0.03;
	yield[0] = bVals[0] * gVals[0];
	yieldHighErrs[0] = ((bVals[0]+bHighErrs[0]) * (gVals[0]+gHighErrs[0])) - yield[0];
	yieldLowErrs[0] = yield[0] - ((bVals[0]-bLowErrs[0]) * (gVals[0]-gLowErrs[0]));
	qfVals[0] = yield[0] / eVals[0] * 100.;
	qfHighErrs[0] = yieldHighErrs[0] / eVals[0] * 100.;
	qfLowErrs[0] = yieldLowErrs[0] / eVals[0] * 100.;
	br_qfVals[0] = yield[0] / brVals[0] * 100.;
	br_qfHighErrs[0] = yieldHighErrs[0] / brVals[0] * 100.;
	br_qfLowErrs[0] = yieldLowErrs[0] / brVals[0] * 100.;
	ce_qfVals[0] = yield[0] / ceVals[0] * 100.;
	ce_qfHighErrs[0] = yieldHighErrs[0] / ceVals[0] * 100.;
	ce_qfLowErrs[0] = yieldLowErrs[0] / ceVals[0] * 100.;
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
	brVals[1] = 65.49;
	brLowErrs[1] = 2.80;
	brHighErrs[1] = 2.80;
	ceVals[1] = 36.82;
	ceLowErrs[1] = 1.81;
	ceHighErrs[1] = 1.81;
	gVals[1] = 6.40;
	gHighErrs[1] = 0.53;
	gLowErrs[1] = 0.50;
	bVals[1] = 0.46;
	bHighErrs[1] = 0.04;
	bLowErrs[1] = 0.03;
	yield[1] = bVals[1] * gVals[1];
	yieldHighErrs[1] = ((bVals[1]+bHighErrs[1]) * (gVals[1]+gHighErrs[1])) - yield[1];
	yieldLowErrs[1] = yield[1] - ((bVals[1]-bLowErrs[1]) * (gVals[1]-gLowErrs[1]));
	qfVals[1] = yield[1] / eVals[1] * 100.;
	qfHighErrs[1] = yieldHighErrs[1] / eVals[1] * 100.;
	qfLowErrs[1] = yieldLowErrs[1] / eVals[1] * 100.;
	br_qfVals[1] = yield[1] / brVals[1] * 100.;
	br_qfHighErrs[1] = yieldHighErrs[1] / brVals[1] * 100.;
	br_qfLowErrs[1] = yieldLowErrs[1] / brVals[1] * 100.;
	ce_qfVals[1] = yield[1] / ceVals[1] * 100.;
	ce_qfHighErrs[1] = yieldHighErrs[1] / ceVals[1] * 100.;
	ce_qfLowErrs[1] = yieldLowErrs[1] / ceVals[1] * 100.;
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
	brVals[2] = 52.87;
	brLowErrs[2] = 2.24;
	brHighErrs[2] = 2.24;
	ceVals[2] = 30.04;
	ceLowErrs[2] = 1.44;
	ceHighErrs[2] = 1.44;
	gVals[2] = 2.37;
	gHighErrs[2] = 0.27;
	gLowErrs[2] = 0.20;
	bVals[2] = 0.88;
	bHighErrs[2] = 0.07;
	bLowErrs[2] = 0.08;
	yield[2] = bVals[2] * gVals[2];
	yieldHighErrs[2] = ((bVals[2]+bHighErrs[2]) * (gVals[2]+gHighErrs[2])) - yield[2];
	yieldLowErrs[2] = yield[2] - ((bVals[2]-bLowErrs[2]) * (gVals[2]-gLowErrs[2]));
	qfVals[2] = yield[2] / eVals[2] * 100.;
	qfHighErrs[2] = yieldHighErrs[2] / eVals[2] * 100.;
	qfLowErrs[2] = yieldLowErrs[2] / eVals[2] * 100.;
	br_qfVals[2] = yield[2] / brVals[2] * 100.;
	br_qfHighErrs[2] = yieldHighErrs[2] / brVals[2] * 100.;
	br_qfLowErrs[2] = yieldLowErrs[2] / brVals[2] * 100.;
	ce_qfVals[2] = yield[2] / ceVals[2] * 100.;
	ce_qfHighErrs[2] = yieldHighErrs[2] / ceVals[2] * 100.;
	ce_qfLowErrs[2] = yieldLowErrs[2] / ceVals[2] * 100.;
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
	brVals[3] = 39.50;
	brLowErrs[3] = 2.43;
	brHighErrs[3] = 2.43;
	ceVals[3] = 22.34;
	ceLowErrs[3] = 1.34;
	ceHighErrs[3] = 1.34;
	gVals[3] = 1.51;
	gHighErrs[3] = 0.18;
	gLowErrs[3] = 0.14;
	bVals[3] = 0.65;
	bHighErrs[3] = 0.06;
	bLowErrs[3] = 0.06;
	yield[3] = bVals[3] * gVals[3];
	yieldHighErrs[3] = ((bVals[3]+bHighErrs[3]) * (gVals[3]+gHighErrs[3])) - yield[3];
	yieldLowErrs[3] = yield[3] - ((bVals[3]-bLowErrs[3]) * (gVals[3]-gLowErrs[3]));
	qfVals[3] = yield[3] / eVals[3] * 100.;
	qfHighErrs[3] = yieldHighErrs[3] / eVals[3] * 100.;
	qfLowErrs[3] = yieldLowErrs[3] / eVals[3] * 100.;
	br_qfVals[3] = yield[3] / brVals[3] * 100.;
	br_qfHighErrs[3] = yieldHighErrs[3] / brVals[3] * 100.;
	br_qfLowErrs[3] = yieldLowErrs[3] / brVals[3] * 100.;
	ce_qfVals[3] = yield[3] / ceVals[3] * 100.;
	ce_qfHighErrs[3] = yieldHighErrs[3] / ceVals[3] * 100.;
	ce_qfLowErrs[3] = yieldLowErrs[3] / ceVals[3] * 100.;
	sigVals[3] = 736.94;
	sigHighErrs[3] = 55.65;
	sigLowErrs[3] = 53.73;
	noiseVals[3] = 21316.77;
	noiseHighErrs[3] = 153.10;
	noiseLowErrs[3] = 147.78;
	//cout << "QF Vals are " << qfVals[3] << endl;
	//cout << " + " << qfHighErrs[3] << endl;
	//cout << " - " << qfLowErrs[3] << endl;
	
	eVals[4] = 16.32;
	eLowErrs[4] = 1.25;
	eHighErrs[4] = 1.25;
	brVals[4] = 23.27;
	brLowErrs[4] = 1.98;
	brHighErrs[4] = 1.98;
	ceVals[4] = 13.31;
	ceLowErrs[4] = 0.94;
	ceHighErrs[4] = 0.94;
	gVals[4] = 3.35;
	gHighErrs[4] = 0.38;
	gLowErrs[4] = 0.44;
	bVals[4] = 0.28;
	bHighErrs[4] = 0.03;
	bLowErrs[4] = 0.03;
	yield[4] = bVals[4] * gVals[4];
	yieldHighErrs[4] = ((bVals[4]+bHighErrs[4]) * (gVals[4]+gHighErrs[4])) - yield[4];
	yieldLowErrs[4] = yield[4] - ((bVals[4]-bLowErrs[4]) * (gVals[4]-gLowErrs[4]));
	qfVals[4] = yield[4] / eVals[4] * 100.;
	qfHighErrs[4] = yieldHighErrs[4] / eVals[4] * 100.;
	qfLowErrs[4] = yieldLowErrs[4] / eVals[4] * 100.;
	br_qfVals[4] = yield[4] / brVals[4] * 100.;
	br_qfHighErrs[4] = yieldHighErrs[4] / brVals[4] * 100.;
	br_qfLowErrs[4] = yieldLowErrs[4] / brVals[4] * 100.;
	ce_qfVals[4] = yield[4] / ceVals[4] * 100.;
	ce_qfHighErrs[4] = yieldHighErrs[4] / ceVals[4] * 100.;
	ce_qfLowErrs[4] = yieldLowErrs[4] / ceVals[4] * 100.;
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
	brVals[5] = 16.32;
	brLowErrs[5] = 1.41;
	brHighErrs[5] = 1.41;
	ceVals[5] = 9.20;
	ceLowErrs[5] = 0.84;
	ceHighErrs[5] = 0.84;
	gVals[5] = 1.82;
	gHighErrs[5] = 0.13;
	gLowErrs[5] = 0.10;
	bVals[5] = 0.36;
	bHighErrs[5] = 0.02;
	bLowErrs[5] = 0.02;
	yield[5] = bVals[5] * gVals[5];
	yieldHighErrs[5] = ((bVals[5]+bHighErrs[5]) * (gVals[5]+gHighErrs[5])) - yield[5];
	yieldLowErrs[5] = yield[5] - ((bVals[5]-bLowErrs[5]) * (gVals[5]-gLowErrs[5]));
	qfVals[5] = yield[5] / eVals[5] * 100.;
	qfHighErrs[5] = yieldHighErrs[5] / eVals[5] * 100.;
	qfLowErrs[5] = yieldLowErrs[5] / eVals[5] * 100.;
	br_qfVals[5] = yield[5] / brVals[5] * 100.;
	br_qfHighErrs[5] = yieldHighErrs[5] / brVals[5] * 100.;
	br_qfLowErrs[5] = yieldLowErrs[5] / brVals[5] * 100.;
	ce_qfVals[5] = yield[5] / ceVals[5] * 100.;
	ce_qfHighErrs[5] = yieldHighErrs[5] / ceVals[5] * 100.;
	ce_qfLowErrs[5] = yieldLowErrs[5] / ceVals[5] * 100.;
	sigVals[5] = 3140.23;
	sigHighErrs[5] = 118.91;
	sigLowErrs[5] = 122.32;
	noiseVals[5] = 30710.24;
	noiseHighErrs[5] = 197.74;
	noiseLowErrs[5] = 197.64;
	//cout << "QF Vals are" << qfVals[5] << endl;
	//cout << " + " << qfHighErrs[5] << endl;
	//cout << " - " << qfLowErrs[5] << endl;
	
	eVals[6] = 5.44;
	eLowErrs[6] = 0.63;
	eHighErrs[6] = 0.63;
	brVals[6] = 8.08;
	brLowErrs[6] = 0.93;
	brHighErrs[6] = 0.93;
	ceVals[6] = 4.66;
	ceLowErrs[6] = 0.54;
	ceHighErrs[6] = 0.54;
	gVals[6] = 1.66;
	gHighErrs[6] = 0.03;
	gLowErrs[6] = 0.03;
	bVals[6] = 0.21;
	bHighErrs[6] = 0.01;
	bLowErrs[6] = 0.01;
	yield[6] = bVals[6] * gVals[6];
	yieldHighErrs[6] = ((bVals[6]+bHighErrs[6]) * (gVals[6]+gHighErrs[6])) - yield[6];
	yieldLowErrs[6] = yield[6] - ((bVals[6]-bLowErrs[6]) * (gVals[6]-gLowErrs[6]));
	qfVals[6] = yield[6] / eVals[6] * 100.;
	qfHighErrs[6] = yieldHighErrs[6] / eVals[6] * 100.;
	qfLowErrs[6] = yieldLowErrs[6] / eVals[6] * 100.;
	br_qfVals[6] = yield[6] / brVals[6] * 100.;
	br_qfHighErrs[6] = yieldHighErrs[6] / brVals[6] * 100.;
	br_qfLowErrs[6] = yieldLowErrs[6] / brVals[6] * 100.;
	ce_qfVals[6] = yield[6] / ceVals[6] * 100.;
	ce_qfHighErrs[6] = yieldHighErrs[6] / ceVals[6] * 100.;
	ce_qfLowErrs[6] = yieldLowErrs[6] / ceVals[6] * 100.;
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
	outlier_e[0] = 1.55;
	double out_e_HighErrs[1];
	out_e_HighErrs[0] = 0.38;
	double out_e_LowErrs[1];
	out_e_LowErrs[0] = 0.38;
	double outlier_br[1];
	outlier_br[0] = 1.94;
	double out_br_HighErrs[1];
	out_br_HighErrs[0] = 0.42;
	double out_br_LowErrs[1];
	out_br_LowErrs[0] = 0.42;
	double outlier_ce[1];
	outlier_ce[0] = 1.07;
	double out_ce_HighErrs[1];
	out_ce_HighErrs[0] = 0.32;
	double out_ce_LowErrs[1];
	out_ce_LowErrs[0] = 0.32;
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
	double out_yield[1];
	double out_yield_LowErrs[1];
	double out_yield_HighErrs[1];
	out_yield[0] = outlier_b[0] * outlier_g[0];
	out_yield_HighErrs[0] = ((outlier_b[0]+out_b_HighErr[0]) * (outlier_g[0]+out_g_HighErr[0])) - out_yield[0];
	out_yield_LowErrs[0] = out_yield[0] - (outlier_b[0]-out_b_LowErr[0]) * (outlier_g[0]-out_g_LowErr[0]);
	double outlier_qfVal[1];
	outlier_qfVal[0] = out_yield[0] / outlier_e[0] * 100.;
	double out_qfHighErr[1];
	out_qfHighErr[0] = out_yield_HighErrs[0] / outlier_e[0] * 100.;
	double out_qfLowErr[1];
	out_qfLowErr[0] = out_yield_LowErrs[0] / outlier_e[0] * 100.;
	double outlier_br_qfVal[1];
	outlier_br_qfVal[0] = out_yield[0] / outlier_br[0] * 100.;
	double out_br_qfHighErr[1];
	out_br_qfHighErr[0] = out_yield_HighErrs[0] / outlier_br[0] * 100.;
	double out_br_qfLowErr[1];
	out_br_qfLowErr[0] = out_yield_LowErrs[0] / outlier_br[0] * 100.;
	double outlier_ce_qfVal[1];
	outlier_ce_qfVal[0] = out_yield[0] / outlier_ce[0] * 100.;
	double out_ce_qfHighErr[1];
	out_ce_qfHighErr[0] = out_yield_HighErrs[0] / outlier_ce[0] * 100.;
	double out_ce_qfLowErr[1];
	out_ce_qfLowErr[0] = out_yield_LowErrs[0] / outlier_ce[0] * 100.;
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
	/*
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
	*/
	//Make our TGraphs.

	///////////////////////////////////
	///// Quenching Factor Graphs /////
	///////////////////////////////////

	TGraph* qf_graph = new TGraphAsymmErrors(numPoints, eVals, qfVals, eLowErrs, eHighErrs, qfLowErrs, qfHighErrs);
	qf_graph->SetTitle("Average Recoil Energy");
	TGraph* br_qf_graph = new TGraphAsymmErrors(numPoints, brVals, br_qfVals, brLowErrs, brHighErrs, br_qfLowErrs, br_qfHighErrs);
	br_qf_graph->SetTitle("Assuming Bromine Recoils");
	TGraph* ce_qf_graph = new TGraphAsymmErrors(numPoints, ceVals, ce_qfVals, ceLowErrs, ceHighErrs, ce_qfLowErrs, ce_qfHighErrs);
	ce_qf_graph->SetTitle("Assuming Cerium Recoils");
	TGraph* outlier = new TGraphAsymmErrors(1, outlier_e, outlier_qfVal, out_e_LowErrs, out_e_HighErrs, out_qfLowErr, out_qfHighErr);
	TGraph* br_outlier = new TGraphAsymmErrors(1, outlier_br, outlier_br_qfVal, out_br_LowErrs, out_br_HighErrs, out_br_qfLowErr, out_br_qfHighErr);
	TGraph* ce_outlier = new TGraphAsymmErrors(1, outlier_ce, outlier_ce_qfVal, out_ce_LowErrs, out_ce_HighErrs, out_ce_qfLowErr, out_ce_qfHighErr);
	//TGraph* qf_graph_pr2 = new TGraphAsymmErrors(4, eVals_PR2, qfVals_PR2, eLowErrs_PR2, eHighErrs_PR2, qfLowErrs_PR2, qfHighErrs_PR2);
	//qf_graph_pr2->SetTitle("Production Run 2");
	TMultiGraph* qf_MultiGraph = new TMultiGraph();
	qf_MultiGraph->Add(qf_graph,"AP");
	qf_MultiGraph->Add(br_qf_graph,"AP");
	qf_MultiGraph->Add(ce_qf_graph,"AP");
	qf_MultiGraph->Add(outlier,"AP");
	qf_MultiGraph->Add(br_outlier,"AP");
	qf_MultiGraph->Add(ce_outlier,"AP");

	////////////////////////
	///// Yield Graphs /////
	////////////////////////

	TGraph* yield_graph = new TGraphAsymmErrors(numPoints, eVals, yield, eLowErrs, eHighErrs, yieldLowErrs, yieldHighErrs);
	yield_graph->SetTitle("Average Recoil Energy");
	TGraph* br_yield_graph = new TGraphAsymmErrors(numPoints, brVals, yield, brLowErrs, brHighErrs, yieldLowErrs, yieldHighErrs);
	br_yield_graph->SetTitle("Assuming Bromine Recoils");
	TGraph* ce_yield_graph = new TGraphAsymmErrors(numPoints, ceVals, yield, ceLowErrs, ceHighErrs, yieldLowErrs, yieldHighErrs);
	ce_yield_graph->SetTitle("Assuming Cerium Recoils");
	TGraph* yield_outlier = new TGraphAsymmErrors(1, outlier_e, out_yield, out_e_LowErrs, out_e_HighErrs, out_yield_LowErrs, out_yield_HighErrs);
	TGraph* br_yield_outlier = new TGraphAsymmErrors(1, outlier_br, out_yield, out_br_LowErrs, out_br_HighErrs, out_yield_LowErrs, out_yield_HighErrs);
	TGraph* ce_yield_outlier = new TGraphAsymmErrors(1, outlier_ce, out_yield, out_ce_LowErrs, out_ce_HighErrs, out_yield_LowErrs, out_yield_HighErrs);
	TMultiGraph* yield_MultiGraph = new TMultiGraph();
	yield_MultiGraph->Add(yield_graph,"AP");
	yield_MultiGraph->Add(br_yield_graph,"AP");
	yield_MultiGraph->Add(ce_yield_graph,"AP");
	yield_MultiGraph->Add(yield_outlier,"AP");
	yield_MultiGraph->Add(br_yield_outlier,"AP");
	yield_MultiGraph->Add(ce_yield_outlier,"AP");

	////////////////////////////
	///// Parameter Graphs /////
	////////////////////////////

	TGraph* gamma_graph = new TGraphAsymmErrors(numPoints, eVals, gVals, eLowErrs, eHighErrs, gLowErrs, gHighErrs);
	gamma_graph->SetTitle("Production Run 1");
	TGraph* g_outlier = new TGraphAsymmErrors(1, outlier_e, outlier_g, out_e_LowErrs, out_e_HighErrs, out_g_LowErr, out_g_HighErr);
	g_outlier->SetTitle("Production Run 1 Outlier");
	//TGraph* gamma_graph_pr2 = new TGraphAsymmErrors(4, eVals_PR2, gVals_PR2, eLowErrs_PR2, eHighErrs_PR2, gLowErrs_PR2, gHighErrs_PR2);
	//gamma_graph_pr2->SetTitle("Production Run 2");
	TMultiGraph* g_MultiGraph = new TMultiGraph();
	g_MultiGraph->Add(gamma_graph,"AP");
	g_MultiGraph->Add(g_outlier,"AP");
	//g_MultiGraph->Add(gamma_graph_pr2,"AP");
	TGraph* beta_graph = new TGraphAsymmErrors(numPoints, eVals, bVals, eLowErrs, eHighErrs, bLowErrs, bHighErrs);
	beta_graph->SetTitle("Production Run 1");
	TGraph* b_outlier = new TGraphAsymmErrors(1, outlier_e, outlier_b, out_e_LowErrs, out_e_HighErrs, out_b_LowErr, out_b_HighErr);
	b_outlier->SetTitle("Production Run 1 Outlier");
	//TGraph* beta_graph_pr2 = new TGraphAsymmErrors(4, eVals_PR2, bVals_PR2, eLowErrs_PR2, eHighErrs_PR2, bLowErrs_PR2, bHighErrs_PR2);
	//beta_graph_pr2->SetTitle("Production Run 2");
	TMultiGraph* b_MultiGraph = new TMultiGraph();
	b_MultiGraph->Add(beta_graph,"AP");
	b_MultiGraph->Add(b_outlier,"AP");
	//b_MultiGraph->Add(beta_graph_pr2,"AP");
	TGraph* signal_graph = new TGraphAsymmErrors(numPoints, eVals, sigVals, eLowErrs, eHighErrs, sigLowErrs, sigHighErrs);
	signal_graph->SetTitle("Production Run 1");
	TGraph* sig_outlier = new TGraphAsymmErrors(1, outlier_e, outlier_sigVal, out_e_LowErrs, out_e_HighErrs, out_sigLowErr, out_sigHighErr);
	sig_outlier->SetTitle("Production Run 1 Outlier");
	//TGraph* signal_graph_pr2 = new TGraphAsymmErrors(4, eVals_PR2, sigVals_PR2, eLowErrs, eHighErrs, sigLowErrs_PR2, sigHighErrs_PR2);
	//signal_graph_pr2->SetTitle("Production Run 2");
	TMultiGraph* sig_MultiGraph = new TMultiGraph();
	sig_MultiGraph->Add(signal_graph,"AP");
	sig_MultiGraph->Add(sig_outlier,"AP");
	//sig_MultiGraph->Add(signal_graph_pr2,"AP");
	TGraph* noise_graph = new TGraphAsymmErrors(numPoints, eVals, noiseVals, eLowErrs, eHighErrs, noiseLowErrs, noiseHighErrs);
	noise_graph->SetTitle("Production Run 1");
	TGraph* noise_outlier = new TGraphAsymmErrors(1, outlier_e, outlier_noiseVal, out_e_LowErrs, out_e_HighErrs, out_noiseLowErr, out_noiseHighErr);
	noise_outlier->SetTitle("Production Run 1 Outlier");
	//TGraph* noise_graph_pr2 = new TGraphAsymmErrors(4, eVals_PR2, noiseVals_PR2, eLowErrs, eHighErrs, noiseLowErrs_PR2, noiseHighErrs_PR2);
	//noise_graph_pr2->SetTitle("Production Run 2");
	TMultiGraph* noise_MultiGraph = new TMultiGraph();
	noise_MultiGraph->Add(noise_graph,"AP");
	noise_MultiGraph->Add(noise_outlier,"AP");
	//noise_MultiGraph->Add(noise_graph_pr2,"AP");

	////////////////////////
	///// Style Points /////
	////////////////////////

	qf_graph->SetMarkerColor(46);
	qf_graph->SetMarkerStyle(21);
	br_qf_graph->SetMarkerColor(30);
	br_qf_graph->SetMarkerStyle(20);
	ce_qf_graph->SetMarkerColor(40);
	ce_qf_graph->SetMarkerStyle(20);
	yield_graph->SetMarkerColor(46);
	yield_graph->SetMarkerStyle(21);
	br_yield_graph->SetMarkerColor(9);
	br_yield_graph->SetMarkerStyle(20);
	ce_yield_graph->SetMarkerColor(8);
	ce_yield_graph->SetMarkerStyle(20);
	//qf_graph_pr2->SetMarkerColor(40);
	//qf_graph_pr2->SetMarkerStyle(21);
	qf_MultiGraph->SetTitle("Nuclear Quenching Factors in CeBr3");
	qf_MultiGraph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	qf_MultiGraph->GetXaxis()->CenterTitle();
	qf_MultiGraph->GetYaxis()->SetTitle("Nuclear Quenching Factor (%)");
	qf_MultiGraph->GetYaxis()->CenterTitle();
	qf_MultiGraph->GetYaxis()->SetRangeUser(0,20);
	qf_MultiGraph->GetXaxis()->SetRangeUser(0,85);
	yield_MultiGraph->SetTitle("Nuclear Recoil Signal Yield in CeBr3");
	yield_MultiGraph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	yield_MultiGraph->GetXaxis()->CenterTitle();
	yield_MultiGraph->GetYaxis()->SetTitle("Detectable Signal [keVee]");
	yield_MultiGraph->GetYaxis()->CenterTitle();
	yield_MultiGraph->GetYaxis()->SetRangeUser(0,5);
	yield_MultiGraph->GetXaxis()->SetRangeUser(0,85);
	outlier->SetMarkerColor(46);
	outlier->SetMarkerStyle(25);
	br_outlier->SetMarkerColor(30);
	br_outlier->SetMarkerStyle(24);
	ce_outlier->SetMarkerColor(40);
	ce_outlier->SetMarkerStyle(24);
	yield_outlier->SetMarkerColor(kGreen-10);
	yield_outlier->SetMarkerStyle(25);
	br_yield_outlier->SetMarkerColor(9);
	br_yield_outlier->SetMarkerStyle(24);
	ce_yield_outlier->SetMarkerColor(8);
	ce_yield_outlier->SetMarkerStyle(24);
	gamma_graph->SetMarkerColor(9);
	gamma_graph->SetMarkerStyle(21);
	//gamma_graph_pr2->SetMarkerColor(28);
	//gamma_graph_pr2->SetMarkerStyle(21);
	g_MultiGraph->SetTitle("Gamma vs. Recoil Energy");
	g_MultiGraph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	g_MultiGraph->GetXaxis()->CenterTitle();
	g_MultiGraph->GetYaxis()->SetTitle("Gamma");
	g_MultiGraph->GetYaxis()->CenterTitle();
	g_outlier->SetMarkerColor(9);
	g_outlier->SetMarkerStyle(24);
	beta_graph->SetMarkerColor(41);
	beta_graph->SetMarkerStyle(21);
	//beta_graph_pr2->SetMarkerColor(45);
	//beta_graph_pr2->SetMarkerStyle(21);
	b_MultiGraph->SetTitle("Beta vs. Recoil Energy");
	b_MultiGraph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	b_MultiGraph->GetXaxis()->CenterTitle();
	b_MultiGraph->GetYaxis()->SetTitle("Beta");
	b_MultiGraph->GetYaxis()->CenterTitle();
	b_outlier->SetMarkerColor(41);
	b_outlier->SetMarkerStyle(24);
	signal_graph->SetMarkerStyle(21);
	signal_graph->SetMarkerColor(29);
	//signal_graph_pr2->SetMarkerStyle(21);
	//signal_graph_pr2->SetMarkerColor(49);
	sig_MultiGraph->SetTitle("Signal Counts vs. Recoil Energy");
	sig_MultiGraph->GetXaxis()->SetTitle("Recoil Energy [keVnr]");
	sig_MultiGraph->GetXaxis()->CenterTitle();
	sig_MultiGraph->GetYaxis()->SetTitle("Signal Counts");
	sig_MultiGraph->GetYaxis()->CenterTitle();
	sig_outlier->SetMarkerColor(29);
	sig_outlier->SetMarkerStyle(24);
	noise_graph->SetMarkerStyle(21);
	noise_graph->SetMarkerColor(13);
	//noise_graph_pr2->SetMarkerStyle(21);
	//noise_graph_pr2->SetMarkerColor(42);
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
	//c1->BuildLegend(0.62,0.67,0.87,0.85);
	c1->Modified();
	c1->Update();
	c1->cd(2);
	b_MultiGraph->Draw("a");
	//c1->BuildLegend(0.62,0.67,0.87,0.85);
	c1->Modified();
	c1->Update();
	c1->cd(3);
	sig_MultiGraph->Draw("a");
	//c1->BuildLegend(0.62,0.67,0.87,0.85);
	c1->Modified();
	c1->Update();
	c1->cd(4);
	noise_MultiGraph->Draw("a");
	//c1->BuildLegend(0.62,0.67,0.87,0.85);
	c1->Modified();
	c1->Update();
	
	TCanvas* c2 = new TCanvas("c2","Quenching Factor in CeBr3");
	c2->cd();
	qf_MultiGraph->Draw("a");
	c2->BuildLegend(0.62,0.67,0.87,0.85);
	c2->Modified();
	c2->Update();

	TCanvas* c3 = new TCanvas("c3","Signal Yield in CeBr3");
	c3->cd();
	yield_MultiGraph->Draw("a");
	TLine* highQF = new TLine(0,0,50,5);
	highQF->SetLineColor(33);
	highQF->SetLineWidth(3);
	TLine* medQF = new TLine(0,0,83.33,5);
	medQF->SetLineColor(33);
	medQF->SetLineWidth(3);
	medQF->SetLineStyle(9);
	TLine* lowQF = new TLine(0,0,85,1.7);
	lowQF->SetLineColor(33);
	lowQF->SetLineWidth(3);
	highQF->Draw("same");
	lowQF->Draw("same");
	medQF->Draw("same");
	c3->BuildLegend(0.62,0.67,0.87,0.85);
	c3->Modified();
	c3->Update();
	
	TString imagePath = "/var/phy/project/phil/cma46/CeBr3/cebr3-qf-analysis/Plots/";
	c1->Print(imagePath + "Parameter_Plots.png");
	c2->Print(imagePath + "QF_Plot.png");
	c3->Print(imagePath + "Yield_Plot.png");
}
	
	
	
	
	
	
	
