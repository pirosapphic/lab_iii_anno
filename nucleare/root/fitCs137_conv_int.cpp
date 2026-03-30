#include <iostream>
#include <fstream>
#include <math.h>   

#include "TH1F.h"
#include "TF1.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"

void fitCs137_conv_int(string input = "../data_sorg/A8_Cs137_total.txt") 
{

    ifstream parInput(input.c_str());

    // in ADC counts
    float xMin = 532.5; float xMax = 35567.5;  float binWidth = 65.;

    TCanvas* c1 = new TCanvas("c1","c1",20,20,800,600);
    c1->SetGrid();
    float x;  float y; string parName; 
    TH1F* theHisto = new TH1F("Energy Spectrum", "Electron energy ^{137}Cs", (int)((xMax-xMin)/binWidth), xMin, xMax);
    //theHisto->Sumw2();

    if (parInput.is_open()) {
	while ( parInput.good() ) {
	    parInput >> x >> y;   
	    int theBin = theHisto->FindBin(x);
	    theHisto->SetBinContent(theBin,float(y));
	    theHisto->SetBinError(theBin,sqrt(float(y)));
	}
	parInput.close();
    }

    // rebinning  
    //theHisto->Rebin(2);	//constant
    gStyle->SetOptStat(0);	//toglie tabellina di merda
    // fit

    //estremi di fit
    double min = 11000;
    double max = 30000;

    TF1* ffit = new TF1("ffit","gaus(0)+gaus(3)+[6]",0.,50000.);
    ffit->SetParameter(0,700); //norm
    ffit->SetParameter(1,20000);//media
    ffit->SetParameter(2,2000);//sigma
    
    ffit->SetParameter(3,3000);
    ffit->SetParameter(4,9000);
    ffit->SetParameter(5,4000);


    ffit->SetParameter(6,0);
    
    theHisto->Draw("e1");
    theHisto->Fit("ffit","","e1",min,max);
    std::cout << "Chi^2:" << ffit->GetChisquare() << ", number of DoF: " << ffit->GetNDF();
    std::cout << " (Probability: " << ffit->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;


    //draw components
    TF1* lilgaus = (TF1*)ffit->Clone();
    lilgaus->SetRange(min,max);
    lilgaus->SetLineColor(kBlue);   lilgaus->SetLineStyle(kDashed);
    lilgaus->SetParameter(3,0.);
    TF1* biggaus = (TF1*)ffit->Clone();
    biggaus->SetRange(min,max);
    biggaus->SetLineColor(kRed);   biggaus->SetLineStyle(kDashed);
    biggaus->SetParameter(0,0.);
    biggaus->Draw("same");
    lilgaus->Draw("same");
    
    //getting the gaussian peak
    double peak = ffit->GetParameter(1);  //mean
    double s_peak = ffit->GetParError(1);
    double sigma = ffit->GetParameter(2);//std dev
    double s_sigma = ffit->GetParError(2);
    std::cout << "Picco gamma = ("<<peak<<" +/- "<<s_peak<<")CHN\n";
    double k = 2.4e-5;		double s_k = 5e-6; //MeV/CHN
    double T_gamma = k*peak;
    double s_T_gamma = T_gamma * sqrt(pow(s_k/k,2)+pow(s_peak/peak,2));
    std::cout << "T_gamma = ("<<T_gamma*1000<<" +/- "<<s_T_gamma*1000<<")keV\n";
    double T_teo = 0.629; //MeV
    double s_T_teo = 0.016;
    std::cout << "E_teorica = ("<<T_teo*1000<<" +/- "<<s_T_teo*1000<<")keV\n";
    double z = (T_teo - T_gamma)/(sqrt(pow(s_T_teo,2)+pow(s_T_gamma,2)));
    std::cout << "z-score = "<<z<<std::endl;
    double R_exp = sigma/peak; //risoluzione sperimentiale
    double s_R_exp = R_exp*sqrt(pow(s_sigma/sigma,2)+pow(s_peak/peak,2));
    std::cout << "Risoluzione sperimentale = "<<R_exp<<" +/- "<<s_R_exp<<std::endl;
}
