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
    ffit->SetParameter(2,2000);//std dev
    
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

}
