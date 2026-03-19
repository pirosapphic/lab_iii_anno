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

void fitCs137(string input = "../data_sorg/A8_Cs137_total.txt") 
{

    ifstream parInput(input.c_str());

    // in ADC counts
    float xMin = 525.; float xMax = 30525.;  float binWidth = 50.;

    TCanvas* c1 = new TCanvas("c1","c1",20,20,800,600);
    float x;  float y; string parName; 
    TH1F* theHisto = new TH1F("theHisto", "Electron energy", (int)((xMax-xMin)/binWidth), xMin, xMax);
    //theHisto->Sumw2();

    if (parInput.is_open()) {
	while ( parInput.good() ) {
	    parInput >> x >> y;   
	    cout << x << " " << y << endl;
	    int theBin = theHisto->FindBin(x);
	    theHisto->SetBinContent(theBin,float(y));
	    //theHisto->SetBinError(theBin,sqrt(float(y)));
	}
	parInput.close();
    }

    // rebinning  
    //theHisto->Rebin(2);    // constant

    // fit

    //estremi
    double min = 11000;
    double max = 30000;
    TF1* ffit = new TF1("ffit","gaus(0)+gaus(3)+[6]",0.,50000.);
    //ffit->SetParLimits(0,0.,10000.);
    //ffit->SetParLimits(1,11000.,13000.);
    //ffit->SetParLimits(2,1000.,4000.);
    ffit->SetParameter(0,700); //norm
    ffit->SetParameter(1,20000);//media
    ffit->SetParameter(2,2000);//std dev
    
    ffit->SetParameter(3,3000);
    ffit->SetParameter(4,9000);
    ffit->SetParameter(5,4000);


    ffit->SetParameter(6,0);
    /*[3]*pow(x-[4],2)+[5]
    ffit->SetParameter(3,1.2e-5); 
    ffit->SetParameter(4,25e3); //x vertice
    ffit->SetParameter(5,0);	//y vertice
*/
    theHisto->Draw("e1");
    theHisto->Fit("ffit","","e1",min,max);
    std::cout << "Chi^2:" << ffit->GetChisquare() << ", number of DoF: " << ffit->GetNDF() << " (Probability: " << ffit->GetProb() << ")." << std::endl;
  cout << "--------------------------------------------------------------------------------------------------------" << std::endl;


    //draw components
    TF1* gauscomp = (TF1*)ffit->Clone();
    gauscomp->SetRange(min,max);
    gauscomp->SetLineColor(kBlue);   gauscomp->SetLineStyle(kDashed);
    gauscomp->SetParameter(3,0.);
    gauscomp->SetParameter(6,0.);
    //gauscomp->SetParameter(5,0.);
    TF1* parabcomp = (TF1*)ffit->Clone();
    parabcomp->SetRange(min,max);
    parabcomp->SetLineColor(kRed);   parabcomp->SetLineStyle(kDashed);
    parabcomp->SetParameter(0,0.);
    parabcomp->Draw("same");
    gauscomp->Draw("same");

    //c1.SaveAs("fitCs137.gif");
}
