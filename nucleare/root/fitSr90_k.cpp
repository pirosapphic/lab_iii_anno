#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>   

#include "TH1F.h"
#include "TF1.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"

void fitSr90_k(string input = "../data_sorg/A8_Sr90_centrato.txt"){    

    ifstream parInput(input.c_str());

    // in ADC counts
    double xMin = 532.5;
    double xMax = 67027.5; 
    double binWidth = 65.;
    int nbins = (int)((xMax-xMin)/binWidth);

    TCanvas* ctot = new TCanvas("c1","c1",20,20,800,600);
    double x;
    double y;
    string parName; 
    TH1D* histoTot = new TH1D("Energy Spectrum", "Electron energy ^{90}Sr",nbins, xMin, xMax);
    //theHisto->Sumw2();

    if (parInput.is_open()) {
	while ( parInput.good() ) {
	    parInput >> x >> y;   
	    int theBin = histoTot->FindBin(x);
	    histoTot->SetBinContent(theBin,double(y));
	    histoTot->SetBinError(theBin,sqrt(double(y)));
	}
	parInput.close();
    }
    histoTot->Draw("e1");
    
    //endpoint
    TCanvas* cep = new TCanvas("c2","c2",20,20,800,600);
    cep->cd();
    double xEndpoint = 50025.5; //first x value in endpoint fit
    int binEnd = (int)((xEndpoint - xMin)/binWidth); //first bin for the end-point fit
    TH1D* histoEp = new TH1D("Endpoint", "Kurie Plot",nbins-binEnd,xEndpoint,xMax);
    for(int i = binEnd; i < nbins; i++){
	int theBin = histoTot->FindBin(i);
	double binContent = histoTot->GetBinContent(theBin);
	histoEp->SetBinContent(i,binContent);
    }
    histoEp->Draw("e1");
}
