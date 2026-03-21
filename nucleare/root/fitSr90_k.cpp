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

void fitSr90_k(string input = "../data_sorg/A8_Sr90_centrato.txt"){    

    ifstream parInput(input.c_str());

    // in ADC counts
    double xMin = 532.5;
    double xMax = 67027.5;
    double xMinep = 57017.5;
    double xMaxep = xMax;

    double binWidth = 65.;
    int nbinstot = (int)((xMax-xMin)/binWidth)+1;
    int nbinsep = (int)((xMaxep-xMinep)/binWidth)+1;
    
    double x;
    double y;
    string parName; 
   // TH1D* histoTot = new TH1D("Energy Spectrum", "Electron energy ^{90}Sr;ADC [CHN];Counts [#]",nbinstot, xMin, xMax);
    TH1D* histoEp = new TH1D("Endopint", "Kurie fit, ^{90}Sr;ADC [CHN];#frac{#sqrt{counts}}{CHN} [#frac{#}{CHN}]",nbinsep, xMinep, xMaxep);

    if (parInput.is_open()) {
	while ( parInput.good() ) {
	    parInput >> x >> y;   
	    double s_y = sqrt(double(y));
	/*
	    int theBin = histoTot->FindBin(x);
	    histoTot->SetBinContent(theBin,double(y));
	    histoTot->SetBinError(theBin,s_y);
	*/
	    if(x >= xMinep and x<=xMaxep){
		int theBinEp = histoEp->FindBin(x);
		histoEp->SetBinContent(theBinEp,double(y));
	        histoEp->SetBinError(theBinEp,s_y);
	    }
	}
	parInput.close();
    }
   // TCanvas* ctot = new TCanvas("c1","c1",20,20,800,600);
    //histoTot->Rebin(2);
   // histoTot->Draw("e1");
    
    //endpoint
    TCanvas* cep = new TCanvas("c2","c2",20,20,800,600);
    cep->cd();
    histoEp->Rebin(31);
    for (double x = xMinep; x<=xMaxep; x+=binWidth){
	int theBin = histoEp->FindBin(x);
	double s_x = histoEp->GetBinWidth(theBin);
	double y = histoEp->GetBinContent(theBin);
	double s_y = histoEp->GetBinError(theBin);
	double new_y = sqrt(y)/x;
	double s_new_y = sqrt(pow(new_y/x *s_x,2)+pow(0.5*new_y/y * s_y,2));
	histoEp->SetBinContent(theBin,new_y);
	histoEp->SetBinError(theBin,s_new_y);
    }
    histoEp->Draw("e1");
    TF1* f1 = new TF1("f1","pol1",57000,66200);
    f1->SetParameter(0,0);
    f1->SetParameter(1,8e-6);
    histoEp->Fit("f1","","e1",xMinep,xMaxep);
        std::cout << "Chi^2:" <<f1->GetChisquare() << ", number of DoF: " << f1->GetNDF();
    std::cout << " (Probability: " << f1->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    double m = f1->GetParameter(1);
    double s_m = f1->GetParError(1);
    double q = f1->GetParameter(0);
    double s_q = f1->GetParError(0);
    double CH0 = -q/m;
    double s_CH0 = sqrt(pow(CH0/m * s_m,2)+pow(CH0/q * s_q,2));
    std::cout<<"CHN_0 = ("<<CH0<<" +/- "<<s_CH0<<")CHN\n";
    //T_0 = (2.28 pm 0.04) MeV
    double k = 2.28/CH0;
    double s_k = sqrt(pow(k/CH0 * s_CH0,2)+pow(k/2.28 * 0.04,2));
   std::cout<<"k = ("<< k<<" +/- "<<s_k<<")MeV/CHN\n";
}
