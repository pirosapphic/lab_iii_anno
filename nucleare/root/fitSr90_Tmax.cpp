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

void fitSr90_Tmax(string input = "../data_sorg/A8_Sr90_centrato.txt"){    

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
    TH1D* histoTot = new TH1D("Energy Spectrum", "Electron energy ^{90}Sr;ADC [CHN];Counts [#]",nbinstot, xMin, xMax);

    if (parInput.is_open()) {
	while ( parInput.good() ) {
	    parInput >> x >> y;   
	    double s_y = sqrt(double(y));
	    int theBin = histoTot->FindBin(x);
	    histoTot->SetBinContent(theBin,double(y));
	    histoTot->SetBinError(theBin,s_y);
	}
	parInput.close();
    }
    TCanvas* ctot = new TCanvas("c1","c1",20,20,800,600);
    ctot->SetGrid();
    gStyle->SetOptStat(0);
    histoTot->Rebin(4);
    histoTot->Draw("e1");
   
    double minfit = 4000;
    double maxfit = 10500;
    TF1* f1 = new TF1("f1","pol2",minfit,maxfit);
    histoTot->Fit("f1","","e1",minfit,maxfit);
    std::cout << "Chi^2:" <<f1->GetChisquare() << ", number of DoF: " << f1->GetNDF();
    std::cout << " (Probability: " << f1->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    double a = f1->GetParameter(2);
    double s_a = f1->GetParError(2);
    double b = f1->GetParameter(1);
    double s_b = f1->GetParError(1);
    
    double max = -b/a*0.5;
    double s_max = max *sqrt(pow(s_a/a,2)+pow(s_b/b,2));
    std::cout << "CHN_max = ("<<max<<" +/- "<<s_max<<")CHN\n";
    double  k = 2.4e-5; //MeV/CHN
    double s_k = 5e-6;
    double T_max = k*max; //MeV
    double s_T_max = T_max * sqrt(pow(s_k/k,2)+pow(s_max/max,2));
    std::cout << "T_max = ("<<T_max<<" +/- "<<s_T_max<<")MeV\n";
}
