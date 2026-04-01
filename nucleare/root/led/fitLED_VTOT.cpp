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

//Gsipm(Vbias)
//
//I = 2.5 cost
//Gain = 29 dB cost
//Vbias = 53-57V a step di 1

void fitLED_VTOT(){ 
    int n = 5;
    std::vector<double> Gsipm = {0.71e6,1.39e6,2.1e6,2.71e6,3.36e6};
    std::vector<double> s_Gsipm = {0.09e6,0.08e6,0.1e6,0.11e6,0.13e6};

    std::vector<double> V = {53.07,53.97,55.00,56.03,57.06}; //V
    std::vector<double> s_V ={0.03,0.03,0.03,0.03,0.03};
    
    TCanvas* c1 = new TCanvas("c1","c1",20,20,1098,732);
    c1->SetGrid();
    TGraphErrors * g1 = new TGraphErrors(n,V.data(),Gsipm.data(),s_V.data(),s_Gsipm.data());
    g1->SetTitle("G_{SiPM} vs V_{bias}, con G_{amp}=29dB e intensita' 2.5;V_{bias} [V];G_{SiPM} [#]");
    g1->Draw("AP");
    
    TF1* f1 = new TF1("f1","pol1",52.5,57.5);
    g1->Fit("f1","R+","e1");
        std::cout << "Chi^2:" <<f1->GetChisquare();
    std::cout<< ", number of DoF: " << f1->GetNDF();
    std::cout << " (Probability: " << f1->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
}
