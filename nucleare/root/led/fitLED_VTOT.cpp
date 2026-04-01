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
    std::vector<double> Gsipm = {710.e3,1.39e6,2.05e6,2.71e6,3.36e6};
    std::vector<double> s_Gsipm = {90.e3,8.e4,10.e4,11.e4,13e4};

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
/*
    std::cout<<"Il termine noto è compatibile con 0\n";
    std::cout<<"Il coeff angolare va confrontato con Gsipm*e/k, con k = 40fC/CHN\n";
    std::cout<<"uso la stima già fatta di Gsipm (con Vbias = 55V) per il coeff 'teorico'";
    double Gsipm = 2050000;
    double s_Gsipm = 20000;
    double e = 1.6022e-4; //fC
    double k = 40; //fC/CHN
    double mteo = Gsipm*e/k;
    double s_mteo = s_Gsipm*e/k;
    double msper = f1->GetParameter(1);
    double s_msper = f1->GetParError(1);
    std::cout<<"coeff sper = ("<<f1->GetParameter(1)<<" +/- "<<f1->GetParError(1)<<")CHN"<<std::endl;
    std::cout<<"coeff teo = ("<<mteo<<" +/- "<<s_mteo<<")CHN\n";
    std::cout<<"Z score = "<< (mteo-msper)/sqrt(s_mteo*s_mteo + s_msper*s_msper)<<std::endl;
*/
}
