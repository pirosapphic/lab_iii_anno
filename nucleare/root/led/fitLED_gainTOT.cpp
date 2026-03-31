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

//LINEARITÀ di Deltapp e Gampl
//deltapp = a*Gampl + b, con b=0
			//   a=Gsipm*e/k	
//Vbias = 55 V cost
//I = 2.5 cost
//Gain = 25-35 dB a step di 2

void fitLED_gainTOT(){
    int n = 6;
    std::vector<double> G = {25.,27.,29.,31.,33.,35.}; //dB
    std::vector<double> s_G(6); //1 dB
    for(int i = 0; i<n; i++){
	s_G[i] = 20/log(10)/G[i]*1;
	G[i] = pow(10,G[i]/20);
    }
    std::vector<double> delta = {151.28,185.90,231.68,293.53,369.89,469.58}; //CHN
    std::vector<double> s_delta = {0.04,0.04,0.04,0.05,0.05,0.07};
    
    TCanvas* c1 = new TCanvas("c1","c1",20,20,1098,732);
    c1->SetGrid();
    TGraphErrors * g1 = new TGraphErrors(n,G.data(),delta.data(),s_G.data(),s_delta.data());
    g1->SetTitle("<#Delta_{pp}> vs G_{ampl}, con V_{bias}=55V e intensita' 2.5;G_{ampl} [#];<#Delta_{pp}> [CHN]");
    g1->Draw("AP");
    TF1* f1 = new TF1("f1","pol1",15,60);
    g1->Fit("f1","R+","e1");
        std::cout << "Chi^2:" <<f1->GetChisquare();
    std::cout<< ", number of DoF: " << f1->GetNDF();
    std::cout << " (Probability: " << f1->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;

    std::cout<<"Il termine noto è compatibile con 0\n";
    std::cout<<"Il coeff angolare va confrontato con Gsipm*e/k, con k = 40fC/CHN\n";
    std::cout<<"uso la stima già fatta di Gsipm (con Vbias = 55V) per il coeff 'teorico'\n";
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

}
