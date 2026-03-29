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


void fit_DCRT(){
    std::vector<double> T = {20.4,22.1,22.7,23.0,23.4,23.7,23.9};
    std::vector<double> s_T = {0.5,0.3,0.2,0.1,0.1,0.1,0.1}; //RIVEDI CHIEDI A MAGDA
    std::vector<double> dcrl = {16600.,17600.,18400.,19100.,19700.,20200.,20750.};	//-10mV
    std::vector<double> s_dcrl = {200.,200,300.,200.,200.,190.,160.};
    std::vector<double> dcrh = {270.,290.,276.,290.,290.,300.,293.};		//-24mV
    std::vector<double> s_dcrh = {20.,50.,12.,20.,20.,30.,19.};
    std::vector<double> oct(7);
    std::vector<double> s_oct(7);
    for(int i = 0; i < 7; i++){
	oct[i] = dcrh[i]/dcrl[i]*100; //%
	s_oct[i] = oct[i]*sqrt(pow(s_dcrh[i]/dcrh[i],2)+pow(s_dcrl[i]/dcrl[i],2));
    }
    TCanvas* cl = new TCanvas("cl","cl",20,20,800,600);
    cl->SetGrid();
    cl->cd();
    TGraphErrors* gl = new TGraphErrors(dcrl.size(),T.data(),dcrl.data(),s_T.data(),s_dcrl.data());
    gl->SetTitle("DCR(T), V_{thr} = -10mV, V_{BIAS} = 55V, G_{ampl} = 28dB;T[^{o}C];DCR[Hz]");
    gl->Draw("AP");
    TF1* f1 = new TF1("f1","pol1",18.5,24.5);
    gl->Fit("f1","R+");
    std::cout << "Chi^2:" <<f1->GetChisquare();
    std::cout<< ", number of DoF: " << f1->GetNDF();
    std::cout << " (Probability: " << f1->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;


    TCanvas* ch = new TCanvas("ch","ch",20,20,800,600);
    ch->SetGrid();
    ch->cd();
    TGraphErrors* gh = new TGraphErrors(dcrh.size(),T.data(),dcrh.data(),s_T.data(),s_dcrh.data());
    gh->SetTitle("DCR(T), V_{thr} = -24mV, V_{BIAS} = 55V, G_{ampl} = 28dB;T[^{o}C];DCR[Hz]");
    gh->Draw("AP");
    TF1* f2 = new TF1("f2","pol1",18.5,24.5);
    gh->Fit("f2","R+");
    std::cout << "Chi^2:" <<f2->GetChisquare();
    std::cout<< ", number of DoF: " << f2->GetNDF();
    std::cout << " (Probability: " << f2->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    
    TCanvas* co = new TCanvas("co","co",20,20,800,600);
    co->SetGrid();
    co->cd();
    TGraphErrors* go = new TGraphErrors(oct.size(),T.data(),oct.data(),s_T.data(),s_oct.data());
    go->SetTitle("OCT(T), V_{BIAS} = 55V, G_{ampl} = 28dB;T[^{o}C];OCT[%]");
    go->Draw("AP");
    TF1* f3 = new TF1("f3","pol1",18.5,24.5);
    go->Fit("f3","R+");
    std::cout << "Chi^2:" <<f3->GetChisquare();
    std::cout<< ", number of DoF: " << f3->GetNDF();
    std::cout << " (Probability: " << f3->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
}
