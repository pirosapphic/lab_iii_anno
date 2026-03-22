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

//prima stima di GSiPM,
//Gain = 29dB
//Vbias = 55V
//LED I = 2.8

TH1D* histo_filler(string name, string title, string path); //general purpose

void fitLED_GSiPM(string input = "../../data_SiPM/LED/I/A8_LED55292-8.txt"){    
 
    TH1D* histoTot = histo_filler("histotot","Energy Spectrum, G=29dB, V_{bias}=55V, LED 2.8",input);
    TCanvas* ctot = new TCanvas("c1","c1",20,20,800,600);
    ctot->SetGrid();
    //histoTot->Rebin(8);
    gStyle->SetOptStat(0);
    histoTot->SetLineColor(kBlack);
    histoTot->Draw("e1");
    
    std::vector<double> peak(6);
    std::vector<double> s_peak(6);
    std::vector<double> sigma(6);
    std::vector<double> s_sigma(6);
    
    //PICCO 0
    TF1* gaus0 = new TF1("gaus0","gaus",-15,15);
    gaus0->SetParameter(0,4000);
    gaus0->SetParameter(1,2);
    gaus0->SetParameter(2,6);
    //gaus0->Draw("same");
    histoTot->Fit("gaus0","R+","e1",-15,15);
    std::cout << "Chi^2:" <<gaus0->GetChisquare();
    std::cout<< ", number of DoF: " << gaus0->GetNDF();
    std::cout << " (Probability: " << gaus0->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    peak[0] = gaus0->GetParameter(1); 
    s_peak[0] = gaus0->GetParError(1);
    sigma[0] = gaus0->GetParameter(2); 
    s_sigma[0] = gaus0->GetParError(2);
    
    //PICCO 1
    TF1* gaus1 = new TF1("gaus1","gaus",205,260);
    gaus1->SetLineColor(kYellow);
    histoTot->Fit("gaus1","R+","e1");
    std::cout << "Chi^2:" <<gaus1->GetChisquare();
    std::cout<< ", number of DoF: " << gaus1->GetNDF();
    std::cout << " (Probability: " << gaus1->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    peak[1] = gaus1->GetParameter(1); 
    s_peak[1] = gaus1->GetParError(1);
    sigma[1] = gaus1->GetParameter(2); 
    s_sigma[1] = gaus1->GetParError(2);
    
    //PICCO 2
    TF1* gaus2 = new TF1("gaus2","gaus",426,500);
    gaus2->SetLineColor(kGreen);
    histoTot->Fit("gaus2","R+","e1");
    std::cout << "Chi^2:" <<gaus2->GetChisquare();
    std::cout<< ", number of DoF: " << gaus2->GetNDF();
    std::cout << " (Probability: " << gaus2->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    peak[2] = gaus2->GetParameter(1); 
    s_peak[2] = gaus2->GetParError(1);
    sigma[2] = gaus2->GetParameter(2); 
    s_sigma[2] = gaus2->GetParError(2);

    //PICCO 3
    TF1* gaus3 = new TF1("gaus3","gaus",650,735);
    gaus3->SetLineColor(kCyan);
    histoTot->Fit("gaus3","R+","e1");
    std::cout << "Chi^2:" <<gaus3->GetChisquare();
    std::cout<< ", number of DoF: " << gaus3->GetNDF();
    std::cout << " (Probability: " << gaus3->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    peak[3] = gaus3->GetParameter(1); 
    s_peak[3] = gaus3->GetParError(1);
    sigma[3] = gaus3->GetParameter(2); 
    s_sigma[3] = gaus3->GetParError(2);

    //PICCO 4
    TF1* gaus4 = new TF1("gaus4","gaus",875,975);
    gaus4->SetLineColor(kBlue);
    histoTot->Fit("gaus4","R+","e1");
    std::cout << "Chi^2:" <<gaus4->GetChisquare();
    std::cout<< ", number of DoF: " << gaus4->GetNDF();
    std::cout << " (Probability: " << gaus4->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    peak[4] = gaus4->GetParameter(1); 
    s_peak[4] = gaus4->GetParError(1);
    sigma[4] = gaus4->GetParameter(2); 
    s_sigma[4] = gaus4->GetParError(2);

    //PICCO 5
    TF1* gaus5 = new TF1("gaus5","gaus",1090,1210);
    gaus5->SetLineColor(kMagenta);
    histoTot->Fit("gaus5","R+","e1");
    std::cout << "Chi^2:" <<gaus5->GetChisquare();
    std::cout<< ", number of DoF: " << gaus5->GetNDF();
    std::cout << " (Probability: " << gaus5->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    peak[5] = gaus5->GetParameter(1); 
    s_peak[5] = gaus5->GetParError(1);
    sigma[5] = gaus5->GetParameter(2); 
    s_sigma[5] = gaus5->GetParError(2);
    

    std::vector<double>deltapp(5);
    std::vector<double>s_deltapp(5);
    std::cout<<"Stime di Delta_pp:\n";
    for(int i = 0; i < deltapp.size(); i++){
	deltapp[i] = peak[i+1]-peak[i];
	s_deltapp[i] = sqrt(pow(s_peak[i+1],2)+pow(s_peak[i],2));
	std::cout<<"("<<deltapp[i]<<" +/- "<<s_deltapp[i]<<")CHN\n";
    }
}

TH1D* histo_filler(string name, string title, string path){ //general purpose
    //fills histo with name name and title title with data from path.
    ifstream in_file(path.c_str());
    double x = 666;
    double y = 777;
    std::vector<double> xvec;
    std::vector<double> yvec;
    //in_file.open(path,ios::in);
    while(in_file.good()){
	in_file >> x >> y;
	xvec.push_back(double(x));
	yvec.push_back(double(y));
    }
    xvec.shrink_to_fit();
    yvec.shrink_to_fit();
    int nbins = xvec.size()-1; //dark magic!
    double binsize = xvec[2]-xvec[1];
    double xMin = xvec[0]-(binsize/2);
    double xMax = xvec[nbins-1]+(binsize/2);
    TH1D* theHisto = new TH1D(name.c_str(),title.c_str(),nbins,xMin,xMax);
    for(int i = 0; i < nbins; i++){
	double s_y = sqrt(yvec[i]);
	int theBin = theHisto->FindBin(xvec[i]);
        theHisto->SetBinContent(theBin,yvec[i]);
        theHisto->SetBinError(theBin,s_y);
	//std::cout<<i<<" "<<xvec[i]<<" "<<yvec[i]<<std::endl;
    }
    return theHisto;
}
