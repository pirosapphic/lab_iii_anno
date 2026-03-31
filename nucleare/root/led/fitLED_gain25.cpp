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


TH1D* histo_filler(string name, string title, string path); //general purpose
std::vector<double> w_mean(std::vector<double> val, std::vector<double> s_val);

void fitLED_gain25(string input = "../../data_SiPM/LED/gain/A8_LED55"){    
    //firstly, we draw and fit the histograms, then we calculate <deltapp>;
    std::vector<double> y_deltapp(6);
    std::vector<double> s_y_deltapp(6);
//--------------------------------------------------------------------------------
    //25 dB
    string gain = "25";
    string path = input+gain+"2-5.txt";
    TH1D* histo25 = histo_filler("histo25","Spettro Energetico, G=25dB, V_{bias}=55V, LED 2.5;ADC channels [CHN];Conteggi[#]",path);
    TCanvas* c25 = new TCanvas("c25","c25",20,20,1098,732);
    c25->SetGrid();
    gStyle->SetOptStat(0);
    histo25->SetLineColor(kBlack);
    histo25->Draw("e1");
    
    std::vector<double> peak(6);
    std::vector<double> s_peak(6);
    peak = {0,0,0,0,0,0};
    s_peak = {0,0,0,0,0,0};
    std::vector<double> sigma(6);
    std::vector<double> s_sigma(6);
    sigma = {0,0,0,0,0,0};
    s_sigma = {0,0,0,0,0,0};
    //PICCO 0
    TF1* gaus0 = new TF1("gaus0","gaus",-15,15);
    histo25->Fit("gaus0","R+","e1",-15,15);
    std::cout << "Chi^2:" <<gaus0->GetChisquare();
    std::cout<< ", number of DoF: " << gaus0->GetNDF();
    std::cout << " (Probability: " << gaus0->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    peak[0] = gaus0->GetParameter(1); 
    s_peak[0] = gaus0->GetParError(1);
    sigma[0] = gaus0->GetParameter(2); 
    s_sigma[0] = gaus0->GetParError(2);
    
    //PICCO 1
    TF1* gaus1 = new TF1("gaus1","gaus",132,175);
    gaus1->SetLineColor(kYellow);
    histo25->Fit("gaus1","R+","e1");
    std::cout << "Chi^2:" <<gaus1->GetChisquare();
    std::cout<< ", number of DoF: " << gaus1->GetNDF();
    std::cout << " (Probability: " << gaus1->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    peak[1] = gaus1->GetParameter(1); 
    s_peak[1] = gaus1->GetParError(1);
    sigma[1] = gaus1->GetParameter(2); 
    s_sigma[1] = gaus1->GetParError(2);
    
    //PICCO 2
    TF1* gaus2 = new TF1("gaus2","gaus",279,325);
    gaus2->SetLineColor(kGreen);
    histo25->Fit("gaus2","R+","e1");
    std::cout << "Chi^2:" <<gaus2->GetChisquare();
    std::cout<< ", number of DoF: " << gaus2->GetNDF();
    std::cout << " (Probability: " << gaus2->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    peak[2] = gaus2->GetParameter(1); 
    s_peak[2] = gaus2->GetParError(1);
    sigma[2] = gaus2->GetParameter(2); 
    s_sigma[2] = gaus2->GetParError(2);

    //PICCO 3
    TF1* gaus3 = new TF1("gaus3","gaus",425,480);
    gaus3->SetLineColor(kCyan);
    histo25->Fit("gaus3","R+","e1");
    std::cout << "Chi^2:" <<gaus3->GetChisquare();
    std::cout<< ", number of DoF: " << gaus3->GetNDF();
    std::cout << " (Probability: " << gaus3->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    peak[3] = gaus3->GetParameter(1); 
    s_peak[3] = gaus3->GetParError(1);
    sigma[3] = gaus3->GetParameter(2); 
    s_sigma[3] = gaus3->GetParError(2);

    //PICCO 4
    TF1* gaus4 = new TF1("gaus4","gaus",575,625);
    gaus4->SetLineColor(kBlue);
    histo25->Fit("gaus4","R+","e1");
    std::cout << "Chi^2:" <<gaus4->GetChisquare();
    std::cout<< ", number of DoF: " << gaus4->GetNDF();
    std::cout << " (Probability: " << gaus4->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    peak[4] = gaus4->GetParameter(1); 
    s_peak[4] = gaus4->GetParError(1);
    sigma[4] = gaus4->GetParameter(2); 
    s_sigma[4] = gaus4->GetParError(2);

    //PICCO 5
    TF1* gaus5 = new TF1("gaus5","gaus",720,780);
    gaus5->SetLineColor(kMagenta);
    histo25->Fit("gaus5","R+","e1");
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
    //these are not compatible upon a Gauss test... but we take the mean anyways
    std::cout<<"OCCHIO, non sono compatibili con test Z!!\nMa faccio lo stesso la\n";
    std::vector<double> respp = w_mean(deltapp,s_deltapp);
    double meanpp = respp[0];
    double s_meanpp = respp[1];
    std::cout<<"Media pesata Deltapp = ("<<meanpp<<" +/- "<<s_meanpp<<")CHN\n";
    y_deltapp[0] = meanpp;
    s_y_deltapp[0] = s_meanpp;
}

TH1D* histo_filler(string name, string title, string path){ //general purpose
    //fills histo with name name and title title with data from path.
    ifstream in_file(path.c_str());
    double x = 666;
    double y = 777;
    std::vector<double> xvec;
    std::vector<double> yvec;
    int count = 0;
    while(in_file.good()){
	in_file >> x >> y;
	xvec.push_back(double(x));
	yvec.push_back(double(y));
	count++;
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
	//std::cout<<i<<" "<<xvec[i]<<" "<<yvec[i]<<std::endl; //debug!
    }
    return theHisto;
}

std::vector<double> w_mean(std::vector<double> val, std::vector<double> s_val){
    std::vector<double> res(2);
    double sum=0;
    double sumw=0;
    int n = val.size();
    for(int i = 0; i < n; i++){
	sum += val[i]*pow(s_val[i],-2);
	sumw += pow(s_val[i],-2);
    }
    res[0] = sum/sumw;
    res[1] = sqrt(1/sumw);
    return res;
}
