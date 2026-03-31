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


TH1D* histo_filler(string name, string title, string path); //general purpose
TGraph* graph_filler(string title, string path); //general purpose
std::vector<double> w_mean(std::vector<double> val, std::vector<double> s_val);

void fit_DCR23_3(){
    string temp = "23_3";
    string T = "(23.4#pm0.1)";
    string stairpath = "../../data_SiPM/DCR_temp/stair";
    string DCR = "../../data_SiPM/DCR_temp/dcr";

    TCanvas* cst = new TCanvas("cst","cst",20,20,1098,732);
    cst->SetGrid();
    TGraph* gst = graph_filler("Staircase plot, G_{ampl}=28dB, V_{BIAS}=55V, T="+T+"^{o}C;V_{thr}[mV];#nu[kHz]",stairpath+temp+"c.txt");

    cst->SetLogy();
    gst->SetMarkerStyle(21);
    gst->Draw("APRX");

//------------------------------------------------------------------------------- 
//			LOW
    TCanvas* c17 = new TCanvas("c17","c17",20,20,1098,732);
    c17->SetGrid();
    c17->cd();
    TH1D* h17 = histo_filler("h17","Impulsi ricevuti a V_{thr}=-10mV con gate = 500ms, T="+T+"^{o}C;N. impulsi/gate [x500 ms^{-1}];Conteggi [#]",DCR+"low_"+temp+"C.txt");

    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
    h17->Draw("HIST E1");
    TF1* gaus1 = new TF1("gaus1","gaus",9650,10100);
    h17->Fit("gaus1","R+","e1");
    gaus1->Draw("same");
    std::cout << "Chi^2:" <<gaus1->GetChisquare();
    std::cout<< ", number of DoF: " << gaus1->GetNDF();
    std::cout << " (Probability: " << gaus1->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;


//------------------------------------------------------------------------------- 
//			HIGH
    TCanvas* c44 = new TCanvas("c44","c4",20,20,1098,732);
    c44->SetGrid();
    c44->cd();
    TH1D* h44 = histo_filler("h44","Impulsi ricevuti a V_{thr}=-24mV con gate = 500ms, T="+T+"^{o}C;N. impulsi/gate [x500 ms^{-1}];Conteggi [#]",DCR+"high_"+temp+"C.txt");

    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
    h44->Draw("HIST E1");
    TF1* gaus2 = new TF1("gaus2","gaus",120,170);
    h44->Fit("gaus2","R+","e1");
    gaus2->Draw("same");
    std::cout << "Chi^2:" <<gaus2->GetChisquare();
    std::cout<< ", number of DoF: " << gaus2->GetNDF();
    std::cout << " (Probability: " << gaus2->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    double nu1 = gaus1->GetParameter(1)*2; //Hz
    double s_nu1 = gaus1->GetParameter(2)*2;
    double nu2 = gaus2->GetParameter(1)*2;
    double s_nu2 = gaus2->GetParameter(2)*2;
    std::cout<<"Dark Count rate @ -17mV = ("<<nu1<<" +/- "<<s_nu1<<")Hz\n";
    std::cout<<"Dark Count rate @ -44mV = ("<<nu2<<" +/- "<<s_nu2<<")Hz\n";
    double oct = nu2/nu1 *100;//%
    double s_oct = oct*sqrt(pow(s_nu1/nu1,2)+pow(s_nu2/nu2,2));
    std::cout<<"Optical cross talk = nu2/nu1 = ("<<oct<<" +/- "<<s_oct<<")%"<<std::endl;


}

TH1D* histo_filler(string name, string title, string path){ //general purpose
    //fills histo with name name and title title with data from path.
    ifstream in_file(path.c_str());
    double x = 666;
    double y = 777;
    std::vector<double> xvec;
    std::vector<double> yvec;
    while(in_file.good()){
	in_file >> x >> y;
//	std::cout<<x<<" "<<y<<std::endl;
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
    }
    return theHisto;
}

TGraph* graph_filler(string title, string path){ //general purpose
    //fills histo with name name and title title with data from path.
    ifstream in_file(path.c_str());
    double x = 666;
    double y = 777;
    std::vector<double> xvec;
    std::vector<double> yvec;
    while(in_file.good()){
	in_file >> x >> y;
//	std::cout<<x<<" "<<y<<std::endl;
	xvec.push_back(double(x));
	yvec.push_back(double(y));
    }
    xvec.pop_back();
    yvec.pop_back(); //CHECK IF IT IS WRONG!!!!!!!!!!
    xvec.shrink_to_fit();
    yvec.shrink_to_fit();
    TGraph* theGraph = new TGraph(xvec.size(),xvec.data(),yvec.data());
    theGraph->SetTitle(title.c_str());
    return theGraph;
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
