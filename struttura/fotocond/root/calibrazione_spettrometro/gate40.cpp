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
#include "TFitResult.h"

TH1D* histo_filler(string name, string title, string path); //funzione che cre oggetto istogramma da graficare
TGraphErrors* graph_errors_filler(string title, string path); //general purpose
std::vector<double> w_mean(std::vector<double> val, std::vector<double> s_val);

void gate40(){
    string path = "../../data_fc/calib_spettr/gate40.txt";
    TH1D* h1 = histo_filler("h1","titolo;assex;assey",path); // oggetto istogramma da graficare
    TFitResultPtr fitres; 
    TCanvas* c1 = new TCanvas("c1","c1",20,20,1098,732);
    c1->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h1->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    //nono picco
    fitres = h1->Fit("gaus","RS+","e1",693.,697.);
    std::cout<<"p-value "<<fitres->Prob()<<std::endl;
    cout<< "stima picco 9 (venuto): " << fitres->Parameter(1) << " pm " <<  fitres->Parameter(2)*2.355 <<endl;
    //sesto picco
    fitres = h1->Fit("gaus","RS+","e1",498.,502.);
    std::cout<<"p-value "<<fitres->Prob()<<std::endl;
    cout<< "stima picco 6 (venuto): " << fitres->Parameter(1) << " pm " <<  fitres->Parameter(2)*2.355 <<endl;
    //decimo picco
    fitres = h1->Fit("gaus","RS+","e1",1007.,1010.);
    std::cout<<"p-value "<<fitres->Prob()<<std::endl;
    cout<< "stima picco 10 (venuto): " << fitres->Parameter(1) << " pm " <<  fitres->Parameter(2)*2.355 <<endl;
    //19-esimo picco
    fitres = h1->Fit("gaus","RS+","e1",1666.,1669.);
    std::cout<<"p-value "<<fitres->Prob()<<std::endl;
    cout<< "stima picco 19 (non viene): " << fitres->Parameter(1) << " pm " <<  fitres->Parameter(2)*2.355 <<endl;
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
	double s_y = 0;
	if (yvec[i]>0) s_y=sqrt(yvec[i]);
	int theBin = theHisto->FindBin(xvec[i]);
        theHisto->SetBinContent(theBin,yvec[i]);
        theHisto->SetBinError(theBin,s_y);
	//std::cout<<i<<" "<<xvec[i]<<" "<<yvec[i]<<std::endl; //debug!
    }
    return theHisto;
}

TGraphErrors* graph_errors_filler(string title, string path){ //general purpose
    //fills grapherrors with title title with data from path.
    ifstream in_file(path.c_str());
    double x = 666;
    double y = 777;
    double s_x = 888;
    double s_y = 999;
    std::vector<double> xvec;
    std::vector<double> s_xvec;
    std::vector<double> yvec;
    std::vector<double> s_yvec;
    while(in_file.good()){
	in_file >> x >> s_x >> y >> s_y;
	xvec.push_back(double(x));
	yvec.push_back(double(y));
	s_xvec.push_back(double(s_x));
	s_yvec.push_back(double(s_y));
	std::cout<<x<<"pm"<<s_x<<" "<<y<<"pm"<<s_y<<std::endl; //debug!
    }
    std::cout<<"good "<<in_file.good()<<std::endl;

    xvec.pop_back();
    yvec.pop_back(); //CHECK IF IT IS WRONG!!!!!!!!!!
    s_xvec.pop_back();
    s_yvec.pop_back();
    xvec.shrink_to_fit();
    yvec.shrink_to_fit();
    s_xvec.shrink_to_fit();
    s_yvec.shrink_to_fit();
    TGraphErrors* theGraph = new TGraphErrors(xvec.size(),xvec.data(),yvec.data(),s_xvec.data(),s_yvec.data());
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

