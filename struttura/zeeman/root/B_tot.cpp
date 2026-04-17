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
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TFile.h"

TGraphErrors* graph_errors_filler(string title, string path); //general purpose
std::vector<double> w_mean(std::vector<double> val, std::vector<double> s_val);

void B_tot(){
TMultiGraph* gtot = new TMultiGraph();
//salita 1
    string path = "../data/";
    path += "salita1.txt"; 
    TGraphErrors* g1 = graph_errors_filler("titolo1",path);
    g1->SetMarkerStyle(7);
    g1->SetMarkerColor("kMagenta");
    g1->SetLineColor("kMagenta");
    TCanvas* c = new TCanvas("c1","c1",20,20,1098,732);
    c->SetGrid();
    g1->Fit("pol1","R","",-0.5,10.);
    gtot->Add(g1);
//discesa 1
    path = "../data/";
    path += "discesa1.txt"; 
    TGraphErrors* g2 = graph_errors_filler("titolo2",path);
    g2->SetMarkerStyle(7);
    g2->SetMarkerColor("kBlue");
    g2->Fit("pol1","R","",-0.5,10.);
    gtot->Add(g2);
//salita 2
    path = "../data/";
    path += "salita2.txt"; 
    TGraphErrors* g3 = graph_errors_filler("titolo3",path);
    g3->SetMarkerStyle(7);
    g3->SetMarkerColor("kRed");
    g3->Fit("pol1","R","",-0.5,10.);
    gtot->Add(g3);
//discesa 2
    path = "../data/";
    path += "discesa2.txt"; 
    TGraphErrors* g4 = graph_errors_filler("titolo4",path);
    g4->SetMarkerStyle(7);
    g4->SetMarkerColor("kGreen");
    g4->Fit("pol1","R","",-0.5,10.);
    gtot->Add(g4);
    
    gtot->Draw("AP");
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
	//std::cout<<x<<"pm"<<s_x<<" "<<y<<"pm"<<s_y<<std::endl; //debug!
    }
    //std::cout<<"good "<<in_file.good()<<std::endl;

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
