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

TGraphErrors* graph_errors_filler(string title, string path); //general purpose
std::vector<double> w_mean(std::vector<double> val, std::vector<double> s_val);

void B_discesa1(){
    string path = "../../data/";
    path += "discesa1.txt";
    string title = "discesa 1";
    title = "B vs I, "+title+";I[A];B[mT]";
    TGraphErrors* g = graph_errors_filler(title,path);
    g->SetMarkerStyle(7);
    TCanvas* c = new TCanvas("c1","c1",20,20,1098,732);
    c->SetGrid();
    g->Draw("AP");
    g->Fit("pol1","R+","",-0.5,10.);
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
