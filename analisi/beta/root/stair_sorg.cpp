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


TGraph* graph_filler(string title, string path); //general purpose

void stair_sorg(){
    string stairpath = "../data_sorg/A8_staircase_pt2_-6.txt";

    TCanvas* cst = new TCanvas("cs","cs",20,20,1000,800);
    cst->SetGrid();
    TGraph* gst = graph_filler("Staircase plot;V_{thr}[mV];#nu[kHz]",stairpath);
    cst->SetLogy();
    gst->SetMarkerStyle(21);
    gst->Draw("APRX");
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
