 //usa covarianze
 
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

TGraphErrors* graph_errors_fillers(string title, string path); //general purpose
void prova_cazzo(){
    string path_1="../data_fc/data_V_0.txt";
    TGraphErrors*e1 = graph_errors_fillers("Tensione di fondo;lunghezza d'onda [ua];V_{0,pp} [mV]",path_1);
    TCanvas* c1 = new TCanvas("c1","c1",20,20,1098,732);
    c1->SetGrid();
    e1->SetMarkerStyle(20);
    e1->SetMarkerSize(1.2);
    e1->Draw("APE");
    c1->Update();
    string path_2="../data_fc/data_V_trasm.txt";
    TGraphErrors*e2 = graph_errors_fillers("Tensione trasmessa;lunghezza d'onda [ua];V_{pp} [mV]",path_2);
    TCanvas* c2 = new TCanvas("c2","c2",20,20,1098,732);
    c2->SetGrid();
    e2->SetMarkerStyle(20);
    e2->SetMarkerSize(1.2);
    e2->Draw("APE");
    c2->Update();

    //trasmittanza
    

}

TGraphErrors* graph_errors_fillers(string title, string path){

    ifstream in_file(path.c_str());

    if(!in_file.is_open()){
        std::cout << "ERRORE: file non aperto\n";
        return nullptr;
    }

    double x,y,s_x,s_y;

    vector<double> xvec, yvec, s_xvec, s_yvec;

    while(in_file >> x >> s_x >> y >> s_y){

        xvec.push_back(x);
        s_xvec.push_back(s_x);
        yvec.push_back(y);
        s_yvec.push_back(s_y);
    }

    TGraphErrors* g =
        new TGraphErrors(
            xvec.size(),
            xvec.data(),
            yvec.data(),
            s_xvec.data(),
            s_yvec.data()
        );

    g->SetTitle(title.c_str());

    return g;
}

