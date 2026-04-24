/*Non ho capito cosè I? Devo dividere vpp per R?
 nel dubbio io faccio come il gruppo D che i grafici li ha fati con V0 */

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

//Vbias = 55V
//Gampl = 29dB
//I = 2.8




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



void trasmittanza(string title1= "Trasmittanza",string path1="../data_fc/data_V_0.txt",string title2="Tensione",string path2="../data_fc/data_V_trasm.txt"){

    TGraphErrors*g1 = graph_errors_fillers("Trasmittanza","../data_fc/data_V_0.txt");
    TGraphErrors*g2 = graph_errors_fillers("Tensione", "../data_fc/data_V_trasm.txt");
    TCanvas* c = new TCanvas("c","c",900,600);

    // --- pad principale (sinistra: trasmittanza)
    TPad* p1 = new TPad("p1","p1",0,0,1,1);
    p1->SetTicks();
    p1->Draw();
    p1->cd();

    g1->SetMarkerStyle(20);
    g1->SetLineColor(kBlue);
    g1->Draw("AP");

    g1->GetXaxis()->SetTitle("Lunghezza d’onda (nm)");
    g1->GetYaxis()->SetTitle("Trasmittanza");

    c->cd();

    // --- secondo pad trasparente (destra: fotocorrente)
    TPad* p2 = new TPad("p2","p2",0,0,1,1);
    p2->SetFillStyle(4000); // trasparente
    p2->SetFrameFillStyle(0);
    p2->Draw();
    p2->cd();

    g2->SetMarkerStyle(21);
    g2->SetLineColor(kRed);
    g2->Draw("P SAME");

    // asse destro manuale
    Double_t ymin = g2->GetYaxis()->GetXmin();
    Double_t ymax = g2->GetYaxis()->GetXmax();

    TGaxis* axis = new TGaxis(
        1,0,1,1,
        ymin,ymax,
        510,"+L"
    );

    axis->SetTitle("Fotocorrente");
    axis->SetLineColor(kRed);
    axis->SetLabelColor(kRed);
    axis->Draw();
}

