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



TGraphErrors* graph_errors_fillers(string title, string path);
tuple<
vector<double>,
vector<double>,
vector<double>,
vector<double>
>lettura(string path);

void trasmittanza(){

    string path1= "../data_fc/dati_V0_calib.txt";
    string path2= "../data_fc/dati_Vt_calib.txt";
    auto [a,s_a, V_0,sV_0]= lettura(path1);
    auto [b, s_b,V,sV]= lettura(path2);
    vector<double> trasm, s_trasm;
    double x,s_x;
    for (int i = 0; i < a.size(); i++)
    {
        x = V[i]/V_0[i];
        s_x= x*sqrt(pow(sV_0[i]/V_0[i] ,2)+pow(sV[i]/V_0[i],2));
        trasm.push_back(x);
        s_trasm.push_back(s_x);
        cout<<trasm[i] << "pm" <<s_trasm[i] << endl;
        
    }

    
    //auto c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
    //c1->SetGrid();

    /*const int n= a.size();
    auto gr = new TGraphErrors(n,a.data(),trasm.data(),s_a.data(),s_trasm.data());
    gr->SetTitle("Trasmittanza;lunghezza d'onda [nm],T[#]");
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);*/
    
    const int n= a.size();
    string path_3="../data_fc/dati_I_calib.txt";
    TGraphErrors*g2 = graph_errors_fillers("Fotocorrente;lunghezza d'onda [nm];I_f [mA]",path_3);
    TGraphErrors *g1 = new TGraphErrors(n,a.data(),trasm.data(),s_a.data(),s_trasm.data());


    TCanvas *c = new TCanvas("c","comp",200,10,700,500);

    TPad *p1 = new TPad("p1","",0,0,1,1);
    p1->SetGrid();
    p1->Draw();
    p1->cd();

    g1->SetMarkerStyle(7);
    g1->SetMarkerColor(kBlack);
    g1->Draw("ALP");

    gPad->Update();

    // range asse X
    double xmin = g1->GetXaxis()->GetXmin();
    double xmax = g1->GetXaxis()->GetXmax();

    // range Y seconda curva
    double ymin = TMath::MinElement(n, trasm.data());
    double ymax = TMath::MaxElement(n, trasm.data());

    TPad *p2 = new TPad("p2","",0,0,1,1);
    p2->SetFillStyle(4000);
    p2->Draw();
    p2->cd();

    g2->SetMarkerStyle(7);
    g2->SetMarkerColor(kRed);
    g2->Draw("LP");

    TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax,
                          ymin, ymax, 510, "+L");

    axis->SetLineColor(kRed);
    axis->SetLabelColor(kRed);
    axis->Draw();

    TLegend *leg = new TLegend(0.75,0.75,0.9,0.9);
    leg->AddEntry(g1,"gr1","lp");
    leg->AddEntry(g2,"gr2","lp");
    leg->Draw();

    c->cd();
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



tuple<
vector<double>,
vector<double>,
vector<double>,
vector<double>
>
lettura(string path){

    ifstream in_file(path);

    double x,y,s_x,s_y;

    vector<double> xvec,yvec,s_xvec,s_yvec;

    while(in_file >> x >> s_x >> y >> s_y){
        xvec.push_back(x);
        s_xvec.push_back(s_x);
        yvec.push_back(y);
        s_yvec.push_back(s_y);
    }

    return make_tuple(xvec,s_xvec,yvec,s_yvec);
}
