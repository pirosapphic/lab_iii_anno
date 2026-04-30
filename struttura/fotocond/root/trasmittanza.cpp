#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <cmath>

#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TMath.h"

using namespace std;

tuple<
vector<double>,
vector<double>,
vector<double>,
vector<double>
> lettura(string path);


TGraphErrors* graph_errors_fillers(string title,string path);


void trasmittanza(){

    //-----------------------------------
    // LETTURA DATI TRASMITTANZA
    //-----------------------------------

    string path1="../data_fc/dati_V0_calib.txt";
    string path2="../data_fc/dati_Vt_calib.txt";

    auto [lambda,s_lambda,V0,sV0]=lettura(path1);
    auto [dummy1,dummy2,Vt,sVt]=lettura(path2);

    vector<double> T;
    vector<double> sT;

    for(int i=0;i<lambda.size();i++){

        double t = Vt[i]/V0[i];

        // propagazione errori corretta
        double st =
        t*sqrt(
          pow(sVt[i]/Vt[i],2)
        + pow(sV0[i]/V0[i],2)
        );

        T.push_back(t);
        sT.push_back(st);
    }


    //-----------------------------------
    // GRAFICI
    //-----------------------------------

    int n=lambda.size();

    // g1 = trasmittanza
    TGraphErrors *g1 =
    new TGraphErrors(
      n,
      lambda.data(),
      T.data(),
      s_lambda.data(),
      sT.data()
    );

    // g2 = fotocorrente
    TGraphErrors *g2=
    graph_errors_fillers(
      "Fotocorrente",
      "../data_fc/dati_I_calib.txt"
    );


    //-----------------------------------
    // CANVAS + PAD
    //-----------------------------------

    TCanvas *c =
    new TCanvas("c","Trasmittanza",900,600);

    TPad *p1=
    new TPad("p1","",0,0,1,1);

    TPad *p2=
    new TPad("p2","",0,0,1,1);

    p1->SetGrid();
    p2->SetFillStyle(4000); // trasparente


    //-----------------------------------
    // PRIMO GRAFICO
    //-----------------------------------

    p1->Draw();
    p1->cd();

    g1->SetTitle(
    "Trasmittanza e Fotocorrente;#lambda [nm];T [#]"
    );

    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(kBlue);
    g1->SetLineColor(kBlue);
    TF1 *fT = new TF1("fT","[0] + [1]/(1+exp((x-[2])/[3]))",600,650);

    fT->SetParameters(0.0,0.18,630,5); //0.0,0.18,630,5
    fT->SetLineColor(kBlue);
    // b, A, lambda0, DeltaLambda

    g1->Fit(fT,"R"); // R = usa solo il range dato
    g1->Draw("AP");
    fT->Draw("same");
    double probT= fT->GetProb();
    cout<< probT<<endl;
    //calcolo energy gap 

    double lambda0T= fT-> GetParameter(2);
    double s_lambda0T= fT-> GetParameter(3);
    double EgT  = 1240./lambda0T;
    double s_EgT = 1240.*s_lambda0T/pow(lambda0T,2);
    cout << EgT << " +/- " << s_EgT << " eV" << endl;
   
    gPad->Update();


    //-----------------------------------
    // FONT ASSI
    //-----------------------------------

    Style_t tfont=
    g1->GetHistogram()->GetYaxis()->GetTitleFont();

    Float_t tsize=
    g1->GetHistogram()->GetYaxis()->GetTitleSize();

    Style_t lfont=
    g1->GetHistogram()->GetYaxis()->GetLabelFont();

    Float_t lsize=
    g1->GetHistogram()->GetYaxis()->GetLabelSize();


    //-----------------------------------
    // RANGE PER ASSE DESTRO
    //-----------------------------------

    Double_t xmin=p1->GetUxmin();
    Double_t xmax=p1->GetUxmax();

    Double_t dx=(xmax-xmin)/0.8;

    Double_t ymin=
    g2->GetHistogram()->GetMinimum();

    Double_t ymax=
    g2->GetHistogram()->GetMaximum();

    Double_t dy=(ymax-ymin)/0.8;


    //-----------------------------------
    // SECONDO PAD
    //-----------------------------------

    p2->Range(
      xmin-0.1*dx,
      ymin-0.1*dy,
      xmax+0.1*dx,
      ymax+0.1*dy
    );

    p2->Draw();
    p2->cd();

    g2->SetMarkerStyle(21);
    g2->SetMarkerColor(kRed);
    g2->SetLineColor(kRed);
    
    TF1 *fI = new TF1("fI","[0] + [1]/(1+exp((x-[2])/[3]))",625,680);
    fI->SetParameters(0.1,6,650,10); //0.1,5,650,10

    g2->Fit(fI,"R");
    g2->Draw("P");
    fI->SetLineColor(kRed);
    fI->Draw("same");
    double probI= fI->GetProb();
    cout<< probI<<endl;

    //calcolo energy gap
    double lambda0I= fI-> GetParameter(2);
    double s_lambda0I= fI-> GetParameter(3);
    double EgI  = 1240./lambda0I;
    double s_EgI = 1240.*s_lambda0I/pow(lambda0I,2);
    
    cout << EgI << " +/- " << s_EgI << " eV" << endl;
    


    //-----------------------------------
    // ASSE DESTRO
    //-----------------------------------

    TGaxis *axis=
    new TGaxis(
      xmax,ymin,
      xmax,ymax,
      ymin,ymax,
      510,"+L"
    );

    axis->SetTitle("Fotocorrente [mA]");

    axis->SetTitleFont(tfont);
    axis->SetTitleSize(tsize);
    axis->SetLabelFont(lfont);
    axis->SetLabelSize(lsize);

    axis->SetTitleColor(kRed);
    axis->SetLabelColor(kRed);
    axis->SetLineColor(kRed);

    axis->Draw();


    //-----------------------------------
    // LEGENDA
    //-----------------------------------

    p2->cd();

    TLegend *leg=
    new TLegend(0.55,0.75,0.85,0.88);

    leg->AddEntry(g1,"Trasmittanza","lp");
    leg->AddEntry(g2,"Fotocorrente","lp");

    leg->Draw("same");


    //test Z

    double Z= (EgI-EgT)/sqrt(pow(s_EgI,2)+pow(s_EgT,2));
    cout<<"Z: "<< Z<<endl;

    //media pesata

    double Eg=(EgI/pow(s_EgI,2)+EgT/pow(s_EgT,2))/(1./pow(s_EgT,2)+1./pow(s_EgI,2));
    double s_Eg=pow(1./pow(s_EgT,2)+1./pow(s_EgI,2),-0.5);
    cout<<"Eg medio: "<< Eg<< "pm "<< s_Eg<<endl;
}



TGraphErrors*
graph_errors_fillers(string title,string path){

    ifstream in(path);

    double x,sx,y,sy;

    vector<double> xv,sxv,yv,syv;

    while(in>>x>>sx>>y>>sy){

        xv.push_back(x);
        sxv.push_back(sx);

        yv.push_back(y);
        syv.push_back(sy);
    }

    auto g=
    new TGraphErrors(
       xv.size(),
       xv.data(),
       yv.data(),
       sxv.data(),
       syv.data()
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

    ifstream in(path);

    double x,sx,y,sy;

    vector<double> xv,sxv,yv,syv;

    while(in>>x>>sx>>y>>sy){

        xv.push_back(x);
        sxv.push_back(sx);

        yv.push_back(y);
        syv.push_back(sy);
    }

    return make_tuple(
      xv,sxv,yv,syv
    );
}