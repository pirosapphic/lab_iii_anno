
#include <fstream>
#include <math.h>   
#include <iostream>

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"

using namespace std;

TH1F* istogramma(string input,string name) {
  
  // file di test contiene:
  //
  // 1    4
  // 2    5
  // 3    6
  static const int nBins = 20;
  float x[nBins];  int y[nBins]; 
  ifstream parInput(input.c_str());
  int i = 0;

  if (parInput.is_open()) {
    while ( parInput.good() ) {
      parInput >> x[i] >> y[i];   
      i++;   if (i == nBins) break;
    }
    parInput.close();
  }

  // come si riempie un TH1
  TH1F *myHisto = new TH1F("myHisto","myHisto",nBins,x[0],x[nBins-1]);  
  for (int j=0;j<nBins;j++) {
    myHisto->SetBinContent(j+1,y[j]);
    myHisto->SetBinError(j+1, sqrt(y[j]));
  }

  // come si accorpano gli eventi in classi adiacenti
  //myHisto->Rebin(4);
  // ora myHisto ha nBins/4 classi e ognuna contiene gli eventi di 4 adiacenti
  return myHisto;
}

TF1* FitGauss(TH1F* h, double minFit, double maxFit) {
    if (!h) return nullptr;

    // funzione gaussiana
    TF1* f = new TF1("f","gaus",minFit,maxFit);

    // parametri iniziali
    double maxBin = h->GetMaximum();
    int binMax = h->GetMaximumBin();
    double meanGuess = h->GetBinCenter(binMax);
    double sigmaGuess = (maxFit - minFit)/6.; // stima sigma

    f->SetParameters(maxBin, meanGuess, sigmaGuess);
    h->Fit(f,"R"); // R = range limit
    cout<< endl<< "p-value"<< f->GetProb() << endl;
    return f;
}

void ass_all(){
     // Array dei file
    string files[] = {
        "../data_sorg/A8_count_vuoto.txt",
        "../data_sorg/A8_count_all3.txt",
        "../data_sorg/A8_count_all6.txt",
        "../data_sorg/A8_count_all9.txt",
        "../data_sorg/A8_count_all12.txt",
        "../data_sorg/A8_count_all15.txt",
        "../data_sorg/A8_count_all18.txt",
        "../data_sorg/A8_count_all20.txt"
    };

    float xMin[] = {647.275,588.475,526.975,511.15,478.625,425.45,428.025,399.85};
    float xMax[] = {771.725,720.525,678.025,630.85,578.375,556.55,542.975,508.15};

    int colors[] = {kCyan,kGreen,kBlue,kRed,kYellow,kMagenta,kOrange,kBlack};

    TH1F* histos[8];
    TF1* fits[8];

    for (int i=0; i<8; i++) {
        histos[i] = istogramma( files[i], Form("h%d",i));

        if (i == 2) {   //terzo grafico
            histos[i]->Rebin(2);  
        }

        fits[i] = FitGauss(histos[i], xMin[i], xMax[i]);
    }
    double media[8]= {};
    double erroremedia[8]= {};
    // Creazione canvas separati e disegno
    for (int i=0; i<8; i++) {
        TCanvas* c = new TCanvas(Form("c%d",i), Form("c%d",i), 800, 600);

        histos[i]->SetTitle(Form("Spettro Sr90 con %d fogli ", i*3)); 
        if (i == 7) {   
            histos[i]->SetTitle("Count Rate Sr90 20"); 
        }
          // titolo
        histos[i]->GetXaxis()->SetTitle("Number of pulses");    // asse X
        histos[i]->GetYaxis()->SetTitle("Conteggi");         // asse Y
        histos[i]->Draw("e1");
        if (fits[i]) {
            fits[i]->SetLineColor(colors[i]);
            fits[i]->Draw("same");
            //recupera parametri per grafico
            media[i]=fits[i]->GetParameter(1);
            erroremedia[i]= fits[i]->GetParError(1);// l'errore sulla media è troppo piccolo, ma con sigma il x^2 è troppo piccolo
            
        }
    }
    //grafico 

    //rates

    float rates[8] = {};
    float errrates[8] = {};

    for (int i = 0; i < 8; i++)
    {
        rates[i]= media[i]*2;
        errrates[i]= erroremedia[i]*2;
    }

    auto c8 = new TCanvas("c8","Grafico assorbimento alluminio",200,10,700,500);
    c8->SetFillColor(42);
    c8->SetGrid();
    c8->GetFrame()->SetFillColor(21);
    c8->GetFrame()->SetBorderSize(12);
    
    const Int_t n = 8;
   
    //x del grafico
    float ascissa[8] = {0,3,6,9,12,15,18,20};
    float errorex[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //non credo ci sia errore


    auto gr = new TGraphErrors(n,ascissa,rates,errorex,errrates);
    gr->SetTitle("Assorbimento alluminio");
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->GetXaxis()->SetTitle("Numero di fogli");    // asse X
    gr->GetYaxis()->SetTitle("Rates"); 
    gr->Draw("AP");

   //fit esponenziale

    TF1 *f = new TF1("f","[0] +[1]* exp(-x*[2])", 0, 20); 

    f->SetParameters(900,500,1); 


    gr->Fit(f,"R");
    cout<< endl<< "p-value"<< f->GetProb() << endl;

}