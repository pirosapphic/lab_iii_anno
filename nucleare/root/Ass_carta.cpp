
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

// Funzione generica per fare fit gaussiano
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


// Funzione principale
void Ass_carta() {

    // Array dei file
    string files[] = {
        "../data_sorg/A8_count_vuoto.txt",
        "../data_sorg/A8_count_carta3.txt",
        "../data_sorg/A8_count_carta6.txt",
        "../data_sorg/A8_count_carta9.txt",
        "../data_sorg/A8_count_carta12.txt",
        "../data_sorg/A8_count_carta15.txt",
        "../data_sorg/A8_count_carta18.txt",
        "../data_sorg/A8_count_carta20.txt"
    };

    // Range dei fit e dei bin per ogni file
    float xMin[] = {647.275,518.925,444.1,365.525,317.075,289.6,253.325,230.325};
    float xMax[] = {771.725,668.075,561.9,499.475,433.925,388.4,341.675,318.675};
   // float binWidth[] = {6.55,7.85,6.2,7.05,6.15,5.2,4.65,4.65};

    // Colori dei grafici
    int colors[] = {kCyan,kGreen,kBlue,kRed,kYellow,kMagenta,kOrange,kBlack};

    // Creazione degli istogrammi e dei fit
    TH1F* histos[8];
    TF1* fits[8];

    for (int i=0; i<8; i++) {
        histos[i] = istogramma( files[i], Form("h%d",i));
        fits[i] = FitGauss(histos[i], xMin[i], xMax[i]);
    }
    double media[8]= {};
    double erroremedia[8]= {};
    // Creazione canvas separati e disegno

    for (int i=0; i<8; i++) {
            TCanvas* c = new TCanvas(Form("c%d",i), Form("c%d",i),20,20,1098,732);
            c->SetGrid();
            histos[i]->SetTitle(Form("Spettro ^{90}Sr con %d fogli di carta ", i*3)); 
            if (i == 7) {   
                histos[i]->SetTitle("Count Rate ^{90}Sr 20"); 
            }
          // titolo
            histos[i]->GetXaxis()->SetTitle("N. impulsi/gate [x(500 ms)^{-1}]");    // asse X
            histos[i]->GetYaxis()->SetTitle("Conteggi [#]");         // asse Y
            histos[i]->Draw("e1");
            if (fits[i]) {
                fits[i]->SetLineColor(colors[i]);
                fits[i]->Draw("same");
                //recupera parametri per grafico
                media[i]=fits[i]->GetParameter(1);
                erroremedia[i]= fits[i]->GetParameter(2);
            
            }
    }

    

    //spessore carta 2pm0.05 mm per 20 fogli 
    double Sp1= 2.0/20.0 *1./10.;//cm
    double Sper1= 0.05/20.0 *1./10;//mm
    double Sp[8] = {0, Sp1*3.,Sp1*6.,Sp1*9.,Sp1*12.,Sp1*15.,Sp1*18.,Sp1*20.};
    double Sper[8]= {Sper1,Sper1*3.,Sper1*6.,Sper1*9.,Sper1*12.,Sper1*15.,Sper1*18.,Sper1*20.};
    //l'errore sullo spessore è variabile

    //rates

    double rates[8] = {};
    double errrates[8] = {};

    for (int i = 0; i < 8; i++)
    {
        rates[i]= media[i]*2;
        errrates[i]= erroremedia[i]*2;
        cout<< "rates"<< rates[i] <<"+- "<<errrates[i]<<endl;
        cout<< "spessore"<< Sp[i] <<"+- "<<Sper[i]<<endl;
    }
    
   auto c8 = new TCanvas("c8","Grafico assorbimento carta",20,20,1098,732);

   c8->SetGrid();
   c8->GetFrame()->SetFillColor(21);
   c8->GetFrame()->SetBorderSize(12);
   const Int_t n = 8;
   
   auto gr = new TGraphErrors(n,Sp,rates,Sper,errrates);
   gr->SetTitle("Assorbimento carta");
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->GetXaxis()->SetTitle("Spessore [mm]");    // asse X
   gr->GetYaxis()->SetTitle("Rates [Hz]"); 
   gr->Draw("AP");

   //fit esponenziale

   TF1 *f = new TF1("f","[0] +[1]* exp(-[2] * x)", 0, 2);

   f->SetParameters(200,1000,7); 
   gr->Fit(f,"R");
   cout<< endl<< "p-value"<< f->GetProb() << endl;
    //il coefficente di assorbimento viene super basso, circa 0,8
    //da test Z è però compatibile con quello teorico, il X2 va bene, per me si può fare
    
    

    
}