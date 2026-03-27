
#include <fstream>
#include <math.h>   
#include <iostream>

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"

using namespace std;

// Funzione per leggere file e creare istogramma
TH1F* ReadHisto(float xMin, float xMax, float binWidth, string input, string name) {
    ifstream parInput(input.c_str());
    if (!parInput.is_open()) {
        cout << "Errore: file " << input << " non aperto!" << endl;
        return nullptr;
    }

    float x, y;
    TH1F* h = new TH1F(name.c_str(), name.c_str(), (xMax - xMin)/binWidth, xMin, xMax);

    while (parInput >> x >> y) {
        int bin = h->FindBin(x);
        h->SetBinContent(bin, y);
        h->SetBinError(bin, sqrt(y));
    }

    parInput.close();
    //h->Rebin();
    return h;
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
    float binWidth[] = {6.55,7.85,6.2,7.05,6.15,5.2,4.65,4.65};

    // Colori dei grafici
    int colors[] = {kCyan,kGreen,kBlue,kRed,kYellow,kMagenta,kOrange,kBlack};

    // Creazione degli istogrammi e dei fit
    TH1F* histos[8];
    TF1* fits[8];

    for (int i=0; i<8; i++) {
        histos[i] = ReadHisto(xMin[i], xMax[i], binWidth[i], files[i], Form("h%d",i));
        fits[i] = FitGauss(histos[i], xMin[i], xMax[i]);
    }
    double media[8]= {};
    double erroremedia[8]= {};
    // Creazione canvas separati e disegno
    for (int i=0; i<8; i++) {
        TCanvas* c = new TCanvas(Form("c%d",i), Form("numero fogli c%d",i), 800, 600);
        histos[i]->Draw("e1");
        if (fits[i]) {
            fits[i]->SetLineColor(colors[i]);
            fits[i]->Draw("same");
            //recupera parametri per grafico
            media[i]=fits[i]->GetParameter(1);
            erroremedia[i]= fits[i]->GetParError(1);
            
        }
    }

    /*il X2 è un po' basso, controllalo, dovrei accorpare i bin con meno di 5 elementi?
    come faccio con gli errori?
    */

    //spessore carta 2pm0.05 mm per 20 fogli 
    float Sp1= 2.0/20.0 ;//mm
    float Sper1= 0.05/20.0 ;//mm
    float Sp[8] = {0, Sp1*3,Sp1*6,Sp1*9,Sp1*12,Sp1*15,Sp1*18,Sp1*20};
    float Sper[8]= {Sper1,Sper1*3,Sper1*6,Sper1*9,Sper1*12,Sper1*15,Sper1*18,Sper1*20};
    //l'errore sullo spessore è variabile

    //rates

    float rates[8] = {};
    float errrates[8] = {};

    for (int i = 0; i < 8; i++)
    {
        rates[i]= media[i]*2;
        errrates[i]= erroremedia[i]*2;
    }
    
   auto c8 = new TCanvas("c8","Grafico assorbimento carta",200,10,700,500);
   c8->SetFillColor(42);
   c8->SetGrid();
   c8->GetFrame()->SetFillColor(21);
   c8->GetFrame()->SetBorderSize(12);
   const Int_t n = 8;
   
   auto gr = new TGraphErrors(n,Sp,rates,Sper,errrates);
   gr->SetTitle("Assorbimento carta");
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("AP");

   //fit esponenziale

   TF1 *f = new TF1("f","[0] +[1]* exp(-[2] * x)", 0, 2);

   f->SetParameters(1000,400,1); 
   gr->Fit(f,"R");
   cout<< endl<< "p-value"<< f->GetProb() << endl;
    //il coefficente di assorbimento viene super basso, circa 0,8
    //da test Z è però compatibile con quello teorico, il X2 va bene, per me si può fare
    
    

    
}