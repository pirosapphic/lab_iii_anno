
#include <fstream>
#include <math.h>   
#include <iostream>

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"

using namespace std;


TH1F* istogramma(string input,string name,const int nBins) {
  
 
  //static const int nBins = 128;
  double x[nBins];  int y[nBins]; 
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
  // imposta titoli degli assi
    myHisto->GetXaxis()->SetTitle("ADC channels");
    myHisto->GetYaxis()->SetTitle("Conteggi");
    myHisto->SetTitle("Spettro raggi cosmici");
    return myHisto;

}

void raggi_cosmici(){

    TH1F*h= istogramma("../data_sorg/A8_cosmici.txt","isto",128);
    TCanvas* c1 = new TCanvas("c1","h",800,600);
    h->Draw("hist");

    /*TH1F*hz= istogramma("../data_sorg/A8_cosmici_zoom.txt","isto",124);
    hz->Rebin(4);
    TCanvas* c2 = new TCanvas("c2","hz",800,600);
    hz->Draw("hist");
*/
    TH1F*hs= istogramma("../data_sorg/A8_cosmici_zoom.txt","isto",124);
    hs->Rebin(4);
    TCanvas* c3 = new TCanvas("c3","hs",800,600);
    hs->Draw("HIST E1");
    TF1*f1= new TF1("f1", "gaus", 15250., 47250.);
    f1->SetParameters(9,30750,20000);
    hs->Fit(f1,"R+","e1");
    f1->Draw("same");
    double media= f1->GetParameter(1);
    double emedia= f1->GetParError(1);
    cout << "p-value del fit: " << f1->GetProb() << endl;

    double energia= media*4.96766e-05;
    double eenergia = energia*sqrt(pow(emedia/media,2)+pow(1./5.,2));
    cout<<energia<<" +/- "<<eenergia<<endl;

    double penergia= 2.06*1.; //spessore dello scintillatore piccolo
    // 2.06 MeV/cm e 1 sono dati senza errore, allora prendo come costante senza errore penergia

    //test Z

    double Z= (penergia-energia)/eenergia;
    cout<< Z<< endl;
  /* uso keq invece che k perchè i raggi cosmici non cadono sempre perpenicolari alla 
  superficie e nel suo centro, ma possono passare per vari punti dello scintillatore,
  siccome k aumenta linearmente si usa keq che rappresenta il valore medio di k  per 
  segnali che colpiscono con eguale probabilità qualsiasi punto della lastra ovvero 
  quello misurato per un segnale concentrato a req = 0.3826 l dal centro.
  */

}