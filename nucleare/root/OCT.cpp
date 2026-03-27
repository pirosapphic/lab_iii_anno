
#include <fstream>
#include <math.h>   
#include <iostream>

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"

using namespace std;



TH1F*low(string input = "../data_SiPM/DCR/A8_DCR_17mV.txt") {
  
  
  static const int nBins = 30;
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
  TH1F *h = new TH1F("h","h",nBins,x[0],x[nBins-1]);  
  for (int j=0;j<nBins;j++) {
    h->SetBinContent(j+1,y[j]);
    h->SetBinError(j+1, sqrt(y[j]));
  }
  // imposta titoli degli assi
    h->GetXaxis()->SetTitle("Frequenza [Hz]");
    h->GetYaxis()->SetTitle("Conteggi");
    h->SetTitle("Frequenza alla soglia 17mV");

  return h;
  
}

TH1F*high(string input = "../data_SiPM/DCR/A8_DCR_44mV.txt") {
  
  
  static const int nBins = 30;
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
  TH1F *h = new TH1F("h","h",nBins,x[0],x[nBins-1]);  
  for (int j=0;j<nBins;j++) {
    h->SetBinContent(j+1,y[j]);
    h->SetBinError(j+1, sqrt(y[j]));
  }
  // imposta titoli degli assi
    h->GetXaxis()->SetTitle("Frequenza [Hz]");
    h->GetYaxis()->SetTitle("Conteggi");
    h->SetTitle("Frequenza alla soglia 44mV");

  return h;
  
}

void OCT(){
    
    TH1F*hr1 = low();
    TCanvas* c1 = new TCanvas("c1","R1",800,600);
    hr1->Draw();
    TF1*f1= new TF1("f1", "gaus", 25111.5, 26009.5);
    hr1->Fit(f1,"R");
    float R1= f1->GetParameter(1);
    float er1= f1->GetParError(1);
    cout << "p-value del fit: " << f1->GetProb() << endl;

    TH1F*hr2 = high();
    TCanvas* c2 = new TCanvas("c2","R2",800,600);
    hr2->Draw();
    TF1*f2= new TF1("f1", "gaus",  549.1, 670.9);
    hr2->Fit(f2,"R");
    float R2= f2->GetParameter(1);
    float er2= f2->GetParError(1);
    cout << "p-value del fit: " << f2->GetProb() << endl;

    //OCT

    float oct= R1/R2;
    double eoct= pow(pow(er1/R2,2)+pow(R1*(float)er2/pow(R2,2),2),1/2);
    cout << "OCT: " <<oct << " +/- " << eoct << endl;
    cout<< R1 <<R2 << er1<< er2 <<endl;
}
