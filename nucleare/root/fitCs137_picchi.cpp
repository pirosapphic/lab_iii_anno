#include <fstream>
#include <math.h>   

#include "TH1F.h"
#include "TF1.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"

void fitCs137_picchi(string input = "../data_sorg/A8_Cs137_picchi.txt") 
{
 
  ifstream parInput(input.c_str());

  // in ADC counts
  float xMin = 505.; float xMax = 5615.;  float binWidth = 10.;

  TCanvas* c1 = new TCanvas("c1","c1",20,20,800,600);
  float x;  float y; string parName; 
  TH1F* theHisto = new TH1F("theHisto", "Electron energy", (int)((xMax-xMin)/binWidth), xMin, xMax);

  
  if (parInput.is_open()) {
    while ( parInput.good() ) {
      parInput >> x >> y;   
      //cout << x << " " << y << endl;
      int theBin = theHisto->FindBin(x);
      theHisto->SetBinContent(theBin,float(y));
      theHisto->SetBinError(theBin,sqrt(float(y)));
      //std::cout << theHisto->GetBinError(theBin) << std::endl;
    }
    parInput.close();
  }

  double A[10] = {};

  //prima gaussiana
  TF1* f1 = new TF1("f1","gaus",1370,1470);
  f1->SetParameters(1500, 1425, 90);
  f1->SetLineColor(kRed);

  theHisto->Draw("e1");
  theHisto->Fit(f1,"R");

 A[0]   = f1->GetParameter(0);
  double mean1  = f1->GetParameter(1);
  double sigma1 = f1->GetParameter(2);
  double errM1     = f1->GetParError(1);

 //seconda gaussiana
  TF1* f2 = new TF1("f2","gaus",1570,1640);
  f2->SetParameters(1600, 1600, 70);
  f2->SetLineColor(kBlue);

  theHisto->Fit(f2,"R+");
  double prob = f2->GetProb();
  cout << "p-value = " << prob << endl;

  A[1]    = f2->GetParameter(0);
  double mean2  = f2->GetParameter(1);
  double sigma2 = f2->GetParameter(2);
  double errM2     = f2->GetParError(1);
   //terza gaussiana
  TF1* f3 = new TF1("f3","gaus",1740,1820);
  f2->SetParameters(1700, 1770, 60);
  f3->SetLineColor(kGreen);

  theHisto->Fit(f3,"R+");

  A[2]    = f3->GetParameter(0);
  double mean3  = f3->GetParameter(1);
  double sigma3 = f3->GetParameter(2);
  double errM3     = f3->GetParError(1);
   
   //quarto gaussiana
  TF1* f4 = new TF1("f4","gaus",1900,2000);
  f4->SetParameters(1600, 1950, 80);
  f4->SetLineColor(kBlack);
  theHisto->Fit(f4,"R+");

   A[3]    = f4->GetParameter(0);
  double mean4  = f4->GetParameter(1);
  double sigma4 = f4->GetParameter(2);
  double errM4     = f4->GetParError(1);
   //quinta gaussiana
  TF1* f5 = new TF1("f5","gaus", 2070,2180);
  f5->SetParameters(1400, 2140, 75);
  f5->SetLineColor(kMagenta);
  theHisto->Fit(f5,"R+");

  A[4]    = f5->GetParameter(0);
  double mean5  = f5->GetParameter(1);
  double sigma5 = f5->GetParameter(2);
  double errM5     = f5->GetParError(1);

   //sesta gaussiana
  TF1* f6 = new TF1("f6","gaus", 2240,2350);
  f6->SetParameters(1300, 2310, 90);
  f6->SetLineColor(kCyan);
  theHisto->Fit(f6,"R+");

   A[5]    = f6->GetParameter(0);
  double mean6  = f6->GetParameter(1);
  double sigma6 = f6->GetParameter(2);
  double errM6     = f6->GetParError(1);

  //settima gaussiana
  TF1* f7 = new TF1("f7","gaus", 2425,2525);
  f7->SetParameters(1300, 2480, 80);
  f7->SetLineColor(kYellow);
  theHisto->Fit(f7,"R+");

   A[6]    = f7->GetParameter(0);
  double mean7  = f7->GetParameter(1);
  double sigma7 = f7->GetParameter(2);
  double errM7     = f7->GetParError(1);

  //ottavo gaussiana
  TF1* f8 = new TF1("f8","gaus", 2610,2700);
  f8->SetParameters(1200, 2650, 80);
  f8->SetLineColor(kOrange);

  theHisto->Fit(f8,"R+");

   A[7]    = f8->GetParameter(0);
  double mean8  = f8->GetParameter(1);
  double sigma8 = f8->GetParameter(2);
  double errM8     = f8->GetParError(1);

   //nono gaussiana
  TF1* f9 = new TF1("f9","gaus", 2780,2880);
  f9->SetParameters(1200, 2830, 80);
  f9->SetLineColor(kViolet);
  theHisto->Fit(f9,"R+");

  A[8]   = f9->GetParameter(0);
  double mean9  = f9->GetParameter(1);
  double sigma9 = f9->GetParameter(2);
  double errM9     = f9->GetParError(1);

   //dieci gaussiana
  TF1* f10 = new TF1("f10","gaus", 2960,3050);
  f10->SetParameters(1300, 3000, 70);
  f10->SetLineColor(kPink);
  theHisto->Fit(f10,"R+");

  A[9]   = f10->GetParameter(0);
  double mean10  = f10->GetParameter(1);
  double sigma10 = f10->GetParameter(2);
  double errM10     = f10->GetParError(1);

  //array 

  double media[10]= {mean1, mean2, mean3, mean4, mean5, mean6, mean7, mean8, mean9, mean10};
  double errM[10] ={errM1,errM2,errM3,errM4,errM5,errM6,errM7,errM8,errM9,errM10};
  double sigma[10] = {sigma1,sigma2,sigma3,sigma4,sigma5,sigma6,sigma7,sigma8,sigma9,sigma10};
  //delta p-p
  double deltapp[9] = {};
  double errpp[9] = {};
  for(int i=0;i < 9; i++){
    deltapp[i]= media[i+1]-media[i];
    errpp[i] = sqrt(pow(errM[i],2)+ pow(errM[i+1],2)); //misure indipendenti
    std::cout<<"Deltapp "<<i+1<<"("<<deltapp[i]<<" +/- "<<errpp[i]<<")CHN\n";
  }
  //media pesata
  double a[9] = {};
  double n= 0;
  double d=0;
  double meanpp=0;
  double error=0;
  for(int i=0; i<9; i++){
    a[i] = pow(1/errpp[i],2);
    n += a[i]*deltapp[i];
    d += a[i];

  }

  meanpp= n/d;
  error= pow(1/d,1/2);
  cout << "<deltapp> vale: " << meanpp << " +/- " << error << endl;

  // Nci = CHNci/meanpp

  double Nci[10] = {};
  for (int i = 0; i < 10; i++)
  {
    Nci[i]= media[i]/meanpp; 
  }
  
  //risoluzione sperimentale Rexp
   double Rexp[10] = {};
  for ( int i = 0; i < 10; i++)
  {
    Rexp[i]= sigma[i]/media[i]; 
  }
  
  //risoluzione attesa Rat

  double Rat[10] = {};
  for (int i = 0; i < 10; i++)
  {
    Rat[i]= 1/sqrt(Nci[i])*sqrt(1+pow(sigma[i]/meanpp,2));
  }
  
  for (int i = 0; i < 10; i++) {
        cout << Rexp[i]<< "  " << Rat[i]<< endl;
  }
 
}

 /*cose da fare:chiedere cos'è sigmap, media sigma o diversa per ogni fit?
 cercare errore R e verificare formula Rat
  */
 /*cose da fare:chiedere cos'è sigmap, media sigma o diversa per ogni fit?
 cercare errore R e verificare formula Rat
  */