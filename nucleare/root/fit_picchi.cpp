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

void fit_picchi(string input = "../data_sorg/A8_Cs137_picchi.txt") 
{
 
  ifstream parInput(input.c_str());

  // in ADC counts
  float xMin = 505.; float xMax = 5615.;  float binWidth = 10.;

  TCanvas* c1 = new TCanvas("c1","c1",20,20,800,600);
  float x;  float y; string parName; 
  TH1F* theHisto = new TH1F("theHisto", "Electron energy", (int)((xMax-xMin)/binWidth), xMin, xMax);
  //theHisto->Sumw2();
  
  if (parInput.is_open()) {
    while ( parInput.good() ) {
      parInput >> x >> y;   
      cout << x << " " << y << endl;
      int theBin = theHisto->FindBin(x);
      theHisto->SetBinContent(theBin,float(y));
      // theHisto->SetBinError(theBin,sqrt(float(y)));
    }
    parInput.close();
  }

  double A[10];

  //prima gaussiana
  TF1* f1 = new TF1("f1","gaus",1370,1470);
  f1->SetParameters(1500, 1425, 90);
  f1->SetLineColor(kRed);

  theHisto->Draw("e1");
  theHisto->Fit(f1,"R");

  double A[0]   = f1->GetParameter(0);
  double mean1  = f1->GetParameter(1);
  double sigma1 = f1->GetParameter(2);
  double errM1     = f1->GetParError(1);

 //seconda gaussiana
  TF1* f2 = new TF1("f2","gaus",1550,1650);
  f2->SetParameters(1600, 1600, 70);
  f2->SetLineColor(kBlue);

  theHisto->Fit(f2,"R+");

  double [1]    = f2->GetParameter(0);
  double mean2  = f2->GetParameter(1);
  double sigma2 = f2->GetParameter(2);
  double errM2     = f2->GetParError(1);
   //terza gaussiana
  TF1* f3 = new TF1("f3","gaus",1720,1820);
  f2->SetParameters(1600, 1775, 90);
  f3->SetLineColor(kGreen);

  theHisto->Fit(f3,"R+");

  double A[2]    = f3->GetParameter(0);
  double mean3  = f3->GetParameter(1);
  double sigma3 = f3->GetParameter(2);
  double errM3     = f3->GetParError(1);
   //quarto gaussiana
  TF1* f4 = new TF1("f4","gaus",1900,2000);
  f4->SetParameters(1600, 1960, 80);
  f4->SetLineColor(kBlack);
  theHisto->Fit(f4,"R+");

  double A[3]    = f4->GetParameter(0);
  double mean4  = f4->GetParameter(1);
  double sigma4 = f4->GetParameter(2);
  double errM4     = f4->GetParError(1);
   //quinta gaussiana
  TF1* f5 = new TF1("f5","gaus", 2070,2190);
  f5->SetParameters(1400, 2140, 75);
  f5->SetLineColor(kMagenta);
  theHisto->Fit(f5,"R+");

  double A[4]    = f5->GetParameter(0);
  double mean5  = f5->GetParameter(1);
  double sigma5 = f5->GetParameter(2);
  double errM5     = f5->GetParError(1);

     //sesta gaussiana
  TF1* f6 = new TF1("f6","gaus", 2240,2350);
  f6->SetParameters(1300, 2310, 90);
  f6->SetLineColor(kCyan);
  theHisto->Fit(f6,"R+");

  double A[5]    = f6->GetParameter(0);
  double mean6  = f6->GetParameter(1);
  double sigma6 = f6->GetParameter(2);
  double errM6     = f6->GetParError(1);

  //settima gaussiana
  TF1* f7 = new TF1("f7","gaus", 2425,2525);
  f7->SetParameters(1300, 2480, 80);
  f7->SetLineColor(kYellow);
  theHisto->Fit(f7,"R+");

  double A[6]    = f7->GetParameter(0);
  double mean7  = f7->GetParameter(1);
  double sigma7 = f7->GetParameter(2);
  double errM7     = f7->GetParError(1);

  //ottavo gaussiana
  TF1* f8 = new TF1("f8","gaus", 2610,2700);
  f8->SetParameters(1200, 2650, 80);
  f8->SetLineColor(kOrange);

  theHisto->Fit(f8,"R+");

  double A[7]    = f8->GetParameter(0);
  double mean8  = f8->GetParameter(1);
  double sigma8 = f8->GetParameter(2);
  double errM8     = f8->GetParError(1);

   //nono gaussiana
  TF1* f9 = new TF1("f9","gaus", 2780,2880);
  f9->SetParameters(1200, 2830, 80);
  f9->SetLineColor(kViolet);
  theHisto->Fit(f9,"R+");

  double A[8]   = f9->GetParameter(0);
  double mean9  = f9->GetParameter(1);
  double sigma9 = f9->GetParameter(2);
  double errM9     = f9->GetParError(1);

   //dieci gaussiana
  TF1* f10 = new TF1("f10","gaus", 2940,3050);
  f10->SetParameters(1000, 3010, 60);
  f10->SetLineColor(kPink);
  theHisto->Fit(f10,"R+");

  double A[9]   = f10->GetParameter(0);
  double mean10  = f10->GetParameter(1);
  double sigma10 = f10->GetParameter(2);
  double errM10     = f10->GetParError(1);
 

}


//void fitCs137(string input = "cs137-gruppob11.txt") 