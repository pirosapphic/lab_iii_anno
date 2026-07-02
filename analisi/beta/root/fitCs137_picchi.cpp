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

  TCanvas* c1 = new TCanvas("c1","c1",20,20,1098,732);
  c1->SetGrid();
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

    // imposta titoli degli assi
   theHisto->GetXaxis()->SetTitle("ADC channels [CHN]");
   theHisto->GetYaxis()->SetTitle("Conteggi [#]");
   theHisto->SetTitle("Spettro a bassi canali del ^{137}Cs");
  double A[10] = {};

  //prima gaussiana
  TF1* f1 = new TF1("f1","gaus",1370,1490);
  f1->SetParameters(1500, 1425, 90);
  f1->SetLineColor(kRed);

  gStyle->SetOptStat(0);
  theHisto->Draw("e1");
  theHisto->Fit(f1,"R");
  double prob = f1->GetProb();
  cout << "p-value = " << prob << endl;

 A[0]   = f1->GetParameter(0);
  double mean1  = f1->GetParameter(1);
  double sigma1 = f1->GetParameter(2);
  double errM1     = f1->GetParError(1);
  double ersig1 = f1->GetParError(2);

 //seconda gaussiana
  TF1* f2 = new TF1("f2","gaus",1570,1640);
  f2->SetParameters(1600, 1600, 70);
  f2->SetLineColor(kBlue);

  theHisto->Fit(f2,"R+");
  prob = f2->GetProb();
  cout << "p-value = " << prob << endl;

  A[1]    = f2->GetParameter(0);
  double mean2  = f2->GetParameter(1);
  double sigma2 = f2->GetParameter(2);
  double errM2     = f2->GetParError(1);
  double ersig2 = f2->GetParError(2);

   //terza gaussiana
  TF1* f3 = new TF1("f3","gaus",1740,1820);
  f2->SetParameters(1700, 1770, 60);
  f3->SetLineColor(kGreen);

  theHisto->Fit(f3,"R+");

  A[2]    = f3->GetParameter(0);
  double mean3  = f3->GetParameter(1);
  double sigma3 = f3->GetParameter(2);
  double errM3     = f3->GetParError(1);
  double ersig3 = f3->GetParError(2);  

   //quarto gaussiana
  TF1* f4 = new TF1("f4","gaus",1900,2000);
  f4->SetParameters(1600, 1950, 80);
  f4->SetLineColor(kBlack);
  theHisto->Fit(f4,"R+");

   A[3]    = f4->GetParameter(0);
  double mean4  = f4->GetParameter(1);
  double sigma4 = f4->GetParameter(2);
  double errM4     = f4->GetParError(1);
  double ersig4 = f4->GetParError(2);

   //quinta gaussiana
  TF1* f5 = new TF1("f5","gaus", 2070,2180);
  f5->SetParameters(1400, 2140, 75);
  f5->SetLineColor(kMagenta);
  theHisto->Fit(f5,"R+");

  A[4]    = f5->GetParameter(0);
  double mean5  = f5->GetParameter(1);
  double sigma5 = f5->GetParameter(2);
  double errM5     = f5->GetParError(1);
  double ersig5 = f5->GetParError(2);

   //sesta gaussiana
  TF1* f6 = new TF1("f6","gaus", 2240,2350);
  f6->SetParameters(1300, 2310, 90);
  f6->SetLineColor(kCyan);
  theHisto->Fit(f6,"R+");

   A[5]    = f6->GetParameter(0);
  double mean6  = f6->GetParameter(1);
  double sigma6 = f6->GetParameter(2);
  double errM6     = f6->GetParError(1);
  double ersig6 = f6->GetParError(2);

  //settima gaussiana
  TF1* f7 = new TF1("f7","gaus", 2425,2525);
  f7->SetParameters(1300, 2480, 80);
  f7->SetLineColor(kYellow);
  theHisto->Fit(f7,"R+");

   A[6]    = f7->GetParameter(0);
  double mean7  = f7->GetParameter(1);
  double sigma7 = f7->GetParameter(2);
  double errM7     = f7->GetParError(1);
  double ersig7 = f7->GetParError(2);
  //ottavo gaussiana
  TF1* f8 = new TF1("f8","gaus", 2610,2700);
  f8->SetParameters(1200, 2650, 80);
  f8->SetLineColor(kOrange);

  theHisto->Fit(f8,"R+");

   A[7]    = f8->GetParameter(0);
  double mean8  = f8->GetParameter(1);
  double sigma8 = f8->GetParameter(2);
  double errM8     = f8->GetParError(1);
  double ersig8 = f8->GetParError(2);

   //nono gaussiana
  TF1* f9 = new TF1("f9","gaus", 2780,2880);
  f9->SetParameters(1200, 2830, 80);
  f9->SetLineColor(kViolet);
  theHisto->Fit(f9,"R+");

  A[8]   = f9->GetParameter(0);
  double mean9  = f9->GetParameter(1);
  double sigma9 = f9->GetParameter(2);
  double errM9     = f9->GetParError(1); 
  double ersig9 = f9->GetParError(2);
   //dieci gaussiana
  TF1* f10 = new TF1("f10","gaus", 2960,3050);
  f10->SetParameters(1300, 3000, 70);
  f10->SetLineColor(kPink);
  theHisto->Fit(f10,"R+");

  A[9]   = f10->GetParameter(0);
  double mean10  = f10->GetParameter(1);
  double sigma10 = f10->GetParameter(2);
  double errM10     = f10->GetParError(1);
  double ersig10 = f10->GetParError(2);
  //array 

  double media[10]= {mean1, mean2, mean3, mean4, mean5, mean6, mean7, mean8, mean9, mean10};
  //uso sigma come errore
  double s_media[10] ={errM1,errM2,errM3,errM4,errM5,errM6,errM7,errM8,errM9,errM10};
  double sigma[10] = {sigma1,sigma2,sigma3,sigma4,sigma5,sigma6,sigma7,sigma8,sigma9,sigma10};
  double ersig[10] = {ersig1,ersig2,ersig3,ersig4,ersig5,ersig6,ersig7,ersig8,ersig9,ersig10};


  //delta p-p
  double deltapp[9] = {};
  double errpp[9] = {};
  for(int i=0;i < 9; i++){
    deltapp[i]= media[i+1]-media[i];
    errpp[i] = sqrt(pow(sigma[i],2)+ pow(sigma[i+1],2)); //misure indipendenti
    cout<<"deltapp "<< i<<" "<<deltapp[i]<<"+/-"<<errpp[i]<<endl;
  }
  //media pesata
  double a[9] = {};
  double n= 0.0;
  double d=0.0;
  double meanpp=0.0;
  double s_meanpp=0.0;
  for(int i=0; i<9; i++){
    a[i] = pow(1./errpp[i],2);
    n += a[i]*deltapp[i];
    d += a[i];

  }

  meanpp= n/d;
  s_meanpp= pow(1./d,0.5); //errore su deltapp
  cout << "<deltapp> vale: " << meanpp << " +/- " << s_meanpp << endl;
/*
//sigma pesata
   n= 0;
   d=0;
  double sigmap=0;
  double ersigmap=0;
  for(int i=0; i<9; i++){
    a[i] = pow(1/ersig[i],2);
    n += a[i]*sigma[i];
    d += a[i];

  }
 sigmap = n/d;
  ersigmap = pow(1/d,0.5);
  cout << "<sigmap vale: " <<sigmap  << " +/- " << ersigmap  << endl;
  // Nci = CHNci/meanpp
  //sti valori li ho presi dai fit di paolo spero siano giusti

//CHNci = (19655.1 +/- 2580.22) //CHN
//k = (2.4287e-05 +/- 5.44346e-06) //MeV/CHN
  float CHNci= 19655.1;
  float eCHNci= 2580.22;
  float Nci= CHNci/meanpp;

  float eNci= pow(pow(eCHNci/meanpp,2)+pow(CHNci*error/pow(meanpp,2),2),0.5);
  cout << "Nci: " <<Nci << " +/- " << eNci << endl;
  //risoluzione attesa Rat

  double Rat;
  double erat;//devo calcolare l'errore

  Rat= 1/sqrt(Nci)*sqrt(1+pow(sigmap/meanpp,2));
  float parteB= pow(Nci,-0.5);
  float parteA = 1 + pow(sigmap/meanpp,2);


  double dR_dN = -0.5 * pow(Nci,-1.5) * sqrt(parteA);
  double dR_dsigma = parteB * (sigmap / (pow(meanpp,2)*sqrt(parteA)));
  double dR_dmean = -parteB * (pow(sigmap,2) / (pow(meanpp,3)*sqrt(parteA)));

  erat = sqrt(
    pow(dR_dN * eNci, 2) +
    pow(dR_dsigma * ersigmap, 2) +
    pow(dR_dmean * error, 2)
  );

  cout <<  " La risoluzione attesa è: " << Rat<< " +/- " << erat<< endl;

  //test Z picchi
 double Zpp;
  for (int i = 0; i < 8; i++)
  {
    for (int j= i+1; j < 8; j++)
    {
      Zpp= (deltapp[i]-deltapp[j])/sqrt(pow(sigma[i],2)+pow(sigma[j],2));
      cout<< i << " " <<j <<"Z "<< Zpp <<endl;
      if (1.96 >Zpp && -1.96 < Zpp)
      {
       cout<< "buono"<<endl;
      }
      
    }
    
  }
  */
 // //test Z risoluzione
 // double Zr= (Rat-0.128)/sqrt(pow(erat,2)+pow(0.016,2));
  //cout<< "Zr: "<< Zr<< endl;
//ora facciamo il fit delle sigma^2:
//s^2_n = s_0^2 + ns_p^2
std::vector<double>var(10);
std::vector<double>s_var(10);
std::vector<double>nci(10);
std::vector<double>s_nci(10);
std::vector<double>sqnci(10);
std::vector<double>s_sqnci(10);
for(int i = 0;i<10;i++){
    var[i] = pow(sigma[i],2);
    s_var[i] = var[i]*2.*ersig[i]/sigma[i];
    nci[i] = media[i]/meanpp; //calcolo ogni volta quanto vale il n. di celle accese
    s_nci[i] = nci[i]*sqrt(sigma[i]/media[i]+s_meanpp/meanpp);
    nci[i] = 8.+i;
    s_nci[i] = 0.5;
    sqnci[i] = sqrt(nci[i]);
    s_sqnci[i] = s_nci[i]/2./sqnci[i];
}
TGraphErrors* g1 = new TGraphErrors(var.size(),sqnci.data(),sigma,s_sqnci.data(),ersig); //sigma vs sqrt(N)
//TGraphErrors* g2 = new TGraphErrors(var.size(),nci.data(),var.data(),s_nci.data(),s_var.data()); //var vs N
g1->SetMarkerStyle(7);
g1->SetTitle("#sigma_{p,i} vs #sqrt{N};#sqrt{N}[#];#sigma_{p} [CHN]");
//g2->SetMarkerStyle(7);
//g2->SetTitle("#sigma^{2}_{p,N} vs N;N[#];#sigma^{2}_{p} [CHN^2]");
TCanvas* c2 = new TCanvas("c2","c2",20,20,1098,732);
c2->cd();
c2->SetGrid();
g1->Draw("ALP");
TF1* lin1 = new TF1("lin1","pol1");
g1->Fit(lin1,"R+","",2.9,4.1);
std::cout<<"p="<<lin1->GetProb()<<std::endl;
double sigmap = lin1->GetParameter(1);
double s_sigmap = lin1->GetParError(1);
for(int i=0;i<var.size();i++){
    std::cout<<"N celle "<<nci[i]<<"pm"<<s_nci[i]<<" s^2 "<<var[i]<<"pm"<<s_var[i]<<std::endl;
}
//ORA calcolo la risoluzione attesa del picco di conversione interna


    double CHNci= 19704.8;
    double eCHNci= 2512.8;
    double Nci= CHNci/meanpp;

  double s_Nci= pow(pow(eCHNci/meanpp,2)+pow(CHNci*s_meanpp/pow(meanpp,2),2),0.5);
  cout << "Nci: " <<Nci << " +/- " << s_Nci << endl;



    double R_att = pow(Nci,-0.5)*sqrt(1+sigmap*sigmap/meanpp/meanpp);
    double s_R_att = pow(sqrt(1+sigmap*sigmap/meanpp/meanpp)*1.5*pow(Nci,-1.5)*s_Nci,2);
    s_R_att+= pow(sigmap/Nci/meanpp/meanpp/sqrt(1./Nci+sigmap*sigmap/meanpp/meanpp/Nci)*s_sigmap,2);
    s_R_att+=pow(0.5*sigmap*sigmap*3*pow(meanpp,-3)/Nci/sqrt(1./Nci+sigmap*sigmap/meanpp/meanpp/Nci)*s_meanpp,2);
    s_R_att = sqrt(s_R_att);


//    double s_R_att=pow(3./2.*R_att*s_Nci/Nci,2);
 //   s_R_att+=pow(pow(Nci,-0.5)*sigmap*pow(1+sigmap*sigmap/meanpp/meanpp,-0.5)*s_sigmap,2);
 //   s_R_att+=pow(pow(Nci,-0.5)*pow(1+sigmap*sigmap/meanpp/meanpp,-0.5)*pow(meanpp,-3)*s_meanpp,2);
  //  s_R_att = sqrt(s_R_att); //i REALLY hope this is correct.
std::cout<<"R_att = "<<R_att<<"pm"<<s_R_att<<"\n";

}
