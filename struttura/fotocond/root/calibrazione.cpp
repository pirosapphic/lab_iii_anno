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

tuple<vector<double>,vector<double>,vector<double>,vector<double>>lettura(string path);

void calibrazione(){ 
    string path= "../data_fc/data_fotocorr.txt";
    auto [ua,s_ua,V_0,s_V_0]= lettura(path);
    vector<double> pixel,s_pixel;
    double x,s_x, cov;
    double p0=-417;
    double s_p0=400;
    double p1=2.2;
    double s_p1=1.2;
    double p2= 0.0006;
    double s_p2= 0.0009;
    double s_p01=-418.8; //bisogna tenere conto del segno nella covarianza
    double s_p02=-0.3287;
    double s_p21=0.001044;
    for (int i = 0; i < ua.size(); i++)
    {
        x= p0 + ua[i]*p1+pow(ua[i],2)*p2; //passo da ua a pixel relazione polinomiale calib_monocromatore
        cov= 2*ua[i]*s_p01 +2*pow(ua[i],2)*s_p02+ 2*pow(ua[i],3)*s_p21;
        s_x= sqrt(pow(s_p0,2)+pow(ua[i]*s_p1,2)+pow(pow(ua[i],2)*s_p2,2)+pow(p1+2.*ua[i]*p2,2)*pow(s_ua[i],2)+cov);
        pixel.push_back(x);
        s_pixel.push_back(s_x);
       // cout<< pixel[i] <<" pm "<< s_x <<endl;
    }



    vector<double> lambda,s_lambda;
    double f0=178.;
    double s_f0=3.;
    double f1=0.388;
    double s_f1=0.006;
    double f2= -2.17142e-05;
    double s_f2= 3.e-06;
    double s_f01=-0.01515 ; //bisogna tenere conto del segno nella covarianza
    double s_f02=7.e-06;
    double s_f21=-1.8e-08;
    for (int i = 0; i < ua.size(); i++)
    {
        x= f0 + pixel[i]*f1+pow(pixel[i],2)*f2; //passo da pixel a nm relazione polinomiale calib_monocromatore
        cov= 2*pixel[i]*s_f01 +2*pow(pixel[i],2)*s_f02+ 2*pow(pixel[i],3)*s_f21;
        s_x= sqrt(pow(s_f0,2)+pow(pixel[i]*s_f1,2)+pow(pow(pixel[i],2)*s_f2,2)+pow(f1+2.*pixel[i]*f2,2)*pow(s_pixel[i],2)+cov);
        lambda.push_back(x);
        s_lambda.push_back(s_x);
        cout<< lambda[i] <<" "<< s_lambda[i] <<" "<< V_0[i] <<" "<<s_V_0[i]<<endl;
    }
    


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
