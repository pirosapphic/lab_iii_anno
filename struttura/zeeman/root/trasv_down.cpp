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



std::vector<std::vector<double>> tsvreader(string path, const int ncols);
std::vector<double> w_mean(std::vector<double> val, std::vector<double> s_val);

void trasv_down(){
    std::vector<double> I = {0.,5.66,6.30,7.04,7.60,9.02,9.58};
    int n = I.size();
    std::vector<double> s_I(n); 
    for(int i = 0; i<n; i++){
	s_I[i] = 0.025*I[i]+0.10;
    }
    std::vector<double> B(n);
    std::vector<double> s_B(n);
//B CALIBRATION PARAMETERS B = mI+q
    double m = 63.87;
    double s_m = 1.1875;
    double q = 11.91;
    double s_q = 5.24303;
    for(int i = 0; i<n; i++){
	B[i] = (m*I[i]+q)/1000.; //convert to Tesla
	s_B[i] = sqrt(pow(I[i]*s_m,2)+pow(m*s_I[i],2)+pow(s_q,2))/1000.;
    }
    std::vector<double> delta;
    std::vector<double> s_delta;
    std::vector<double> Delta;
    std::vector<double> s_Delta;
    std::vector<double> d_D(n);
    std::vector<double> s_d_D(n);
//I=0----------------------------------------------------------------------
    string path = "../data/trasv_down/dt0.txt";
    std::vector<std::vector<double>>dati= tsvreader(path,4);
    int ncols = dati[0].size();
    int nrows = dati.size();
    for(int i=0; i<nrows; i++){} //if i remove this, dati ceases to exist!! WTF?!
    
    double delta1	= dati[0][2]-dati[0][0];//doppietto n=1
    double s_delta1	= dati[0][1]+dati[0][3];
    double delta2	= dati[1][2]-dati[1][0];//doppietto n=2
    double s_delta2	= dati[1][1]+dati[1][3];
    double delta3	= dati[2][2]-dati[2][0];//doppietto n=3
    double s_delta3	= dati[2][1]+dati[2][3];
    double Delta1	= dati[1][0]-dati[0][0];//interno-interno 1 e 2 
    double s_Delta1	= dati[1][1]+dati[0][1];
    double Delta2	= dati[2][0]-dati[1][0];//interno-interno 2 e 3
    double s_Delta2	= dati[2][1]+dati[1][1];
    double Delta3	= dati[1][2]-dati[0][2];//esterno-esterno 1 e 2 
    double s_Delta3	= dati[1][3]+dati[0][3];
    double Delta4	= dati[2][2]-dati[1][2];//esterno-esterno 2 e 3
    double s_Delta4	= dati[2][3]+dati[1][3];

    double mdelta = w_mean({delta1,delta2,delta3},{s_delta1,s_delta2,s_delta3})[0];
    double s_mdelta = w_mean({delta1,delta2,delta3},{s_delta1,s_delta2,s_delta3})[1];
    
    double mDelta = w_mean({Delta1,Delta2,Delta3},{s_Delta1,s_Delta2,s_Delta3})[0];
    double s_mDelta = w_mean({Delta1,Delta2,Delta3,Delta4},{s_Delta1,s_Delta2,s_Delta3,s_Delta4})[1];
    std::cout<<"File "+path+": delta = "<<mdelta<<" +/- "<<s_mdelta<<"\t\tDelta = "<<mDelta<<" +/- "<<s_mDelta<<std::endl;
    delta.push_back(mdelta);
    s_delta.push_back(s_mdelta);
    Delta.push_back(mDelta);
    s_Delta.push_back(s_mDelta);
   // dati.clear();

//I=2----------------------------------------------------------------------
    path = "../data/trasv_down/dt5.txt";
    dati = tsvreader(path,4);
    ncols = dati[0].size();
    nrows = dati.size();
    for(int i=0; i<nrows; i++){} //if i remove this, dati ceases to exist!! WTF?!
    delta1	= dati[0][2]-dati[0][0];//doppietto n=1
    s_delta1	= dati[0][1]+dati[0][3];
    delta2	= dati[1][2]-dati[1][0];//doppietto n=2
    s_delta2	= dati[1][1]+dati[1][3];
    delta3	= dati[2][2]-dati[2][0];//doppietto n=3
    s_delta3	= dati[2][1]+dati[2][3];
    Delta1	= dati[1][0]-dati[0][0];//interno-interno 1 e 2 
    s_Delta1	= dati[1][1]+dati[0][1];
    Delta2	= dati[2][0]-dati[1][0];//interno-interno 2 e 3
    s_Delta2	= dati[2][1]+dati[1][1];
    Delta3	= dati[1][2]-dati[0][2];//esterno-esterno 1 e 2 
    s_Delta3	= dati[1][3]+dati[0][3];
    Delta4	= dati[2][2]-dati[1][2];//esterno-esterno 2 e 3
    s_Delta4	= dati[2][3]+dati[1][3];
    mdelta = w_mean({delta1,delta2,delta3},{s_delta1,s_delta2,s_delta3})[0];
    s_mdelta = w_mean({delta1,delta2,delta3},{s_delta1,s_delta2,s_delta3})[1];
    
    mDelta = w_mean({Delta1,Delta2,Delta3},{s_Delta1,s_Delta2,s_Delta3})[0];
    s_mDelta = w_mean({Delta1,Delta2,Delta3,Delta4},{s_Delta1,s_Delta2,s_Delta3,s_Delta4})[1];
    std::cout<<"File "+path+": delta = "<<mdelta<<" +/- "<<s_mdelta<<"\tDelta = "<<mDelta<<" +/- "<<s_mDelta<<std::endl;
    delta.push_back(mdelta);
    s_delta.push_back(s_mdelta);
    Delta.push_back(mDelta);
    s_Delta.push_back(s_mDelta);
    
//I=3----------------------------------------------------------------------
    path = "../data/trasv_down/dt6.txt";
    dati = tsvreader(path,4);
    ncols = dati[0].size();
    nrows = dati.size();
    for(int i=0; i<nrows; i++){} //if i remove this, dati ceases to exist!! WTF?!
    delta1	= dati[0][2]-dati[0][0];//doppietto n=1
    s_delta1	= dati[0][1]+dati[0][3];
    delta2	= dati[1][2]-dati[1][0];//doppietto n=2
    s_delta2	= dati[1][1]+dati[1][3];
    delta3	= dati[2][2]-dati[2][0];//doppietto n=3
    s_delta3	= dati[2][1]+dati[2][3];
    Delta1	= dati[1][0]-dati[0][0];//interno-interno 1 e 2 
    s_Delta1	= dati[1][1]+dati[0][1];
    Delta2	= dati[2][0]-dati[1][0];//interno-interno 2 e 3
    s_Delta2	= dati[2][1]+dati[1][1];
    Delta3	= dati[1][2]-dati[0][2];//esterno-esterno 1 e 2 
    s_Delta3	= dati[1][3]+dati[0][3];
    Delta4	= dati[2][2]-dati[1][2];//esterno-esterno 2 e 3
    s_Delta4	= dati[2][3]+dati[1][3];
    mdelta = w_mean({delta1,delta2,delta3},{s_delta1,s_delta2,s_delta3})[0];
    s_mdelta = w_mean({delta1,delta2,delta3},{s_delta1,s_delta2,s_delta3})[1];
    
    mDelta = w_mean({Delta1,Delta2,Delta3},{s_Delta1,s_Delta2,s_Delta3})[0];
    s_mDelta = w_mean({Delta1,Delta2,Delta3,Delta4},{s_Delta1,s_Delta2,s_Delta3,s_Delta4})[1];
    std::cout<<"File "+path+": delta = "<<mdelta<<" +/- "<<s_mdelta<<"\tDelta = "<<mDelta<<" +/- "<<s_mDelta<<std::endl;
    delta.push_back(mdelta);
    s_delta.push_back(s_mdelta);
    Delta.push_back(mDelta);
    s_Delta.push_back(s_mDelta);
//I=4----------------------------------------------------------------------
    path = "../data/trasv_down/dt7.txt";
    dati = tsvreader(path,4);
    ncols = dati[0].size();
    nrows = dati.size();
    for(int i=0; i<nrows; i++){} //if i remove this, dati ceases to exist!! WTF?!
    delta1	= dati[0][2]-dati[0][0];//doppietto n=1
    s_delta1	= dati[0][1]+dati[0][3];
    delta2	= dati[1][2]-dati[1][0];//doppietto n=2
    s_delta2	= dati[1][1]+dati[1][3];
    delta3	= dati[2][2]-dati[2][0];//doppietto n=3
    s_delta3	= dati[2][1]+dati[2][3];
    Delta1	= dati[1][0]-dati[0][0];//interno-interno 1 e 2 
    s_Delta1	= dati[1][1]+dati[0][1];
    Delta2	= dati[2][0]-dati[1][0];//interno-interno 2 e 3
    s_Delta2	= dati[2][1]+dati[1][1];
    Delta3	= dati[1][2]-dati[0][2];//esterno-esterno 1 e 2 
    s_Delta3	= dati[1][3]+dati[0][3];
    Delta4	= dati[2][2]-dati[1][2];//esterno-esterno 2 e 3
    s_Delta4	= dati[2][3]+dati[1][3];
    mdelta = w_mean({delta1,delta2,delta3},{s_delta1,s_delta2,s_delta3})[0];
    s_mdelta = w_mean({delta1,delta2,delta3},{s_delta1,s_delta2,s_delta3})[1];
    
    mDelta = w_mean({Delta1,Delta2,Delta3},{s_Delta1,s_Delta2,s_Delta3})[0];
    s_mDelta = w_mean({Delta1,Delta2,Delta3,Delta4},{s_Delta1,s_Delta2,s_Delta3,s_Delta4})[1];
    std::cout<<"File "+path+": delta = "<<mdelta<<" +/- "<<s_mdelta<<"\tDelta = "<<mDelta<<" +/- "<<s_mDelta<<std::endl;
    delta.push_back(mdelta);
    s_delta.push_back(s_mdelta);
    Delta.push_back(mDelta);
    s_Delta.push_back(s_mDelta);
//I=5----------------------------------------------------------------------
    path = "../data/trasv_down/dt7_6.txt";
    dati = tsvreader(path,4);
    ncols = dati[0].size();
    nrows = dati.size();
    for(int i=0; i<nrows; i++){} //if i remove this, dati ceases to exist!! WTF?!
    delta1	= dati[0][2]-dati[0][0];//doppietto n=1
    s_delta1	= dati[0][1]+dati[0][3];
    delta2	= dati[1][2]-dati[1][0];//doppietto n=2
    s_delta2	= dati[1][1]+dati[1][3];
    delta3	= dati[2][2]-dati[2][0];//doppietto n=3
    s_delta3	= dati[2][1]+dati[2][3];
    Delta1	= dati[1][0]-dati[0][0];//interno-interno 1 e 2 
    s_Delta1	= dati[1][1]+dati[0][1];
    Delta2	= dati[2][0]-dati[1][0];//interno-interno 2 e 3
    s_Delta2	= dati[2][1]+dati[1][1];
    Delta3	= dati[1][2]-dati[0][2];//esterno-esterno 1 e 2 
    s_Delta3	= dati[1][3]+dati[0][3];
    Delta4	= dati[2][2]-dati[1][2];//esterno-esterno 2 e 3
    s_Delta4	= dati[2][3]+dati[1][3];
    mdelta = w_mean({delta1,delta2,delta3},{s_delta1,s_delta2,s_delta3})[0];
    s_mdelta = w_mean({delta1,delta2,delta3},{s_delta1,s_delta2,s_delta3})[1];
    
    mDelta = w_mean({Delta1,Delta2,Delta3},{s_Delta1,s_Delta2,s_Delta3})[0];
    s_mDelta = w_mean({Delta1,Delta2,Delta3,Delta4},{s_Delta1,s_Delta2,s_Delta3,s_Delta4})[1];
    std::cout<<"File "+path+": delta = "<<mdelta<<" +/- "<<s_mdelta<<"\tDelta = "<<mDelta<<" +/- "<<s_mDelta<<std::endl;
    delta.push_back(mdelta);
    s_delta.push_back(s_mdelta);
    Delta.push_back(mDelta);
    s_Delta.push_back(s_mDelta);
//I=6----------------------------------------------------------------------
    path = "../data/trasv_down/dt9.txt";
    dati = tsvreader(path,4);
    ncols = dati[0].size();
    nrows = dati.size();
    for(int i=0; i<nrows; i++){} //if i remove this, dati ceases to exist!! WTF?!
    delta1	= dati[0][2]-dati[0][0];//doppietto n=1
    s_delta1	= dati[0][1]+dati[0][3];
    delta2	= dati[1][2]-dati[1][0];//doppietto n=2
    s_delta2	= dati[1][1]+dati[1][3];
    delta3	= dati[2][2]-dati[2][0];//doppietto n=3
    s_delta3	= dati[2][1]+dati[2][3];
    Delta1	= dati[1][0]-dati[0][0];//interno-interno 1 e 2 
    s_Delta1	= dati[1][1]+dati[0][1];
    Delta2	= dati[2][0]-dati[1][0];//interno-interno 2 e 3
    s_Delta2	= dati[2][1]+dati[1][1];
    Delta3	= dati[1][2]-dati[0][2];//esterno-esterno 1 e 2 
    s_Delta3	= dati[1][3]+dati[0][3];
    Delta4	= dati[2][2]-dati[1][2];//esterno-esterno 2 e 3
    s_Delta4	= dati[2][3]+dati[1][3];
    mdelta = w_mean({delta1,delta2,delta3},{s_delta1,s_delta2,s_delta3})[0];
    s_mdelta = w_mean({delta1,delta2,delta3},{s_delta1,s_delta2,s_delta3})[1];
    
    mDelta = w_mean({Delta1,Delta2,Delta3},{s_Delta1,s_Delta2,s_Delta3})[0];
    s_mDelta = w_mean({Delta1,Delta2,Delta3,Delta4},{s_Delta1,s_Delta2,s_Delta3,s_Delta4})[1];
    std::cout<<"File "+path+": delta = "<<mdelta<<" +/- "<<s_mdelta<<"\tDelta = "<<mDelta<<" +/- "<<s_mDelta<<std::endl;
    delta.push_back(mdelta);
    s_delta.push_back(s_mdelta);
    Delta.push_back(mDelta);
    s_Delta.push_back(s_mDelta);
//I=7.2----------------------------------------------------------------------
    path = "../data/trasv_down/dt9_5.txt";
    dati = tsvreader(path,4);
    ncols = dati[0].size();
    nrows = dati.size();
    for(int i=0; i<nrows; i++){} //if i remove this, dati ceases to exist!! WTF?!
    delta1	= dati[0][2]-dati[0][0];//doppietto n=1
    s_delta1	= dati[0][1]+dati[0][3];
    delta2	= dati[1][2]-dati[1][0];//doppietto n=2
    s_delta2	= dati[1][1]+dati[1][3];
    delta3	= dati[2][2]-dati[2][0];//doppietto n=3
    s_delta3	= dati[2][1]+dati[2][3];
    Delta1	= dati[1][0]-dati[0][0];//interno-interno 1 e 2 
    s_Delta1	= dati[1][1]+dati[0][1];
    Delta2	= dati[2][0]-dati[1][0];//interno-interno 2 e 3
    s_Delta2	= dati[2][1]+dati[1][1];
    Delta3	= dati[1][2]-dati[0][2];//esterno-esterno 1 e 2 
    s_Delta3	= dati[1][3]+dati[0][3];
    Delta4	= dati[2][2]-dati[1][2];//esterno-esterno 2 e 3
    s_Delta4	= dati[2][3]+dati[1][3];
    mdelta = w_mean({delta1,delta2,delta3},{s_delta1,s_delta2,s_delta3})[0];
    s_mdelta = w_mean({delta1,delta2,delta3},{s_delta1,s_delta2,s_delta3})[1];
    
    mDelta = w_mean({Delta1,Delta2,Delta3},{s_Delta1,s_Delta2,s_Delta3})[0];
    s_mDelta = w_mean({Delta1,Delta2,Delta3,Delta4},{s_Delta1,s_Delta2,s_Delta3,s_Delta4})[1];
    std::cout<<"File "+path+": delta = "<<mdelta<<" +/- "<<s_mdelta<<"\tDelta = "<<mDelta<<" +/- "<<s_mDelta<<std::endl;
    delta.push_back(mdelta);
    s_delta.push_back(s_mdelta);
    Delta.push_back(mDelta);
    s_Delta.push_back(s_mDelta);
//------------------------------FINE CALCOLO DELTE---------------------------------
    delta.shrink_to_fit();
    s_delta.shrink_to_fit();
    Delta.shrink_to_fit();
    s_Delta.shrink_to_fit();
    delta[0]=1.e-12; //to avoid dividing by 0
    std::cout<<"Punti del fit\n";
    for(int i = 0;i<n;i++){
	d_D[i] = -delta[i]/Delta[i];
	s_d_D[i] = abs(d_D[i])*sqrt(pow(s_delta[i]/delta[i],2)+pow(s_Delta[i]/Delta[i],2));
	d_D[0]=0.;
        std::cout<<"B = ("<<B[i]*1000.<<"+/-"<<s_B[i]*1000.<<")mT,\t d/D = "<<d_D[i]<<"+/-"<<s_d_D[i]<<std::endl;
    }
    TGraphErrors* g = new TGraphErrors(n,B.data(),d_D.data(),s_B.data(),s_d_D.data());
    g->SetMarkerStyle(7);
    g->SetTitle("<#delta>/<#Delta> vs B, #Delta m_{L} = -1;B[T];#frac{<#delta>}{<#Delta>}[#]");
    TCanvas* c1 = new TCanvas("c1","c1",20,20,1098,732);
    c1->SetGrid();
    g->Draw("AP");
    TF1* f1=new TF1("f1","pol1",0,0.7);
    f1->SetParameter(0,0);
    f1->SetParameter(1,-0.4);
    //f1->Draw("same");    
    g->Fit("f1","R+","",0,0.7);
    std::cout<<"p-value = "<<f1->GetProb()<<std::endl;
    double mu = 1.4560; //indice di rifr
    double hc = 6.62607015e-34*299792458; //SI units
    double t = 0.003; //m, spessore
    double muB = -f1->GetParameter(1)*hc/2/mu/t;
    double s_muB = f1->GetParError(1)*hc/2/mu/t;
    std::cout<<"Stima magnetone mu_B = ("<<muB<<" +/- "<<s_muB<<")J/T\n";
    double truemuB = 9.2740100657e-24;
    std::cout<<"Migliore stima mu_B = "<<truemuB<<", Z = "<<(muB-truemuB)/s_muB<<"\n";
    std::cout<<"Compatibile !\n";
    double muup = 9.6e-24;
    double s_muup = 4.e-25;
    double mudown = muB;
    double s_mudown = s_muB;
    double mul=10.e-24;
    double s_mul = 3.e-25;
    std::vector<double> muBfinale = w_mean({muup,mudown,mul},{s_muup,s_mudown,s_mul});
    double s_muBfinale = w_mean({muup,mudown,mul},{s_muup,s_mudown,s_mul})[1];
    std::cout<<"Stima finale mu_B = ("<<muBfinale[0]<<" +/- "<<muBfinale[1]<<")J/T, Z = "<<(muBfinale[0]-truemuB)/muBfinale[1]<<"\n";
;
    
}

std::vector<std::vector<double>> tsvreader(string path, const int ncols){
    ifstream in_file(path.c_str());
    std::vector<std::vector<double>>res;
    while(in_file.good()){
	double x;
	std::vector<double>tmp(ncols);
	for(int j = 0; j<ncols; j++){
	    in_file >> x;
	    tmp[j] = x;
	}
	res.push_back(tmp);
    }
    res.pop_back();
    res.shrink_to_fit();
    return res;
}

std::vector<double> w_mean(std::vector<double> val, std::vector<double> s_val){
    std::vector<double> res(2);
    double sum=0;
    double sumw=0;
    int n = val.size();
    for(int i = 0; i < n; i++){
	sum += val[i]*pow(s_val[i],-2);
	sumw += pow(s_val[i],-2);
    }
    res[0] = sum/sumw;
    res[1] = sqrt(1/sumw);
    return res;
}
