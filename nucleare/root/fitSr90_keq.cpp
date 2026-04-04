#include <iostream>
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

TH1D* histo_filler(string name, string title, string path); //general purpose
std::vector<double> w_mean(std::vector<double> val, std::vector<double> s_val);

void fitSr90_keq(string input = "../data_sorg/A8_Sr90_decentrato.txt"){    
    TH1D* histoTot = histo_filler("histoTot","Spettro energetico ^{90}Sr, decentrato;ADC channels [CHN];Conteggi[#]",input);
    TCanvas* ctot = new TCanvas("c1","c1",20,20,1098,732);
    gStyle->SetOptStat(0);
    ctot->SetGrid();
    histoTot->Draw("e1");
    TH1D *histoEp = (TH1D*)histoTot->Clone("histoEp");
    histoEp->Rebin(16);
    histoEp->SetTitle("Plot di Kurie, ^{90}Sr decentrato;ADC channels [CHN];#frac{#sqrt{counts}}{CHN} [#frac{#}{CHN}]");
    double xMinep = 37000;
    double xMaxep = 43000;
    int nbins = histoEp->GetNbinsX();
    double binwidth = histoEp->GetBinWidth(1);
    double xlow = histoEp->GetBinLowEdge(1);
    double xmax = histoEp->GetBinLowEdge(nbins)+binwidth;
    for( double x = xlow; x<xmax; x+=binwidth ){
	int theBin = histoEp->FindBin(x);
	if(x < xMaxep and x > xMinep){
	    //continue;
	    double y = histoEp->GetBinContent(theBin);
	    double s_y = histoEp->GetBinError(theBin);
	    double newy = sqrt(y)/x;
	    double s_newy = sqrt(pow(newy/x *binwidth,2)+pow(0.5*newy/y * s_y,2));
	    histoEp->SetBinContent(theBin,newy);
	    histoEp->SetBinError(theBin,s_newy);
	}
	else{
	    histoEp->SetBinContent(theBin,0);
	    histoEp->SetBinError(theBin,0);
	}
    }
    histoEp->SetAxisRange(xMinep,xMaxep+1000);
    TCanvas* cep = new TCanvas("c2","c2",20,20,1098,732);
    cep->cd();
    gStyle->SetOptStat(0);
    cep->SetGrid();
    histoEp->Draw("e1");
    

    TF1* f1 = new TF1("f1","pol1",38000,xMaxep);
    f1->SetParameter(0,0);
    f1->SetParameter(1,8e-6);
    histoEp->Fit("f1","R+","e1");
    std::cout << "Chi^2:" <<f1->GetChisquare() << ", number of DoF: " << f1->GetNDF();
    std::cout << " (Probability: " << f1->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    double m = f1->GetParameter(1);
    double s_m = f1->GetParError(1);
    double q = f1->GetParameter(0);
    double s_q = f1->GetParError(0);
    double CH0 = -q/m;
    double s_CH0 = sqrt(pow(CH0/m * s_m,2)+pow(CH0/q * s_q,2));
    std::cout<<"CHN_0 = ("<<CH0<<" +/- "<<s_CH0<<")CHN\n";
    //T_0 = (2.28 pm 0.04) MeV
    double k = 2.28/CH0;
    double s_k = sqrt(pow(k/CH0 * s_CH0,2)+pow(k/2.28 * 0.04,2));
   std::cout<<"keq = ("<< k<<" +/- "<<s_k<<")MeV/CHN\n";

    /*
    //endpoint
    TCanvas* cep = new TCanvas("c2","c2",20,20,800,600);
    cep->SetGrid();
    cep->cd();
    histoEp->Rebin(31);
    for (double x = xMinep; x<=xMaxep; x+=binWidth){
	int theBin = histoEp->FindBin(x);
	double s_x = histoEp->GetBinWidth(theBin);
	double y = histoEp->GetBinContent(theBin);
	double s_y = histoEp->GetBinError(theBin);
	double new_y = sqrt(y)/x;
	double s_new_y = sqrt(pow(new_y/x *s_x,2)+pow(0.5*new_y/y * s_y,2));
	histoEp->SetBinContent(theBin,new_y);
	histoEp->SetBinError(theBin,s_new_y);
    }
    gStyle->SetOptStat(0);
    histoEp->Draw("e1");
    TF1* f1 = new TF1("f1","pol1",57000,66200);
    f1->SetParameter(0,0);
    f1->SetParameter(1,8e-6);
    histoEp->Fit("f1","","e1",xMinep,xMaxep);
        std::cout << "Chi^2:" <<f1->GetChisquare() << ", number of DoF: " << f1->GetNDF();
    std::cout << " (Probability: " << f1->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    double m = f1->GetParameter(1);
    double s_m = f1->GetParError(1);
    double q = f1->GetParameter(0);
    double s_q = f1->GetParError(0);
    double CH0 = -q/m;
    double s_CH0 = sqrt(pow(CH0/m * s_m,2)+pow(CH0/q * s_q,2));
    std::cout<<"CHN_0 = ("<<CH0<<" +/- "<<s_CH0<<")CHN\n";
    //T_0 = (2.28 pm 0.04) MeV
    double k = 2.28/CH0;
    double s_k = sqrt(pow(k/CH0 * s_CH0,2)+pow(k/2.28 * 0.04,2));
   std::cout<<"k = ("<< k<<" +/- "<<s_k<<")MeV/CHN\n";
    */
}

TH1D* histo_filler(string name, string title, string path){ //general purpose
    //fills histo with name name and title title with data from path.
    ifstream in_file(path.c_str());
    double x = 666;
    double y = 777;
    std::vector<double> xvec;
    std::vector<double> yvec;
    int count = 0;
    while(in_file.good()){
	in_file >> x >> y;
	xvec.push_back(double(x));
	yvec.push_back(double(y));
	count++;
    }
    xvec.shrink_to_fit();
    yvec.shrink_to_fit();
    int nbins = xvec.size()-1; //dark magic!
    double binsize = xvec[2]-xvec[1];
    double xMin = xvec[0]-(binsize/2);
    double xMax = xvec[nbins-1]+(binsize/2);
    TH1D* theHisto = new TH1D(name.c_str(),title.c_str(),nbins,xMin,xMax);
    for(int i = 0; i < nbins; i++){
	double s_y = sqrt(yvec[i]);
	int theBin = theHisto->FindBin(xvec[i]);
        theHisto->SetBinContent(theBin,yvec[i]);
        theHisto->SetBinError(theBin,s_y);
	//std::cout<<i<<" "<<xvec[i]<<" "<<yvec[i]<<std::endl; //debug!
    }
    return theHisto;
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

