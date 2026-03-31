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

//Vbias = 55V
//Gampl = 29dB
//I = 2.8

TH1D* histo_filler(string name, string title, string path); //general purpose
std::vector<double> w_mean(std::vector<double> val, std::vector<double> s_val);

void fitLED_I2_8(string input = "../../data_SiPM/LED/I/A8_LED5529"){    
    //firstly, we draw and fit the histograms, then we calculate <deltapp>;
//--------------------------------------------------------------------------------
    
    string Int = "2-8.txt";
    string path = input+Int;
    //path = "../../data_SiPM/LED/gain/A8_LED55252-5.txt";
    TH1D* histo27 = histo_filler("histo27","Spettro Energetico, G=29dB, V_{bias}=55V, LED 2.8;ADC channels [CHN];Conteggi[#]",path);
    TCanvas* c27 = new TCanvas("c27","c27",20,20,1098,732);
    c27->SetGrid();
    gStyle->SetOptStat(0);
    histo27->SetLineColor(kBlack);
    //histo27->Rebin(2);
    histo27->Draw("e1");
    int npeaks = 6;
    double xMin = -49.5;
    double xMax = 1000;
    
    //histo27->SetAxisRange(xMin,xMax);
    std::vector<double> peak(npeaks);
    std::vector<double> s_peak(npeaks);
    std::vector<double> sigma(npeaks);
    std::vector<double> s_sigma(npeaks);
    std::vector<double> norm(npeaks);
    std::vector<double> s_norm(npeaks);

    //PICCO 0
    TF1* gaus0 = new TF1("gaus0","gaus",-20,20);
    histo27->Fit("gaus0","R+","e1",-20,20);
    std::cout << "Chi^2:" <<gaus0->GetChisquare();
    std::cout<< ", number of DoF: " << gaus0->GetNDF();
    std::cout << " (Probability: " << gaus0->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    norm[0] = gaus0->GetParameter(0); 
    s_norm[0] = gaus0->GetParError(0);
    peak[0] = gaus0->GetParameter(1); 
    s_peak[0] = gaus0->GetParError(1);
    sigma[0] = gaus0->GetParameter(2); 
    s_sigma[0] = gaus0->GetParError(2);
    
    //PICCO 1
    TF1* gaus1 = new TF1("gaus1","gaus",205,260);
    gaus1->SetLineColor(kYellow);
    histo27->Fit("gaus1","R+","e1");
    std::cout << "Chi^2:" <<gaus1->GetChisquare();
    std::cout<< ", number of DoF: " << gaus1->GetNDF();
    std::cout << " (Probability: " << gaus1->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    norm[1] = gaus1->GetParameter(0); 
    s_norm[1] = gaus1->GetParError(0);
    peak[1] = gaus1->GetParameter(1); 
    s_peak[1] = gaus1->GetParError(1);
    sigma[1] = gaus1->GetParameter(2); 
    s_sigma[1] = gaus1->GetParError(2);
    
    //PICCO 2
    TF1* gaus2 = new TF1("gaus2","gaus",425,500);
    gaus2->SetLineColor(kGreen);
    histo27->Fit("gaus2","R+","e1");
    std::cout << "Chi^2:" <<gaus2->GetChisquare();
    std::cout<< ", number of DoF: " << gaus2->GetNDF();
    std::cout << " (Probability: " << gaus2->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    norm[2] = gaus2->GetParameter(0); 
    s_norm[2] = gaus2->GetParError(0);
    peak[2] = gaus2->GetParameter(1); 
    s_peak[2] = gaus2->GetParError(1);
    sigma[2] = gaus2->GetParameter(2); 
    s_sigma[2] = gaus2->GetParError(2);

    //PICCO 3
    TF1* gaus3 = new TF1("gaus3","gaus",650,730);
    gaus3->SetLineColor(kCyan);
    histo27->Fit("gaus3","R+","e1");
    std::cout << "Chi^2:" <<gaus3->GetChisquare();
    std::cout<< ", number of DoF: " << gaus3->GetNDF();
    std::cout << " (Probability: " << gaus3->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    norm[3] = gaus3->GetParameter(0); 
    s_norm[3] = gaus3->GetParError(0);
    peak[3] = gaus3->GetParameter(1); 
    s_peak[3] = gaus3->GetParError(1);
    sigma[3] = gaus3->GetParameter(2); 
    s_sigma[3] = gaus3->GetParError(2);
  
    //PICCO 4
    TF1* gaus4 = new TF1("gaus4","gaus",870,973);
    gaus4->SetLineColor(kBlue);
    histo27->Fit("gaus4","R+","e1");
    std::cout << "Chi^2:" <<gaus4->GetChisquare();
    std::cout<< ", number of DoF: " << gaus4->GetNDF();
    std::cout << " (Probability: " << gaus4->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    norm[4] = gaus4->GetParameter(0); 
    s_norm[4] = gaus4->GetParError(0);
    peak[4] = gaus4->GetParameter(1); 
    s_peak[4] = gaus4->GetParError(1);
    sigma[4] = gaus4->GetParameter(2); 
    s_sigma[4] = gaus4->GetParError(2);


    //PICCO 5
    TF1* gaus5 = new TF1("gaus5","gaus",1090,1210);
    gaus5->SetLineColor(kMagenta);
    histo27->Fit("gaus5","R+","e1");
    std::cout << "Chi^2:" <<gaus5->GetChisquare();
    std::cout<< ", number of DoF: " << gaus5->GetNDF();
    std::cout << " (Probability: " << gaus5->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    norm[5] = gaus5->GetParameter(0); 
    s_norm[5] = gaus5->GetParError(0);
    peak[5] = gaus5->GetParameter(1); 
    s_peak[5] = gaus5->GetParError(1);
    sigma[5] = gaus5->GetParameter(2); 
    s_sigma[5] = gaus5->GetParError(2);

    std::vector<double>deltapp(npeaks-1);
    std::vector<double>s_deltapp(npeaks-1);
    std::cout<<"Stime di Delta_pp:\n";
    for(int i = 0; i < deltapp.size(); i++){
	deltapp[i] = peak[i+1]-peak[i];
	s_deltapp[i] = sqrt(pow(s_peak[i+1],2)+pow(s_peak[i],2));
	std::cout<<"("<<deltapp[i]<<" +/- "<<s_deltapp[i]<<")CHN\n";
    }
    //these are not compatible upon a Gauss test... but we take the mean anyways
    std::cout<<"OCCHIO, non sono compatibili con test Z!!\nMa faccio lo stesso la\n";
    std::vector<double> respp = w_mean(deltapp,s_deltapp);
    double meanpp = respp[0];
    double s_meanpp = respp[1];
    std::cout<<"Media pesata Deltapp = ("<<meanpp<<" +/- "<<s_meanpp<<")CHN\n";
    std::cout<<"Calcolo del numero di fotoni incidenti\n";

    std::cout<<"I) Metodo approssimativo \\mu_ph = <CH>/<deltapp>\n";
    double meanh = histo27->GetMean();
    double s_meanh = histo27->GetMeanError();
    std::cout<<"Media histo <CH> = ("<<meanh<<" +/- "<<s_meanh<<")CHN\n";
    double mu_app = meanh/meanpp;
    double s_mu_app = mu_app*sqrt(pow(s_meanh/meanh,2)+pow(s_meanpp/meanpp,2));
    std::cout<<"Valore approx \\mu_ph = "<<mu_app<<" +/- "<<s_mu_app<<"\n";
    
    std::cout<<"II) Metodo accurato\n";
    std::cout<<"si prendono gli integrali dei picchi gaussiani (=#eventi) e se ne fa un fit poissoniano\n";
    std::cout<<"QUESTO METODO NON VIENE PERCHÉ CI SONO POCHI PICCHI.\n\n";

    std::cout<<"verifica dell'andamento lineare tra Npicco e la sua varianza\n";
    std::cout<<"non viene, ma non viene nemmeno agli altri: si puo dire che i primi dati sono circa lineari\n";
    std::vector<double> var(npeaks-1);
    std::vector<double> s_var(npeaks-1);
    std::vector<double> N(npeaks-1);
    std::vector<double> s_N(npeaks-1);
    for(int i = 1; i < npeaks; i++){
	int j = i-1;
	var[j] = pow(sigma[i],2);
	s_var[j] = 2*sigma[i]*s_sigma[i];
	N[j] = i;
	s_N[j] = 0;
    }
    TCanvas *c1 = new TCanvas("c1","c1",20,20,800,600);
    c1->cd();
    c1->SetGrid();
    TGraphErrors* g1 = new TGraphErrors(npeaks,N.data(),var.data(),s_N.data(),s_var.data());
    g1->SetTitle("Varianze vs N. di picco;N[#];#sigma_{N}^{2} [CHN^{2}]");
    g1->Draw("AP");
    //g1->GetXaxis()->SetRange(-1,100);
    TF1* f1 = new TF1("f1","pol1",1,5);
    g1->Fit("f1","R+");
    std::cout << "Chi^2:" <<f1->GetChisquare();
    std::cout<< ", number of DoF: " << f1->GetNDF();
    std::cout << " (Probability: " << f1->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    
/* //qui c'è metodo integrale che non funziona
    std::vector<double> I(npeaks);
    std::vector<double> s_I(npeaks);
    for(int i = 0; i<npeaks; i++){
	I[i] = norm[i]*sigma[i]*sqrt(2*3.1415926);
	s_I[i] = I[i]*sqrt(pow(s_norm[i]/norm[i],2)+pow(s_sigma[i]/sigma[i],2));
    }
    
    TCanvas *c1 = new TCanvas("c1","c1",20,20,800,600);
    c1->cd();
    c1->SetGrid();

    TH1D *h1 = new TH1D("h1", "h1 title", 8, -0.5, 7.5);
    for(int i = 0; i<npeaks; i++){
	int theBin = h1->FindBin(i);
	h1->SetBinContent(theBin,norm[i]);
	h1->SetBinError(theBin,s_norm[i]);
    }
    h1->Draw("e1");
    
    TF1* f1 = new TF1("f1","[0]*TMath::Poisson(x,[1])",0,4.5);
    f1->SetNpx(10000);
    f1->SetParameter(0,40000);
    f1->SetParameter(1,0.8);
    //f1->Draw("same");
    h1->Fit("f1","R+","e1");
    std::cout << "Chi^2:" <<f1->GetChisquare();
    std::cout<< ", number of DoF: " << f1->GetNDF();
    std::cout << " (Probability: " << f1->GetProb() << ")." << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
*/
}

TH1D* histo_filler(string name, string title, string path){ //general purpose
    //fills histo with name name and title title with data from path.
    ifstream in_file(path.c_str());
    double x = 666;
    double y = 777;
    std::vector<double> xvec;
    std::vector<double> yvec;
    while(in_file.good()){
	in_file >> x >> y;
	xvec.push_back(double(x));
	yvec.push_back(double(y));
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
