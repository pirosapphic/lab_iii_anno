#include <iostream> // scrivere a video
#include <fstream> //legge i file
#include <vector> 
#include <math.h> //funzioni matematiche   

#include "TH1F.h" //istogramma
#include "TF1.h" //funzioni per i fit
#include "TLine.h" 
#include "TRandom3.h"
#include "TMath.h"
#include "TGraphErrors.h" //grafici con barre d'errore
#include "TMultiGraph.h" // + rafici nella stessa canvas
#include "TCanvas.h" // tela per i grafici
#include "TFile.h" 
#include "TFitResult.h"

TGraphErrors* graph_errors_filler(string title, string path); //general purpose
std::vector<double> w_mean(std::vector<double> val, std::vector<double> s_val);

void function (){
    string path = "../data/discesa1.txt";  
    string title = "Y vs X;X axis;Y axis";
    //title = "B vs I, "+title+";I[A];B[mT]";
    TGraphErrors* g = graph_errors_filler(title,path); //creazione del grafico
    g->SetMarkerStyle(7);
    TCanvas* c = new TCanvas("c","c",20,20,1098,732); //oggetto su cui disegno
    // g->GetXaxis()->SetTitle("X axis");
    // g->GetYaxis()->SetTitle("Y axis");
    c->SetGrid();
    g->Draw("AP"); //A -> assi P->punti
    TFitResultPtr fitres = g->Fit("pol1","RS+","",-0.,9.8);
    std::cout<<"p-value "<<fitres->Prob()<<std::endl;
}
// costruisce l'oggetto da graficare. file text : x  \\  s_x  \\  y  \\   y_s
TGraphErrors* graph_errors_filler(string title, string path){ 
    //fills grapherrors with title title with data from path.
    ifstream in_file(path.c_str());
    double x = 666;
    double y = 777;
    double s_x = 888;
    double s_y = 999;
    std::vector<double> xvec;
    std::vector<double> s_xvec;
    std::vector<double> yvec;
    std::vector<double> s_yvec;
    while(in_file.good()){
	in_file >> x >> s_x >> y >> s_y;
	xvec.push_back(double(x));
	yvec.push_back(double(y));
	s_xvec.push_back(double(s_x));
	s_yvec.push_back(double(s_y));
	//std::cout<<x<<"pm"<<s_x<<" "<<y<<"pm"<<s_y<<std::endl; //debug!
    }
    //std::cout<<"good "<<in_file.good()<<std::endl;

    xvec.pop_back();
    yvec.pop_back(); //CHECK IF IT IS WRONG!!!!!!!!!!
    s_xvec.pop_back();
    s_yvec.pop_back();
    xvec.shrink_to_fit();
    yvec.shrink_to_fit();
    s_xvec.shrink_to_fit();
    s_yvec.shrink_to_fit();
    TGraphErrors* theGraph = new TGraphErrors(xvec.size(),xvec.data(),yvec.data(),s_xvec.data(),s_yvec.data());
    theGraph->SetTitle(title.c_str());
    return theGraph;
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