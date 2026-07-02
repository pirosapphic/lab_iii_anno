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
#include "TFitResult.h"

TH1D* histo_filler(string name, string title, string path); //funzione che cre oggetto istogramma da graficare
TGraphErrors* graph_errors_filler(string title, string path); //general purpose
std::vector<double> w_mean(std::vector<double> val, std::vector<double> s_val);

void picchi_monoc(){

    //picco a 540ua
    string path540 = "../../data_fc/calib_monoc/540.txt";

    TH1D* h540 = histo_filler("h540","540ua;assex;assey",path540); // oggetto istogramma da graficare 
    TCanvas* c540 = new TCanvas("c540","c540",20,20,1098,732);
    c540->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h540->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres540;

    fitres540 = h540->Fit("gaus","RS+","e1",950.,968.);
    std::cout<<"p-value "<<fitres540->Prob()<<std::endl;
    cout<< "stima picco a 540ua (viene!!): " << fitres540->Parameter(1) << " pm " <<  fitres540->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 550ua
    string path550 = "../../data_fc/calib_monoc/550.txt";

    TH1D* h550 = histo_filler("h550","550ua;assex;assey",path550); // oggetto istogramma da graficare  
    TCanvas* c550 = new TCanvas("c550","c550",20,20,1098,732);
    c550->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h550->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres550;

    fitres550 = h550->Fit("gaus","RS+","e1",979.,997.);
    std::cout<<"p-value "<<fitres550->Prob()<<std::endl;
    cout<< "stima picco a 550ua (viene!): " << fitres550->Parameter(1) << " pm " <<  fitres550->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 560ua
    string path560 = "../../data_fc/calib_monoc/560.txt";

    TH1D* h560 = histo_filler("h560","560ua;assex;assey",path560); // oggetto istogramma da graficare  
    TCanvas* c560 = new TCanvas("c560","c560",20,20,1098,732);
    c560->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h560->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres560;

    fitres560 = h560->Fit("gaus","RS+","e1",1009.,1025.);
    std::cout<<"p-value "<<fitres560->Prob()<<std::endl;
    cout<< "stima picco a 560ua (viene!!): " << fitres560->Parameter(1) << " pm " <<  fitres560->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 570ua
    string path570 = "../../data_fc/calib_monoc/570.txt";

    TH1D* h570 = histo_filler("h570","570ua;assex;assey",path570); // oggetto istogramma da graficare  
    TCanvas* c570 = new TCanvas("c570","c570",20,20,1098,732);
    c570->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h570->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres570;

    fitres570 = h570->Fit("gaus","RS+","e1",1038.,1054.);
    std::cout<<"p-value "<<fitres570->Prob()<<std::endl;
    cout<< "stima picco a 570ua (viene!): " << fitres570->Parameter(1) << " pm " <<  fitres570->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 580ua
    string path580 = "../../data_fc/calib_monoc/580.txt";

    TH1D* h580 = histo_filler("h580","580ua;assex;assey",path580); // oggetto istogramma da graficare  
    TCanvas* c580 = new TCanvas("c580","c580",20,20,1098,732);
    c580->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h580->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres580;

    fitres580 = h580->Fit("gaus","RS+","e1",1067.,1084.);
    std::cout<<"p-value "<<fitres580->Prob()<<std::endl;
    cout<< "stima picco a 580ua (viene!): " << fitres580->Parameter(1) << " pm " <<  fitres580->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza


//picco a 590ua
    string path590 = "../../data_fc/calib_monoc/590.txt";

    TH1D* h590 = histo_filler("h590","590ua;assex;assey",path590); // oggetto istogramma da graficare  
    TCanvas* c590 = new TCanvas("c590","c590",20,20,1098,732);
    c590->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h590->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres590;

    fitres590 = h590->Fit("gaus","RS+","e1",1096.,1113.);
    std::cout<<"p-value "<<fitres590->Prob()<<std::endl;
    cout<< "stima picco a 590ua (viene!!): " << fitres590->Parameter(1) << " pm " <<  fitres590->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 600ua
    string path600 = "../../data_fc/calib_monoc/600.txt";

    TH1D* h600 = histo_filler("h600","600ua;assex;assey",path600); // oggetto istogramma da graficare  
    TCanvas* c600 = new TCanvas("c600","c600",20,20,1098,732);
    c600->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h600->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres600;

    fitres600 = h600->Fit("gaus","RS+","e1",1126.,1142.);
    std::cout<<"p-value "<<fitres600->Prob()<<std::endl;
    cout<< "stima picco a 600ua (viene!): " << fitres600->Parameter(1) << " pm " <<  fitres600->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 610ua
    string path610 = "../../data_fc/calib_monoc/610.txt";

    TH1D* h610 = histo_filler("h610","610ua;assex;assey",path610); // oggetto istogramma da graficare  
    TCanvas* c610 = new TCanvas("c610","c610",20,20,1098,732);
    c610->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h610->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres610;

    fitres610 = h610->Fit("gaus","RS+","e1",1157.,1172.);
    std::cout<<"p-value "<<fitres610->Prob()<<std::endl;
    cout<< "stima picco a 610ua (viene!): " << fitres610->Parameter(1) << " pm " <<  fitres610->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 620ua
    string path620 = "../../data_fc/calib_monoc/620.txt";

    TH1D* h620 = histo_filler("h620","620ua;assex;assey",path620); // oggetto istogramma da graficare  
    TCanvas* c620 = new TCanvas("c620","c620",20,20,1098,732);
    c620->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h620->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres620;

    fitres620 = h620->Fit("gaus","RS+","e1",1185.,1201.);
    std::cout<<"p-value "<<fitres620->Prob()<<std::endl;
    cout<< "stima picco a 620ua (viene!): " << fitres620->Parameter(1) << " pm " <<  fitres620->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 630ua
    string path630 = "../../data_fc/calib_monoc/630.txt";

    TH1D* h630 = histo_filler("h630","630ua;assex;assey",path630); // oggetto istogramma da graficare  
    TCanvas* c630 = new TCanvas("c630","c630",20,20,1098,732);
    c630->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h630->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres630;

    fitres630 = h630->Fit("gaus","RS+","e1",1215.,1232.);
    std::cout<<"p-value "<<fitres630->Prob()<<std::endl;
    cout<< "stima picco a 630ua (viene!!): " << fitres630->Parameter(1) << " pm " <<  fitres630->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza




    //picco a 640ua
    string path640 = "../../data_fc/calib_monoc/640.txt";

    TH1D* h640 = histo_filler("h640","640ua;assex;assey",path640); // oggetto istogramma da graficare  
    TCanvas* c640 = new TCanvas("c640","c640",20,20,1098,732);
    c640->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h640->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres640;

    fitres640 = h640->Fit("gaus","RS+","e1",1245.,1261.);
    std::cout<<"p-value "<<fitres640->Prob()<<std::endl;
    cout<< "stima picco a 640ua (viene!!): " << fitres640->Parameter(1) << " pm " <<  fitres640->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 650ua
    string path650 = "../../data_fc/calib_monoc/650.txt";

    TH1D* h650 = histo_filler("h650","650ua;assex;assey",path650); // oggetto istogramma da graficare  
    TCanvas* c650 = new TCanvas("c650","c650",20,20,1098,732);
    c650->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h650->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres650;

    fitres650 = h650->Fit("gaus","RS+","e1",1275.,1291.);
    std::cout<<"p-value "<<fitres650->Prob()<<std::endl;
    cout<< "stima picco a 650ua (viene!): " << fitres650->Parameter(1) << " pm " <<  fitres650->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 660ua
    string path660 = "../../data_fc/calib_monoc/660.txt";

    TH1D* h660 = histo_filler("h660","660ua;assex;assey",path660); // oggetto istogramma da graficare  
    TCanvas* c660 = new TCanvas("c660","c660",20,20,1098,732);
    c660->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h660->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres660;

    fitres660 = h660->Fit("gaus","RS+","e1",1305.,1321.);
    std::cout<<"p-value "<<fitres660->Prob()<<std::endl;
    cout<< "stima picco a 660ua (viene!!): " << fitres660->Parameter(1) << " pm " <<  fitres660->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 670ua
    string path670 = "../../data_fc/calib_monoc/670.txt";

    TH1D* h670 = histo_filler("h670","670ua;assex;assey",path670); // oggetto istogramma da graficare  
    TCanvas* c670 = new TCanvas("c670","c670",20,20,1098,732);
    c670->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h670->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres670;

    fitres670 = h670->Fit("gaus","RS+","e1",1335.,1351.);
    std::cout<<"p-value "<<fitres670->Prob()<<std::endl;
    cout<< "stima picco a 670ua (viene!): " << fitres670->Parameter(1) << " pm " <<  fitres670->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 680ua
    string path680 = "../../data_fc/calib_monoc/680.txt";

    TH1D* h680 = histo_filler("h680","680ua;assex;assey",path680); // oggetto istogramma da graficare  
    TCanvas* c680 = new TCanvas("c680","c680",20,20,1098,732);
    c680->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h680->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres680;

    fitres680 = h680->Fit("gaus","RS+","e1",1366.,1383.);
    std::cout<<"p-value "<<fitres680->Prob()<<std::endl;
    cout<< "stima picco a 680ua (viene!!): " << fitres680->Parameter(1) << " pm " <<  fitres680->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 690ua
    string path690 = "../../data_fc/calib_monoc/690.txt";

    TH1D* h690 = histo_filler("h690","690ua;assex;assey",path690); // oggetto istogramma da graficare  
    TCanvas* c690 = new TCanvas("c690","c690",20,20,1098,732);
    c690->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h690->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres690;

    fitres690 = h690->Fit("gaus","RS+","e1",1397.,1412.);
    std::cout<<"p-value "<<fitres690->Prob()<<std::endl;
    cout<< "stima picco a 690ua (viene!): " << fitres690->Parameter(1) << " pm " <<  fitres690->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 700ua
    string path700 = "../../data_fc/calib_monoc/700.txt";

    TH1D* h700 = histo_filler("h700","700ua;assex;assey",path700); // oggetto istogramma da graficare  
    TCanvas* c700 = new TCanvas("c700","c700",20,20,1098,732);
    c700->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h700->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres700;

    fitres700 = h700->Fit("gaus","RS+","e1",1428.,1443.);
    std::cout<<"p-value "<<fitres700->Prob()<<std::endl;
    cout<< "stima picco a 700ua (viene!): " << fitres700->Parameter(1) << " pm " <<  fitres700->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza




    //picco a 710ua
    string path710 = "../../data_fc/calib_monoc/710.txt";

    TH1D* h710 = histo_filler("h710","710ua;assex;assey",path710); // oggetto istogramma da graficare  
    TCanvas* c710 = new TCanvas("c710","c710",20,20,1098,732);
    c710->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h710->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres710;

    fitres710 = h710->Fit("gaus","RS+","e1",1459.,1474.);
    std::cout<<"p-value "<<fitres710->Prob()<<std::endl;
    cout<< "stima picco a 710ua (viene!): " << fitres710->Parameter(1) << " pm " <<  fitres710->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 720ua
    string path720 = "../../data_fc/calib_monoc/720.txt";

    TH1D* h720 = histo_filler("h720","720ua;assex;assey",path720); // oggetto istogramma da graficare  
    TCanvas* c720 = new TCanvas("c720","c720",20,20,1098,732);
    c720->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h720->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres720;

    fitres720 = h720->Fit("gaus","RS+","e1",1489.,1505.);
    std::cout<<"p-value "<<fitres720->Prob()<<std::endl;
    cout<< "stima picco a 720ua (viene!!): " << fitres720->Parameter(1) << " pm " <<  fitres720->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 730ua
    string path730 = "../../data_fc/calib_monoc/730.txt";

    TH1D* h730 = histo_filler("h730","730ua;assex;assey",path730); // oggetto istogramma da graficare  
    TCanvas* c730 = new TCanvas("c730","c730",20,20,1098,732);
    c730->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h730->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres730;

    fitres730 = h730->Fit("gaus","RS+","e1",1520.,1535.);
    std::cout<<"p-value "<<fitres730->Prob()<<std::endl;
    cout<< "stima picco a 730ua (viene!): " << fitres730->Parameter(1) << " pm " <<  fitres730->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza



    //picco a 740ua
    string path740 = "../../data_fc/calib_monoc/740.txt";

    TH1D* h740 = histo_filler("h740","740ua;assex;assey",path740); // oggetto istogramma da graficare  
    TCanvas* c740 = new TCanvas("c740","c740",20,20,1098,732);
    c740->SetGrid();
    gStyle->SetOptStat(""); //toglie statistiche dell'istogramma dal grafico
    //gStyle->SetOptFit(1);
    h740->Draw("e1"); //e1-> errori sui punti;   histo-> disegna come un istogramma
    TFitResultPtr fitres740;

    fitres740 = h740->Fit("gaus","RS+","e1",1552.,1568.);
    std::cout<<"p-value "<<fitres740->Prob()<<std::endl;
    cout<< "stima picco a 730ua (viene!): " << fitres740->Parameter(1) << " pm " <<  fitres740->Parameter(2)*2.355 <<endl; //errore ampiezza a metà altezza
    

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
	double s_y = 0;
	if (yvec[i]>0) s_y=sqrt(yvec[i]);
	int theBin = theHisto->FindBin(xvec[i]);
        theHisto->SetBinContent(theBin,yvec[i]);
        theHisto->SetBinError(theBin,s_y);
	//std::cout<<i<<" "<<xvec[i]<<" "<<yvec[i]<<std::endl; //debug!
    }
    return theHisto;
}

TGraphErrors* graph_errors_filler(string title, string path){ //general purpose
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
	std::cout<<x<<"pm"<<s_x<<" "<<y<<"pm"<<s_y<<std::endl; //debug!
    }
    std::cout<<"good "<<in_file.good()<<std::endl;

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
