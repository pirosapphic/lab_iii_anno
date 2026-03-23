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


TH1F* ReadHisto(float xMin,float xMax,float binWidth, string input, string name) {

    ifstream parInput(input.c_str());

    float x, y;


    TH1F* h = new TH1F(name.c_str(), name.c_str(),
                       (xMax - xMin)/binWidth, xMin, xMax);

    while (parInput >> x >> y) {
        int bin = h->FindBin(x);
        h->SetBinContent(bin, y);
        h->SetBinError(bin, sqrt(y)); 
    }

    return h;
}

void Ass_carta() {


    // avevo 20 bin, ho raddoppiato la larghezza e ora ne ho 10
    TH1F* h0 = ReadHisto(647.275,771.725,13.1,"../data_sorg/A8_count_vuoto.txt", "h0");
    TH1F* h3 = ReadHisto(518.925,668.075,15.7,"../data_sorg/A8_count_carta3.txt", "h3");
    TH1F* h6 = ReadHisto(444.1,561.9,12.4,"../data_sorg/A8_count_carta6.txt", "h6");
    TH1F* h9 = ReadHisto(365.525,499.475,14.1,"../data_sorg/A8_count_carta9.txt", "h9");
    TH1F* h12 = ReadHisto(317.075,433.925,12.3,"../data_sorg/A8_count_carta12.txt", "h12");
    TH1F* h15 = ReadHisto(289.6, 388.4,10.4,"../data_sorg/A8_count_carta15.txt", "h15");
    TH1F* h18 = ReadHisto(253.325,341.675,9.3,"../data_sorg/A8_count_carta18.txt", "h18");
    TH1F* h20 = ReadHisto(230.325,318.675,9.3,"../data_sorg/A8_count_carta20.txt", "h20");   

    h0->SetLineColor(kCyan);
    h3->SetLineColor(kGreen);
    h6->SetLineColor(kBlue);
    h9->SetLineColor(kRed);
    h9->SetLineColor(kYellow);
    h9->SetLineColor(kMagenta);
    h9->SetLineColor(kOrange);
    h9->SetLineColor(kBlack);
    
    TCanvas* c0 = new TCanvas("c0","c0",800,600);
    h0->Draw("hist");

    TCanvas* c3 = new TCanvas("c3","c3",800,600);
    h3->Draw("hist");

    TCanvas* c6 = new TCanvas("c6","c6",800,600);
    h6->Draw("hist");

    TCanvas* c9 = new TCanvas("c9","c9",800,600);
    h9->Draw("hist");

    TCanvas* c12 = new TCanvas("c12","c12",800,600);
    h12->Draw("hist");

    TCanvas* c15 = new TCanvas("c15","c15",800,600);
    h15->Draw("hist");

    TCanvas* c18 = new TCanvas("c18","c18",800,600);
    h18->Draw("hist");

    TCanvas* c20 = new TCanvas("c20","c20",800,600);
    h20->Draw("hist");

    /* in generale, il rate è definito come numero di impulsi che superano una certa
soglia (threshold, in mV e sempre intesa come negativa, dato che si tratta di elettroni)
in un intervallo di tempo definito (gate, in ms).*/
}