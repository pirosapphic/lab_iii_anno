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

void normalizzazione(){ 
    string path= "../data_fc/dati_I_calib.txt";
    auto [x,s_x,I,s_I]= lettura(path);
    string path1= "../data_fc/dati_V0_calib.txt";
    auto [l,s_l,V0,s_V0]= lettura(path1);

   
    vector<double> In,s_In;
    double y,s_y;
    

    for (int i = 0; i < x.size(); i++)
    {   
        y= I[i]/V0[i];
        s_y =y*sqrt(pow(s_I[i]/I[i],2)+ pow(s_V0[i]/V0[i],2));
        
        In.push_back(y);
        s_In.push_back(s_y);
        cout<< x[i] <<" "<< s_x[i] <<" "<< In[i] <<" "<<s_In[i]<<endl;
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
