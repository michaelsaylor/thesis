#include "/home/mike/thesis/root/include/Riostream.h"
void dNdW(TString textfile) {
//Read data from a .txt or .dat file and plots cross section as a function of transverse
//momentum transfer with constraints on energy and x
//Author: Rene Brun


// read file $ROOTSYS/tutorials/tree/basic.dat
// this file has 3 columns of float data
    TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
    dir.ReplaceAll("basic.C","");
    dir.ReplaceAll("/./","/");
    ifstream infile(textfile);

    double W;
    double xxx;
    double tt;
    double dpy;
    double dsigmadt;
    double gluoncont;
    double quarkcont;
    double dNdW;
    double dNdWquark;
    double dNdWgluon;
    double xtintegral;


    TGraph *g1 = new TGraph();
    g1->SetLineWidth(3);
    g1->SetTitle("Quark contribution");
    g1->SetDrawOption("AP");
    g1->SetLineColor(1);

    TGraph *g2 = new TGraph();
    g2->SetLineWidth(3);
    g2->SetTitle("Gluon contribution");
    g2->SetDrawOption("AP");
    g2->SetLineColor(4);

    TGraph *g3 = new TGraph();
    g3->SetLineWidth(3);
    g3->SetTitle("Total cross section");
    g3->SetDrawOption("AP");
    g3->SetLineColor(2);

    TCanvas *c1 = new TCanvas();
    TMultiGraph *mg = new TMultiGraph();
    mg->SetName("Name for graph");
    mg->SetTitle("dN/dW vs. W     x#in(0.1,0.3) ;Energy W [GeV] ;Number of events dN/dW [GeV^{-1}]");

    Int_t pt = 0;

    while(infile.good())
    {
        //infile >> W >> dNdWquark >> dNdWgluon >> dNdW;
        infile >> W >> dNdW >> dNdWgluon >> dNdWquark;

        if(dNdW > 0 /*& W != 50 & W!= 140*/) //0.150153 and 1.5
        {
            g1->SetPoint(pt, W, dNdWquark);
            g2->SetPoint(pt, W, dNdWgluon);
            g3->SetPoint(pt, W, dNdW);
            cout << W << "     " << "     " << dNdWquark << "     " << dNdWgluon <<"     " << dNdW << "     " << "\n";
            pt++;
        }
    }
    //cout << "x= " << xfix << "      " << "W= " << Wfix << "\n";

    infile.close();
    mg->Add(g1);
    mg->Add(g2);
    mg->Add(g3);
    c1->SetLogy();


    //g1->GetXaxis()->SetLimits(0.5,12.5);
    mg->Draw("AL");
    c1->BuildLegend();
}