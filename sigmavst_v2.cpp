#include "/home/mike/thesis/root/include/Riostream.h"
void sigmavst_v2(TString textfile) {
//  Read data from an ascii file and create a root file with an histogram and an ntuple.
//   see a variant of this macro in basic2.C
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
    double xfix[4] = {0.01,0.05,0.1,0.3};
    double Wfix = 130; //492.858, 991.673


    TCanvas *c1 = new TCanvas();
    TGraph *graph1[4];
    TGraph *graph2[4];
    TGraph *graph3[4];

    TAxis *axis[4];


    TMultiGraph *mg[4];
    //mg->SetName("Name for graph");


    c1->Divide(2,2);

    for(Int_t i=0; i<4;)
    {
        c1->cd(i+1);
        c1->cd(i+1)->SetTicky(1);
        mg[i] = new TMultiGraph();
        axis[i] = mg[i]->GetXaxis();
        axis[i]->SetLimits(0,30);

        mg[i]->SetTitle(Form("d#sigma/dt vs. t  x=%g",xfix[i]));
        mg[i]->GetYaxis()->SetTitle("Differential cross section d#sigma/dx [nb]");
        mg[i]->GetXaxis()->SetTitle("-t");


        graph1[i] = new TGraph();
        graph1[i]->SetLineWidth(3);
        graph1[i]->SetTitle("Quark contribution");
        graph1[i]->SetDrawOption("AP");
        graph1[i]->SetLineColor(1);

        graph2[i] = new TGraph();
        graph2[i]->SetLineWidth(3);
        graph2[i]->SetTitle("Gluon contribution");
        graph2[i]->SetDrawOption("AP");
        graph2[i]->SetLineColor(4);


        graph3[i] = new TGraph();
        graph3[i]->SetLineWidth(3);
        graph3[i]->SetTitle("Total contribution");
        graph3[i]->SetDrawOption("AP");
        graph3[i]->SetLineColor(2);




        Int_t pt = 0;
        if(i>0)
        {
            infile.open(textfile);
        }

        while(infile.good())
        {
            infile >> xxx >> tt >> W >> quarkcont >> gluoncont >> dsigmadt;

            if(W==Wfix & xxx==xfix[i] & dsigmadt!=0)
            {
                graph1[i]->SetPoint(pt, tt, quarkcont);
                graph2[i]->SetPoint(pt, tt, gluoncont);
                graph3[i]->SetPoint(pt, tt, dsigmadt);
                cout << tt << "     " << quarkcont << "     " << gluoncont <<"     " << dsigmadt << "\n";

                pt++;

            }

        }
        pt = 0;
        infile.close();


        cout << "t= " << xfix[i] << "      " << "W= " << Wfix << "\n";




        mg[i]->Add(graph1[i]);
        mg[i]->Add(graph2[i]);
        mg[i]->Add(graph3[i]);

        mg[i]->GetXaxis()->SetRangeUser(0.,30);
        mg[i]->Draw("AL");
        //mg[i]->SetMaximum(1);
        //mg[i]->GetXaxis()->SetAxisRange(1,30);

        gPad->SetLogy();
        if(i==0)
        {gPad->BuildLegend();}

        c1->Modified();
        c1->Update();
        i++;
    }

}
