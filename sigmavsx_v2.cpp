#include "/home/mike/thesis/root/include/Riostream.h"
void sigmavsx_v2(TString textfile) {
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
    double Wfix = 1000; //991.673, 492.858
    double tfix[4] = {5, 10, 20, 30}; //5.06, 9.99, 19.995, 29.855


    TCanvas *c1 = new TCanvas();
    TGraph *graph1[4];
    TGraph *graph2[4];
    TGraph *graph3[4];

    TMultiGraph *mg[4];
    //mg->SetName("Name for graph");


    c1->Divide(2,2);

    for(Int_t i=0; i<4;)
    {
        c1->cd(i+1);
        c1->cd(i+1)->SetTicky(1);
        mg[i] = new TMultiGraph();

        mg[i]->SetTitle(Form("d#sigma/dx vs. x  t=%g",tfix[i]));
        mg[i]->GetYaxis()->SetTitle("Differential cross section d#sigma/dx [nb]");
        mg[i]->GetXaxis()->SetTitle("x");


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

            if(W==Wfix & tt==tfix[i] /*& dsigmadt!=0*/)
            {
                graph1[i]->SetPoint(pt, xxx, quarkcont);
                graph2[i]->SetPoint(pt, xxx, gluoncont);
                graph3[i]->SetPoint(pt, xxx, dsigmadt);
                cout << xxx << "     " << quarkcont << "     " << gluoncont <<"     " << dsigmadt << "\n";

                pt++;

            }

        }
        pt = 0;
        infile.close();


        cout << "t= " << tfix[i] << "      " << "W= " << Wfix << "\n";




        mg[i]->Add(graph1[i]);
        mg[i]->Add(graph2[i]);
        mg[i]->Add(graph3[i]);

        mg[i]->Draw("AL");
        //mg[i]->GetXaxis()->SetLimits(0.,1.);
        //mg[i]->SetMaximum(1);

        gPad->SetLogy();
        gPad->SetLogx();
        if(i==0)
        {gPad->BuildLegend();}

        //c1->Modified();
        //c1->Update();
        i++;
    }

}
