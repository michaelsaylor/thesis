#include "/home/mike/thesis/root/include/Riostream.h"
void sigmavsw_v2(TString textfile) {
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
    double tfix = 1; //5.06, 9.99, 19.995, 29.855


    TCanvas *c1 = new TCanvas();
    c1->SetTitle(Form("J/psi_sigmavsW_t=%g", tfix));
    TGraph *graph1[4];
    TGraph *graph2[4];
    TGraph *graph3[4];

    TMultiGraph *mg[4];
    TAxis *axis[4];
    //mg->SetName("Name for graph");

    c1->Divide(2,2);

    for(Int_t i=0; i<4;)
    {
        c1->cd(i+1);
        c1->cd(i+1)->SetTicky(1);
        mg[i] = new TMultiGraph();
        axis[i] = mg[i]->GetXaxis();
        axis[i]->SetLimits(0,1000);

        mg[i]->SetTitle(Form("d#sigma/dW vs. W  x=%g",xfix[i]));
        mg[i]->GetYaxis()->SetTitle("Differential cross section d#sigma/dW [nb/GeV]");
        mg[i]->GetXaxis()->SetTitle("W [GeV]");


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



        Int_t pt = -1;
        double prevsigma = 0;
        double prevW = 0;
        double logprevsigma = 0;
        double logprevW = 0;
        double logderivative = 0;

        if(i>0)
        {
            infile.open(textfile);
        }

        while(infile.good())
        {
            infile >> xxx >> tt >> W >> quarkcont >> gluoncont >> dsigmadt;

            if(tt==tfix & xxx==xfix[i] & dsigmadt > 0)
            {
                graph1[i]->SetPoint(pt, W, quarkcont);
                graph2[i]->SetPoint(pt, W, gluoncont);
                graph3[i]->SetPoint(pt, W, dsigmadt);

                /*if(pt>-1)
                {
                    logderivative = (log(dsigmadt) - logprevsigma)/(log(W)-logprevW);
                    graph3[i]->SetPoint(pt, W, logderivative);      //plots the logarithmic derivative d(ln sigma)/d(ln W)
                }
                prevsigma = dsigmadt;
                prevW = W;
                logprevsigma = log(prevsigma);
                logprevW = log(prevW);
                */
                cout << W << "     " << quarkcont << "     " << gluoncont <<"     " << dsigmadt << "     " << logderivative << "\n";

                pt++;

            }

        }
        pt = -1;
        infile.close();


        cout << "x= " << xfix[i] << "      " << "t= " << tfix << "\n";




        mg[i]->Add(graph1[i]);
        mg[i]->Add(graph2[i]);
        mg[i]->Add(graph3[i]);

        mg[i]->Draw("AL");
        mg[i]->GetXaxis()->SetRangeUser(50.,1000.);
        //mg[i]->SetMaximum(1);

        //gPad->SetLogx();
        gPad->SetLogy();
        if(i==0)
        {gPad->BuildLegend();}

        c1->Modified();
        c1->Update();
        i++;
    }

}

