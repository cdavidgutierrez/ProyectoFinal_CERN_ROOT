#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TChain.h>

#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooArgusBG.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "RooVoigtian.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"

#include "RooCBShape.h"

#include "TPaveText.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TLatex.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooNumIntConfig.h"
#include "RooRandom.h"

using namespace RooFit;

void ToyMCHiggs()
{
    

    Int_t seed=3;
    Int_t n_total=20;

  
    // ---------------------------------------------
    // READ WORKSPACE.
    // ---------------------------------------------
    
    TFile *f = new TFile("HiggsWs.root");
    
    // Read workspace from file
    RooWorkspace *w = (RooWorkspace *)f->Get("workspace");
    
    // Model and data from workspace
    RooRealVar *hgg_mass = w->var("CMS_hgg_mass");
    RooAbsPdf *model = w->pdf("model");

    // Params from workspace
    RooRealVar *MH = w->var("MH");
    RooRealVar *meanMass = new RooRealVar("meanMass", "meanMass", MH->getVal());


    model->Print("t");


    RooRandom::randomGenerator()->SetSeed(seed);
    
    // Verbose 0:
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    RooMsgService::instance().setSilentMode(kTRUE);
    RooMsgService::instance().setStreamStatus(1,false);
    RooMsgService::instance().getStream(1).removeTopic(Integration);  
    RooMsgService::instance().getStream(1).removeTopic(Minimization);  
    RooMsgService::instance().getStream(1).removeTopic(Fitting);  
    RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
    RooMsgService::instance().getStream(1).removeTopic(Optimization);
    RooMsgService::instance().getStream(1).removeTopic(ObjectHandling);
    RooMsgService::instance().getStream(1).removeTopic(Eval);
    RooMsgService::instance().Print();


    // Params to fit 
    RooRealVar Mmean("mean"," Mass mean", 120, 125.3, 127.0);
    RooDataSet* MupullData = new RooDataSet("MupullData", "Mupull dataset", RooArgSet(Mmean));


    // ---------------------------------------------
    // FITTING & TOY MONTE CARLO 
    // ---------------------------------------------

    RooFitResult* fitResult;
    Double_t Mupull;


    for (int i = 0; i < n_total; i++) {

        // Generate a toy dataset
        RooDataSet* toyData = model->generate(RooArgSet(*hgg_mass), Extended(kTRUE));

        // Fit the model to the toy dataset
        fitResult = model->fitTo(*toyData,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));


        // Nspull = (Ns.getVal()-nSignal)/Ns.getError();
        // Nbpull = (Nb.getVal()-nBkg)/Nb.getError();   

        cout << "meanMass -------------- " << meanMass->getVal() << "MH -------------- " << MH->getVal() << endl;

        // Params and errors 
        Mupull = (MH->getVal() - meanMass->getVal()) / MH->getError();

        // Add Mupull value to the dataset
        Mmean.setVal(Mupull);
        MupullData->add(RooArgSet(Mmean));

            
        delete toyData;
        delete fitResult;

    }

    // Save the Mupull dataset to the workspace
    w->import(*MupullData);
    TFile outputFile("MC_Toy_Results.root", "RECREATE");

    w->Write();
    outputFile.Close();


}