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
    

    Int_t seed = 4;
    Int_t n_total = 300;
    Double_t MminP = -6.0;
    Double_t MmaxP = 6.0; 

  
    // ---------------------------------------------
    // READ WORKSPACE.
    // ---------------------------------------------
    
    TFile *f = new TFile("../HiggsWs.root");
    
    // Read workspace from file
    RooWorkspace *w = (RooWorkspace *)f->Get("workspace");
    
    // Get data and model from workspace
    RooRealVar *hgg_mass = w->var("CMS_hgg_mass");
    RooAbsPdf *model = w->pdf("model");
    model->Print("t");

    RooRealVar *meanMass = new RooRealVar("meanMass", "meanMass", w->var("MH")->getVal());
    
    // -------------------------------------------
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

    // -------------------------------------------

    
    // Params to fit 
    RooRealVar Mu("Mu","Mu", 0.0, MminP, MmaxP);
    RooDataSet* MupullData = new RooDataSet("MupullData", "Mupull dataset", RooArgSet(Mu));


    // ---------------------------------------------
    // FITTING - TOY MONTE CARLO 
    // ---------------------------------------------

    RooFitResult* fitResult;
    Double_t Mupull;
    Int_t status,covQual;


    for (int i = 0; i < n_total; i++) {

        // Generate a toy dataset
        RooDataSet* toyData = model->generate(RooArgSet(*hgg_mass), Extended(kTRUE));

        // Fit the model to the toy dataset
        fitResult = model->fitTo(*toyData,Extended(),Minos(kFALSE),Save(kTRUE));


        if(i==0){
            TCanvas *c1 = new TCanvas();
            c1->SetLeftMargin(0.12);
            c1->SetRightMargin(0.07);
            c1->SetTopMargin(0.09);
            c1->SetBottomMargin(0.14); 

            RooPlot* Mframe = hgg_mass->frame();
            toyData->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.5)); 
            w->pdf("model")->plotOn(Mframe);
            w->pdf("model")->plotOn(Mframe, RooFit::Components("signal"), 
                                    RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
            w->pdf("model")->plotOn(Mframe, RooFit::Components("exp"), 
                                    RooFit::LineColor(kOrange), RooFit::LineStyle(kDashed));

            Mframe->Draw();
            Mframe->SetXTitle("m_{#gamma #gamma} [GeV]");
            gStyle->SetOptTitle(0/1);

            TLegend *legend1 = new TLegend(0.5,0.7,0.8,0.88);
            legend1->SetTextSize(0.04); //text size in pixels                                 
            legend1->SetFillColor(0);
            legend1->SetBorderSize(0);
            legend1->SetFillStyle(0); 
            legend1->AddEntry((TObject*)nullptr, TString::Format("m_H = %.3f #pm %.3f GeV", 
                                w->var("MH")->getVal(), w->var("MH")->getError()), ""); 
            legend1->AddEntry((TObject*)nullptr, TString::Format("#alpha = %.3f #pm %.3f", 
                                w->var("alpha")->getVal(), w->var("alpha")->getError()), ""); 
            legend1->AddEntry((TObject*)nullptr, TString::Format("N_b = %.3f #pm %.3f", 
                                w->var("norm_b")->getVal(), w->var("norm_b")->getError()), ""); 
            legend1->AddEntry((TObject*)nullptr, TString::Format("N_s = %.3f #pm %.3f", 
                                w->var("norm_s")->getVal(), w->var("norm_s")->getError()), ""); 

            legend1->Draw();

            c1->Print(Form("../plots/Fit_Higgs_ToyMC_%1i.png",seed));


        }

        // -----------------
        
        status = fitResult->status();
        covQual = fitResult->covQual();

        if(status!=0)continue;
        if(covQual!=3)continue;

        // -----------------

        // Pull - Mean mass
        // Higgs mass from workspace MH
        Mupull = (w->var("MH")->getVal() - meanMass->getVal()) / w->var("MH")->getError();   


        cout << i << " meanMass -------------- " << meanMass->getVal() << "MH -------------- " << w->var("MH")->getVal() << endl;
        cout << "Mupull -------------- " << Mupull << "Error" << w->var("MH")->getError() << endl;

        // Add Mupull, Nspull, Nspull value to the dataset
        Mu.setVal(Mupull);
        MupullData->add(RooArgSet(Mu));

            
        delete toyData;
        delete fitResult;

    }

    // Save the Mupull dataset to the workspace
    w->import(*MupullData);
    //w->writeToFile("MC_Toy_Results.root");
    TFile outputFile("MC_Toy_Results.root", "RECREATE");

    w->Write();
    outputFile.Close();


}