#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"
using namespace RooFit;

void AnalysisMCHiggs()
{

    Int_t seed = 5;

    Double_t MminP = -3.0;
    Double_t MmaxP = 3.0; 


    // ---------------------------------------------
    // READ WORKSPACE.
    // ---------------------------------------------
    
    TFile *f = new TFile("MC_Toy_Results.root");

    // Read workspace from file
    RooWorkspace *w = (RooWorkspace *)f->Get("workspace");
    
    // Model and data from workspace
    RooRealVar *Mu = w->var("mean");
    //RooRealVar *meanMass = new RooRealVar("meanMass", "meanMass", Mu->getVal());


    cout << "meanMass -------------- " << Mu->getVal() << endl;

    // ---------------------------------------------

    RooRealVar Mmass("Mu","Mu", MminP,MmaxP); 
    RooDataSet dataMu("dataMu","dataMu",RooArgSet(Mmass));




    // Model Gaussian Mass Mean  
    RooRealVar meanMu("meanMu","meanMu",0.0,MminP,MmaxP);
    RooRealVar sigmaMu("sigmaMu","sigmaMu",1.0,0.0,5.0);
    RooGaussian SigMu("SigMu","SignalMu", Mmass, meanMu, sigmaMu); 
    

    // -------------------------------------------
    // FITTING AND CREATE CANVAS
    // -------------------------------------------

    // Fit Mass mean 
    RooFitResult* fitMu = SigMu.fitTo(dataMu, Minos(kFALSE),Save(kTRUE), NumCPU(4));
    fitMu->Print("v");

    // TCanvas* canv_Mupull = CreateCanvas("canv_Mu", fitMu, dataMu, Mu, MmaxP, MminP, SigMu, sigmaMu, meanMu);
    // canv_Mupull->Print(Form("plots/Pull_MuBc_ToyMC_%1i.png",seed));


}