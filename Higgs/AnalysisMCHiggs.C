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

    Double_t MminP = -20.0;
    Double_t MmaxP = 20.0; 


    // ---------------------------------------------
    // READ WORKSPACE.
    // ---------------------------------------------
    
    TFile *f = new TFile("MC_Toy_Results.root");

    // Read workspace from file
    RooWorkspace *w = (RooWorkspace *)f->Get("workspace");

    // Access the MupullData RooDataSet from the workspace
    RooDataSet *MupullData = (RooDataSet *)w->data("MupullData");



    RooRealVar Mmass("Mu","Mu", MminP,MmaxP); 
    //RooDataSet dataMu("dataMu","dataMu",RooArgSet(Mmass));


    // Model Gaussian Mass Mean  
    RooRealVar meanMu("meanMu","meanMu",0.0,MminP,MmaxP);
    RooRealVar sigmaMu("sigmaMu","sigmaMu",1.0,0.0,5.0);
    RooGaussian SigMu("SigMu","SignalMu", Mmass, meanMu, sigmaMu); 
    

    // -------------------------------------------
    // FITTING AND CREATE CANVAS
    // -------------------------------------------

    // Fit Mass mean 
    RooFitResult* fitMu = SigMu.fitTo(*MupullData, Minos(kFALSE),Save(kTRUE), NumCPU(4));
    fitMu->Print("v");

    // Create a canvas for plotting
    TCanvas* c = new TCanvas("c", "Mupull Fit", 800, 600);

    // Plot the fitted result
    RooPlot* Mframe = Mmass.frame();
    MupullData->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.5)); 
    SigMu.plotOn(Mframe);
    Mframe->Draw();
                                                                                                                      
    Mframe->SetLabelSize(0.04,"XY");
    Mframe->SetTitleSize(0.05,"XY");
    Mframe->GetYaxis()->CenterTitle();   
    Mframe->GetXaxis()->CenterTitle();   
    Mframe->SetTitleOffset(1.0,"X");
    Mframe->SetTitleOffset(0.9,"Y");
    Mframe->SetTitleSize(0.06,"XY");
    Mframe->SetMaximum(7.0);
    Mframe->Draw();  
    gStyle->SetOptTitle(0/1);


}