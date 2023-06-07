#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include <fstream>
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooArgusBG.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
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
#include "TObject.h"

using namespace RooFit;

void ToyMC_Class()
{

    // ---------------------------------------------
    // READ WORKSPACE.
    // ---------------------------------------------
    
    TFile *f = new TFile("../HiggsWs.root");
    
    // Read workspace from file
    RooWorkspace *w = (RooWorkspace *)f->Get("workspace");
    
    // Get data and model from workspace
    RooRealVar *hgg_mass = w->var("CMS_hgg_mass");
    RooAbsPdf *model = w->pdf("model");

    // ---------------------------------------------
    // TOY MONTE CARLO 
    // ---------------------------------------------

    RooMCStudy *ToyMC = new RooMCStudy(*model, *hgg_mass, Binned(false), 
                                      Silence(true), Extended(true), FitOptions(Save(true), PrintEvalErrors(0)));            
    ToyMC->generateAndFit(1000);


    // ---------------------------------------------
    // PLOT PULL 
    // ---------------------------------------------
    TCanvas *MC_canvas = new TCanvas("Estudio MC", "Estudio MC", 600, 600);
    RooPlot *meanFrame = ToyMC->plotPull(*(w->var("MH")), Bins(40), FitGauss(true));
    meanFrame->Draw();

    MC_canvas->Print(Form("../plots/ToyMC_Class.png"));
    
}