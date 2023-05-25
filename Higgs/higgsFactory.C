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

void higgsFactory()
{
    TFile *file = TFile::Open("tutorial.root");
    //Dentro del archivo hay algo llamado RooWorkspace.
    file->ls();

    RooWorkspace *wspace = (RooWorkspace*) file->Get("workspace");
    //El workspace contiene una variable y un dataset.
    wspace->Print("v");

    RooDataSet *hgg_data = (RooDataSet*) wspace->data("dataset");
    RooRealVar *hgg_mass = (RooRealVar*) wspace->var("CMS_hgg_mass");

    plot = hgg_mass->frame();

    hgg_data->plotOn(plot,RooFit::Binning(160)); 

    TCanvas *hggcan = new TCanvas();
    plot->Draw();
    hggcan->Update();
    hggcan->Draw();

    //Fitting the data.
    RooRealVar alpha("alpha","#alpha",-0.05,-0.2,0.01);
    RooExponential expo("exp","exponential function",*hgg_mass,alpha);

    RooNLLVar *nll = (RooNLLVar*) expo.createNLL(*hgg_data);
    nll->Print("v");

    RooMinimizer minim(*nll);
    minim.minimize("Minuit2","migrad"); 

    expo.plotOn(plot);
    expo.paramOn(plot);
    plot->Draw();
    hggcan->Update();
    hggcan->Draw();

}