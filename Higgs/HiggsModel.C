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

void HiggsModel() {
    TFile *file = TFile::Open("tutorial.root");
    //Dentro del archivo hay algo llamado RooWorkspace.
    file->ls();

    RooRealVar MH("MH","mass of the Hypothetical Boson (H-boson) in GeV",125,120,130);
    RooRealVar mass("m","m (GeV)",100,80,200);
    RooRealVar sigma("resolution","#sigma",10,0,20);

    RooWorkspace *wspace = (RooWorkspace*) file->Get("workspace");
    //El workspace contiene una variable y un dataset.
    wspace->Print("v");

    //Leyendo los datos y el observable para la masa del Higgs.
    RooDataSet *hgg_data = (RooDataSet*) wspace->data("dataset");
    RooRealVar *hgg_mass = (RooRealVar*) wspace->var("CMS_hgg_mass");

    RooPlot *plot = hgg_mass->frame();
    hgg_data->plotOn(plot,RooFit::Binning(160)); 

    TCanvas *hggcan = new TCanvas();
    plot->Draw();
    hggcan->Update();
    hggcan->Draw();

    //------------------- Modelo del Background -----------------------
    RooRealVar alpha("alpha","#alpha",-0.05,-0.2,0.01);
    RooExponential expo("exp","exponential function", *hgg_mass, alpha);

    //------------------- Modelo de la señal --------------------------
    sigma.setVal(1.);
    MH.setVal(125);
    sigma.setConstant();
    MH.setConstant(false);
    RooGaussian hgg_signal("signal","Gaussian PDF", *hgg_mass, MH, sigma);

    //------------------- Modelo combinado ----------------------------
    RooRealVar norm_s("norm_s","N_{s}",10,100);
    RooRealVar norm_b("norm_b","N_{b}",0,1000);
    const RooArgList components(hgg_signal,expo);
    const RooArgList coeffs(norm_s,norm_b);
    // Modelo Exp + Gauss.
    RooAddPdf model("model","f_{s+b}",components,coeffs);

    //------------------ Ajuste de los datos ---------------------------
    model.fitTo(*hgg_data,RooFit::Extended());
    

    model.plotOn(plot,RooFit::Components("exp"),RooFit::LineColor(kGreen));
    model.plotOn(plot,RooFit::LineColor(kRed));
    model.paramOn(plot);

    hggcan->Clear();
    plot->Draw();
    hggcan->Update();
    hggcan->Draw();

    //------------------- Guardando el modelo en un workspace ------------
    wspace->import(model);  
    RooArgSet *params = model.getParameters(*hgg_data);
    wspace->saveSnapshot("nominal_values",*params); //<----------- Guardar el estado actual de los parámetros.
    wspace->writeToFile("HiggsWs.root");
    wspace->Print("V");
}