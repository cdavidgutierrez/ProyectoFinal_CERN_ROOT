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

void ToyClass()
{

    // ---------------------------------------------
    // READ WORKSPACE.
    // ---------------------------------------------
    
    TFile *f = new TFile("HiggsWs.root");
    
    // Read workspace from file
    RooWorkspace *w = (RooWorkspace *)f->Get("workspace");
    
    // Get data and model from workspace
    RooRealVar *hgg_mass = w->var("CMS_hgg_mass");
    RooAbsPdf *model = w->pdf("model");

   // --------------
    int const nbin = 1000;
    RooDataSet* DataSet = model->generate(*hgg_mass , 10000); 
    RooFitResult* FitResult = model->fitTo(*DataSet , Extended(true) , Save(true));

    // ---------


    const RooArgSet Par = RooArgSet(w->var("MH")->getVal(), w->var("norm_b")->getVal(), w->var("norm_s")->getVal());

    int NumPar = Par.getSize();
    TIterator* ParIter = Par.createIterator();

    RooMCStudy *MC = new RooMCStudy(*model,*hgg_mass, Binned(false), Silence(true), Extended(true), FitOptions(Save(true), PrintEvalErrors(0))); 
    MC->generateAndFit(20); //---Fit "Rapido"---


    //-------------------------------------------
    //-------------------------------------------


    // Dibujando Canvas
    TCanvas *c = new TCanvas("c", "c", x, y);
    c->Divide(1, 2 , 0 ,0);

    c->cd(1); gPad->SetRightMargin(0.01);

    // Primer Frame
    RooPlot* frame = hgg_mass->frame();
    DataSet->plotOn(frame);
    model->plotOn(frame);
    frame->Draw("A");
    frame->GetYaxis()->CenterTitle(); 
    frame->GetYaxis()->SetTitleSize(0.07); 
    frame->GetYaxis()->SetTitleOffset(0.5);
    frame->GetYaxis()->SetLabelSize(0.045);

    c->cd(2);gPad->SetRightMargin(0.01); gPad->SetBottomMargin(0.3);

    // Segundo Frame (Pull)
    RooHist* pullHist = frame->pullHist();
    RooPlot* pullFrame = hgg_mass->frame();
    pullFrame->addPlotable(pullHist, "P");
    pullFrame->GetYaxis()->SetNdivisions(505);
    pullFrame->GetYaxis()->SetTitle("(Data - Fit) /#sigma");
    pullFrame->GetYaxis()->CenterTitle(); 
    pullFrame->GetYaxis()->SetTitleSize(0.07); 
    pullFrame->GetYaxis()->SetTitleOffset(0.5);
    pullFrame->GetYaxis()->SetLabelSize(0.045);
    pullFrame->GetXaxis()->CenterTitle(); 
    pullFrame->GetXaxis()->SetTitleSize(0.07); 
    pullFrame->GetXaxis()->SetTitleOffset(0.9);
    pullFrame->GetXaxis()->SetLabelSize(0.045);
    pullFrame->Draw("A");
    
    TLine* zeroLine = new TLine(6.05, 0, 6.5, 0);
    zeroLine->SetLineStyle(2);
    zeroLine->Draw("same");
    c->Update();

    gStyle->SetOptStat(0);
    TCanvas *MC_canvas = new TCanvas();
    MC_canvas->Divide(NumPar , 2);


    //-------------------------------------------
    //-------------------------------------------

    // Iteraciones por cada parámetro
    auto var = ParIter->Next();
    int i = 1;
    while (var) {
      
      RooPlot *ParMeanFrame ;
      RooPlot *ParMeanPullFrame ;

      if (i>NumPar){std::cout<<"Pare en la iteracion "<<i<<endl;break;} ;

      ParMeanFrame = MC->plotParam(*(RooRealVar*)(var), Bins(nbin)); //desreferrenciar el puntero
      ParMeanPullFrame = MC->plotPull(*(RooRealVar*)(var), Bins(nbin), FitGauss(true));

      MC_canvas->cd(i); gPad->SetLeftMargin(0.1);
      ParMeanFrame->GetYaxis()->SetTitleOffset(1.4);
      ParMeanFrame->Draw();
      MC_canvas->Update();

      MC_canvas->cd(i+NumPar); gPad->SetLeftMargin(0.1);
      ParMeanPullFrame->GetYaxis()->SetTitleOffset(1.4);
      ParMeanPullFrame->Draw();
      MC_canvas->Update();

      var = ParIter->Next();
      i += 1;

    }
}