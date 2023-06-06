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

void ToyClass()
{
    int NumPar = Par.getSize();
    TIterator* ParIter = Par.createIterator();

    RooMCStudy *MC = new RooMCStudy(*Model,*Obs, Binned(false), Silence(true), Extended(true), FitOptions(Save(true), PrintEvalErrors(0))); 
    //MC->generateAndFit(10000); //---- FIT DE 10 HORAS!!----
    MC->generateAndFit(1000); //---Fit "Rapido"---

    gStyle->SetOptStat(0);
    TCanvas *MC_canvas = new TCanvas("Estudio MC", "Estudio MC", x, y);
    MC_canvas->Divide(NumPar , 2);

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