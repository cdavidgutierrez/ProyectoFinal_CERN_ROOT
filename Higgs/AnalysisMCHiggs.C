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


// -------------------------------------------
// PLOT DATA & FIT
// -------------------------------------------
TCanvas* CreateCanvas(TString cname, //RooFitResult* result, 
RooDataSet* data,  RooRealVar M, Double_t supM, Double_t infM,  
RooGaussian MassModel,  RooRealVar bwg1,  RooRealVar bwm1)  
{


    int H = 800;
    int W = 1000;

    TCanvas *c1 = new TCanvas(cname,cname,50,50,W,H);

    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.07);
    c1->SetTopMargin(0.09);
    c1->SetBottomMargin(0.14); 

    // Plot the fitted result
    RooPlot* Mframe = M.frame();
    data->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.5)); 
    MassModel.plotOn(Mframe);
    Mframe->Draw();

    Mframe->SetYTitle("Events / 0.2 ");                                                                                                                   
    Mframe->SetLabelSize(0.04,"XY");
    Mframe->SetTitleSize(0.05,"XY");
    Mframe->GetYaxis()->CenterTitle();   
    Mframe->GetXaxis()->CenterTitle();   
    Mframe->SetTitleOffset(1.0,"X");
    Mframe->SetTitleOffset(0.9,"Y");
    Mframe->SetTitleSize(0.06,"XY");
    Mframe->SetMaximum(30.0);
    Mframe->Draw();  
    gStyle->SetOptTitle(0/1);


    TLegend *legpar = new TLegend(0.5,0.7,0.8,0.88);
    legpar->SetTextSize(0.04); //text size in pixels                                 
    legpar->SetFillColor(0);
    legpar->SetBorderSize(0);
    legpar->SetFillStyle(0); 
    legpar->AddEntry("",Form("mean = %1.4f #pm %1.4f GeV ",bwm1.getVal(), bwm1.getError()),"");
    //legpar->AddEntry("",Form("#sigma_{1} = %1.4f #pm %1.4f GeV",G, Ge),"");
    legpar->Draw();

    TLegend *legMass = new TLegend(0.64,0.57,0.83,0.65);
    legMass->SetTextFont(43); 
    legMass->SetTextSize(20);  
    legMass->SetFillColor(0); 
    legMass->SetBorderSize(0);
    legMass->SetFillStyle(0); 
    legMass->Draw(); 

    TLatex *   tex1 = new TLatex(0.88,0.926,"ToyMC");
    tex1->SetNDC();
    tex1->SetTextAlign(31);
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.04); 
    tex1->SetLineWidth(2);
    
    TLatex *tex2 = new TLatex(0.15,0.926,"CMS");
    tex2->SetNDC();
    tex2->SetTextFont(61);
    tex2->SetTextSize(0.04); 
    tex2->SetLineWidth(2);
    
    TLatex *tex3 = new TLatex(0.24,0.926,"Preliminary");
    tex3->SetNDC();
    tex3->SetTextFont(52);
    tex3->SetTextSize(0.04); 
    tex3->SetLineWidth(2);
    
    tex1->Draw();  
    tex2->Draw();
    tex3->Draw();
    
    c1->Modified();
    return c1;


}


void AnalysisMCHiggs()
{

    Int_t seed = 5;

    Double_t MminP = -6.0;
    Double_t MmaxP = 6.0; 


    // ---------------------------------------------
    // READ WORKSPACE.
    // ---------------------------------------------
    
    TFile *f = new TFile("MC_Toy_Results.root");

    // Read workspace from file
    RooWorkspace *w = (RooWorkspace *)f->Get("workspace");
    

    // Access the MupullData RooDataSet from the workspace
    RooDataSet *MupullData = (RooDataSet *)w->data("MupullData");
    RooRealVar Mu("Mu","Mu", MminP,MmaxP); 


    // Model Gaussian Mass Mean  
    RooRealVar meanMu("meanMu","meanMu",0.0,MminP,MmaxP);
    RooRealVar sigmaMu("sigmaMu","sigmaMu",1.0,0.0,5.0);
    RooGaussian SigMu("SigMu","SignalMu", Mu, meanMu, sigmaMu); 
    

    // -------------------------------------------
    // FITTING AND CREATE CANVAS
    // -------------------------------------------

    // Fit Mass mean 
    RooFitResult* fitMu = SigMu.fitTo(*MupullData, Minos(kFALSE),Save(kTRUE));
    fitMu->Print("v");

    TCanvas* canv_Mupull = CreateCanvas("canv_Mu", MupullData, Mu, MminP, MmaxP, SigMu, sigmaMu, meanMu);
    canv_Mupull->Print(Form("Pull_MuHiggs_ToyMC_%1i.png",seed));


}