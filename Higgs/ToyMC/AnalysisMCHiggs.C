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
// PLOT MC - PULL
// -------------------------------------------
TCanvas* CreateCanvas(TString cname, RooAbsPdf* model,  
RooDataSet* data,  RooRealVar M, RooRealVar* mean,  RooRealVar* sigma)  
{

    int H = 800;
    int W = 1000;

    TCanvas *c1 = new TCanvas(cname);

    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.07);
    c1->SetTopMargin(0.09);
    c1->SetBottomMargin(0.14); 

    // Plot the fitted model and data
    //RooPlot* Mframe = Mu.frame();
    //MupullData->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.5));
    RooPlot* Mframe = M.frame();
    data->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.5));
    
    // Plot the model 
    //w->pdf("modelPull")->plotOn(Mframe);
    model->plotOn(Mframe);
    Mframe->Draw();
                                                                                                                   
    Mframe->SetLabelSize(0.04,"XY");
    Mframe->SetTitleSize(0.05,"XY");
    Mframe->GetYaxis()->CenterTitle();   
    Mframe->GetXaxis()->CenterTitle();   
    Mframe->SetTitleOffset(1.0,"X");
    Mframe->SetTitleOffset(0.9,"Y");
    Mframe->SetTitleSize(0.06,"XY");
    Mframe->SetMaximum(20.0);
    Mframe->Draw();  
    gStyle->SetOptTitle(0/1);

    TLegend *legend1 = new TLegend(0.7,0.75,0.8,0.88);
    legend1->SetTextSize(0.04); //text size in pixels                                 
    legend1->SetFillColor(0);
    legend1->SetBorderSize(0);
    legend1->SetFillStyle(0); 
    // legend1->AddEntry((TObject*)nullptr, TString::Format("#mu = %.3f #pm %.3f", w->var("meanPull")->getVal(), w->var("meanPull")->getError()), ""); // Add a blank entry with the legend text
    // legend1->AddEntry((TObject*)nullptr, TString::Format("#sigma = %.3f #pm %.3f", w->var("sigmaPull")->getVal(), w->var("sigmaPull")->getError()), ""); 
    legend1->AddEntry((TObject*)nullptr, TString::Format("#mu = %.3f #pm %.3f", mean->getVal(), mean->getError()), ""); // Add a blank entry with the legend text
    legend1->AddEntry((TObject*)nullptr, TString::Format("#sigma = %.3f #pm %.3f", sigma->getVal(), sigma->getError()), ""); 
    legend1->Draw(); 
    

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

    Int_t seed = 3;

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

    // Model Gaussian Mass Mean  
    RooRealVar Mu("Mu","m_H mass Pull", 0, MminP,MmaxP); 
    w->factory("Gaussian::modelPull(Mu, meanPull[0.0,MminP,MmaxP], sigmaPull[1.0,0.0,5.0])");

    // -------------------------------------------
    // FITTING AND CREATE CANVAS
    // -------------------------------------------

    // Fit the model to the data and obtain the fit result
    RooFitResult* fitMu = w->pdf("modelPull")->fitTo(*MupullData, RooFit::Extended(), RooFit::Save(true));

    // Print the fit result
    if (fitMu) fitMu->Print();


    TCanvas* canv_Mupull = CreateCanvas("canv_Mu", w->pdf("modelPull"), MupullData, Mu, 
                                        w->var("meanPull"), w->var("sigmaPull"));
    //canv_Mupull->Print(Form("../plots/Pull_MuHiggs_ToyMC_%1i.png",seed));


}