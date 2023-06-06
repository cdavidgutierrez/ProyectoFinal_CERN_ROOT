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

void DataFit()
{
    TFile *file = TFile::Open("../tutorial.root");

    //Dentro del archivo hay algo llamado RooWorkspace.
    file->ls();

    RooWorkspace *wspace = (RooWorkspace*) file->Get("workspace");

    //El workspace contiene una variable y un dataset.
    wspace->Print("v");

    RooDataSet *hgg_data = (RooDataSet*) wspace->data("dataset");
    RooRealVar *hgg_mass = (RooRealVar*) wspace->var("CMS_hgg_mass");


    cout << "Number of entries: " << hgg_data->numEntries() << endl;
  

    // -------------------------------------------
    // FIT DATA
    // -------------------------------------------

    RooRealVar alpha("alpha","#alpha",-0.05,-0.2,0.01);
    RooExponential expo("exp","exponential function",*hgg_mass,alpha);

    RooNLLVar *nll = (RooNLLVar*) expo.createNLL(*hgg_data); //???
    nll->Print("v");

    RooMinimizer minim(*nll);
    minim.minimize("Minuit2","migrad"); 


    // -------------------------------------------
    // PLOT DATA
    // -------------------------------------------

    int H = 800;
    int W = 1000;
    TCanvas *hggcan = new TCanvas();
    hggcan->cd();

    hggcan->SetLeftMargin(0.005);
    hggcan->SetRightMargin(0.01);
    hggcan->SetTopMargin(0.09);
    hggcan->SetBottomMargin(0.1);

    // Two Pads = fit + pull 
    TPad *pad1 = new TPad("pad1", "pad1", 0.01, 0.41, 0.9903769, 0.99);
    pad1->SetLeftMargin(0.09);   
    pad1->SetRightMargin(0.019);
    pad1->SetTopMargin(0.09);
    pad1->SetBottomMargin(0.0);  

    TPad *pad2 = new TPad("pad2", "pad2", 0.01, 0.01, 0.9903769, 0.41);
    pad2->SetLeftMargin(0.09);
    pad2->SetRightMargin(0.019);  
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.25);
    pad2->SetFillColor(0);
    pad2->SetGridx(0);
    pad2->SetGridy(0);

    pad1->Draw();
    pad2->Draw();
    
    // --- PAD 1 ---
    pad1->cd();

    Double_t supH = 110;
    Double_t infH = 150;
    Double_t nbin = 100;//((supH - infH) / 0.010);


    RooPlot *Hframe = hgg_mass->frame();

    // Data
    hgg_data->plotOn(Hframe,RooFit::Binning(130), Name("data")); 

    // Model
    expo.plotOn(Hframe, Name("fittotal"));

    Hframe->SetYTitle("Events / (0.25)");
    Hframe->SetLabelSize(0.06, "XY");
    Hframe->SetTitleSize(0.07, "XY");
    Hframe->GetYaxis()->CenterTitle();
    Hframe->GetXaxis()->CenterTitle();
    Hframe->GetYaxis()->SetNdivisions(506, 1);
    Hframe->GetXaxis()->SetNdivisions(510, 1);
    Hframe->GetXaxis()->SetDecimals(1);
    Hframe->SetTitleOffset(0.9, "X");
    Hframe->SetTitleOffset(0.8, "Y");
    Hframe->SetTitleSize(0.06, "XY");
    Hframe->SetMinimum(-3);
    Hframe->SetMaximum(24);
    Hframe->Draw();
    gStyle->SetOptTitle(0/1);


    // Legend 
    TLegend *legend = new TLegend(0.60,0.65,0.95,0.88); 
    legend->SetTextSize(0.07);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(Hframe->findObject("data"), "Data", "ep"); 
    legend->AddEntry(Hframe->findObject("fittotal"), "Fit", "l");
    legend->AddEntry("",Form("#alpha = %1.2f #pm %1.2f", alpha.getVal(), alpha.getError()),"");
    legend->Draw();



    // -------------------------------------------
    // PLOT PULL
    // -------------------------------------------
    
    RooHist* pullHiggs = Hframe->pullHist();

    // --- PAD 2 ---
    pad2->cd();
    
    // Frame with the pull distribution 
    RooPlot* framePull = hgg_mass->frame();
    framePull->addPlotable(pullHiggs,"P");

    framePull->SetYTitle("(Data - Fit)/#sigma");
    framePull->SetXTitle("m_{#gamma #gamma} [GeV]");
    framePull->SetLabelSize(0.1,"XY");
    framePull->SetTitleSize(0.13,"X");
    framePull->SetTitleSize(0.11,"Y");
    framePull->GetYaxis()->CenterTitle();
    framePull->GetXaxis()->CenterTitle();
    framePull->GetYaxis()->SetNdivisions(505,1);
    framePull->GetXaxis()->SetNdivisions(510,1);
    framePull->GetXaxis()->SetTickLength(0.07);
    framePull->SetTitleOffset(0.9,"X");
    framePull->SetTitleOffset(0.4,"Y");
    framePull->SetMaximum(3.5);
    framePull->SetMinimum(-3.5);
    
    framePull->Draw();

    // line at 0 
    TLine *line1 = new TLine(infH,0.0,supH,0.0);
    line1->SetLineColor(1);
    line1->SetLineWidth(1);
    line1->Draw();


    hggcan->Modified();
    gPad->Update();
    gPad->RedrawAxis();


    hggcan->Print(Form("plots/HigssDataFit.png"));


}