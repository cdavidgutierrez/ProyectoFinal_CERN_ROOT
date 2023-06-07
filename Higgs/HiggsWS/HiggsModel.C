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

using namespace std;
using namespace RooFit;

void HiggsModel() {

    TFile *file = TFile::Open("../tutorial.root");
    
    //Dentro del archivo hay algo llamado RooWorkspace.
    file->ls();

    RooRealVar MH("MH","mass of the Hypothetical Boson (H-boson) in GeV",124.79,124.1,125.1);
    RooRealVar mass("m","m (GeV)",100,80,200);
    RooRealVar sigma("resolution","#sigma",10,0,20);

    RooWorkspace *wspace = (RooWorkspace*) file->Get("workspace");

    //El workspace contiene una variable y un dataset.
    wspace->Print("v");

    //Leyendo los datos y el observable para la masa del Higgs.
    RooDataSet *hgg_data = (RooDataSet*) wspace->data("dataset");
    RooRealVar *hgg_mass = (RooRealVar*) wspace->var("CMS_hgg_mass");

    //------------------- Modelo del Background -----------------------
    RooRealVar alpha("alpha","#alpha",-0.05,-0.2,0.01);
    RooExponential expo("exp","exponential function", *hgg_mass, alpha);

    //------------------- Modelo de la señal --------------------------
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

    //------------------- Guardando el modelo en un workspace ------------
    wspace->import(model);  
    RooArgSet *params = model.getParameters(*hgg_data);

    cout<<"Importando el modelo completo en un WS........"<<endl;
    // Guardar el estado actual de los parámetros.
    wspace->saveSnapshot("nominal_values",*params); 
    wspace->writeToFile("../HiggsWs.root");
    wspace->Print("V");



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
    hgg_data->plotOn(Hframe,RooFit::Binning(40), Name("data")); 

    // Model
    model.plotOn(Hframe,RooFit::Components("exp"),RooFit::LineColor(kGreen), Name("bkg"));
    model.plotOn(Hframe,RooFit::LineColor(kRed),Name("fittotal"));

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
    Hframe->SetMinimum(-2);
    //Hframe->SetMaximum(33);
    Hframe->Draw();
    gStyle->SetOptTitle(0/1);


    // Legend 1
    TLegend *legend1 = new TLegend(0.50,0.70,0.90,0.88); 
    legend1->SetTextSize(0.07);
    legend1->SetFillColor(0);
    legend1->SetBorderSize(0);
    legend1->SetFillStyle(0);
    legend1->AddEntry(Hframe->findObject("data"), "Data", "ep"); 
    legend1->AddEntry(Hframe->findObject("fittotal"), "Sig+Bkg Fit", "l");
    legend1->AddEntry(Hframe->findObject("bkg"), "Bkg.", "l");
    legend1->Draw();

    // Legend 2
    TLegend *legend2 = new TLegend(0.50,0.45,0.90,0.70); 
    legend2->SetTextSize(0.07);
    legend2->SetFillColor(0);
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->AddEntry("",Form("#alpha = %1.2f #pm %1.2f", alpha.getVal(), alpha.getError()),"");
    legend2->AddEntry("",Form("N_{s} = %1.2f #pm %1.2f", norm_s.getVal(), norm_s.getError()),"");
    legend2->AddEntry("",Form("N_{b} = %1.2f #pm %1.2f", norm_b.getVal(), norm_b.getError()),"");
    legend2->AddEntry("",Form("m_{H} = %1.2f #pm %1.2f GeV", MH.getVal(), MH.getError()),"");
    legend2->Draw();



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


    hggcan->Print(Form("../plots/HigssModel.png"));
}
