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
 
void readWorkspace()
{
    // LECTURA DE UN WORKSPACE.
    // -----------------------------------------------
    
    TFile *f = new TFile("myWorkspace.root");
    
    // Recuperar el workspace del archivo.
    RooWorkspace *w = (RooWorkspace *)f->Get("w");
    
    // Recuperar x, el modelo y los datos del workspace
    RooRealVar *x = w->var("x");
    RooAbsPdf *model = w->pdf("model");
    RooAbsData *data = w->data("modelData");
    
    //-------------------------------------------------------------
    cout<<"Contenido del MODELO........"<<endl;
    model->Print("t");
    model->fitTo(*data);
    
    RooPlot *xframe = x->frame(Title("Model and data read from workspace"));
    data->plotOn(xframe);
    model->plotOn(xframe);
    
    model->plotOn(xframe, Components("bkg"), LineStyle(kDashed));
    
    model->plotOn(xframe, Components("bkg,sig2"), LineStyle(kDotted));
    
    new TCanvas("rf503_wspaceread", "rf503_wspaceread", 600, 600);
    gPad->SetLeftMargin(0.15);
    xframe->GetYaxis()->SetTitleOffset(1.4);
    xframe->Draw();
}