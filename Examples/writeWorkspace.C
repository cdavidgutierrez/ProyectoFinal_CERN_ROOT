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
 
void writeWorkspace()
{
 
    RooRealVar x("x", "x", 0, 10);
    
    RooRealVar mean("mean", "mean of gaussians", 5, 0, 10);
    RooRealVar sigma1("sigma1", "width of gaussians", 0.5);
    RooRealVar sigma2("sigma2", "width of gaussians", 1);
    
    RooGaussian sig1("sig1", "Signal component 1", x, mean, sigma1);
    RooGaussian sig2("sig2", "Signal component 2", x, mean, sigma2);
    

    RooRealVar a0("a0", "a0", 0.5, 0., 1.);
    RooRealVar a1("a1", "a1", 0.2, 0, 1.);
    RooChebychev bkg("bkg", "Background", x, RooArgSet(a0, a1));
    

    RooRealVar sig1frac("sig1frac", "fraction of component 1 in signal", 0.8, 0., 1.);
    RooAddPdf sig("sig", "Signal", RooArgList(sig1, sig2), sig1frac);
    

    RooRealVar bkgfrac("bkgfrac", "fraction of background", 0.5, 0., 1.);
    RooAddPdf model("model", "g1+g2+a", RooArgList(bkg, sig), bkgfrac);
    

    RooDataSet *data = model.generate(x, 1000);

    // USO DE WORKSPACES.
    
    // Crear un nuevo workspace.
    RooWorkspace *w = new RooWorkspace("w", "workspace");
    
    // Importar el modelo y todos sus componentes en el workspace.
    w->import(model);
    
    // Importar los datos en el workspace.
    w->import(*data);
    
    // Imprimir el contenido del workspace.
    w->Print();
    
    // Guardar el workspace en un archivo de root.
    w->writeToFile("myWorkspace.root");

}
