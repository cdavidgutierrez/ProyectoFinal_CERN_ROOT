#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH2.h"
#include "RooFitResult.h"
#include "TStyle.h"
#include "TDirectory.h"
 
using namespace RooFit;
 
void myToy()
{
   // C r e a t e   m o d e l
   // -----------------------

    TFile *file = TFile::Open("../tutorial.root");
    RooWorkspace *wspace = (RooWorkspace*) file->Get("workspace");
    RooRealVar *hgg_mass = wspace->var("CMS_hgg_mass");
   // Create variables and PDFs directly in the workspace factory
    wspace->factory("Exponential::expo(CMS_hgg_mass, alpha[-0.05, -0.2, 0.01])");
    wspace->factory("Gaussian::hgg_signal(CMS_hgg_mass, mean[125, 120, 130], sigma[10, 0.001, 20])");
    wspace->factory("SUM::modelHiggs(norm_s[10, 100]*hgg_signal, norm_b[0, 1000]*expo)");
 
   // Declare observable x
   RooRealVar x("x", "x", 0, 10);
   x.setBins(40);
 
   // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
   RooRealVar mean("mean", "mean of gaussians", 5, 0, 10);
   RooRealVar sigma1("sigma1", "width of gaussians", 0.5);
   RooRealVar sigma2("sigma2", "width of gaussians", 1);
 
   RooGaussian sig1("sig1", "Signal component 1", x, mean, sigma1);
   RooGaussian sig2("sig2", "Signal component 2", x, mean, sigma2);
 
   // Build Chebychev polynomial pdf
   RooRealVar a0("a0", "a0", 0.5, 0., 1.);
   RooRealVar a1("a1", "a1", -0.2, -1, 1.);
   RooChebychev bkg("bkg", "Background", x, RooArgSet(a0, a1));
 
   // Sum the signal components into a composite signal pdf
   RooRealVar sig1frac("sig1frac", "fraction of component 1 in signal", 0.8, 0., 1.);
   RooAddPdf sig("sig", "Signal", RooArgList(sig1, sig2), sig1frac);
 
   // Sum the composite signal and background
   RooRealVar nbkg("nbkg", "number of background events,", 150, 0, 1000);
   RooRealVar nsig("nsig", "number of signal events", 150, 0, 1000);
   RooAddPdf model("model", "g1+g2+a", RooArgList(bkg, sig), RooArgList(nbkg, nsig));

   RooAbsPdf *pdf = wspace->pdf("modelHiggs");
   pdf->Print();

   const RooArgSet *observables = pdf->getVariables();
observables->Print();
 
   // C r e a t e   m a n a g e r
   // ---------------------------
 
   // Instantiate RooMCStudy manager on model with x as observable and given choice of fit options
   //
   // The Silence() option kills all messages below the PROGRESS level, leaving only a single message
   // per sample executed, and any error message that occur during fitting
   //
   // The Extended() option has two effects:
   //    1) The extended ML term is included in the likelihood and
   //    2) A poisson fluctuation is introduced on the number of generated events
   //
   // The FitOptions() given here are passed to the fitting stage of each toy experiment.
   // If Save() is specified, the fit result of each experiment is saved by the manager
   //
   // A Binned() option is added in this example to bin the data between generation and fitting
   // to speed up the study at the expense of some precision
 
   RooMCStudy *mcstudy =
      new RooMCStudy(model, x, Binned(true), Silence(), Extended(), FitOptions(Save(true), PrintEvalErrors(0)));
 
   // G e n e r a t e   a n d   f i t   e v e n t s
   // ---------------------------------------------
 
   // Generate and fit 1000 samples of Poisson(nExpected) events
   mcstudy->generateAndFit(1000);
 
   // E x p l o r e   r e s u l t s   o f   s t u d y
   // ------------------------------------------------
 
   // Make plots of the distributions of mean, the error on mean and the pull of mean
   RooPlot *frame1 = mcstudy->plotParam(mean, Bins(40));
   RooPlot *frame2 = mcstudy->plotError(mean, Bins(40));
   RooPlot *frame3 = mcstudy->plotPull(mean, Bins(40), FitGauss(true));
 
   // Plot distribution of minimized likelihood
   RooPlot *frame4 = mcstudy->plotNLL(Bins(40));
 
   // Make some histograms from the parameter dataset
   TH1 *hh_cor_a0_s1f = mcstudy->fitParDataSet().createHistogram("hh", a1, YVar(sig1frac));
   TH1 *hh_cor_a0_a1 = mcstudy->fitParDataSet().createHistogram("hh", a0, YVar(a1));
 
   // Access some of the saved fit results from individual toys
   TH2 *corrHist000 = mcstudy->fitResult(0)->correlationHist("c000");
   TH2 *corrHist127 = mcstudy->fitResult(127)->correlationHist("c127");
   TH2 *corrHist953 = mcstudy->fitResult(953)->correlationHist("c953");
 
   // Draw all plots on a canvas
   gStyle->SetOptStat(0);
   TCanvas *c = new TCanvas("rf801_mcstudy", "rf801_mcstudy", 900, 900);
   
   gPad->SetLeftMargin(0.15);
   frame3->GetYaxis()->SetTitleOffset(1.4);
   frame3->Draw();
 
   // Make RooMCStudy object available on command line after
   // macro finishes
   gDirectory->Add(mcstudy);
}