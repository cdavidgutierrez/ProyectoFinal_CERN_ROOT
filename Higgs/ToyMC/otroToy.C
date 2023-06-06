using namespace RooFit;

void otroToy() {
    TFile *file = TFile::Open("../HiggsWs.root");
    RooWorkspace *wspace = (RooWorkspace*) file->Get("workspace");
    RooAddPdf *modelHiggs = dynamic_cast<RooAddPdf*>(wspace->pdf("model"));
    RooRealVar *CMS_hgg_mass = wspace->var("CMS_hgg_mass");
    RooRealVar *mean = dynamic_cast<RooRealVar*>(modelHiggs->getVariables()->find("MH"));

    modelHiggs->Print();
    mean->Print();
    wspace->getSnapshots()->Print();

    RooMCStudy *mcStudy = new RooMCStudy(*modelHiggs, *CMS_hgg_mass, Binned(false), Silence(), Extended());

    mcStudy->generateAndFit(500); // Generate and fit 1000 toy experiments

    RooPlot *frame = mcStudy->plotPull(*mean, Bins(50), FitGauss(true));

    // Draw all plots on a canvas
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("rf801_mcstudy", "rf801_mcstudy", 600, 600);
    
    gPad->SetLeftMargin(0.15);
    frame->GetYaxis()->SetTitleOffset(1.4);
    frame->Draw();
}