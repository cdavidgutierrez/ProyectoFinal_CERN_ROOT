using namespace RooFit;

void HiggsModelFactory() {
    TFile *file = TFile::Open("../tutorial.root");
    RooWorkspace *wspace = (RooWorkspace*) file->Get("workspace");
    wspace->Print("v");


    RooDataSet *hgg_data = (RooDataSet*) wspace->data("dataset");


    // Access the variable directly from the workspace
    RooRealVar *hgg_mass = wspace->var("CMS_hgg_mass");


    // Create variables and PDFs directly in the workspace factory
    wspace->factory("Exponential::expo(CMS_hgg_mass, alpha[-0.05, -0.2, 0.01])");
    wspace->factory("Gaussian::hgg_signal(CMS_hgg_mass, mean[125, 120, 130], sigma[10, 0.001, 20])");
    wspace->factory("SUM::model(norm_s[10, 100]*hgg_signal, norm_b[0, 1000]*expo)");


    // Fit the model to the data and obtain the fit result
    RooFitResult* fitResult = wspace->pdf("model")->fitTo(*hgg_data, RooFit::Extended(), RooFit::Save(true));


    // Print the fit result
    if (fitResult)
        fitResult->Print();


    float massValue = wspace->var("mean")->getVal();
    float massError = wspace->var("mean")->getError();


    TCanvas* canvas = new TCanvas("canvas", "Higgs Fit Result", 800, 600);

    canvas->cd(1);
    canvas->SetLeftMargin(0.005);
    canvas->SetRightMargin(0.01);
    canvas->SetTopMargin(0.09);
    canvas->SetBottomMargin(0.1);

    TPad *pad1 = new TPad("pad1", "padi",0.01,0.311,0.9903769, 0.99 );
    pad1->SetLeftMargin(0.09);   
    pad1->SetRightMargin(0.019);
    pad1->SetTopMargin(0.09);
    pad1->SetBottomMargin(0.08);  

    TPad *pad2 = new TPad("pad2", "pad2",0.01,0.01,0.9903769,0.31);
    pad2->SetLeftMargin(0.09);
    pad2->SetRightMargin(0.019);  
    pad2->SetTopMargin(0.04);
    pad2->SetBottomMargin(0.3);
    //pad2->SetTickx(0);
    pad2->SetFillColor(0);
    pad2->SetGridx(0);
    pad2->SetGridy(0);

    pad1->Draw();
    pad2->Draw();
    pad1->cd();

    RooPlot* frame = hgg_mass->frame();
    frame->GetXaxis()->SetTitle("Mass (GeV)");
    hgg_data->plotOn(frame);
    wspace->pdf("model")->plotOn(frame, RooFit::Name("model"));
    wspace->pdf("model")->plotOn(frame, RooFit::Components("hgg_signal"), RooFit::LineColor(kRed),
                                   RooFit::LineStyle(kDashed), RooFit::Name("signal"));
    wspace->pdf("model")->plotOn(frame, RooFit::Components("expo"), RooFit::LineColor(kGreen),
                                   RooFit::LineStyle(kDashed), RooFit::Name("background"));

    TLegend* legend = new TLegend(0.7, 0.7, 0.95, 0.85);
    legend->AddEntry(frame->findObject("model"), "Total Fit", "l");
    legend->AddEntry(frame->findObject("signal"), "Signal Component", "l");
    legend->AddEntry(frame->findObject("background"), "Background Component", "l");
    legend->SetBorderSize(0);


    TLegend* legpar = new TLegend(0.6, 0.5, 0.95, 0.69);
    legpar->AddEntry("", Form("M(H#rightarrow#gamma#gamma) = %1.4f #pm %1.4f (GeV)",massValue, massError), "");
    legpar->AddEntry("", Form("N_{sgn} = %1.0f #pm %1.0f",wspace->var("norm_s")->getVal(), wspace->var("norm_s")->getError()), "");
    legpar->AddEntry("", Form("N_{bkg} = %1.0f #pm %1.0f",wspace->var("norm_b")->getVal(), wspace->var("norm_b")->getError()), "");
    legpar->SetBorderSize(0);                                   

    frame->Draw(); 
    legend->Draw();
    legpar->Draw();

    RooHist* pullHist = frame->pullHist(); //Calculando el pull del MODELO.

    pad2->cd();       

    RooPlot* pullFrame = hgg_mass->frame();    
    pullFrame->addPlotable(pullHist, "P");
    pullFrame->GetXaxis()->SetTitle("Mass (GeV)"); 
    pullFrame->GetXaxis()->SetTitleSize(0.08); 
    pullFrame->GetXaxis()->SetLabelSize(0.08); 
    pullFrame->GetYaxis()->SetTitle(Form("#frac{Data-Fit}{#sigma}")); 
    pullFrame->SetTitleOffset(0.3,"Y");
    pullFrame->GetYaxis()->SetTitleSize(0.08); 
    pullFrame->GetYaxis()->SetLabelSize(0.08);
    pullFrame->Draw(); 

    canvas->Print(Form("HigssFactory.png")); 

    TCanvas *MC_canvas = new TCanvas("Estudio MC", "Estudio MC", 600, 600);

    RooMCStudy *ToyMC = new RooMCStudy(*(wspace->pdf("model")), *hgg_mass, Binned(false), Silence(true), Extended(true), FitOptions(Save(true), PrintEvalErrors(0)));            
    ToyMC->generateAndFit(500);

    RooPlot *meanFrame = ToyMC->plotPull(*(wspace->var("mean")), Bins(40), FitGauss(true));
    meanFrame->Draw();

    MC_canvas->Print(Form("toyMC_fromFactory.png"));

}