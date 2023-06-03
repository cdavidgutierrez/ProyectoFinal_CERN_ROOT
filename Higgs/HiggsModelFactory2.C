void HiggsModelFactory2() {
    TFile *file = TFile::Open("tutorial.root");
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

    // Plot the fitted model and data
    RooPlot* frame = hgg_mass->frame();
    hgg_data->plotOn(frame);
    wspace->pdf("model")->plotOn(frame);
    wspace->pdf("model")->plotOn(frame, RooFit::Components("hgg_signal"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
    wspace->pdf("model")->plotOn(frame, RooFit::Components("expo"), RooFit::LineColor(kOrange), RooFit::LineStyle(kDashed));

    // TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    // legend->AddEntry(frame->findObject("model_Norm[ model ]"), "Total Fit", "l");
    // legend->AddEntry(frame->findObject("model_Norm[ hgg_signal ]"), "Signal Component", "l");
    // legend->AddEntry(frame->findObject("model_Norm[ expo ]"), "Background Component", "l");
    // legend->SetBorderSize(0);

    frame->Draw();
    //legend->Draw();
}
