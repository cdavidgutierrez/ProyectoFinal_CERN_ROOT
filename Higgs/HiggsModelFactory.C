void HiggsModelFactory() {
    TFile *file = TFile::Open("tutorial.root");
    RooWorkspace *wspace = (RooWorkspace*) file->Get("workspace");
    wspace->Print("v");

    RooDataSet *hgg_data = (RooDataSet*) wspace->data("dataset");
    RooRealVar *hgg_mass = (RooRealVar*) wspace->var("CMS_hgg_mass");

    //------------------- Modelo del Background -----------------------
    wspace->factory("Exponential::expo(hgg_mass[100,80,200], alpha[-0.05,-0.2,0.01])");
    //------------------- Modelo de la señal --------------------------
    wspace->factory("Gaussian::hgg_signal(hgg_mass[100,80,200],mean[125,120,130],sigma[10, 0.001, 20])");
    //------------------- Modelo combinado ----------------------------
    wspace->factory("SUM::model(norm_s[10,100]*hgg_signal,norm_b[0,1000]*expo)");

    wspace->pdf("model")->fitTo(*hgg_data, RooFit::Extended());

    wspace->hgg_signal

    RooPlot* frame = hgg_mass->frame();
    hgg_data->plotOn(frame);
    wspace->pdf("model")->plotOn(frame);
    frame->Draw();
}