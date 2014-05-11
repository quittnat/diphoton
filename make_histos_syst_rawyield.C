void make_histos_syst_rawyield(TString file_syst, TString file_default, TString out_tag){

  TFile *f1 = new TFile(file_syst.Data(),"read");
  TFile *f2 = new TFile(file_default.Data(),"read");
  TDirectoryFile *dir1 = (TDirectoryFile*)(f1->Get("effunf"));
  TDirectoryFile *dir2 = (TDirectoryFile*)(f2->Get("effunf"));
  TList *list = dir1->GetListOfKeys();

  TFile *f = new TFile(Form("plots/ratiosyst_%s.root",out_tag.Data()),"recreate");

  for (int i=0; i<list->GetSize(); i++){
    TString name = dir1->GetListOfKeys()->At(i)->GetName();
    if (!(name.Contains("hreco_"))) continue;
    TObject *obj1 = dir1->Get(name.Data());
    assert(obj1);
    TObject *obj2 = dir2->Get(name.Data());
    assert(obj2);
    TString newname = name;
    newname.Append("_ratiosyst");
    if (name.EndsWith("_0")) newname.ReplaceAll("_0_","_EBEB_");
    if (name.EndsWith("_1")) newname.ReplaceAll("_1_","_EBEE_");
    if (name.EndsWith("_2")) newname.ReplaceAll("_2_","_EEEE_");
    TH1F *h = (TH1F*)(((TH1F*)obj1)->Clone(newname.Data()));
    h->SetTitle(h->GetName());
    h->Divide((TH1F*)obj2);
    for (int j=0; j<h->GetNbinsX(); j++) h->SetBinError(j+1,0);
    for (int j=0; j<h->GetNbinsX(); j++) h->SetBinContent(j+1,1+fabs(1-h->GetBinContent(j+1)));
    f->cd();
    h->Write();
  }

}
