void copyTH1Contents(TH1F *hIN, TH1D *hOUT);
void convertTH1FtoTH1D(TString fINstr) {

  TFile *fIN = TFile::Open(fINstr,"UPDATE");
  TObject *obj;
  TKey *key;
  TIter next(fIN->GetListOfKeys());
  while ((key = (TKey *) next())) {
    obj = fIN->Get(key->GetName());
    TString obj_name     = obj->GetName();
    TString obj_name_new(obj_name(0,obj_name.Length()-4));
    TString obj_type     = obj->ClassName();

    if (obj_type == "TH1F") {
      TH1F *histIN = (TH1F*)obj->Clone(obj_name_new);
      TH1D *hist_ = new TH1D(obj_name_new,obj_name_new,histIN->GetNbinsX(),histIN->GetBinLowEdge(1),(histIN->GetBinLowEdge(histIN->GetNbinsX())+histIN->GetBinWidth(histIN->GetNbinsX())));
      copyTH1Contents(histIN,hist_);
      hist_->Write(hist_->GetName());
      histIN->Delete();
    }
  }

}

void copyTH1Contents(TH1F *hIN, TH1D *hOUT) {
  std::cout << " CopyTH1Contents" << "\n";
  for (unsigned int i0=0; i0<=hIN->GetNbinsX()+1; ++i0) { hOUT->SetBinContent(i0,hIN->GetBinContent(i0)); }
}
