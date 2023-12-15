void yz() {
     //TFile *f1 = TFile::Open("Q2_level_-20mv_selected.root");
     TFile *f1 = TFile::Open("Q2_level_-20mv_all.root");
     TCanvas *c1 = (TCanvas*)f1->Get("c1");
     c1->Draw();
     //TFile *f2 = TFile::Open("Q2_level_-20mv_all.root");
     TFile *f2 = TFile::Open("Q2_level_-20mv_selected.root");
     TCanvas *c2 = (TCanvas*)f2->Get("c1");
     //TH1 *h2 = (TH1*)c2->GetPrimitive("Q_C2_all");
     TH1 *h2 = (TH1*)c2->GetPrimitive("Q_C2_selected");
     c1->cd();
     h2->Draw("same");
}
