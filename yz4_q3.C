void yz3_q3() {
     //TFile *f1 = TFile::Open("Q2_level_-20mv_selected.root");
     //TFile *f1 = TFile::Open("Q2_level_-20mv_all.root");
     TFile *f1 = TFile::Open("Q3_level_-20mv_all_11mv_header.root");
     TCanvas *c1 = (TCanvas*)f1->Get("c1");
     c1->Draw();
     //TFile *f2 = TFile::Open("Q2_level_-20mv_all.root");
     //TFile *f2 = TFile::Open("Q2_level_-20mv_selected.root");
     TFile *f2 = TFile::Open("Q3_level_-20mv_selected_11mv_header.root");
     TCanvas *c2 = (TCanvas*)f2->Get("c1");
     //TH1 *h2 = (TH1*)c2->GetPrimitive("Q_C2_all");
     TH1 *h2 = (TH1*)c2->GetPrimitive("Q_C3_selected");
     c1->cd();
     h2->Draw("same");
     
     //TFile *f3 = TFile::Open("Q2_level_-20mv_amplitude_cut.root");
     TFile *f3 = TFile::Open("Q3_level_-20mv_Min_C3_cut_11mv_header.root");
     TCanvas *c3 = (TCanvas*)f3->Get("c1");
     TH1 *h3 = (TH1*)c3->GetPrimitive("Q_C3_Min_C3");
     c1->cd();
     h3->Draw("same");

     
     TFile *f4 = TFile::Open("Q3_level_-20mv_T3_CL<1_selected_11mv_header.root");
     TCanvas *c4 = (TCanvas*)f4->Get("c1");
     TH1 *h4 = (TH1*)c4->GetPrimitive("htemp");
     c1->cd();
     h4->Draw("same");
}
