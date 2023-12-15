if (fChain == 0) return;

#define Nmax 1100
//
//#define NpointsCh1 3204
//#define NpointsCh2 3204
//#define NpointsCh1 1282
//#define NpointsCh2 1282
#define NpointsCh2 4002
#define NpointsCh3 4002
#define NpointsCh4 4002
//#define NpointsCh3 2002
// 
#define One_ns 1e-9
//
   TRandom rand;
   int iCount(0);
//
   cout<<"Getting_entries ... Be patient ..."<<endl;
   Long64_t nentries;
   ifstream fcin("nentries");
   if (!fcin) {
     cout << " Unable to find file 'nentries', creating it ..." << endl;
     fstream fcout("nentries",ios::out);
     nentries=fChain->GetEntries();
     fcout << nentries;
     fcout.close();
   }
   else {
     fcin >> nentries;
     cout << " WARNING: Reading # of entries from file 'nentries' ..." << endl;
   }
   fcin.close();
   cout<<"N = "<<nentries<<endl;
   //
   Long64_t nbytes = 0, nb = 0;

   TMultiGraph *mC1 = new TMultiGraph();
   TMultiGraph *mC2 = new TMultiGraph();
   TMultiGraph *mC3 = new TMultiGraph();
   TMultiGraph *mC4 = new TMultiGraph();

   TGraph **gC1 = new TGraph*[nentries];
   TGraph **gC2 = new TGraph*[nentries];
   TGraph **gC3 = new TGraph*[nentries];
   TGraph **gC4 = new TGraph*[nentries];


///////////////////////////////////////////////////////////
// for (Long64_t jentry=0; jentry<500;jentry++) { // <<< for debuging
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     //
     if(10*jentry%nentries == 0) cout<<"Done "<<100.*jentry/nentries<<" %"<<endl;
     //
     Long64_t ientry = LoadTree(jentry);   
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;  // after this line DO an analysis     
     // 
     if(rand.Rndm() < 0.05 && iCount<Nmax){ // .05

//     gC1[iCount] = new TGraph(NpointsCh1,C1_time,C1_vol);
       gC2[iCount] = new TGraph(NpointsCh2,C2_time,C2_vol);
       gC3[iCount] = new TGraph(NpointsCh3,C3_time,C3_vol);
       gC4[iCount] = new TGraph(NpointsCh4,C4_time,C4_ampl);

       mC2->Add(gC2[iCount]);
       mC3->Add(gC3[iCount]);
       mC4->Add(gC4[iCount]);
//       mC1->Add(gC1[iCount]);
 
       iCount++;

      }
     if(iCount==Nmax) break;
   } // here is the end of loop over events
/////////////////////////////////////////////////////////////

//
TFile *file = new TFile("outSimple.root","recreate");
TCanvas *cC2 = new TCanvas("cC2","",1200,800);
TCanvas *cC3 = new TCanvas("cC3","",1200,800);
TCanvas *cC1 = new TCanvas("cC1","",1200,800);
TCanvas *cC4 = new TCanvas("cC4","",1200,800);
//
file->mkdir("Signals");
file->cd("Signals");
cC2->cd();
mC2->Draw("apl");
cC2->Write();
cC3->cd();
mC3->Draw("apl");
cC3->Write();
cC1->cd();
//mC1->Draw("al");
cC1->Write();
cC4->cd();
mC4->Draw("apl");
cC4->Write();


file->Close();
//cout<<"just after of flcose in the file of simple.C"<<endl;
   delete cC2; delete cC3;
   delete cC1; delete cC4;
//   delete cC2t; delete cC3t;
//   delete cC2tCFD;
//   delete cC3tCFD;
//   delete cC3tFit; 
//   delete cDtimeCFD; delete cDtimeCFDboth;
//
   delete mC2; delete mC3;
   delete mC1; delete mC4;

   delete [] gC2; delete [] gC3;
   delete [] gC1; delete [] gC4;



