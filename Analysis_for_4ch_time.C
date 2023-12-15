//                             //
// time is being saved in [ps] //
//                             //
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TRandom.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include <iostream>
#include <vector>
using namespace std;
#define pi 3.141592653589793



#define Unit 1.e+12
#define One_ns 1e-9
//
if (fChain == 0) return;
#define Nmax 50
#define Nparam 2
//
// specify channels containing data from CFD ! (value = 0 not from CFD, =1 from CFD)
//
int NpointsCh1;
int NpointsCh2;
int NpointsCh3;
int NpointsCh4;
//
// time window
/*
#define C2min 35.0e-9  // ok
#define C2max 36.4e-9   // ok
#define C3min 34.5e-9
#define C3max 35.0e-9
*/
//
//double C1Level(-0.1);//modify from 0.012 to 0.12
//double C2Level(-0.008);//change from 0.08. if both set to -0.024, 
double C2Level(-0.008);//change from 0.008.  
//double C3Level(-0.020); // "Q_C2_selected" and "Q_C2" could not coincidence well. 
double C3Level(-0.020); // "Q_C2_selected" and "Q_C2" could not coincidence well. 
//double C4Level(-0.022);
//
int ok(0),null(0);
int nevents(0);
PulseExtrema gExtrema_C1,gExtrema_C2,gExtrema_C3,gExtrema_C4;
//
TF1 **fun_C1,**fun_C2,**fun_C3,**fun_C4;
fun_C1 = new  TF1*() ;
fun_C2 = new  TF1*();
fun_C3 = new  TF1*();
fun_C4 = new  TF1*();
*fun_C1 =NULL;*fun_C2 =NULL;*fun_C3 =NULL;*fun_C4 =NULL;
//
TFile *fTime = new TFile("outTime.root","recreate");
//
   double time_C1,time_C2,time_C3,time_C4;
   double tSimple_C1,tSimple_C2,tSimple_C3,tSimple_C4;
   double timeCFD_C1,timeCFD_C2,timeCFD_C3,timeCFD_C4;
   double Min_C1,Min_C2,Min_C3,Min_C4;
   double Max_C1,Max_C2,Max_C3,Max_C4;
   double Area_C1,Area_C2,Area_C3,Area_C4;
   int    Npeak_C1,Npeak_C2,Npeak_C3,Npeak_C4;
   double Min2_C1,Min2_C2,Min2_C3,Min2_C4;
//
   TTree *tdT = new TTree ("tdT","");
   // CL branches fit
   tdT->Branch("T1_fit_CL",&time_C1,"time_C1/D");
   tdT->Branch("T2_fit_CL",&time_C2,"time_C2/D");
   tdT->Branch("T3_fit_CL",&time_C3,"time_C3/D");
   tdT->Branch("T4_fit_CL",&time_C4,"time_C4/D");
   // CL branches NO fit
   tdT->Branch("T1_CL",&tSimple_C1,"tSimple_C1/D");
   tdT->Branch("T2_CL",&tSimple_C2,"tSimple_C2/D");
   tdT->Branch("T3_CL",&tSimple_C3,"tSimple_C3/D");
   tdT->Branch("T4_CL",&tSimple_C4,"tSimple_C4/D");
// CF branches
   tdT->Branch("T1_fit_CFD",&timeCFD_C1,"timeCFD_C1/D");
   tdT->Branch("T2_fit_CFD",&timeCFD_C2,"timeCFD_C2/D");
   tdT->Branch("T3_fit_CFD",&timeCFD_C3,"timeCFD_C3/D");
   tdT->Branch("T4_fit_CFD",&timeCFD_C4,"timeCFD_C4/D");
//   
   tdT->Branch("Min_C1",&Min_C1,"Min_C1/D");
   tdT->Branch("Min_C2",&Min_C2,"Min_C2/D");
   tdT->Branch("Min_C3",&Min_C3,"Min_C3/D");
   tdT->Branch("Min_C4",&Min_C4,"Min_C4/D");
//
   tdT->Branch("Max_C1",&Max_C1,"Max_C1/D");
   tdT->Branch("Max_C2",&Max_C2,"Max_C2/D");
   tdT->Branch("Max_C3",&Max_C3,"Max_C3/D");
   tdT->Branch("Max_C4",&Max_C4,"Max_C4/D");
   //
   tdT->Branch("Min2_C1",&Min2_C1,"Min2_C1/D");
   tdT->Branch("Max2_C2",&Min2_C2,"Min2_C2/D");
   tdT->Branch("Min2_C3",&Min2_C3,"Min2_C3/D");
   tdT->Branch("Min2_C4",&Min2_C4,"Min2_C4/D");
   //
   Double_t RiseTime_C1,RiseTime_C2,RiseTime_C3,RiseTime_C4;
   tdT->Branch("RiseTime_C1",&RiseTime_C1,"RiseTime_C1/D");
   tdT->Branch("RiseTime_C2",&RiseTime_C2,"RiseTime_C2/D");
   tdT->Branch("RiseTime_C3",&RiseTime_C3,"RiseTime_C3/D");
   tdT->Branch("RiseTime_C4",&RiseTime_C4,"RiseTime_C4/D");
   //
   tdT->Branch("Area_C1",&Area_C1,"Area_C1/D");
   tdT->Branch("Q_C2",&Area_C2,"Area_C2/D");
   tdT->Branch("Q_C3",&Area_C3,"Area_C3/D");
   tdT->Branch("Area_C4",&Area_C4,"Area_C4/D");
   // 
   tdT->Branch("Npeak_C1",&Npeak_C1,"Npeak_C1/I");
   tdT->Branch("Npeak_C2",&Npeak_C2,"Npeak_C2/I");
   tdT->Branch("Npeak_C3",&Npeak_C3,"Npeak_C3/I");
   tdT->Branch("Npeak_C4",&Npeak_C4,"Npeak_C4/I");
   // 
   Int_t event;
   tdT->Branch("Event",&event,"event/I");
   //
   double Chi2NDF_C1,Chi2NDF_C2,Chi2NDF_C3,Chi2NDF_C4;
   tdT->Branch("C1_NChi2",&Chi2NDF_C1,"Chi2NDF_C1/D");
   tdT->Branch("C2_NChi2",&Chi2NDF_C2,"Chi2NDF_C2/D");
   tdT->Branch("C3_NChi2",&Chi2NDF_C3,"Chi2NDF_C3/D");
   tdT->Branch("C4_NChi2",&Chi2NDF_C4,"Chi2NDF_C4/D");
   //
   double a_C1,a_C2,a_C3,a_C4;
   tdT->Branch("a_C1",&a_C1,"a_C1/D");
   tdT->Branch("a_C2",&a_C2,"a_C2/D");
   tdT->Branch("a_C3",&a_C3,"a_C3/D");
   tdT->Branch("a_C4",&a_C4,"a_C4/D");
   
   double RSCFD_C2,RSCFD_C3;//Raw Signal CFD
   tdT->Branch("RSCFD_C2",&RSCFD_C2,"RSCFD_C2/D");
   tdT->Branch("RSCFD_C3",&RSCFD_C3,"RSCFD_C3/D");


   TRandom rand;
   int iCount(0);
   //
   cout<<"Getting_entries ... Be patient ..."<<endl;
   Long64_t nentries,NNentries;

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
     cout << " Reading # of entries from file 'nentries' ..." << endl;
   }
   fcin.close();
   cout<<"N = "<<nentries<<endl;
   NNentries = (nentries%10==0)?nentries:nentries+10-(nentries%10);
cout<<"NNentries="<<NNentries<<endl;

   //
   Long64_t nbytes = 0, nb = 0;

//   TMultiGraph *mC1 = new TMultiGraph();
   TMultiGraph *mC2 = new TMultiGraph();
   TMultiGraph *mC3 = new TMultiGraph();
   TMultiGraph *mC4 = new TMultiGraph();

//   TGraph **gC1 = new TGraph*[Nmax];
   TGraph **gC2 = new TGraph*[Nmax];
   TGraph **gC3 = new TGraph*[Nmax];
   TGraph **gC4 = new TGraph*[Nmax];

//   TGraph **gfunC1 = new TGraph*[Nmax];
   TGraph **gfunC2 = new TGraph*[Nmax];
   TGraph **gfunC3 = new TGraph*[Nmax];
   TGraph **gfunC4 = new TGraph*[Nmax];
   
// start time
time(&start);
current = start;
double DelTime,sumDelT(0.),sumDelT2(0.);
double AvgTime,SigTime;
//
////////////////////////////////////////////////////////////////////////
//for (Long64_t jentry=300; jentry<1600;jentry++) { // <<< for debuging
//Int_t Nchosen(356); for (Long64_t jentry=Nchosen; jentry==Nchosen;jentry++) {
for (Long64_t jentry=0; jentry<nentries;jentry++) {

     if(10*jentry%NNentries == 0) cout<<"Done "<<100*jentry/NNentries<<" %"<<endl;
     //
     Long64_t ientry = LoadTree(jentry);   
     if (ientry < 0) {cout << "ientry = 0"; break;}
     nb = fChain->GetEntry(jentry);   nbytes += nb;  // after this line DO an analysis
     event = jentry;
     //
    int Npoints(4002);
     NpointsCh1 = Npoints;
     NpointsCh2 = Npoints;
     NpointsCh3 = Npoints;
     NpointsCh4 = Npoints;
     //

bool PositivePulse1=0, PositivePulse2=0,PositivePulse3=0,PositivePulse4=0;
//---
// C1 Analysis
//---
//     Double_t TimeShift(140.e-12); // time shift in [s] for CFD algorithm
//     CFDLEVEL = 3.;
//     if(! *fun_C1) *fun_C1 = new TF1();
//     gExtrema_C1 =Measures(NpointsCh1,C1_time,C1_vol, C1min,C1max,C1Level,TimeShift, fun_C1,
//                           PositivePulse1,
//                           /* results */ Min_C1,Max_C1,Min2_C1, tSimple_C1,time_C1,timeCFD_C1,RiseTime_C1, Area_C1);
//     if(*fun_C1){
//       Chi2NDF_C1 = (*fun_C1)->GetChisquare() / (*fun_C1)->GetNDF();
//       a_C1 = (*fun_C1)->GetParameter(1);
//     }
//     Npeak_C1 = gExtrema_C1.Number;




 //  for (Int_t i=0; i< NpointsCh2;++i){
 //  cout << "C2_time["<<i<< "] = " << C2_time[i]<<endl;
 //  cout << "C2_vol["<<i<< "] = " << C2_vol[i]<<endl;
 //  }

	      
//---
// C2 Analysis
//---

// ############# Find the max vol point. #############//

 double min_ampl_2 = 10.,min_time_2 = 0.,  max_ampl_2 = -10.;
 int m2 = 0;
// double tmin = -0.14e-6, tmax = -0.02e-6;
 for(int k=0;k<NpointsCh2;k++){
    if(C2_time[k] >= C2min && C2_time[k] <= C2max){
       if(C2_vol[k] < min_ampl_2  )
        {
          min_ampl_2 = C2_vol[k];
          min_time_2 = C2_time[k];
          m2 = k;
        }
      }
 }

//  ########  End of find the max point  #########3  //


  // ############# Do interpolate  *********** //

//  int n_side = 10, l_side = 500, r_side = 500, m = i;
  int n_side_2 = 22; // original is 30; 
  const int J2 = 50;// number of points want to interpolate.
  int NIT_C2 = 2*n_side_2*J2;// 2000 = (20+20)*50
  vector <double> tmp2_time_interpolated,tmp2_vol_interpolated;
  double T_itv_2 = 2.5e-11;//T_itv_2 = time_interval for TB data.
    //for(int i=0;i<Nmax;++i){
    for(int i=m2-n_side_2;i < m2+n_side_2;++i){//n_side_2 points in each side of max point.
        for(int j=0;j<J2;++j){
            double sum = 0.; 
            tmp2_time_interpolated.push_back(C2_time[i]+T_itv_2*j/J2);
            if(j == 0) sum = C2_vol[i];
            else
            {
                //for(int k=0;k<Nmax;++k){
                for(int k=m2-n_side_2;k<m2+n_side_2;++k){
                double x = C2_time[i] + T_itv_2*j/J2 - C2_time[k];
                //cout << "x= " << x<< endl;
                sum += C2_vol[k]*sin(pi*x/T_itv_2)/(pi*x/T_itv_2);
                }
            }
            tmp2_vol_interpolated.push_back(sum);
       }
    }
   // for (int i = 0; i < vol_interpolated.size(); ++ i)
   //     cout << "vol_interpolated[" << i <<"]= " << vol_interpolated[i] << endl;

// #########  End of do interpolate ########### //


// #########  convert from vector to array ########### //

    int new_size_2 = tmp2_vol_interpolated.size();
    //cout << "size of vol_interpolated[] is, " << new_size << endl;
    Double_t C2_time_interpolated[new_size_2], C2_vol_interpolated[new_size_2];
    for (Int_t i=0; i<new_size_2; ++i) {
    C2_time_interpolated[i] = tmp2_time_interpolated[i];
    C2_vol_interpolated[i] = tmp2_vol_interpolated[i];
    }

// ##############  End of convertion  ############### //

//   for (Int_t i=0; i< new_size_2;++i){
//   cout << "C2_time_interpolated["<<i<< "]" << C2_time_interpolated[i]<<endl;
//   cout << "C2_vol_interpolated["<<i<< "]" << C2_vol_interpolated[i]<<endl;
//   }


     Double_t TimeShift_C2(100.e-12); // time shift in [s] for CFD algorithm, original is 140.
     CFDLEVEL = 3.;
     double max_ampl_cut_C2 = -0.001;//since cC2's signal is quite smaller than cC3.
     double half_width_C2 = 2.5e-10;
     double base_line_shift_C2 = -0.007;
     if(! *fun_C2) *fun_C2 = new TF1();
     gExtrema_C2 =Measures(new_size_2,C2_time_interpolated,C2_vol_interpolated, C2min,C2max,C2Level,TimeShift_C2, fun_C2,PositivePulse2,Min_C2,Max_C2,Min2_C2, tSimple_C2,time_C2,timeCFD_C2,RiseTime_C2, Area_C2, half_width_C2,base_line_shift_C2,RSCFD_C2,max_ampl_cut_C2,NIT_C2);
     //cout << "RSCFD_C2 = "<< RSCFD_C2 << endl;
     if(*fun_C2){
       Chi2NDF_C2 = (*fun_C2)->GetChisquare() / (*fun_C2)->GetNDF();
       a_C2 = (*fun_C2)->GetParameter(1);
     }
     Npeak_C2 = gExtrema_C2.Number;
    Area_C2 = (1.56e+11)/3.305*6.5*Area_C2; //1.56e+11 and 3.305, in the case of PMTB@3.1K; 
    //According to the log, there exist an attenuator, but data shows no, so delete "/1.4".
    //Area_C2 = (2.36e+11)/3.5*Area_C2; //2.36e+11 and 3.5(3.5 by guess), in the case of PMTB@3.05K.
          
     //Area_C2 = ((3.0e+11)/3.9)*Area_C2;//in the case of PMTB@2.9K.
     //Area_C2 = ((1.56e+11)/3.305)*Area_C2;//in the case of PMTB@3.1K.
     //cout << "tSimple_C2=" << tSimple_C2 <<endl;


//---
// C3 Analysis
//---       
     //

// ############# Find the max vol point. #############//

 double min_ampl_3 = 10.,min_time_3 = 0.,  max_ampl_3 = -10.;
 int m3 = 0;
// double tmin = -0.14e-6, tmax = -0.02e-6;
 for(int k=0;k<NpointsCh3;k++){
    if(C3_time[k] >= C3min && C3_time[k] <= C3max){
       if(C3_vol[k] < min_ampl_3  )
        {
          min_ampl_3 = C3_vol[k];
          min_time_3 = C3_time[k];
          m3 = k;
        }
      }
 }

//  ########  End of find the max point  #########3  //


  // ############# Do interpolate  *********** //

//  int n_side = 10, l_side = 500, r_side = 500, m = i;
  int n_side_3 = 21; //original n_side_2 = 21, modified @ Aug 2rd, 2011.  n_side_2 = 20; 
  const int J3 = J2;// number of points want to interpolate.
  int NIT_C3 = 2*n_side_3*J3;// originally, NIT_C3 = 2800 
  vector <double> tmp3_time_interpolated,tmp3_vol_interpolated;
  double T_itv_3 = T_itv_2;//T_itv_3 = time_interval for TB data.
    //for(int i=0;i<Nmax;++i){
    for(int i=m3-n_side_3;i < m3+n_side_3;++i){//n_side points in each side of max point.
        for(int j=0;j<J3;++j){
            double sum = 0.; 
            tmp3_time_interpolated.push_back(C3_time[i]+T_itv_3*j/J3);
            if(j == 0) sum = C3_vol[i];
            else
            {
                //for(int k=0;k<Nmax;++k){
                for(int k=m3-n_side_3;k<m3+n_side_3;++k){
                double x = C3_time[i] + T_itv_3*j/J3 - C3_time[k];
                //cout << "x= " << x<< endl;
                sum += C3_vol[k]*sin(pi*x/T_itv_3)/(pi*x/T_itv_3);
                }
            }
            tmp3_vol_interpolated.push_back(sum);
       }
    }
   // for (int i = 0; i < vol_interpolated.size(); ++ i)
   //     cout << "vol_interpolated[" << i <<"]= " << vol_interpolated[i] << endl;

// #########  End of do interpolate ########### //


// #########  convert from vector to array ########### //

    int new_size_3 = tmp3_vol_interpolated.size();
  //  cout << "size of vol_interpolated[] is, " << new_size << endl;
    Double_t C3_time_interpolated[new_size_3], C3_vol_interpolated[new_size_3];
    for (Int_t i=0; i<new_size_3; ++i) {
    C3_time_interpolated[i] = tmp3_time_interpolated[i];
    C3_vol_interpolated[i] = tmp3_vol_interpolated[i];
    }

// ##############  End of convertion  ############### //


   cout << "########## begin of one set  #############" << endl <<endl;
   for (Int_t i=0; i< new_size_3;++i){
   cout << C3_time_interpolated[i]<< "	" << C3_vol_interpolated[i]<<endl;
   //cout << "C3_vol_interpolated["<<i<< "]" << C3_vol_interpolated[i]<<endl;
   }
   cout << "########## end of one set  #############" << endl << endl;
     Double_t TimeShift_C3(100.e-12); // time shift in [s] for CFD algorithm, original is 140.
     double max_ampl_cut_C3 = -0.002;
     double half_width_C3 = 2.8e-10;
     double base_line_shift_C3 = -0.006;
     if(! *fun_C3) *fun_C3 = new TF1();
     gExtrema_C3 =Measures(new_size_3,C3_time_interpolated,C3_vol_interpolated, C3min,C3max,C3Level,TimeShift_C3, fun_C3, PositivePulse3, Min_C3,Max_C3,Min2_C3, tSimple_C3,time_C3,timeCFD_C3,RiseTime_C3, Area_C3,half_width_C3,base_line_shift_C3,RSCFD_C3,max_ampl_cut_C3,NIT_C3);
     //cout << "RSCFD_C3 = "<< RSCFD_C3 << endl;
     if(*fun_C3){
       Chi2NDF_C3 = (*fun_C3)->GetChisquare() / (*fun_C3)->GetNDF();
       a_C3 = (*fun_C3)->GetParameter(1);
     }
     Npeak_C3 = gExtrema_C3.Number ;
     //Area_C3 = (9.62e+10)/1.2*Area_C3; //This is former laser test, photek@5.0kv. 
     //Area_C3 = (2.08e+11)/1.3*Area_C3; // photek_0605@4.7kv, 1/1.3 was the expected value of CE. 
       
     //Area_C1 = (1.02e+12)/4.112*Area_C1; //1/(current_gain*(1.6e-19)*50), for PMTA@2.8K.
     Area_C3 = ((1.25e+11)/2.62)*Area_C3*3; //1/(current_gain*(1.6e-19)*50), for PMTA@3.25K;3 is amp+att.
    // cout << "Area_C3=" << Area_C3 <<endl;
     //cout << "tSimple_C3=" << tSimple_C3 <<endl;

//---
// C4 Analysis
//---       
     //   cerr<<"\nC4"<<endl;
/*     
     double half_width_C4 = 0;
     double base_line_shift_C4 = 0;
     if(! *fun_C4) *fun_C4 = new TF1();
     gExtrema_C4 =Measures(NpointsCh3,C4_time,C4_ampl, C4min,C4max,C4Level,TimeShift, fun_C4,
			   PositivePulse4,
      // results :
			   Min_C4,Max_C4,Min2_C4, tSimple_C4,time_C4,timeCFD_C4,RiseTime_C4, Area_C4,half_width_C4,base_line_shift_C4,"cfd");
     if(*fun_C4){
       Chi2NDF_C4 = (*fun_C4)->GetChisquare() / (*fun_C4)->GetNDF();
       a_C4 = (*fun_C4)->GetParameter(1);
     }
     Npeak_C4 = gExtrema_C4.Number;
*/
     //
     tdT->Fill();
     //
     if( timeCFD_C1 < 1. && timeCFD_C2 < 1.) ok++;
     //
     // if(jentry%100 == 1) 
/*
     cout<<jentry<<"  "<<time_C2<<", "<< timeCFD_C3
	 <<" funC2="<<*fun_C2<<" funC1="<<*fun_C1<<endl;
*/
     if( iCount<Nmax  &&  rand.Rndm() < 0.1
	 // && time_C2  <1. && timeCFD_C3 < 1.
	 //	 !(timeCFD_C3  <1. && timeCFD_C4 < 1.)
	 //&& Npeak_C2<3 && Npeak_C3<3&& 
	 //(Npeak_C2!=2|| Min2_C2>-0.4) &&(Npeak_C3!=2 || Min2_C3>-0.4)
	 //&& ( (timeCFD_C2-timeCFD_C3)>-1.e-9 ||  (timeCFD_C2-timeCFD_C3)<-1.4e-9 )
	 ){
       //gC1[iCount] = new TGraph(NpointsCh1,C1_time,C1_ampl);
       /*cout<<"NpointsCh2="<<NpointsCh2<<endl;
       for(int i=0;i<NpointsCh2;i++){
	 cout<<i<<"  "<<C1_time[i]<<",  "<<C2_vol[i]<<endl;
       }
       exit(1);*/
      // gC1[iCount] = new TGraph(NpointsCh1,C1_time,C1_vol);
      // gC2[iCount] = new TGraph(NpointsCh2,C1_time,C2_vol);
       gC2[iCount] = new TGraph(new_size_2,C2_time_interpolated,C2_vol_interpolated);
       gC3[iCount] = new TGraph(new_size_3,C3_time_interpolated,C3_vol_interpolated);
       gC4[iCount] = new TGraph(NpointsCh4,C4_time,C4_ampl);

       //gC1[iCount]->SetLineColor(kBlue);
       gC2[iCount]->SetLineColor(kBlue);
       gC3[iCount]->SetLineColor(kBlue);
       gC4[iCount]->SetLineColor(kBlue);

       //
      // mC1->Add(gC1[iCount]);
       mC2->Add(gC2[iCount]);
       mC3->Add(gC3[iCount]);
       mC4->Add(gC4[iCount]);
       //mC1->Add(gC1[iCount]); 
       //
       
//       if(*fun_C1){
//	 gfunC1[iCount] = new TGraph(*fun_C1);
//	 mC1->Add(gfunC1[iCount]);
//       }
       if(*fun_C2){
	 gfunC2[iCount] = new TGraph(*fun_C2);
	 mC2->Add(gfunC2[iCount]);
       }

       if(*fun_C3){
	 gfunC3[iCount] = new TGraph(*fun_C3);
	 mC3->Add(gfunC3[iCount]);
       }
       
       if(*fun_C4){
	 gfunC4[iCount] = new TGraph(*fun_C4);
	 mC4->Add(gfunC4[iCount]);
	 }
       iCount++;
       /*       cout<<"Entry # "<<jentry<<" iCount="<<iCount<<endl;
       cout<<"  --> time_C1 = " <<time_C1<<endl;
       cout<<"  --> tCFD_C1 = " <<timeCFD_C1<<endl;
       */
     }
     //
     //     if(!*fun_C4) null++;
     //
     //     if(iCount==Nmax) break;

     // time 
     last = current;
     time(&current);

     nevents++;
 
     DelTime = difftime (current,last);
     sumDelT += DelTime;
     sumDelT2 += DelTime*DelTime;
     AvgTime = sumDelT / double(nevents);
/*     
     if(nevents>1){
       SigTime = (sumDelT2 - (nevents)*AvgTime*AvgTime)/double(nevents-1);
       SigTime /= double(nevents);
       if(SigTime>=0.) SigTime = Sqrt(SigTime);
     }else{
       SigTime = 0.;
     }
     cout<<"Avg t="<<AvgTime<<" +/- "<<SigTime<<endl;*/
} // here is the end of loop over events
/////////////////////////////////////////////////////////

cout<<"ok (coincidences) = "<<ok<<endl;
//cout<<" *fun_C4 = NULL in "<<null<<"  times"<<endl;
//---------------------------------------------------------

fTime->Write("tdT");

fTime->mkdir("Pulses");
fTime->cd("Pulses");

//TCanvas *cC1 = new TCanvas("cC1","cC1",1200,800);
TCanvas *cC2 = new TCanvas("cC2","cC2",1200,800);
TCanvas *cC3 = new TCanvas("cC3","cC3",1200,800);
TCanvas *cC4 = new TCanvas("cC4","",1200,800);

//cC1->cd(); mC1->Draw("al"); cC1->Write();
cC2->cd(); mC2->Draw("al"); cC2->Write();
cC3->cd(); mC3->Draw("al"); cC3->Write();
cC4->cd(); mC4->Draw("al"); cC4->Write();

fTime->cd("..");
fTime->Close();

//delete cC1; 
delete cC2; 
delete cC3; 
delete cC4;
//  delete mC1; 
delete mC2;   
delete mC3;
delete mC4;

//delete [] gC1; 
delete [] gC2; 
delete [] gC3;
delete [] gC4;
//delete [] gfunC1; 
delete [] gfunC2; 
delete [] gfunC3;
delete [] gfunC4;

delete fun_C1;
delete fun_C2;
delete fun_C3;
delete fun_C4;
