#include "TMath.h"
using namespace TMath;

#include "TGraph.h"
// for TrackFinder function
struct TrackStat {
  Short_t nZerosX;
  Short_t nZerosY;
  Short_t nBadX;
  Short_t nBadY;
  Short_t nEventsX;
  Short_t nEventsY;
  Int_t eventCut;
};

class Extrema {
 public:
  Int_t Number;
  Double_t Time[10],Vol[10];
  Extrema():Number(0){};
};
class PulseExtrema: public Extrema {
 public:
  Double_t MaxWidth; //maximal width of separate extremas widths 
  PulseExtrema():Extrema(),MaxWidth(0.){}
  //
  void FindWidth(Double_t level,int npoints, Double_t *time, Double_t *vol,Int_t iWhere){
    int i(iWhere-1);
    double tStart(1.),tEnd(1.);
    do{
      if(i<0) break;
      tStart = time[i--];
    }while(Abs(vol[i])>Abs(level));
    i = iWhere+1;
    do{
      if(i>=npoints) break;
      tEnd = time[i++];
    }while(Abs(vol[i])>Abs(level));
    if(tStart < 1. && tEnd < 1.) MaxWidth = Max(MaxWidth,tEnd-tStart);
    return;
  }
};
//
//
double Time(int iData, double *TabTime, double *TabAmpl, double Level);
double Time_(int iData, double *TabTime, double *TabAmpl, double Level);
double TimeCFD(int iData,int Shift, double *TabTime, double *TabAmpl);
Double_t TimeCFD(int Nmax,int Shift,double Div, double Threshold, double *T, double *V);

Double_t myFunc(Double_t *x,Double_t *par);
//void Zeros(Int_t &npar, Double_t *gin, Double_t &f, Double_t *x, Int_t iflag);
//Double_t TimeFromFit(int Ndata,Double_t *time,Double_t *ampl,double,double, TF1 &fun_,
//                   bool PositivePulse_,Double_t ampl_level=1.4, Double_t CFD_Delta_t_=140e-12);
void  TimeFromFit(int Ndata,Double_t *time,Double_t *ampl,Double_t tmin, Double_t tmax,
		  TF1 **fun_, bool PositivePulse_, Double_t ampl_level,
		  /*results*/		  Double_t &timeLv,Double_t &timeCFD, Double_t &RiseTime,
		  Double_t CFD_Delta_t_,const char *type);
PulseExtrema Measures(int Npoints, double* time, double* vol,
	      double tmin, double tmax, double level, double DeltaT,
	      TF1 **fun, bool PositivePulse_,
	      double &min_ampl, double &max_ampl, double &min_ampl2,
	      double &t_CL, double &t_fit_CL, double &t_fit_CFD, Double_t &RiseTime,
	      double &area,double half_width, double base_line_shift, const char type[]="normal");
void TimeForCFD(Int_t Npoints, Double_t *time, Double_t *volt, Int_t pos, Double_t cutNoise, Double_t &t_min, Double_t &t_max);
Double_t GetVolRMS(Int_t Nsize,Double_t *time, Double_t* vol, Double_t tmin, Double_t tmax, TH1D* hist, Double_t &GetMean);
Int_t I_Max(double *tab, int nMax,double Level);
Int_t I_MaxPlus(double *tab, int nMax,double Level);
Int_t I_Min(double *tab, int nMax);
Int_t I_MinPlus(double *tab, int nMax);
Double_t Area(Double_t *time, Double_t *vol, Int_t nMax, Double_t tmin, Double_t tmax, Double_t min_time, Double_t Half_width, Double_t Base_line_shift);

void FindGoodEvents(TTree *hitsX,Short_t &NEventsX,Short_t *eventCutX,vector<Float_t> &resultX);
void TrackFinder(TFile *f,TFile *fout,Short_t *eventCutX,TrackStat &outData,vector<Float_t> &resultX,vector<Float_t> &resultY);
void cdc(int Ndata, Double_t *time, Double_t *ampl);
