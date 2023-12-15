//#include <iostream>
//#include <vector>
//
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
//
//#define CFDLEVEL 3. //2.5 //1.5
//#define Delay 10e-12 // internal delay of Luc cfd module
#define LevelOfNewPulse 1.0
//
// auxiliary
double AmplLevel;
double CFD_Delta_t;
//-----------------

class WhereFit {
public:
  int n2Fit; //how many points taken to a fit
  int i_begin,i_end;
  double Chi2_nof;
  TF1 fun;
  
  WhereFit():n2Fit(0),i_begin(0),i_end(1),Chi2_nof(-1){}
  WhereFit(short N2fit, int i_begin_, int i_end_, double Chi2_nof_, TF1 &fun_):
    n2Fit(N2fit),i_begin(i_begin_),i_end(i_end_),Chi2_nof(Chi2_nof_),fun(fun_){}
  int Get_Npoints() const{
    return (i_end - i_begin + 1);
  }
  inline bool operator<(const WhereFit &in) const;
};
bool WhereFit::operator<(const WhereFit &in) const{
  return (this->Chi2_nof < in.Chi2_nof)?kTRUE:kFALSE;
  }
/*bool WhereFit::operator<(const WhereFit &in) const{
  return ((this->Chi2_nof)/(Power(this->Get_Npoints(),2)) 
	  < in.Chi2_nof/Power(in.Get_Npoints(),2))?kTRUE:kFALSE;
	  }*/

double Time_(int Nmax, double *TabTime, double *TabAmpl, double Level){
  double timeout;
  int i;
  for(i=0;i<Nmax;i++) if(Abs(TabAmpl[i])>Abs(Level)) break;
  if(i == Nmax) return 1.;
  int iData(i);
  //  cerr<<"iData= "<<iData<<endl;
  //
  if(Abs(TabAmpl[iData]-Level) < Abs(Level-TabAmpl[iData-1])){
    timeout = TabTime[iData] -(TabTime[iData]-TabTime[iData-1])/(TabAmpl[iData]-TabAmpl[iData-1])*(TabAmpl[iData]-Level);
  }else{
    timeout = TabTime[iData-1] +(TabTime[iData]-TabTime[iData-1])/(TabAmpl[iData]-TabAmpl[iData-1])*(Level-TabAmpl[iData-1]);
  }
  return timeout;
}

double TimeCFD(int iData,int Shift, double *TabTime, double *TabAmpl){
  double timeout;
  double Div(3.),DeltaT(TabTime[iData]-TabTime[iData-1]);
  timeout = DeltaT*(TabAmpl[iData-1]/3. -TabAmpl[iData-1-Shift])/(TabAmpl[iData-Shift]-TabAmpl[iData-1-Shift]+(TabAmpl[iData-1]-TabAmpl[iData])/3.);
  timeout += TabTime[iData-1];
  return timeout;
}
//------
Double_t Time(int Nmax, double *T, double *V, double Level){
  Int_t i;
  for(i=0;i<Nmax;i++) if(Abs(V[i])>Abs(Level)) break;
  if(i == Nmax) return 1.;
  Double_t a((V[i]-V[i-1])/(T[i]-T[i-1]));//a = slope
  Double_t b(0.5*(V[i]+V[i-1]-a*(T[i]+T[i-1])));//b=intercept
  //
  return (Level - b)/a ;//return = the time of signal crossing Level
}
Double_t TimeCFD(int Nmax,int Shift,double Div, double Threshold, double *T, double *V){
  Int_t i;
  if(Shift>=Nmax-1) return 1.;
  //
  for(i=1+Shift;i<Nmax;i++){
    if(V[i-1]/Div > V[i-1-Shift] && V[i-Shift]> V[i]/Div 
       && Abs(V[i-1-Shift])>Abs(Threshold)){
      Double_t DeltaT(T[i]-T[i-1]);
      Double_t a_s((V[i-Shift]-V[i-1-Shift])/DeltaT);
      Double_t a_div((V[i]-V[i-1])/DeltaT/Div);
      Double_t b_s(0.5*(V[i-Shift]+V[i-1-Shift]-a_s*(T[i]+T[i-1])));
      Double_t b_div(0.5*((V[i]+V[i-1])/Div -a_div*(T[i]+T[i-1])));
      if(a_s!=a_div){
	return (b_div-b_s)/(a_s - a_div);
      }
    }
    /*else{
      cerr<<bool(V[i-1]/Div > V[i-1-Shift])<<"  "<<bool(V[i-Shift]> V[i]/Div)
	  <<"  "<<bool(Abs(V[i-1-Shift])>Abs(Threshold))<<endl;
	  }*/
  }
  return 1.;
}
//------
#define Nparam 2
Double_t PARAM[Nparam];

Double_t myFunc(Double_t *x,Double_t *par){
  // in fact this is a polynominal of power Nparam-1
  Double_t X = x[0]*1.e+9;
  Double_t val(0.);
  //
  for(short i=0;i<Nparam;i++){
    val += par[i]*Power(X,i);
    PARAM[i] = par[i];
  }
  /*
  Double_t val = par[0];
  PARAM[0] = par[0];
  for(short i=1;i<Nparam;i++){
    val *= (X-par[i]);
    PARAM[i] = par[i];
  }
  */
  return val;
}

void Zeros(Int_t &npar, Double_t *gin, Double_t &f, Double_t *x, Int_t iflag)
{
  Double_t delta;
// x is an answer - we look here for (x?)  ; same in CFD function.
  delta = (AmplLevel - myFunc(x,PARAM))*1.e+11;
  f = delta*delta;
  /*
  cout<<endl;
  for(short i=0;i<Nparam;i++){
    cout<<"par["<<i<<"]="<<PARAM[i]<<endl;
  }
  */
}
//#define DeltaT 140e-12
void CFD(Int_t &npar, Double_t *gin, Double_t &f, Double_t *x, Int_t iflag)
{
  Double_t delta;
  Double_t *p = new Double_t (x[0]-CFD_Delta_t);

  delta = (myFunc(p,PARAM) - myFunc(x,PARAM)/CFDLEVEL) *1.e+11; //3.)*1.e+11;
  f = delta*delta;
  delete p;
}
double TimeCFD(){
  return 0.;
}

// convert discrete signal to continuous signal
/*void cd2c(int Ndata, Double_t *time, Double_t *ampl)
{
  Double_t T = 2.5e-11;//25ps
  Double_t t;//variable of t
  Double_t pi = 3.1415926;
  TF1 *cd2c;
  for (i=0; i<=Ndata; i++)
  {
   cd2c += ampl[i]*(TMath::Sin(pi*(t-time[i])/T))/(pi*(t-time[i]/T));//
   
  }
  
}
*/

void  TimeFromFit(int Ndata,Double_t *time,Double_t *ampl,Double_t tmin, Double_t tmax,
		  TF1 **fun_, bool PositivePulse_, Double_t ampl_level,
		  /*results*/		  Double_t &timeLv,Double_t &timeCFD, Double_t &RiseTime,
		  Double_t CFD_Delta_t_,const char *type){
  //
  CFD_Delta_t = CFD_Delta_t_;
  AmplLevel = ampl_level;
  //
  Int_t iMax; Double_t vMax;
  Int_t iLow(-1);
  vMax = (PositivePulse_)? -10. : 10. ;
  //
  //The first time of definition to iMax
  if(strcmp(type,"cfd")==0){
    for(Int_t i=0;i<Ndata;i++){
      if( (ampl[i] > ampl_level && PositivePulse_ 
       ||  ampl[i] < ampl_level && !PositivePulse_  ) && time[i]>=tmin ){
	if(!PositivePulse_){
	  if(ampl[i]<vMax) { 
	  //if(ampl[i]<vMax||ampl[i]>-0.6) { 
	    vMax = ampl[i];
	    iMax = i;
	  }else{
	    break;
	  }
	}else{
	  // Positive pulse
	  if(ampl[i]>vMax) { 
	    vMax = ampl[i];
	    iMax = i;
	  }else{
	    break;
	  }
	}
      }
      if(time[i]>tmax ) break;
    }// end of for loop
    iLow = 0;
    for(Int_t i=iMax-1; i>=0; i--){
      if(Abs(ampl[i])<=Abs(ampl_level)){
	iLow = i;
	break;
      }
    }
  }else if(strcmp(type,"normal")==0){
  //
    for(Int_t i=0;i<Ndata;i++){
      if((ampl[i] > vMax && PositivePulse_  && time[i]>=tmin)
	 || (ampl[i] < vMax && !PositivePulse_  && time[i]>=tmin) ){
	vMax = ampl[i];
	iMax = i;//iMax is the item corresponding to max amplitude.
      }
      //
      if(time[i]>tmax ) break;
    }
    //
    //  iMax--;
  }else{
    cerr<<"Unknown type="<<type<<endl;
    exit(-2);
  }
  //
  
  ////////////////////Beginning of iLow iMax modification ///////////////
  
    
  //cout<<"iLow = "<<iLow<<"    iMax= "<<iMax<<endl; 
  Double_t LowFrac(0.1),HighFrac(0.7);//original is 0.2 and 0.8 separatelly.
  //
  const int nMax(iMax);
  double Abs_AmplLevel(Abs(AmplLevel));
  //for(Int_t i=nMax-1; i>=0; i--){
  for(Int_t i=nMax-1; i>0; i--){
    if(PositivePulse_ && ampl[i] <  Abs_AmplLevel){  iLow = i+1; break; }
    if(!PositivePulse_ && ampl[i]> -Abs_AmplLevel){  iLow = i+1; break; }
    if((!PositivePulse_ && ampl[i]>= LowFrac*ampl[nMax]) || (PositivePulse_ && ampl[i]<= LowFrac*ampl[nMax]) ){
      iLow =(Abs(ampl[i]-LowFrac*ampl[nMax])<Abs(ampl[i+1]-LowFrac*ampl[nMax]))?i:i+1;
      break; 
    }
  }
//  cout<<"first iMax =" <<iMax<<"\tiLow = "<<iLow<<endl;
  if(iLow == -1){//cout<<"amp[iMax]="<<ampl[iMax]<<"\ttime[iMax]="<<time[iMax]<<endl; 
                 //cout<<"amp[iLow]="<<ampl[iLow]<<"\ttime[iLow]="<<time[iLow]<<endl;
                 cerr<<"iLow=-1"<<endl; exit(2);}
  //The second time of definition to iMax
  for(Int_t i=nMax-1; i>=0; i--){
    if((!PositivePulse_ && ampl[i]> HighFrac*ampl[nMax]) || (PositivePulse_ && ampl[i]< HighFrac*ampl[nMax]) ){ 
      iMax = (Abs(ampl[i]-HighFrac*ampl[nMax])<Abs(ampl[i+1]-HighFrac*ampl[nMax]))?i:i+1;
      break;
    }
  }
    
    
  //////////////End of iLow & iMax modification //////////////////////

//  cout<<"iLow = "<<iLow<<"    iMax= "<<iMax<<endl; 
//  cout<<"second iMax =" <<iMax<<"\tiLow = "<<iLow<<endl;
  //  cout<<"\nvol_max="<<ampl[nMax]<<"   vol_0.9="<<ampl[iMax]<<" vol_0.1="<<ampl[iLow]<<endl;
  //  
  //
  RiseTime = (time[iMax]-time[iLow]);///(HighFrac - LowFrac);
  //
  vector<WhereFit> Fits;
  Int_t n2Fit(iMax-iLow+1);
  //  cout<<"NtoFit= "<<n2Fit<<endl;
  //  if(n2Fit<2) cerr<<"n2Fit="<<n2Fit<<endl;
  /*
  if(!strcmp(type,"cfd")){
    cout<<"NtoFit= "<<n2Fit<<endl;
    }*/
  //
  TGraph *gr = new TGraph(Ndata,time,ampl);
  Int_t Low(iLow),High(iMax);
  Double_t t_min(time[Low]), t_max(time[High]);
//cout<<"Line 252"<<endl;
//    cout<<"t_min, t_max="<<t_min<<", "<<t_max<<"\tvol_tmin="<<ampl[Low]<<", vol_tmax="<<ampl[High]<<endl;
  //  cout<<"  n2Fit="<<n2Fit<<endl;
  Int_t Min2Fit(5);//original is 20
  if(n2Fit>100){
    Low  = Int_t((1.7*iLow + 0.3*iMax)*0.5);
    High = Int_t((0.3*iLow + 1.7*iMax)*0.5);
    //cout<<"iLow,iMax="<<iLow<<" "<<iMax<<"  Low,High="<<Low<<" "<<High<<endl;
    //cout<<"Min2Fit="<<Min2Fit;
    Min2Fit = Max(Min2Fit,Int_t(0.9*(High-Low+1)));
    //cout<<"\tMin2Fit="<<Min2Fit<<endl;
    iLow = Low; iMax = High;
    n2Fit = iMax-iLow+1;
  }
  // cout<<"n2Fit\="<<n2Fit<<endl<<endl;

  if( !strcmp(type,"cfd") ) Min2Fit = 5;//original is 8
  if( !strcmp(type,"g1")  ) Min2Fit = 3;
//  cerr<<"Min2Fit, n2Fit  "<<Min2Fit<<"  "<<n2Fit;
  //
  while(n2Fit >= Min2Fit){
    //cout<<"Low="<<Low<<"  High="<<High<<"  n2Fit="<<n2Fit;
    t_min = time[Low]; t_max = time[High];
    //
    TF1 *fun=new TF1("func",&myFunc,t_min,t_max,Nparam);
    //
    fun->SetParameter(0,0.5);  
    if( PositivePulse_ ) { fun->SetParameter(1,2.);}
    else { fun->SetParameter(1,-2.);}
    gr->Fit(fun,"QF","same",t_min,t_max);    
    double Chi2ndf(fun->GetChisquare()/fun->GetNDF());
    //    if(Chi2ndf>2e-4) 
    Fits.push_back(WhereFit(n2Fit,Low,High,Chi2ndf,*fun));
    delete fun;
    //cout<<"  Chi2_nof="<<Fits.back().Chi2_nof<<endl;
      //
    if(Low == iLow){
      Low++;
      n2Fit = (High-Low+1);
      High = iMax; Low = High - n2Fit + 1;
    }else{
      Low--; High--;
    }
    //    cout<<"in do ... while     n2Fit="<<n2Fit<<endl;
  } // number points taken to a fit - not number degrees of freedom
  delete gr;
  //
  if(Fits.size() == 0 ){
     timeLv = 1.; timeCFD = 1.;
     return;
  }
     //     cout<<"Fits.size() = "<<Fits.size()<<endl;

  sort(Fits.begin(),Fits.end());
  //  if(Fits[0].Chi2_nof < 0.035) break;

  **fun_ = Fits[0].fun;
//  cerr<<"   chosen Npoints="<<Fits[0].n2Fit<<endl;
  /*
  fun_ = Fits[Fits.size()-1].fun;
  for(Int_t i;i<Fits.size();i++){
    if(Fits[i].Chi2_nof > 2e-4 && Fits[i].Chi2_nof < 0.035){
      fun_ = Fits[i].fun;
      break;
    }
    }*/
  /* for(Int_t i=0;i<Fits.size();i++){
    cout<<"Chi2_nof="<<Fits[i].Chi2_nof<<"  begin="<<Fits[i].i_begin
	<<"  end="<<Fits[i].i_end<<endl;
	}*/
  
  if(Nparam == 2){
    double b( (*fun_)->GetParameter(0) );
    double a( (*fun_)->GetParameter(1) );
    double tLv((ampl_level-b)/a*1e-9);//ample_level = C#level in "Analyse...time.C"
    double tCFD(CFDLEVEL*CFD_Delta_t_/(CFDLEVEL-1.) - b/a*1e-9);
    timeLv = tLv;//timeLv is the time that linear fit function cross the C#level
    timeCFD= tCFD;
    /*
    if(Abs(a*CFD_Delta_t_*1e+9/(CFDLEVEL-1.))<Abs(ampl_level) ||
       tCFD<time[iLow] || tCFD>time[iMax] ) timeCFD= 1.;
    //
    */
    //
    /*
    if(!strcmp(type,"cfd")){
      cout<<"a="<<a<<"  b="<<b<<endl;
      cout<<"tLv="<<timeLv<<"  tCFD="<<timeCFD<<"  time[iLow-1]="<<  time[iLow-1]<<"  time[iMax]="<<time[iMax]<<endl;
      }*/
    if(Abs(a*1e+9*(tCFD+Delay)+b)<Abs(ampl_level) 
       //  ||  tCFD<time[iLow-1] || tCFD>time[iMax] 
       ) timeCFD= 1.;
    //    if(tLv<time[iLow-1] || tLv>time[iMax])  timeLv = 1.;    
    /*
    if(Fits[0].n2Fit < 6){
      timeCFD= 1.; timeLv = 1.; 
      }*/
    //
  }else{
    //
    cerr<<"STRENGE !!!"<<endl;
    //
    Double_t arglist[10];
    Double_t vstart[1],step[1],ErrX0;
    Int_t ierflg = 0;
    //
    TMinuit *gEvalEq = new TMinuit(1);
    gEvalEq->SetPrintLevel(-1);
    gEvalEq->SetFCN(Zeros);
    arglist[0] = 1;
    gEvalEq->mnexcm("SET ERR", arglist ,1,ierflg);
    vstart[0]=0.5*(t_min+t_max); step[0]   = 1.e-14; 
    gEvalEq->mnparm(0, "x0", vstart[0], step[0], 0,0,ierflg);
    arglist[0] = 500;
    arglist[1] = 1.;
    //cout<<"TimeFromFit :: Zero"<<endl;
    gEvalEq->mnexcm("MIGRAD", arglist ,2,ierflg);
    //cout<<"Zero /end"<<endl;
    gEvalEq->GetParameter(0,timeLv,ErrX0);   
    //
    gEvalEq->SetFCN(CFD);
    arglist[0] = 1;
    gEvalEq->mnexcm("SET ERR", arglist ,1,ierflg);
    vstart[0]=0.5*(t_min+t_max); step[0]   = 1.e-14; 
    gEvalEq->mnparm(0, "x0", vstart[0], step[0], 0,0,ierflg);
    arglist[0] = 500;
    arglist[1] = 1.;
    //cout<<"TimeFromFit :: CFD"<<endl;
    gEvalEq->mnexcm("MIGRAD", arglist ,2,ierflg);
    //cout<<"CFD /end"<<endl;
    gEvalEq->GetParameter(0,timeCFD,ErrX0);   

    //  cout<<" analitical:  Lv="<<tLv<<"  Cfd="<<tCFD<<endl;
    //cout<<" fitted    :  Lv="<<timeLv<<"  Cfd="<<timeCFD<<endl;
    delete gEvalEq;
  }
  //
  return ;
}

PulseExtrema Measures(int Npoints, double* time, double* vol,
	      double tmin, double tmax, double level, double DeltaT,
	      TF1 **fun, bool PositivePulse_,
	      double &min_ampl, double &max_ampl, double &min_ampl2,
	      double &t_CL, double &t_fit_CL, double &t_fit_CFD, Double_t &RiseTime,
		      double &area, double half_width, double base_line_shift, const char *type){
  //
//  cout << "positive ? = " << PositivePulse_ <<endl;
  double level_(level);
  double min_time = 0;
  //  level *= 1.5;
  //
  PulseExtrema outExtrema;
  Int_t nExtrema(0);
  //
  min_ampl = 10.;
  max_ampl = -10.; 
  //
  for(int k=0;k<Npoints;k++){
    if(time[k] >= tmin && time[k] <= tmax){
      //if((vol[k])*(vol[k]) < (max_ampl)*(max_ampl)  ) max_ampl = vol[k]; 
      //if((vol[k])*(vol[k]) > (min_ampl)*(min_ampl)  ) min_ampl = vol[k]; 
      if(vol[k] > max_ampl && vol[k] < -1.0e-5 ) max_ampl = vol[k]; 
      if(vol[k] < min_ampl  ) 
        {
          min_ampl = vol[k]; 
          min_time = time[k];
        }
    }
  }
//    cout << "Npoints =" << Npoints <<endl;
  //
  t_CL = Time(Npoints,time,vol,level);
//cout << "positive ? = " << PositivePulse_ << "  , level=" <<level<<"  ,min="<<min_ampl<<"  ,max="<<max_ampl<<endl;
  //
  if((!PositivePulse_ && level<0. && level > min_ampl ) ||
     (PositivePulse_ && level>0. && level< max_ampl) ){ 
    //
//    cout << "line 424 headers.C" <<endl;
    Double_t localMinimum(10.),localMaximum(-10.);
    Int_t WhereExt(0);
    Bool_t foundMin(0),foundMax(0);
    for(int k=0;k<Npoints-1;k++){
      if(PositivePulse_){
	// new
	if( vol[k] > LevelOfNewPulse && level > 0 && vol[k]>localMaximum ){ 
	  localMaximum = vol[k]; 
	  WhereExt = k;
	}
	else if( vol[k]<=LevelOfNewPulse){
	  if( localMaximum > LevelOfNewPulse)  foundMax = 1;
	}
	//
	if(foundMax){
	  localMaximum = -10.;     foundMax = 0;
	  if(nExtrema<10){
	    outExtrema.Time[nExtrema] = time[WhereExt];
	    outExtrema.Vol[nExtrema]  = vol[WhereExt];
	    //	      WhereExt = 0;
	  }
	  outExtrema.Number = ++nExtrema;
	}
	// new
      }else{
	if( vol[k]< -LevelOfNewPulse  && level < 0 && vol[k]<localMinimum ){ 
	  localMinimum = vol[k]; 
	  WhereExt = k;
	}
	else if( vol[k]>=-LevelOfNewPulse){
	  if( localMinimum < -LevelOfNewPulse)  foundMin = 1;
	}
	//
	if(foundMin){
	  localMinimum = 10.;     foundMin = 0;
	  if(nExtrema<10){
	    outExtrema.Time[nExtrema] = time[WhereExt];
	    outExtrema.Vol[nExtrema]  = vol[WhereExt];
	    //	      WhereExt = 0;
	  }
	  outExtrema.Number = ++nExtrema;
	  // neeeeew :)
	  outExtrema.FindWidth(-LevelOfNewPulse,Npoints,time,vol,WhereExt);
	}
      }// end if
    }
    if(WhereExt>0 && WhereExt<Npoints && (foundMin || foundMax)){
      //      cout<<"Dodaje extra punkt :-)"<<endl;
      if(nExtrema<10){
	outExtrema.Time[nExtrema] = time[WhereExt];
	outExtrema.Vol[nExtrema]  = vol[WhereExt];
      }
      WhereExt = 0;
      outExtrema.Number = ++nExtrema;
    } 
  //  cout << "line 475 headers.C" <<endl;
    //
    TimeFromFit(Npoints,time,vol,tmin,tmax,fun,PositivePulse_,level_,t_fit_CL,t_fit_CFD,RiseTime,DeltaT, type);
  }else{
    t_fit_CL    = 1.;
    t_fit_CFD = 1.;
  }
  if( t_fit_CL  == 1.){
    delete *fun; *fun = NULL;
  }
  area = Area(time,vol,Npoints,tmin,tmax,min_time,half_width,base_line_shift);
  if(outExtrema.Number>0){
    double maxMin(-10.),minMax(10.);
    for(short n=0;n<Min(outExtrema.Number,10);n++){
      if(!PositivePulse_){
	if(outExtrema.Vol[n]>maxMin) maxMin = outExtrema.Vol[n];
      }else{ // PositivePulse_ = true
	if(outExtrema.Vol[n]<minMax) minMax = outExtrema.Vol[n];
      }
    }
    min_ampl2 =(!PositivePulse_)? maxMin : minMax;
  }else{
    min_ampl2 =(!PositivePulse_)? min_ampl : max_ampl;
  }
  return outExtrema;
}

void  TimeForCFD(Int_t Npoints, Double_t *time, Double_t *volt, Int_t pos, Double_t cutNoise, Double_t &t_min, Double_t &t_max) {
  // returns min and max time needed by 'TimeFromFit' function
  Double_t pointsDistance=0.03;
  Int_t i=0;

  if (pos==1) {
    for (i=0; i<Npoints && volt[i]<=cutNoise; i++)
      ;
    if (i<Npoints-1) {
      if ((volt[i]-volt[i-1])<pointsDistance)
        t_min=time[i-1];
      else
        t_min=time[i-2];  // the edge of pulse is very sharp so this should be sufficient
      if ((volt[i+2]-volt[i+1])<pointsDistance)
        t_max=time[i+1];
      else
        t_max=time[i+2];
    }
      //else
	//cout << "WARNING from 'TimeForCFD': Didn't pass cut on noise." << endl;
  }
  else {
    for (i=0; i<Npoints && volt[i]>=cutNoise; i++)
      ;
    if (i<Npoints-1) {
      if ((volt[i]-volt[i-1])>-pointsDistance)
        t_min=time[i-1];
      else
        t_min=time[i-2];  // the edge of pulse is very sharp so this should be sufficient
      if ((volt[i+2]-volt[i+1])>-pointsDistance)
        t_max=time[i+1];
      else
        t_max=time[i+2];
      }
      //else
	//cout << "WARNING from 'TimeForCFD': Didn't pass cut on noise." << endl;
    //}
  }
}
Double_t GetVolRMS(Int_t Nsize,Double_t *time, Double_t* vol,
		   Double_t tmin, Double_t tmax, TH1D* hist, Double_t &GetMean){
  for(int i=0;i<Nsize;i++){
    if(time[i]>tmin && time[i]<tmax){
      hist->Fill(vol[i]);
    }
    if(time[i]>tmax) break;
  }
  GetMean = hist->GetMean();
  return hist->GetRMS();
}

Int_t I_Max(double *tab, int nMax,double Level){
  int out(0); 
  for(int i=0;i<nMax;i++){
    if(tab[i]<Level) return i-1;
  }
  return out;
}

Int_t I_MaxPlus(double *tab, int nMax,double Level){
  int out(0);
  for(int i=0;i<nMax;i++){
    if(tab[i]>Level) return i-1;
  }
  return out;
}

Int_t I_Min(double *tab, int nMax){
  int out(0); 
  double Min_tab(1e+10);
  for(int i=0;i<nMax;i++){
    if(tab[i]<Min_tab){ 
      Min_tab = tab[i];
      out = i;
    }
  }
  
  for(int i=0;i<=out;i++){
    if(tab[i]<0.80*tab[out]) return i-1;
  }
  //cout<<"\nmin_tab["<<out<<"]= "<<Min_tab<<" = "<<tab[out]<<endl;
  return out;
}

Int_t I_MinPlus(double *tab, int nMax){
  int out(0);
  double Min_tab(1e+10);
  for(int i=0;i<nMax;i++){
    if(tab[i]>Min_tab){
      Min_tab = tab[i];
      out = i;
    }
  }
  /*
  for(int i=0;i<=out;i++){
    if(tab[i]<0.85*tab[out]) return i-1;
  }
  */
  //cout<<"\nmin_tab["<<out<<"]= "<<Min_tab<<" = "<<tab[out]<<endl;
  return out;
}

Int_t I_Max(double *tab, int nMax){
  int out(0); 
  double Min_tab(-1e+10);
  for(int i=0;i<nMax;i++){
    if(tab[i]>Min_tab){ 
      Min_tab = tab[i];
      out = i;
    }
  }
  for(int i=0;i<=out;i++){
    if(tab[i]>0.98*tab[out]) return i-1;
  }
  //cout<<"\nmin_tab["<<out<<"]= "<<Min_tab<<" = "<<tab[out]<<endl;
  return out;
}

/*Double_t Area(Double_t *time, Double_t *vol, Int_t nMax, Double_t tmin, Double_t tmax){
  //cout << "nMax = " << nMax <<endl<<endl; //nMax = 2002.
  //cout << "min_time = " << min_time <<endl<<endl; 
  Double_t area(0.),deltaT(time[1]-time[0]);
  for(Int_t i=1; i<nMax; i++){
    ////if( time[i]>= tmin && time[i]<tmax && vol[i-1] <= -1.0e-2 && vol[i] <= -1.0e-2) area += vol[i-1] + vol[i];
    ///if( time[i]>= tmin && time[i]<tmax )  area += vol[i-1] + vol[i];
  ///}
  
    if( time[i]>= tmin && time[i]<tmax )
    { 
     if(vol[i-1] < -1e-4 && vol[i] < -1e-4)  area += vol[i-1] + vol[i];//from -1e-8 to -1e-4.
     if(vol[i-1] < -1e-4 && vol[i] > -1e-4)  area += vol[i-1] -1* vol[i];
     if(vol[i-1] > -1e-4 && vol[i] < -1e-4)  area += -1*vol[i-1] + vol[i];
     if(vol[i-1] > -1e-4 && vol[i] > -1e-4)  area += -1*vol[i-1] -1*vol[i];
   //  if(vol[i-1] < -1e-2 && vol[i] < -1e-2)  area += vol[i-1] + vol[i];//from -1e-4 to -1e-2.
   //  if(vol[i-1] > 1.1e-2 && vol[i] > 1.1e-2)  area += -1*vol[i-1] -1*vol[i];
    }
}
  area *= 0.5* deltaT;
 // cout << "tmin = " << tmin << "\ttmax = " << tmax<<"\tdeltaT ="<<deltaT<<endl;
  area *= -1;
  return area;
}*/

//Double_t Half_width = half_width;
//Double_t Base_line_shift = base_line_shift;
Double_t Area(Double_t *time, Double_t *vol, Int_t nMax, Double_t tmin, Double_t tmax, Double_t min_time, Double_t half_width, Double_t base_line_shift){
  Double_t area(0.),deltaT(time[1]-time[0]);
  for(Int_t i=1; i<nMax; i++){
    if( time[i] >= (min_time - half_width) && time[i]<= (min_time + half_width))  
    {
     // area += vol[i-1] + vol[i] -0.08;//-0.08 = -0.04 for each vol[], for Ham signal.
      area += vol[i-1] + vol[i] + base_line_shift;//-0.016 = -0.003 for each vol[], for Photek signal.
  // if(vol[i-1] < -1e-4 && vol[i] < -1e-4)  area += vol[i-1] + vol[i]-0.04;//-0.04 for Photek signal.
  //   if(vol[i-1] < -1e-4 && vol[i] > -1e-4)  area += vol[i-1] -1* vol[i]-0.04;
  //   if(vol[i-1] > -1e-4 && vol[i] < -1e-4)  area += -1*vol[i-1] + vol[i]-0.04;
  //   if(vol[i-1] > -1e-4 && vol[i] > -1e-4)  area += -1*vol[i-1] -1*vol[i]-0.04;
    }
  }
  area *= 0.5* deltaT;
  area *= -1;
  return area;
}
