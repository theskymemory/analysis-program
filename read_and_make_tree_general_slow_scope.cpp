/* *********************************************************
********************************************************* */

/* Libraries */
#include <iostream>
using namespace std;
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include <stdio.h>

#include <sys/types.h>
#include <sys/dir.h>
#include <sys/param.h>
#include <time.h>
#include <unistd.h>
#include <signal.h>
#include <math.h>

/* Root Libraries */
#include "TTree.h"
#include "TStyle.h"
#include "TFile.h"
#include "TDatime.h"

#define NFiles     5100 //maximum number of files
#define MaxEventNum 500  // max number of events in one root file
#define VERY_BIG_int 9999999 // just a VERY B_I_G integer

int String2Int(char *);

int main(int argc, char **argv) {
  //
  char dumystr[50],dumystr2[50],dumystr3[50],dumystr4[50];
  char dumyline[200];
  int which(0);
  //
  // TRandom RAND;

  int file_select_C2(const struct dirent64 *);
  int file_select_C3(const struct dirent64 *);
  int file_select_C4(const struct dirent64 *);
  int file_select_C1(const struct dirent64 *);
  void readfile_i(ifstream&,int ,char* , double time[],double vol[]);
  void readfile(ifstream &datafile,char* fname,int Nlines,double area[],double counts[]);
  void CorrelatedFiles(struct dirent64 **,struct dirent64 **,
		     struct dirent64 **,struct dirent64 **,  
		       int ,int ,int ,int, int &);

  if (argc!=2) { std::cout << "Usage:\n read_and_make_tree <data-directory-path>" << std::endl; return 1; }

  char dir[256];
  strcpy(dir,argv[1]);
  struct dirent64 **filelist_C2 = new  struct dirent64* [NFiles];
  struct dirent64 **filelist_C3 = new  struct dirent64* [NFiles];
  struct dirent64 **filelist_C4 = new  struct dirent64* [NFiles];
  struct dirent64 **filelist_C1 = new  struct dirent64* [NFiles];
  for(int itab=0;itab<NFiles;itab++){
    filelist_C2[itab] = new struct dirent64;
    filelist_C3[itab] = new struct dirent64;
    filelist_C4[itab] = new struct dirent64;
    filelist_C1[itab] = new struct dirent64;
  }
  int count_C2(0),count_C3(0),count_C4(0),count_C1(0),count_OK(NFiles+1);
  count_C2 = scandir64(dir, &filelist_C2, file_select_C2, alphasort64);
  count_C3 = scandir64(dir, &filelist_C3, file_select_C3, alphasort64);
  count_C4 = scandir64(dir, &filelist_C4, file_select_C4, alphasort64);
  count_C1 = scandir64(dir, &filelist_C1, file_select_C1, alphasort64);
  //
  CorrelatedFiles(filelist_C1,filelist_C2,filelist_C3,filelist_C4,
		  count_C1,count_C2, count_C3, count_C4, count_OK);
  //
  if (  count_C2>NFiles ){ 
    std::cerr<<"\nERROR! file list array for C2 out of range ("<<count_C2<<">"<<NFiles<<")"<<std::endl;
    return 11;
  }
  if (  count_C3>NFiles ){ 
    std::cerr<<"\nERROR! file list array for C3 out of range ("<<count_C3<<">"<<NFiles<<")"<<std::endl;
    return 11;
  }
  if (  count_C4>NFiles ){ 
    std::cerr<<"\nERROR! file list array for C4 out of range ("<<count_C4<<">"<<NFiles<<")"<<std::endl;
    return 11;
  }
  if (  count_C1>NFiles ){ 
    std::cerr<<"\nERROR! file list array for C1 out of range ("<<count_C1<<">"<<NFiles<<")"<<std::endl;
    return 11;
  }

  Int_t Nevents(0),TotalNevents(0);
  Int_t Npoints,Npoints_C2(0),Npoints_C3(0),Npoints_C4(0),Npoints_C1(0); // number of points in one event
  {
    char  fname_C2[256],fname_C3[256],fname_C4[256],fname_C1[256];
    strcpy(fname_C2,dir);
    strcpy(fname_C3,dir);
    strcpy(fname_C4,dir);
    strcpy(fname_C1,dir);
    //
    Double_t dumdouble;
    if(count_C2>0){
      if ( fname_C2[strlen(fname_C2)-1] != '/' ) { strcat(fname_C2,"/"); }
      strcat(fname_C2,filelist_C2[0]->d_name);
      std::ifstream datafile_C2(fname_C2);
      if (!datafile_C2.good() || !datafile_C2.is_open() ){
	std::cout << "\nWARNING! Cannot read from file: " << fname_C2 << std::endl;
	return 1;
      }
      datafile_C2>>dumystr>>dumystr3>>Npoints_C2>>dumystr4>>dumdouble>>dumdouble;
      datafile_C2.close();
    }
    //
    if(count_C3>0){
      if ( fname_C3[strlen(fname_C3)-1] != '/' ) { strcat(fname_C3,"/"); }
      strcat(fname_C3,filelist_C3[0]->d_name);
      std::ifstream datafile_C3(fname_C3);
      if (!datafile_C3.good() || !datafile_C3.is_open() ){
	std::cout << "\nWARNING! Cannot read from file: " << fname_C3 << std::endl;
	return 1;
      }
      datafile_C3>>dumystr>>dumystr3>>Npoints_C3>>dumystr4>>dumdouble>>dumdouble;
      datafile_C3.close();
    }
    //
    //exit(1);// STOP
    if(count_C4>0){
      if ( fname_C4[strlen(fname_C4)-1] != '/' ) { strcat(fname_C4,"/"); }
      strcat(fname_C4,filelist_C4[0]->d_name);
      std::ifstream datafile_C4(fname_C4);
      if (!datafile_C4.good() || !datafile_C4.is_open() ){
	std::cout << "\nWARNING! Cannot read from file: " << fname_C4 << std::endl;
	return 1;
      }
      datafile_C4>>dumystr>>dumystr3>>Npoints_C4>>dumystr4>>dumdouble>>dumdouble;
      datafile_C4.close();
    }
    //
    if(count_C1>0){
      if ( fname_C1[strlen(fname_C1)-1] != '/' ) { strcat(fname_C1,"/"); }
      strcat(fname_C1,filelist_C1[0]->d_name);
      std::ifstream datafile_C1(fname_C1);
      if (!datafile_C1.good() || !datafile_C1.is_open() ){
	std::cout << "\nWARNING! Cannot read from file: " << fname_C1 << std::endl;
	return 1;
      }
      datafile_C1>>dumystr>>dumystr3>>Npoints_C1>>dumystr4>>dumdouble>>dumdouble;
      datafile_C1.close();
    }
  }
  //
  double C2_vol[Npoints_C2],C2_time[Npoints_C2];
  double C3_vol[Npoints_C3],C3_time[Npoints_C3];
  double time_C4[Npoints_C4],vol_C4[Npoints_C4];
  double time_C1[Npoints_C1],vol_C1[Npoints_C1];
  //
  ostringstream str_time1,str_vol1;
  ostringstream str_time2,str_vol2;
  ostringstream str_time3,str_vol3;
  ostringstream str_time4,str_vol4;
  //
  str_time1<<"C2_time["<<Npoints_C2<<"]/D"<<endl;
  str_vol1 <<"C2_vol["<<Npoints_C2<<"]/D"<<endl;
  //
  str_time2<<"C3_time["<<Npoints_C3<<"]/D"<<endl;
  str_vol2 <<"C3_vol["<<Npoints_C3<<"]/D"<<endl;
  //
  str_time3<<"time_C4["<<Npoints_C4<<"]/D"<<endl;
  str_vol3 <<"vol_C4["<<Npoints_C4<<"]/D"<<endl;
  //
  str_time4<<"time_C1["<<Npoints_C1<<"]/D"<<endl;
  str_vol4 <<"vol_C1["<<Npoints_C1<<"]/D"<<endl;
  //  
  TTree *tree = new TTree("scopetree","Data from a scope");
  if(count_C2>0){
    tree->Branch("C2_time",C2_time,str_time1.str().c_str());
    tree->Branch("C2_vol",C2_vol,str_vol1.str().c_str());
  }
  if(count_C3>0){
    tree->Branch("C3_time",C3_time,str_time2.str().c_str());
    tree->Branch("C3_vol",C3_vol,str_vol2.str().c_str());
  }
  if(count_C4>0){
    tree->Branch("C4_time",time_C4,str_time3.str().c_str());
    tree->Branch("C4_ampl",vol_C4,str_vol3.str().c_str());
  }
  if(count_C1>0){
    tree->Branch("C1_time",time_C1,str_time4.str().c_str());
    tree->Branch("C1_ampl",vol_C1,str_vol4.str().c_str());
  }
 // tree->Branch("Npoints",&Npoints,"Npoints/I");//copy from read_and_make_tree.general_cosmic.cpp

  //
  bool EndOfFile(0);
  TFile *f;
  // begin of a loop over files
  for(int i=1; i<count_OK+1; i++) {
    EndOfFile = 0;
    //
    char  fname_C2[256],fname_C3[256],fname_C4[256],fname_C1[256];
    ifstream datafile_C2,datafile_C3,datafile_C4,datafile_C1;
    if(count_C2>0){
      strcpy(fname_C2,dir);
      if ( fname_C2[strlen(fname_C2)-1] != '/' ) { strcat(fname_C2,"/"); }
      strcat(fname_C2,filelist_C2[i-1]->d_name);
      datafile_C2.open(fname_C2);
      if (!datafile_C2.good() || !datafile_C2.is_open() ){
	std::cout << "\nWARNING! Cannot read from file: " << fname_C2 << std::endl;
	return 1;
      }
    }
    //
    if(count_C3>0){
      strcpy(fname_C3,dir);
      if ( fname_C3[strlen(fname_C3)-1] != '/' ) { strcat(fname_C3,"/"); }
      strcat(fname_C3,filelist_C3[i-1]->d_name);
      datafile_C3.open(fname_C3);
      if (!datafile_C3.good() || !datafile_C3.is_open() ){
	std::cout << "\nWARNING! Cannot read from file: " << fname_C3 << std::endl;
	return 1;
      }
    }
    //
    if(count_C4>0){
      strcpy(fname_C4,dir);
      if ( fname_C4[strlen(fname_C4)-1] != '/' ) { strcat(fname_C4,"/"); }
      strcat(fname_C4,filelist_C4[i-1]->d_name);
      datafile_C4.open(fname_C4);
      if (!datafile_C4.good() || !datafile_C4.is_open() ){
	std::cout << "\nWARNING! Cannot read from file: " << fname_C4 << std::endl;
	return 1;
      }
    }
    //
    if(count_C1>0){
      strcpy(fname_C1,dir);
      if ( fname_C1[strlen(fname_C1)-1] != '/' ) { strcat(fname_C1,"/"); }
      strcat(fname_C1,filelist_C1[i-1]->d_name);
      datafile_C1.open(fname_C1);
      if (!datafile_C1.good() || !datafile_C1.is_open() ){
	std::cout << "\nWARNING! Cannot read from file: " << fname_C1 << std::endl;
	return 1;
      }
    }
    //
    //////////////////////////////////////////////////////////////
    do{
      for(int j=0;j<Npoints_C2;j++){
	// for loop over items in one event
	if(count_C2>0 && !datafile_C2.eof() ){
	  //	  cout<<"file="<<fname_C2<<endl;
	  readfile_i(datafile_C2,j,fname_C2,C2_time,C2_vol);
	}
      }
      for(int j=0;j<Npoints_C3;j++){
	// for loop over items in one event
	if(count_C3>0 && !datafile_C3.eof() ){
	  readfile_i(datafile_C3,j,fname_C3,C3_time,C3_vol);
	}
      }
      for(int j=0;j<Npoints_C4;j++){
	// for loop over items in one event
	if(count_C4>0 && !datafile_C4.eof() ){
	  readfile_i(datafile_C4,j,fname_C4,time_C4,vol_C4);
	}
      }
      for(int j=0;j<Npoints_C1;j++){
	// for loop over items in one event
	if(count_C1>0 && !datafile_C1.eof() ){
	  readfile_i(datafile_C1,j,fname_C1,time_C1,vol_C1);
	}
      }
      //
      if((count_C2>0 && datafile_C2.eof())||(count_C3>0 && datafile_C3.eof())
	 ||(count_C1>0 && datafile_C1.eof())||(count_C4>0 && datafile_C4.eof()) ) EndOfFile = 1;
      if(EndOfFile && count_C2 > 0 && datafile_C2.good()) cout<<"\nWARNING! File: "<<fname_C2<<" has not reached the end - contain more data then other files, but nevertheless files number "<<i-1<<" are closed"<<endl;
      if(EndOfFile && count_C3 > 0 && datafile_C3.good()) cout<<"\nWARNING! File: "<<fname_C3<<" has not reached the end - contain more data then other files, but nevertheless files number "<<i-1<<" are closed"<<endl;
      //
      if(EndOfFile && count_C1 > 0 && datafile_C1.good()) cout<<"\nWARNING! File: "<<fname_C1<<" has not reached the end - contain more data then other files, but nevertheless files number "<<i-1<<" are closed"<<endl;
      if(EndOfFile && count_C4 > 0 && datafile_C4.good()) cout<<"\nWARNING! File: "<<fname_C4<<" has not reached the end - contain more data then other files, but nevertheless files number "<<i-1<<" are closed"<<endl;
      //
      if( ! EndOfFile ){
	Nevents++;
	tree->Fill();
      }	
      //  ***  // end of for loop over items in particular event
      if( (Nevents >= MaxEventNum || (i == count_OK && EndOfFile)) ){
	which++;
	cout<<"\nTree number="<<which<<" created with events of data="<<Nevents<<" read";
	TotalNevents += Nevents;
	Nevents = 0;
	stringstream ss ( stringstream::out);
	ss<<"scope"<<which<<".root";
	f = new TFile(ss.str().c_str(),"RECREATE");
	tree->Write();
	f->Close();
	delete tree; tree = NULL;
	delete f; f = NULL;
	std::cout<<" and saved"<<std::endl;
	if(!(i == count_OK && EndOfFile)){
	  tree = new TTree("scopetree","Data from a scope");
	  if(count_C2>0){
	    tree->Branch("C2_time",C2_time,str_time1.str().c_str());
	    tree->Branch("C2_vol",C2_vol,str_vol1.str().c_str());
	  }
	  if(count_C3>0){
	    tree->Branch("C3_time",C3_time,str_time2.str().c_str());
	    tree->Branch("C3_vol",C3_vol,str_vol2.str().c_str());
	  }
	  if(count_C4>0){
	    tree->Branch("C4_time",time_C4,str_time3.str().c_str());
	    tree->Branch("C4_ampl",vol_C4,str_vol3.str().c_str());
	  }
	  if(count_C1>0){
	    tree->Branch("C1_time",time_C1,str_time4.str().c_str());
	    tree->Branch("C1_ampl",vol_C1,str_vol4.str().c_str());
	  }
	}else{
	  cout<<"\nTotal number of data events saved="<<TotalNevents <<endl<<endl;
	}
      }
    }while(! EndOfFile ); // end of loop over one file ?
    ///////////////////////////////////////////////////////////
  }//end loop over files
  // create a new ROOT file
   //   tree->Print();

  return 0;
}

int file_select_C2(const struct dirent64   *entry) {
  if ( strstr(entry->d_name, ".txt") && strstr(entry->d_name,"C2") )
    return 1;
  else 
    return 0;
}
int file_select_C3(const struct dirent64   *entry) {
  if ( strstr(entry->d_name, ".txt") && strstr(entry->d_name,"C3") )
    return 1;
  else 
    return 0;
}
int file_select_C4(const struct dirent64   *entry) {
  if ( strstr(entry->d_name, ".txt") && strstr(entry->d_name,"C4") )
    return 1;
  else 
    return 0;
}
int file_select_C1(const struct dirent64   *entry) {
  if ( strstr(entry->d_name, ".txt") && strstr(entry->d_name,"C1") )
    return 1;
  else 
    return 0;
}
int file_any(const struct dirent64   *entry) {
  if ( strstr(entry->d_name, ".txt") )
    return 1;
  else 
    return 0;
}

 void readfile_i(ifstream &datafile,int j,char* fname,double time[],double vol[]){ 
   char dumystr[50];
   int Npoints;
   double time_,vol_;
   if(j==0){
     datafile>>dumystr>>dumystr>>Npoints>>dumystr>>time_>>vol_;
   }else if( j==5 || j<=3 ){
     datafile>>dumystr>>dumystr>>dumystr>>dumystr>>time_>>vol_;
     //	cout<<"1 "<<dumystr<<",2 "<<dumystr2<<",3 "<<dumystr3<<",4 "<<dumystr4<<" - "<<time_C4[j]<<" "<<vol_C4[j]<<endl;
   }else{
     datafile>>time_>>vol_;
     //     cout<<time_<<vol_<<endl;
   }
   //   
   if( datafile.eof()){
     if( j>0 ){
       std::cerr<<"\nERROR: End of File "<<fname<<" has been reached before read all pieces of data\n"
		<<"- it read "<<j<<std::endl;
       datafile.close();
       exit(50);
     }else if(j==0){
       datafile.close();
     }
   }else{
     time[j] = time_; vol[j] = vol_;
   }
   //
   return;
 }
 void readfile(ifstream &datafile,char* fname,int Nlines,double area[],double counts[]){ 
   if(datafile.eof()){
     datafile.close();
     return;
   }
   char dumyline[200];
   int _nlines(0);
   double Area,Counts;
   datafile.getline(dumyline,200);
   cout<<"Nlines="<<Nlines<<" in file="<<fname<<endl;
   do{
     datafile>>Area>>Counts;
     if(!datafile.eof()){
       area  [_nlines] = Area;
       counts[_nlines] = Counts;
       _nlines++;
     }
   }while(!datafile.eof());
   //  
   if(_nlines!=Nlines && datafile.eof()){
     std::cerr<<"\nERROR: End of File "<<fname<<" has been reached before read all pieces of data\n"
	      <<"- it read "<<_nlines<<" of "<<Nlines<<std::endl;
     exit(50);
   }
   //
   return;
 }

void CorrelatedFiles(struct dirent64 **file1,struct dirent64 **file2,
		     struct dirent64 **file3,struct dirent64 **file4,  
		     int N1, int N2, int N3, int N4,
  /*results*/	     int &Ncorr){
 
  if(N1 == N2 && N2 == N3 && N3 == N4 && N4 == 0){
    Ncorr = 0;
    return;
  }
  int Min1,Min2,Min3,Min4;
  int Max1,Max2,Max3,Max4;
  char *pointer,name[6];
  
  if(N1>0){
    pointer = strstr(file1[0]->d_name,".txt");
    pointer -= 5;
    strcpy(name,pointer);
    Min1 = String2Int(name);
    //
    pointer = strstr(file1[N1-1]->d_name,".txt");
    pointer -= 5;
    strcpy(name,pointer);
    Max1 = String2Int(name);
  }else{
    Min1 = -1;
    Max1 = NFiles;
  }
  //
  if(N2>0){
    pointer = strstr(file2[0]->d_name,".txt");
    pointer -= 5;
    strcpy(name,pointer);
    Min2 = String2Int(name);
    //
    pointer = strstr(file2[N2-1]->d_name,".txt");
    pointer -= 5;
    strcpy(name,pointer);
    Max2 = String2Int(name);
  }else{
    Min2 = -1;
    Max2 = NFiles;
  }
  //
  if(N3>0){
    pointer = strstr(file3[0]->d_name,".txt");
    pointer -= 5;
    strcpy(name,pointer);
    Min3 = String2Int(name);
    //
    pointer = strstr(file3[N3-1]->d_name,".txt");
    pointer -= 5;
    strcpy(name,pointer);
    Max3 = String2Int(name);
  }else{
    Min3 = -1;
    Max3 = NFiles;
  }
  //
  if(N4>0){
    pointer = strstr(file4[0]->d_name,".txt");
    pointer -= 5;
    strcpy(name,pointer);
    Min4 = String2Int(name);
    //
    pointer = strstr(file4[N4-1]->d_name,".txt");
    pointer -= 5;
    strcpy(name,pointer);
    Max4 = String2Int(name);
  }else{
    Min4 = -1;
    Max4 = NFiles;
  }
  //
  //--- --- --- --- --- --- --- --- --- --- ---
  //
  int start(max(Min1,max(Min2,max(Min3,Min4))));
  int stop(min(Max1,min(Max2,min(Max3,Max4))));
  cout<<"Min1="<<Min1<<" Min2="<<Min2<<" Min3="<<Min3
      <<" Min4="<<Min4<<endl;
  cout<<"Max1="<<Max1<<" Max2="<<Max2<<" Max3="<<Max3
      <<" Max4="<<Max4<<endl;
  cout<<"Common Min="<<start<<", Max="<<stop<<endl; 
  N1 += (Min1 - start);
  N2 += (Min2 - start);
  N3 += (Min3 - start);
  N4 += (Min4 - start);
  //
  N1 += (stop - Max1);
  N2 += (stop - Max2);
  N3 += (stop - Max3);
  N4 += (stop - Max4);
  cout<<"n1="<<N1<<"  n2="<<N2<<"  n3="<<N3<<"  n4="<<N4<<endl;
  //
  Ncorr = VERY_BIG_int;
  if(N1>0)  Ncorr = N1;
  if(N2>0)  Ncorr = min(N2,Ncorr);
  if(N3>0)  Ncorr = min(N3,Ncorr);
  if(N4>0)  Ncorr = min(N4,Ncorr);

  if(Ncorr == VERY_BIG_int){
    cerr<<"There is no C1*.txt nor C2*.txt nor C3*.txt and even C4*.txt files\nSTOP"<<endl;
    exit(9999);
  }
  //
  struct dirent64 *file2cpy_1,*file2cpy_2,*file2cpy_3,*file2cpy_4;
  for(int i=0;i<Ncorr;i++){
    if(start==Min1 && start == Min2 && start == Min3
       && start == Min4) break;
    if(start > Min1 && Min1 >= 0 ){
      file2cpy_1 = file1[start-Min1 +i];
      file1[i] = file2cpy_1;
    }
    //
    if(start > Min2 && Min2 >= 0 ){
      file2cpy_2 = file2[start-Min2 +i];
      file2[i] = file2cpy_2;
    }
    //
    if(start > Min3 && Min3 >= 0 ){
      file2cpy_3 = file3[start-Min3 +i];
      file3[i] = file2cpy_3;
    }
    //
    if(start > Min4 && Min4 >= 0 ){
      file2cpy_4 = file4[start-Min4 +i];
      file4[i] = file2cpy_4;
    }
  }
//cout << "POINT 1" << endl;
  //
  char *Name1,*Name2,*Name3,*Name4;
  Name1 = new char[15]; Name2 = new char[15]; Name3 = new char[15]; Name4 = new char[15];
  int nowInt;
  for(int i=0;i<Ncorr;i++){
    if(Min1 >= 0){
      pointer = strstr(file1[i]->d_name,"Trace"); 
      strcpy(Name1,pointer);
    }else{ delete [] Name1; Name1 = NULL;}
    if(Min2 >= 0){
      pointer = strstr(file2[i]->d_name,"Trace"); 
      strcpy(Name2,pointer);
    }else{ delete [] Name2; Name2 = NULL;}
    if(Min3 >= 0){
      pointer = strstr(file3[i]->d_name,"Trace"); 
      strcpy(Name3,pointer);
    }else{ delete [] Name3; Name3 = NULL;}
    if(Min4 >= 0){
      pointer = strstr(file4[i]->d_name,"Trace"); 
      strcpy(Name4,pointer);
    }else{ delete [] Name4; Name4 = NULL;}
    //
//cout <<
//cout << "POINT 2" << endl;

    if( Name1 && Name2 && strcmp(Name1,Name2)!=0){
      cerr<<" Correlation lost !  C1= "<<file1[i]->d_name<<" but C2= "<<file2[i]->d_name<<endl;
      cerr<<" Move file with smaller number to some subdiractory i.e. \"NotNeeded\"\nSTOP"<<endl;
      exit(13);
    }
    if( Name1 && Name3 && strcmp(Name1,Name3)!=0){
      cerr<<" Correlation lost !  C1= "<<file1[i]->d_name<<" but C3= "<<file3[i]->d_name<<endl;
      cerr<<" Move file with smaller number to some subdiractory i.e. \"NotNeeded\"\nSTOP"<<endl;
      exit(13);
    }
    if( Name1 && Name4 && strcmp(Name1,Name4)!=0){
      cerr<<" Correlation lost !  C1= "<<file1[i]->d_name<<" but C4= "<<file4[i]->d_name<<endl;
      cerr<<" Move file with smaller number to some subdiractory i.e. \"NotNeeded\"\nSTOP"<<endl;
      exit(13);
    }
    if( Name3 && Name2 && strcmp(Name3,Name2)!=0){
      cerr<<" Correlation lost !  C2= "<<file2[i]->d_name<<" but C3= "<<file3[i]->d_name<<endl;
      cerr<<" Move file with smaller number to some subdiractory i.e. \"NotNeeded\"\nSTOP"<<endl;
      exit(13);
    }
    if( Name4 && Name2 && strcmp(Name4,Name2)!=0){
      cerr<<" Correlation lost !  C2= "<<file2[i]->d_name<<" but C4= "<<file4[i]->d_name<<endl;
      cerr<<" Move file with smaller number to some subdiractory i.e. \"NotNeeded\"\nSTOP"<<endl;
      exit(13);
    }
    if( Name3 && Name4 && strcmp(Name3,Name4)!=0){
      cerr<<" Correlation lost !  C3= "<<file3[i]->d_name<<" but C4= "<<file4[i]->d_name<<endl;
      cerr<<" Move file with smaller number to some subdiractory i.e. \"NotNeeded\"\nSTOP"<<endl;
      exit(13);
    }
    //
    stringstream sch;
    pointer = NULL;
    if(N1>0) pointer = strstr(file1[i]->d_name,".txt"); 
    if(N2>0) pointer = strstr(file2[i]->d_name,".txt"); 
    if(N3>0) pointer = strstr(file3[i]->d_name,".txt"); 
    if(N4>0) pointer = strstr(file4[i]->d_name,".txt"); 
    if(!pointer){
      cerr<<"No files at all N1<0 && N2<0 && N3<0 && N4<0\nSTOP"<<endl;
      exit(11);
    }
    //
    for(short k=0;k<5;k++) sch<<*(pointer-5+k);
    sch>>nowInt;
    if(nowInt == stop) break;
  }
  if(!Name1){ delete [] Name1; }else{ delete Name1; Name1 = NULL;}
  if(!Name2){ delete [] Name2; }else{ delete Name2; Name2 = NULL;}
  if(!Name3){ delete [] Name3; }else{ delete Name3; Name3 = NULL;}
  if(!Name4){ delete [] Name4; }else{ delete Name4; Name4 = NULL;}
  //
  cout<<"All Files = "<<Ncorr<<endl;
  return;
}
int String2Int(char *str){
  stringstream ss;
  for(short i=0;i<5;i++){
    ss<<str[i];
  }
  int number;
  ss>>number;
  return number;
}
