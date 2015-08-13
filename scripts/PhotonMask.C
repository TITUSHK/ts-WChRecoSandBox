#include "TRandom3.h"
#include "TMath.h"
#include <iomanip>
#include <iostream>
#include "../include/WCLTreeWriter.hh"
#include "../include/SandBoxPMTcoverage.hh"

using namespace std;

void PhotonMask(TString filename="../../FullEvent.root", TString genfilename="../../generatorcardfile.root",TString outfilename="out.root")
{
  //gSystem->Load("../lib/libWCLAnalysis.so");

  //************************************************************//
  //								//
  // 	Set the value of "opt" to choose between 		//
  //	cylinder or cube detector.				//
  //								//
  // 	Set the value of "LAPPDs" to choose if there is		//
  //	LAPPDs or not.						//
  //								//
  //	And set other values to choose the time resolution,	//
  //	quantum efficiency, size, shape of the PMTs or LAPPDs	//
  //								//
  //************************************************************//

  TRandom3 numberp(42484);

  // Load Data
  // =========
  // WCLTREEREADER
  //**********************************************************
  // first parameter in the construct: 0=direct output of WCLite
  //                                   1=a "smeared" file made from WCL files
  // second parameter in the constructor: 1 means "include gen level tree"
  //                                      0 means no gen level tree
  WCLTreeReader *mTR = new WCLTreeReader(0,1);



  // LoadData is an overloaded function. If you specify only one name,
  // only the main tree is loaded. If you specify two files, it loads
  // both the geant output file and the file containing the gen tree
  // (in that order)
  mTR->LoadData(filename,genfilename);

  cout << "Done loading input data " << endl;
  // WCLTREEWRITER
  // *********************************************************
  // specify the output filename and the 2nd parameter is the  mode
  // mode=2 is for actual PMT hits
  WCLTreeWriter *mTW = new WCLTreeWriter(outfilename,2);

  // set the option: 0=TITUS, 1=ANNIE
  int opt = 0;

  // set the presence of LAPPDs
  bool LAPPDs = false;

  // set the dimensions of the detector.
  double Rdet = 5500.;
  double LdetCylinder = 22000.; 
  double LdetCube = 3000.;
  double Ldet;

  // set the time resolution
  double tResPMT = 2.5;
  double tResLAPPD = 0.1;
  // set the quantum efficiency
  double QE_PMT = 22.;
  double QE_LAPPD = 22.;
  double QE_Wavelengths[20] = {280., 300., 320., 340., 360., 380., 400., 420., 440., 460., 480., 500., 520., 540., 560., 580., 600., 620., 640., 660.};
  double QE_Factor[20] = {0.00, .066, .405, .801, .962, .976, 1., .957, .899, .791, .664, .550, .382, .205, .126, .069, .036, .024, .007, 0.00};
  // double the QE if true because WChSandbox divide by 2 the number of photons
  bool QEset = true;
  // set the size of the detectors
  double sizePMT = 304.8;
  double sizeLAPPD = 203.2;
  double Sdet;
  // set the PMT coverage
  double cover = 40.;
  // set the shape of the detectors: 0=circular, 1=square
  int shapePMT = 0;
  int shapeLAPPD = 1;

  // Number of PMTs in each wall
  int NFront = 0;
  int NBack = 0;
  int NrowCurve = 0;
  int NcolCurve = 0;
  // set the number of PMTs in each wall for the option=1.
  int NbyN = 5;

  if(shapePMT==0) Sdet=3.14*sizePMT*sizePMT/4;
  if(shapePMT==1) Sdet=sizePMT*sizePMT;
  double Ntot = (Rdet*Rdet*3.14+LdetCylinder*2*Rdet*3.14)/Sdet*cover/100;
  if(opt==0) {
    NFront = (int) (sqrt(0.8*Ntot/3.14)-fmod(sqrt(0.8*Ntot/3.14),1));
    NBack = (int) (sqrt(0.8*Ntot/3.14)-fmod(sqrt(0.8*Ntot/3.14),1));
    NcolCurve = (int) (sqrt(0.4*3.14*Ntot)-fmod(sqrt(0.4*3.14*Ntot),1));
    NrowCurve = (int) ((0.8*Ntot/NcolCurve)-fmod(0.8*Ntot/NcolCurve,1));
  }
  if(opt==1) {
    NbyN = (int) (sqrt(Ntot/6)-fmod(sqrt(Ntot/6),1));
  }

  if(QEset){
    QE_PMT = QE_PMT*2.;
    QE_LAPPD = QE_LAPPD*2.;
  }

  // SANDBOXPMTCOVERAGE
  //**********************************************************
  // specify NbyN = the number of PMTs in a single row (or column) 
  // in a square grid.
  SandBoxPMTcoverage* sbPMT = new SandBoxPMTcoverage();

  // set the dimensions of the detector L/2, R
  if( opt==0 ){
    Ldet=LdetCylinder;
  }
  else {
    Ldet=LdetCube;
  }
  sbPMT->SetBoxDimensions(opt,Ldet/2.,Rdet);

  //configure each wall for cylinder option
  if(opt==0){
    //                 front-wall, detectors
    sbPMT->SetWallConfiguration(5,1,shapePMT,sizePMT,NFront,NFront,QE_PMT, QE_Wavelengths, QE_Factor,shapeLAPPD,sizeLAPPD,NFront,NFront,QE_LAPPD);
    cout << "Front wall:  " << NFront << " rows/cols of PMTs" << endl;
    //                back-wall, detectors
    sbPMT->SetWallConfiguration(6,1,shapePMT,sizePMT,NBack,NBack,QE_PMT, QE_Wavelengths, QE_Factor,shapeLAPPD,sizeLAPPD,NBack,NBack,QE_LAPPD);
    cout << "Back wall:  " << NBack << " rows/cols of PMTs" << endl;
    //                 curve-wall, detectors 
    sbPMT->SetWallConfiguration(0,1,shapePMT,sizePMT,NrowCurve,NcolCurve,QE_PMT, QE_Wavelengths, QE_Factor,shapeLAPPD,sizeLAPPD,NrowCurve,NcolCurve,QE_LAPPD);
    cout << "Barrel:  " << NrowCurve << " by " << NcolCurve << " PMTs" << endl;
  }

  //configure each wall for cube option
  if(opt==1){
    // walls are all the same.
    for(int ii=1;ii<7;ii++){
     sbPMT->SetWallConfiguration(ii,1,shapePMT,sizePMT,NbyN,NbyN,QE_PMT,QE_Wavelengths, QE_Factor, shapeLAPPD,sizeLAPPD,NbyN,NbyN,QE_LAPPD);
    }
  }

  int nentries = mTR->GetEntries();
  cout<<nentries<<endl;

  int maxNpmts = opt==0?NFront*NFront+NBack*NBack+NcolCurve*NrowCurve:6*NbyN*NbyN;
  int * pmtIDs = new int[maxNpmts];
  double * pmtXs = new double[maxNpmts];
  double * pmtYs = new double[maxNpmts];
  double * pmtZs = new double[maxNpmts];
  int nPMT = 0;
  double size = LAPPDs ? sizeLAPPD : sizePMT;
  double qe = LAPPDs ? QE_LAPPD : QE_PMT;
  if(QEset) qe/=2.0;
  double tRes = LAPPDs ? tResLAPPD : tResPMT;
  int pmtID;
  //add PMTs for cylinder option
  if(opt==0){
    double pmtX, pmtY, pmtZ, lrow, lcol;
    //Front wall
    lrow = (2.*Rdet)/NFront;
    lcol = (2.*Rdet)/NFront;
    pmtZ = Ldet/2.;
    for(int row = 0; row < NFront; row++){
      pmtX = -Rdet+lrow*(row+0.5);
      for(int col = 0; col < NFront; col++){
        pmtY = -Rdet+lcol*(col+0.5);
        if(TMath::Sqrt(pmtX*pmtX+pmtY*pmtY)+size/2<Rdet){
          pmtID = (row - (NFront/2)) + (col - (NFront/2))*100 + 50000;
          mTW->AddPMT(nPMT,pmtX,pmtY,pmtZ,size,qe,tRes,LAPPDs);
          pmtXs[nPMT]=pmtX;
          pmtYs[nPMT]=pmtY;
          pmtZs[nPMT]=pmtZ;
          pmtIDs[nPMT++]=pmtID;
//          cout << pmtID << " " << row << " " << col << " (" << hitPMTx << ", " << hitPMTy << ", " << hitPMTz << ")" << endl;
        }
      }
    }
    //Back wall
    lrow = (2.*Rdet)/NBack;
    lcol = (2.*Rdet)/NBack;
    pmtZ = -Ldet/2.;
    for(int row = 0; row < NBack; row++){
      pmtX = -Rdet+lrow*(row+0.5);
      for(int col = 0; col < NBack; col++){
        pmtY = -Rdet+lcol*(col+0.5);
        if(TMath::Sqrt(pmtX*pmtX+pmtY*pmtY)+size/2<Rdet){
          pmtID = (row - (NBack/2)) + (col - (NBack/2))*100 + 60000;
          mTW->AddPMT(nPMT,pmtX,pmtY,pmtZ,size,qe,tRes,LAPPDs);
          pmtXs[nPMT]=pmtX;
          pmtYs[nPMT]=pmtY;
          pmtZs[nPMT]=pmtZ;
          pmtIDs[nPMT++]=pmtID;
//          cout << pmtID << " " << row << " " << col << " (" << hitPMTx << ", " << hitPMTy << ", " << hitPMTz << ")" << endl;
        }
      }
    }
    //Barrel
    lrow = Ldet/NrowCurve;
    lcol = (2.*3.1416)/NcolCurve;
    double pmtTheta;
    for(int row = 0; row < NrowCurve; row++){
      pmtZ = -Ldet/2.+lrow*(row+0.5);
      for(int col = 0; col < NcolCurve; col++){
        pmtTheta = -3.1416+lcol*(col+0.5);
        pmtX = Rdet*TMath::Cos(pmtTheta);
        pmtY = Rdet*TMath::Sin(pmtTheta);
        pmtID = (row - (NrowCurve/2)) + (col - (NcolCurve/2))*100;
        mTW->AddPMT(nPMT,pmtX,pmtY,pmtZ,size,qe,tRes,LAPPDs);
          pmtXs[nPMT]=pmtX;
          pmtYs[nPMT]=pmtY;
          pmtZs[nPMT]=pmtZ;
        pmtIDs[nPMT++]=pmtID;
//        cout << pmtID << " " << row << " " << col << " (" << hitPMTx << ", " << hitPMTy << ", " << hitPMTz << ")" << endl;
      }
    }
  }
  //add PMTs for cube option
  else if(opt==1){
    double pmtRow, pmtCol;
    double l = (2.*Ldet)/NbyN;
    for(int row = 0; row < NbyN; row++){
      pmtRow = -Ldet+l*(row+0.5);
      for(int col = 0; col < NbyN; col++){
        pmtCol = -Ldet+l*(col+0.5);
        pmtID = (row - (NbyN/2)) + (col - (NbyN/2))*100;
        mTW->AddPMT(nPMT,Ldet,pmtRow,pmtCol,size,qe,tRes,LAPPDs);
        pmtIDs[nPMT++]=pmtID+10000;
        mTW->AddPMT(nPMT,-Ldet,pmtRow,pmtCol,size,qe,tRes,LAPPDs);
        pmtIDs[nPMT++]=pmtID+20000;
        mTW->AddPMT(nPMT,pmtRow,Ldet,pmtCol,size,qe,tRes,LAPPDs);
        pmtIDs[nPMT++]=pmtID+30000;
        mTW->AddPMT(nPMT,pmtRow,-Ldet,pmtCol,size,qe,tRes,LAPPDs);
        pmtIDs[nPMT++]=pmtID+40000;
        mTW->AddPMT(nPMT,pmtRow,pmtCol,Ldet,size,qe,tRes,LAPPDs);
        pmtIDs[nPMT++]=pmtID+50000;
        mTW->AddPMT(nPMT,pmtRow,pmtCol,-Ldet,size,qe,tRes,LAPPDs);
        pmtIDs[nPMT++]=pmtID+60000;
      }
    }
  }

  // Loop over number of entries
  for(int i=0; i<nentries; i++){
    if(i%10==0) cout<<"event: "<<i<<endl;

    // Load the event variables into TreeReader
    mTR->LoadEvent(i);

    /// grab the relevant variables from the TreeReader event
//    Int_t* phot_isScat = mTR->get_phot_isScat();
//    Int_t* phot_processStart = mTR->get_phot_processStart();
//    Int_t* phot_trackid = mTR->get_phot_trackid(); 
//    Int_t* phot_parentid = mTR->get_phot_parentid();
//    Int_t* phot_capnum = mTR->get_phot_capnum();
//    Double_t* phot_xStart = mTR->get_phot_xStart();
//    Double_t* phot_yStart = mTR->get_phot_yStart();
//    Double_t* phot_zStart = mTR->get_phot_zStart();
//    Double_t* phot_tStart = mTR->get_phot_tStart();
    Double_t* phot_xEnd = mTR->get_phot_xEnd();
    Double_t* phot_yEnd = mTR->get_phot_yEnd();
    Double_t* phot_zEnd = mTR->get_phot_zEnd();
    Double_t* phot_tEnd = mTR->get_phot_tEnd();
    Double_t* phot_wavelength = mTR->get_phot_wavelength();

    //Intialize a new event
    mTW->InitializeEvent();
    //Add all of the information from the event, except for the photon hits
    mTW->AddWholeBranches(mTR,0,0,1,1,1);

    //Loop over photons to add those which hit a PMT or LAPPD
    int nphot = mTR->get_nphot();
    for(int k=0; k<nphot; k++){
      //cut photons outside 300-600nm for now, use proper wavelength dependent QE in future
      //if(phot_wavelength[k] > 600 || phot_wavelength[k] <300) continue;
//      if( (k%10000==0) && (i%10==0) ) cout<<k<<" photons out of "<<nphot<<endl;  
      int PMTid=0;
//      int isHit=0;
//      int isscat = phot_isScat[k];
//      int process = phot_processStart[k];
//      int trackid = phot_trackid[k];
//      int parentid = phot_parentid[k];
//      int capnum = phot_capnum[k];

//      double xS = phot_xStart[k]; double yS = phot_yStart[k]; double zS = phot_zStart[k]; 
//      double tS = phot_tStart[k];
      double xE = phot_xEnd[k];
      double yE = phot_yEnd[k];
      double zE = phot_zEnd[k];
      double tE = phot_tEnd[k];
      double wl = phot_wavelength[k];

/*
      int Nrow;
      int Ncol;
      int mcase;
      if(opt==1){
        Nrow = NbyN;
        Ncol = NbyN;
        if(xE==Ldet/2.) mcase=1;
        if(xE==-Ldet/2.) mcase=2;
        if(yE==Ldet/2.) mcase=3;
        if(yE==-Ldet/2.) mcase=4;
        if(zE==Ldet/2.) mcase=5;
        if(zE==-Ldet/2.) mcase=6;
      }
      else {
        if(zE==Ldet/2.){
          Nrow = NFront;
          Ncol = NFront;
          mcase=5;
        }
        if(zE==-Ldet/2.){
          Nrow = NBack;
          Ncol = NBack;
          mcase=6;
        }
        if((xE*xE + yE*yE)==(Rdet*Rdet)){
          Nrow = NrowCurve;
          Ncol = NcolCurve;
          mcase=0;
        }
      }*/

      int WhichDet = sbPMT->isActiveHit(xE, yE, zE, wl, PMTid, LAPPDs); // 0=none, 1=PMT, 2=LAPPD
      if(WhichDet==0)
        continue;
      double timeRes=WhichDet==1 ? tResPMT : tResLAPPD;
//      isHit=2;
      double hitTime = numberp.Gaus(tE,timeRes);
      int hitPMT = 0;
      while(pmtIDs[hitPMT]!=PMTid)
        hitPMT++;
//      if(fabs(pmtXs[hitPMT]-xE)>size/2.) cout << "Bad X " << PMTid << " " << pmtXs[hitPMT] << " " << xE << endl;
//      if(fabs(pmtYs[hitPMT]-yE)>size/2.) cout << "Bad Y " << PMTid << " " << pmtYs[hitPMT] << " " << yE << endl;
//      if(fabs(pmtZs[hitPMT]-zE)>size/2.) cout << "Bad Z " << PMTid << " " << pmtZs[hitPMT] << " " << zE << endl;
      mTW->AddHit(hitTime, pmtXs[hitPMT], pmtYs[hitPMT], pmtZs[hitPMT], hitPMT);
    }
    mTW->FillEvent();
  }
  mTW->WriteTreeToFile();
  delete mTW;
  return;
}
