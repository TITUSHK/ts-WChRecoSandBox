//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 11 16:57:41 2015 by ROOT version 5.34/26
// from TTree LikelihoodGenerator/LikelihoodGenerator
// found on file: mu_FullEvent.root
//////////////////////////////////////////////////////////

#ifndef LikelihoodGenerator_h
#define LikelihoodGenerator_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH3D.h>
#include <TString.h>
#include <iostream>
#include <TObject.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class LikelihoodGenerator : public TObject {
public :
   TH2D           *hPDF;
   TH2D           *hPhotonCount;
   TH2D           *hTime;
//   TH3D           *h3;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           evt;
   Int_t           nphot;
   Int_t           npart;
   Int_t           ncapturecount;
   Int_t           neutroncount;
   Double_t        phot_xStart[10000000];   //[nphot]
   Double_t        phot_yStart[10000000];   //[nphot]
   Double_t        phot_zStart[10000000];   //[nphot]
   Double_t        phot_tStart[10000000];   //[nphot]
   Double_t        phot_xEnd[10000000];   //[nphot]
   Double_t        phot_yEnd[10000000];   //[nphot]
   Double_t        phot_zEnd[10000000];   //[nphot]
   Double_t        phot_tEnd[10000000];   //[nphot]
   Double_t        phot_wavelength[10000000];   //[nphot]
   Int_t           phot_processStart[10000000];   //[nphot]
   Int_t           phot_isScat[10000000];   //[nphot]
   Int_t           phot_parentid[10000000];   //[nphot]
   Int_t           phot_trackid[10000000];   //[nphot]
   Int_t           phot_hit[10000000];   //[nphot]
   Int_t           phot_capnum[10000000];   //[nphot]
   Double_t        part_xStart[10000];   //[npart]
   Double_t        part_yStart[10000];   //[npart]
   Double_t        part_zStart[10000];   //[npart]
   Double_t        part_tStart[10000];   //[npart]
   Double_t        part_xEnd[10000];   //[npart]
   Double_t        part_yEnd[10000];   //[npart]
   Double_t        part_zEnd[10000];   //[npart]
   Double_t        part_tEnd[10000];   //[npart]
   Double_t        part_pxStart[10000];   //[npart]
   Double_t        part_pyStart[10000];   //[npart]
   Double_t        part_pzStart[10000];   //[npart]
   Double_t        part_pxEnd[10000];   //[npart]
   Double_t        part_pyEnd[10000];   //[npart]
   Double_t        part_pzEnd[10000];   //[npart]
   Double_t        part_KEstart[10000];   //[npart]
   Double_t        part_KEend[10000];   //[npart]
   Int_t           part_processStart[10000];   //[npart]
   Int_t           part_processEnd[10000];   //[npart]
   Int_t           part_parentid[10000];   //[npart]
   Int_t           part_trackid[10000];   //[npart]
   Int_t           part_pid[10000];   //[npart]
   Double_t        capt_x[2];   //[ncapturecount]
   Double_t        capt_y[2];   //[ncapturecount]
   Double_t        capt_z[2];   //[ncapturecount]
   Double_t        capt_t0[2];   //[ncapturecount]
   Double_t        capt_E[2];   //[ncapturecount]
   Int_t           capt_num[2];   //[ncapturecount]
   Int_t           capt_pid[2];   //[ncapturecount]
   Int_t           capt_nucleus[2];   //[ncapturecount]
   Int_t           capt_nphot[2];   //[ncapturecount]
   Int_t           capt_ngamma[2];   //[ncapturecount]

   // List of branches
   TBranch        *b_evt;   //!
   TBranch        *b_nphot;   //!
   TBranch        *b_npart;   //!
   TBranch        *b_ncapturecount;   //!
   TBranch        *b_neutroncount;   //!
   TBranch        *b_phot_xStart;   //!
   TBranch        *b_phot_yStart;   //!
   TBranch        *b_phot_zStart;   //!
   TBranch        *b_phot_tStart;   //!
   TBranch        *b_phot_xEnd;   //!
   TBranch        *b_phot_yEnd;   //!
   TBranch        *b_phot_zEnd;   //!
   TBranch        *b_phot_tEnd;   //!
   TBranch        *b_phot_wavelength;   //!
   TBranch        *b_phot_processStart;   //!
   TBranch        *b_phot_isScat;   //!
   TBranch        *b_phot_parentid;   //!
   TBranch        *b_phot_trackid;   //!
   TBranch        *b_phot_hit;   //!
   TBranch        *b_phot_capnum;   //!
   TBranch        *b_part_xStart;   //!
   TBranch        *b_part_yStart;   //!
   TBranch        *b_part_zStart;   //!
   TBranch        *b_part_tStart;   //!
   TBranch        *b_part_xEnd;   //!
   TBranch        *b_part_yEnd;   //!
   TBranch        *b_part_zEnd;   //!
   TBranch        *b_part_tEnd;   //!
   TBranch        *b_part_pxStart;   //!
   TBranch        *b_part_pyStart;   //!
   TBranch        *b_part_pzStart;   //!
   TBranch        *b_part_pxEnd;   //!
   TBranch        *b_part_pyEnd;   //!
   TBranch        *b_part_pzEnd;   //!
   TBranch        *b_part_KEstart;   //!
   TBranch        *b_part_KEend;   //!
   TBranch        *b_part_processStart;   //!
   TBranch        *b_part_processEnd;   //!
   TBranch        *b_part_parentid;   //!
   TBranch        *b_part_trackid;   //!
   TBranch        *b_part_pid;   //!
   TBranch        *b_capt_x;   //!
   TBranch        *b_capt_y;   //!
   TBranch        *b_capt_z;   //!
   TBranch        *b_capt_t0;   //!
   TBranch        *b_capt_E;   //!
   TBranch        *b_capt_num;   //!
   TBranch        *b_capt_pid;   //!
   TBranch        *b_capt_nucleus;   //!
   TBranch        *b_capt_nphot;   //!
   TBranch        *b_capt_ngamma;   //!

   LikelihoodGenerator(TTree *tree=0);
   LikelihoodGenerator(TString file);
   virtual ~LikelihoodGenerator();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(bool isElectron = false, int start=0, int end=0);
   virtual void     Loop2(int start=0, int end=0);
   virtual void     Loop3(int start=0, int end=0);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   ClassDef(LikelihoodGenerator,0)

};

#endif

