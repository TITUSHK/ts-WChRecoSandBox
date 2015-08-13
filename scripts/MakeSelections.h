//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 20 13:54:03 2015 by ROOT version 5.34/23
// from TTree Debug/Debug
// found on file: nu_numu_1000_reco_12in.root
//////////////////////////////////////////////////////////

#ifndef MakeSelections_h
#define MakeSelections_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class MakeSelections {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TTree          *nRooTracker;
   Int_t           fCurrent; //!current Tree number in a TChain
   Int_t           fCurrent2; //!current Tree number in a TChain

   Int_t           flav_count[4];
   bool            IsAntiNu;

   // Declaration of leaf types
   Int_t           evt;
   Int_t           nClusters;
   Int_t nSubevents;
   Int_t           cluster[20];   //[nSubevents]
   Int_t           ring[20];   //[nSubevents]
   Bool_t          isHighE[20];   //[nSubevents]
   Int_t           ringPEs[20];   //[nSubevents]
   Double_t        trueVtxX;
   Double_t        trueVtxY;
   Double_t        trueVtxZ;
   Double_t        trueTime;
   Double_t        trueDirX;
   Double_t        trueDirY;
   Double_t        trueDirZ;
   Double_t        trueKE;
   Double_t        diffKE;
   Double_t        diffVtxX;
   Double_t        diffVtxY;
   Double_t        diffVtxZ;
   Double_t        diffTime;
   Double_t        diffDirX;
   Double_t        diffDirY;
   Double_t        diffDirZ;
   Double_t        diffVtxAbs;
   Double_t        diffDirAbs;
   Int_t           mode;
   Double_t        neutrinoE;
   Double_t        neutrinoDirX;
   Double_t        neutrinoDirY;
   Double_t        neutrinoDirZ;
   Int_t           neutrinoPID;
   Int_t           nNeutrons;
   Int_t           neutronCount;
   Int_t           nCaptures;

   Double_t        recoVtxXLowE[20];   //[nSubevents]
   Double_t        recoVtxYLowE[20];   //[nSubevents]
   Double_t        recoVtxZLowE[20];   //[nSubevents]
   Double_t        recoTimeLowE[20];   //[nSubevents]
   Double_t        recoDirXLowE[20];   //[nSubevents]
   Double_t        recoDirYLowE[20];   //[nSubevents]
   Double_t        recoDirZLowE[20];   //[nSubevents]
   Double_t        recoChkvAngleLowE[20];   //[nSubevents]
   Double_t        recoEnergyLowE[20];   //[nSubevents]

   Double_t        recoVtxXHighEMuon[20];   //[nSubevents]
   Double_t        recoVtxYHighEMuon[20];   //[nSubevents]
   Double_t        recoVtxZHighEMuon[20];   //[nSubevents]
   Double_t        recoTimeHighEMuon[20];   //[nSubevents]
   Double_t        recoDirXHighEMuon[20];   //[nSubevents]
   Double_t        recoDirYHighEMuon[20];   //[nSubevents]
   Double_t        recoDirZHighEMuon[20];   //[nSubevents]
   Double_t        recoChkvAngleHighEMuon[20];   //[nSubevents]
   Double_t        recoEnergyHighEMuon[20];   //[nSubevents]
   Double_t        recoLnLHighEMuon[20];   //[nSubevents]

   Double_t        recoVtxXHighEElectron[20];   //[nSubevents]
   Double_t        recoVtxYHighEElectron[20];   //[nSubevents]
   Double_t        recoVtxZHighEElectron[20];   //[nSubevents]
   Double_t        recoTimeHighEElectron[20];   //[nSubevents]
   Double_t        recoDirXHighEElectron[20];   //[nSubevents]
   Double_t        recoDirYHighEElectron[20];   //[nSubevents]
   Double_t        recoDirZHighEElectron[20];   //[nSubevents]
   Double_t        recoChkvAngleHighEElectron[20];   //[nSubevents]
   Double_t        recoEnergyHighEElectron[20];   //[nSubevents]
   Double_t        recoLnLHighEElectron[20];   //[nSubevents]

   Double_t        recoVtxX[20];   //[nSubevents]
   Double_t        recoVtxY[20];   //[nSubevents]
   Double_t        recoVtxZ[20];   //[nSubevents]
   Double_t        recoTime[20];   //[nSubevents]
   Double_t        recoDirX[20];   //[nSubevents]
   Double_t        recoDirY[20];   //[nSubevents]
   Double_t        recoDirZ[20];   //[nSubevents]
   Double_t        recoChkvAngle[20];   //[nSubevents]
   Double_t        recoEnergy[20];   //[nSubevents]
   Int_t           recoPID[20];   //[nSubevents]
   Int_t           recoNRings[20];   //[nSubevents]
/*
   Double_t        neutrino_E;
   Int_t           neutrino_id;
   Double_t        neutrino_px;
   Double_t        neutrino_py;
   Double_t        neutrino_pz;
   Int_t           ntrks;
   Int_t           nneutrons;
   Double_t        vtxx;
   Double_t        vtxy;
   Double_t        vtxz;
   Int_t           mpid[100];   //[ntrks]
   Double_t        px[100];   //[ntrks]
   Double_t        py[100];   //[ntrks]
   Double_t        pz[100];   //[ntrks]
   Double_t        KE[100];   //[ntrks]
*/
   Int_t           NEneutmode;
   Int_t           NEnvc;
   Int_t           NEipvc[100];   //[NEnvc]
   Float_t         NEpvc[100][3];   //[NEnvc]
   Int_t           NEiorgvc[100];   //[NEnvc]
   Int_t           NEiflgvc[100];   //[NEnvc]
   Int_t           NEicrnvc[100];   //[NEnvc]
   Float_t         NEcrsx;
   Float_t         NEcrsy;
   Float_t         NEcrsz;
   Float_t         NEcrsphi;
   Int_t           NEnvert;
   Float_t         NEposvert[100][3];   //[NEnvert]
   Int_t           NEiflgvert[100];   //[NEnvert]
   Int_t           NEnvcvert;
   Float_t         NEdirvert[100][3];   //[NEnvcvert]
   Float_t         NEabspvert[100];   //[NEnvcvert]
   Float_t         NEabstpvert[100];   //[NEnvcvert]
   Int_t           NEipvert[100];   //[NEnvcvert]
   Int_t           NEiverti[100];   //[NEnvcvert]
   Int_t           NEivertf[100];   //[NEnvcvert]
   Int_t           StdHepPdg[24];   //[StdHepN]

   // List of branches
   TBranch        *b_evt;   //!
   TBranch        *b_nClusters;   //!
   TBranch        *b_nSubevents;   //!
   TBranch        *b_cluster;   //!
   TBranch        *b_ring;   //!
   TBranch        *b_isHighE;   //!
   TBranch        *b_ringPEs;   //!
   TBranch        *b_trueVtxX;   //!
   TBranch        *b_trueVtxY;   //!
   TBranch        *b_trueVtxZ;   //!
   TBranch        *b_trueTime;   //!
   TBranch        *b_trueDirX;   //!
   TBranch        *b_trueDirY;   //!
   TBranch        *b_trueDirZ;   //!
   TBranch        *b_trueKE;   //!
   TBranch        *b_diffKE;   //!
   TBranch        *b_diffVtxX;   //!
   TBranch        *b_diffVtxY;   //!
   TBranch        *b_diffVtxZ;   //!
   TBranch        *b_diffTime;   //!
   TBranch        *b_diffDirX;   //!
   TBranch        *b_diffDirY;   //!
   TBranch        *b_diffDirZ;   //!
   TBranch        *b_diffVtxAbs;   //!
   TBranch        *b_diffDirAbs;   //!
   TBranch        *b_mode;   //!
   TBranch        *b_neutrinoE;   //!
   TBranch        *b_neutrinoDirX;   //!
   TBranch        *b_neutrinoDirY;   //!
   TBranch        *b_neutrinoDirZ;   //!
   TBranch        *b_neutrinoPID;   //!
   TBranch        *b_nNeutrons;   //!
   TBranch        *b_neutronCount;   //!
   TBranch        *b_nCaptures;   //!

   TBranch        *b_recoVtxXLowE;   //!
   TBranch        *b_recoVtxYLowE;   //!
   TBranch        *b_recoVtxZLowE;   //!
   TBranch        *b_recoTimeLowE;   //!
   TBranch        *b_recoDirXLowE;   //!
   TBranch        *b_recoDirYLowE;   //!
   TBranch        *b_recoDirZLowE;   //!
   TBranch        *b_recoChkvAngleLowE;   //!
   TBranch        *b_recoEnergyLowE;   //!

   TBranch        *b_recoVtxXHighEMuon;   //!
   TBranch        *b_recoVtxYHighEMuon;   //!
   TBranch        *b_recoVtxZHighEMuon;   //!
   TBranch        *b_recoTimeHighEMuon;   //!
   TBranch        *b_recoDirXHighEMuon;   //!
   TBranch        *b_recoDirYHighEMuon;   //!
   TBranch        *b_recoDirZHighEMuon;   //!
   TBranch        *b_recoChkvAngleHighEMuon;   //!
   TBranch        *b_recoEnergyHighEMuon;   //!
   TBranch        *b_recoLnLHighEMuon;   //!

   TBranch        *b_recoVtxXHighEElectron;   //!
   TBranch        *b_recoVtxYHighEElectron;   //!
   TBranch        *b_recoVtxZHighEElectron;   //!
   TBranch        *b_recoTimeHighEElectron;   //!
   TBranch        *b_recoDirXHighEElectron;   //!
   TBranch        *b_recoDirYHighEElectron;   //!
   TBranch        *b_recoDirZHighEElectron;   //!
   TBranch        *b_recoChkvAngleHighEElectron;   //!
   TBranch        *b_recoEnergyHighEElectron;   //!
   TBranch        *b_recoLnLHighEElectron;   //!

   TBranch        *b_recoVtxX;   //!
   TBranch        *b_recoVtxY;   //!
   TBranch        *b_recoVtxZ;   //!
   TBranch        *b_recoTime;   //!
   TBranch        *b_recoDirX;   //!
   TBranch        *b_recoDirY;   //!
   TBranch        *b_recoDirZ;   //!
   TBranch        *b_recoChkvAngle;   //!
   TBranch        *b_recoEnergy;   //!
   TBranch        *b_recoPID;   //!
   TBranch        *b_recoNRings;   //!
/*
   TBranch        *b_neutrino_E;   //!
   TBranch        *b_neutrino_id;   //!
   TBranch        *b_neutrino_px;   //!
   TBranch        *b_neutrino_py;   //!
   TBranch        *b_neutrino_pz;   //!
   TBranch        *b_ntrks;   //!
   TBranch        *b_nneutrons;   //!
   TBranch        *b_vtxx;   //!
   TBranch        *b_vtxy;   //!
   TBranch        *b_vtxz;   //!
   TBranch        *b_mpid;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_KE;   //!
*/
   TBranch        *b_NEneutmode;   //!
   TBranch        *b_NEnvc;   //!
   TBranch        *b_NEipvc;   //!
   TBranch        *b_NEpvc;   //!
   TBranch        *b_NEiorgvc;   //!
   TBranch        *b_NEiflgvc;   //!
   TBranch        *b_NEicrnvc;   //!
   TBranch        *b_NEcrsx;   //!
   TBranch        *b_NEcrsy;   //!
   TBranch        *b_NEcrsz;   //!
   TBranch        *b_NEcrsphi;   //!
   TBranch        *b_NEnvert;   //!
   TBranch        *b_NEposvert;   //!
   TBranch        *b_NEiflgvert;   //!
   TBranch        *b_NEnvcvert;   //!
   TBranch        *b_NEdirvert;   //!
   TBranch        *b_NEabspvert;   //!
   TBranch        *b_NEabstpvert;   //!
   TBranch        *b_NEipvert;   //!
   TBranch        *b_NEiverti;   //!
   TBranch        *b_NEivertf;   //!
   TBranch        *b_StdHepPdg;   //!

   MakeSelections(bool isAntiNu=false, TTree *tree=0, TTree * tracker=0);
   virtual ~MakeSelections();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, TTree *tracker);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   TTree *MakeNtuple(TString name, float Abspvc[100], int &Numatom, double &Erec, double &Etrue, Float_t pnu[100],
                     Float_t dirnu[100][3], Int_t &Numbndn, Int_t &Numbndp, Int_t &Numfrep, Int_t &Ibound, Int_t &flavEvents);
};
#endif

#ifdef MakeSelections_cxx
MakeSelections::MakeSelections(bool isAntiNu, TTree *tree, TTree *tracker) : fChain(0), nRooTracker(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   IsAntiNu = isAntiNu;
   if (tree == 0) {
      TChain *f=new TChain("Final_Reconstruction");
      TChain *le=new TChain("Low_E");
      TChain *hee=new TChain("High_E_Electron");
      TChain *hem=new TChain("High_E_Muon");
      TChain *d=new TChain("Debug");
      TChain *t=new TChain("tcardfile");
      TChain *tr = new TChain("nRooTracker");
      TString flavs[4] = {"numu","nue","antinumu","antinue"};
      const char * a = (isAntiNu) ? "antinu" :"nu";
      for(int j=0; j<4; j++) {
         flav_count[j]=0;
         const char *s = flavs[j].Data();
         for (int i = 1000; i < 1100; i++) {
            //Missing vector files for antinu mode:
            if(isAntiNu && j==0 && (i==1096 || i==1037)) continue;
            if(isAntiNu && j==1 && (i==1017 || i==1053)) continue;
            if(isAntiNu && j==2 && i==1078) continue;
            if(isAntiNu && j==3 && (i==1018 || i==1019 || i==1053)) continue;
            char *file = Form("/data/hyperk/wchsandbox_reco/flav_%s/%s_%s_%i/%s_%s_%i_12in_reco.root", s, a, s, i, a, s, i);
            f->AddFile(file,1000);
            le->AddFile(file,1000);
            hee->AddFile(file,1000);
            hem->AddFile(file,1000);
            d->AddFile(file,1000);
            t->AddFile(Form("/data/hyperk/wchsandbox_reco/flav_%s/%s_%s_%i/%s_%s_%i_generatorcardfile.root",s,a,s,i,a,s,i),1000);
            tr->AddFile(Form("/data/hyperk/wchsandbox_reco/vectors/v00-01/flav_%s/genev_%s_cylinder_r551_z2200_Z_%s_1721827_%i.root",s,a,s,i),1000);
            flav_count[j]+=1000;
         }
      }
      f->AddFriend(le);
      f->AddFriend(hem);
      f->AddFriend(hee);
      f->AddFriend(f);
      f->AddFriend(d);
      tree = f;
      tracker = tr;
   }
   tree->AddFriend(tracker);
   Init(tree,tracker);
}

MakeSelections::~MakeSelections()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   if (!nRooTracker) return;
   delete nRooTracker->GetCurrentFile();
}

Int_t MakeSelections::GetEntry(Long64_t entry)
{
// Read contents of entry.
   int nb = 0;
   if (!fChain) return 0;
   nb += fChain->GetEntry(entry);
   if (!nRooTracker) return 0;
   nb += nRooTracker->GetEntry(entry);
   return nb;
}
Long64_t MakeSelections::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   if (!nRooTracker) return -5;
   Long64_t centry2 = nRooTracker->LoadTree(entry);
   if (centry2 < 0) return centry2;
   if (nRooTracker->GetTreeNumber() != fCurrent2) {
      fCurrent2 = nRooTracker->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MakeSelections::Init(TTree *tree, TTree *tracker)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or fChain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   if (!tracker) return;
   fChain = tree;
   nRooTracker = tracker;
   fCurrent = -1;
   fCurrent2 = -1;
   fChain->SetMakeClass(1);
   nRooTracker->SetMakeClass(1);

   fChain->SetBranchAddress("evt", &evt, &b_evt);
//   fChain->SetBranchAddress("nClusters", &nClusters, &b_nClusters);
//   fChain->SetBranchAddress("nSubevents", &nSubevents, &b_nSubevents);
//   fChain->SetBranchAddress("cluster", cluster, &b_cluster);

   fChain->SetBranchAddress("ring", ring, &b_ring);
   fChain->SetBranchAddress("isHighE", isHighE, &b_isHighE);
//   fChain->SetBranchAddress("ringPEs", ringPEs, &b_ringPEs);
   fChain->SetBranchAddress("trueVtxX", &trueVtxX, &b_trueVtxX);
   fChain->SetBranchAddress("trueVtxY", &trueVtxY, &b_trueVtxY);
   fChain->SetBranchAddress("trueVtxZ", &trueVtxZ, &b_trueVtxZ);
   fChain->SetBranchAddress("trueTime", &trueTime, &b_trueTime);
   fChain->SetBranchAddress("trueDirX", &trueDirX, &b_trueDirX);
   fChain->SetBranchAddress("trueDirY", &trueDirY, &b_trueDirY);
   fChain->SetBranchAddress("trueDirZ", &trueDirZ, &b_trueDirZ);
   fChain->SetBranchAddress("trueKE", &trueKE, &b_trueKE);
   fChain->SetBranchAddress("diffKE", &diffKE, &b_diffKE);
   fChain->SetBranchAddress("diffVtxX", &diffVtxX, &b_diffVtxX);
   fChain->SetBranchAddress("diffVtxY", &diffVtxY, &b_diffVtxY);
   fChain->SetBranchAddress("diffVtxZ", &diffVtxZ, &b_diffVtxZ);
   fChain->SetBranchAddress("diffTime", &diffTime, &b_diffTime);
   fChain->SetBranchAddress("diffDirX", &diffDirX, &b_diffDirX);
   fChain->SetBranchAddress("diffDirY", &diffDirY, &b_diffDirY);
   fChain->SetBranchAddress("diffDirZ", &diffDirZ, &b_diffDirZ);
   fChain->SetBranchAddress("diffVtxAbs", &diffVtxAbs, &b_diffVtxAbs);
   fChain->SetBranchAddress("diffDirAbs", &diffDirAbs, &b_diffDirAbs);
   fChain->SetBranchAddress("mode", &mode, &b_mode);
   fChain->SetBranchAddress("neutrinoE", &neutrinoE, &b_neutrinoE);
   fChain->SetBranchAddress("neutrinoDirX", &neutrinoDirX, &b_neutrinoDirX);
   fChain->SetBranchAddress("neutrinoDirY", &neutrinoDirY, &b_neutrinoDirY);
   fChain->SetBranchAddress("neutrinoDirZ", &neutrinoDirZ, &b_neutrinoDirZ);
   fChain->SetBranchAddress("neutrinoPID", &neutrinoPID, &b_neutrinoPID);
   fChain->SetBranchAddress("nNeutrons", &nNeutrons, &b_nNeutrons);
   fChain->SetBranchAddress("neutronCount", &neutronCount, &b_neutronCount);
   fChain->SetBranchAddress("nCaptures", &nCaptures, &b_nCaptures);

   fChain->SetBranchAddress("recoVtxXLowE", recoVtxXLowE, &b_recoVtxXLowE);
   fChain->SetBranchAddress("recoVtxYLowE", recoVtxYLowE, &b_recoVtxYLowE);
   fChain->SetBranchAddress("recoVtxZLowE", recoVtxZLowE, &b_recoVtxZLowE);
   fChain->SetBranchAddress("recoTimeLowE", recoTimeLowE, &b_recoTimeLowE);
   fChain->SetBranchAddress("recoDirXLowE", recoDirXLowE, &b_recoDirXLowE);
   fChain->SetBranchAddress("recoDirYLowE", recoDirYLowE, &b_recoDirYLowE);
   fChain->SetBranchAddress("recoDirZLowE", recoDirZLowE, &b_recoDirZLowE);
   fChain->SetBranchAddress("recoChkvAngleLowE", recoChkvAngleLowE, &b_recoChkvAngleLowE);
   fChain->SetBranchAddress("recoEnergyLowE", recoEnergyLowE, &b_recoEnergyLowE);

   fChain->SetBranchAddress("recoVtxXHighEMuon", recoVtxXHighEMuon, &b_recoVtxXHighEMuon);
   fChain->SetBranchAddress("recoVtxYHighEMuon", recoVtxYHighEMuon, &b_recoVtxYHighEMuon);
   fChain->SetBranchAddress("recoVtxZHighEMuon", recoVtxZHighEMuon, &b_recoVtxZHighEMuon);
   fChain->SetBranchAddress("recoTimeHighEMuon", recoTimeHighEMuon, &b_recoTimeHighEMuon);
   fChain->SetBranchAddress("recoDirXHighEMuon", recoDirXHighEMuon, &b_recoDirXHighEMuon);
   fChain->SetBranchAddress("recoDirYHighEMuon", recoDirYHighEMuon, &b_recoDirYHighEMuon);
   fChain->SetBranchAddress("recoDirZHighEMuon", recoDirZHighEMuon, &b_recoDirZHighEMuon);
   fChain->SetBranchAddress("recoChkvAngleHighEMuon", recoChkvAngleHighEMuon, &b_recoChkvAngleHighEMuon);
   fChain->SetBranchAddress("recoEnergyHighEMuon", recoEnergyHighEMuon, &b_recoEnergyHighEMuon);
   fChain->SetBranchAddress("recoLnLHighEMuon", recoLnLHighEMuon, &b_recoLnLHighEMuon);

   fChain->SetBranchAddress("recoVtxXHighEElectron", recoVtxXHighEElectron, &b_recoVtxXHighEElectron);
   fChain->SetBranchAddress("recoVtxYHighEElectron", recoVtxYHighEElectron, &b_recoVtxYHighEElectron);
   fChain->SetBranchAddress("recoVtxZHighEElectron", recoVtxZHighEElectron, &b_recoVtxZHighEElectron);
   fChain->SetBranchAddress("recoTimeHighEElectron", recoTimeHighEElectron, &b_recoTimeHighEElectron);
   fChain->SetBranchAddress("recoDirXHighEElectron", recoDirXHighEElectron, &b_recoDirXHighEElectron);
   fChain->SetBranchAddress("recoDirYHighEElectron", recoDirYHighEElectron, &b_recoDirYHighEElectron);
   fChain->SetBranchAddress("recoDirZHighEElectron", recoDirZHighEElectron, &b_recoDirZHighEElectron);
   fChain->SetBranchAddress("recoChkvAngleHighEElectron", recoChkvAngleHighEElectron, &b_recoChkvAngleHighEElectron);
   fChain->SetBranchAddress("recoEnergyHighEElectron", recoEnergyHighEElectron, &b_recoEnergyHighEElectron);
   fChain->SetBranchAddress("recoLnLHighEElectron", recoLnLHighEElectron, &b_recoLnLHighEElectron);

   fChain->SetBranchAddress("recoVtxX", recoVtxX, &b_recoVtxX);
   fChain->SetBranchAddress("recoVtxY", recoVtxY, &b_recoVtxY);
   fChain->SetBranchAddress("recoVtxZ", recoVtxZ, &b_recoVtxZ);
   fChain->SetBranchAddress("recoTime", recoTime, &b_recoTime);
   fChain->SetBranchAddress("recoDirX", recoDirX, &b_recoDirX);
   fChain->SetBranchAddress("recoDirY", recoDirY, &b_recoDirY);
   fChain->SetBranchAddress("recoDirZ", recoDirZ, &b_recoDirZ);
   fChain->SetBranchAddress("recoChkvAngle", recoChkvAngle, &b_recoChkvAngle);
   fChain->SetBranchAddress("recoEnergy", recoEnergy, &b_recoEnergy);
   fChain->SetBranchAddress("recoPID", recoPID, &b_recoPID);
   fChain->SetBranchAddress("recoNRings", recoNRings, &b_recoNRings);
/*
   fCard->SetBranchAddress("neutrino_E", &neutrino_E, &b_neutrino_E);
   fCard->SetBranchAddress("neutrino_id", &neutrino_id, &b_neutrino_id);
   fCard->SetBranchAddress("neutrino_px", &neutrino_px, &b_neutrino_px);
   fCard->SetBranchAddress("neutrino_py", &neutrino_py, &b_neutrino_py);
   fCard->SetBranchAddress("neutrino_pz", &neutrino_pz, &b_neutrino_pz);
   fCard->SetBranchAddress("ntrks", &ntrks, &b_ntrks);
   fCard->SetBranchAddress("nneutrons", &nneutrons, &b_nneutrons);
   fCard->SetBranchAddress("vtxx", &vtxx, &b_vtxx);
   fCard->SetBranchAddress("vtxy", &vtxy, &b_vtxy);
   fCard->SetBranchAddress("vtxz", &vtxz, &b_vtxz);
   fCard->SetBranchAddress("mpid", &mpid, &b_mpid);
   fCard->SetBranchAddress("px", &px, &b_px);
   fCard->SetBranchAddress("py", &py, &b_py);
   fCard->SetBranchAddress("pz", &pz, &b_pz);
   fCard->SetBranchAddress("KE", &KE, &b_KE);
*/
   nRooTracker->SetBranchAddress("NEneutmode", &NEneutmode, &b_NEneutmode);
   nRooTracker->SetBranchAddress("NEnvc", &NEnvc, &b_NEnvc);
   nRooTracker->SetBranchAddress("NEipvc", NEipvc, &b_NEipvc);
   nRooTracker->SetBranchAddress("NEpvc", NEpvc, &b_NEpvc);
   nRooTracker->SetBranchAddress("NEiorgvc", NEiorgvc, &b_NEiorgvc);
   nRooTracker->SetBranchAddress("NEiflgvc", NEiflgvc, &b_NEiflgvc);
   nRooTracker->SetBranchAddress("NEicrnvc", NEicrnvc, &b_NEicrnvc);
   nRooTracker->SetBranchAddress("NEcrsx", &NEcrsx, &b_NEcrsx);
   nRooTracker->SetBranchAddress("NEcrsy", &NEcrsy, &b_NEcrsy);
   nRooTracker->SetBranchAddress("NEcrsz", &NEcrsz, &b_NEcrsz);
   nRooTracker->SetBranchAddress("NEcrsphi", &NEcrsphi, &b_NEcrsphi);
   nRooTracker->SetBranchAddress("NEnvert", &NEnvert, &b_NEnvert);
   nRooTracker->SetBranchAddress("NEposvert", NEposvert, &b_NEposvert);
   nRooTracker->SetBranchAddress("NEiflgvert", NEiflgvert, &b_NEiflgvert);
   nRooTracker->SetBranchAddress("NEnvcvert", &NEnvcvert, &b_NEnvcvert);
   nRooTracker->SetBranchAddress("NEdirvert", NEdirvert, &b_NEdirvert);
   nRooTracker->SetBranchAddress("NEabspvert", NEabspvert, &b_NEabspvert);
   nRooTracker->SetBranchAddress("NEabstpvert", NEabstpvert, &b_NEabstpvert);
   nRooTracker->SetBranchAddress("NEipvert", NEipvert, &b_NEipvert);
   nRooTracker->SetBranchAddress("NEiverti", NEiverti, &b_NEiverti);
   nRooTracker->SetBranchAddress("NEivertf", NEivertf, &b_NEivertf);
   nRooTracker->SetBranchAddress("StdHepPdg", StdHepPdg, &b_StdHepPdg);

   Notify();
}

Bool_t MakeSelections::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MakeSelections::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
   if (!nRooTracker) return;
   nRooTracker->Show(entry);
}
Int_t MakeSelections::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MakeSelections_cxx
