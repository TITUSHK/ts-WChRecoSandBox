//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 20 13:54:03 2015 by ROOT version 5.34/23
// from TTree Debug/Debug
// found on file: nu_numu_1000_reco_12in.root
//////////////////////////////////////////////////////////

#ifndef EnergyLookups_h
#define EnergyLookups_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class EnergyLookups {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TTree          *fChain2;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   Int_t           fCurrent2; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nSubEvts;
   Bool_t          isHighE[20];   //[nSubEvts]
   Int_t           ringPEs[20];   //[nSubEvts]
   Double_t        trueKE;
   Double_t        diffVtxX;
   Double_t        diffVtxY;
   Double_t        diffVtxZ;
   Int_t           mode;
   Double_t        recoToWall[20];   //[nSubEvts]
   Double_t        recoDWall[20];   //[nSubEvts]
   Double_t        recoDirX[20];   //[nSubEvts]
   Double_t        recoDirY[20];   //[nSubEvts]
   Double_t        recoDirZ[20];   //[nSubEvts]
   Int_t           nSubEvts2;
   Bool_t          isHighE2[20];   //[nSubEvts]
   Int_t           ringPEs2[20];   //[nSubEvts]
   Double_t        trueKE2;
   Double_t        diffVtxX2;
   Double_t        diffVtxY2;
   Double_t        diffVtxZ2;
   Int_t           mode2;
   Double_t        recoToWall2[20];   //[nSubEvts]
   Double_t        recoDWall2[20];   //[nSubEvts]
   Double_t        recoDirX2[20];   //[nSubEvts]
   Double_t        recoDirY2[20];   //[nSubEvts]
   Double_t        recoDirZ2[20];   //[nSubEvts]

   // List of branches
   TBranch        *b_nSubEvts;   //!
   TBranch        *b_isHighE;   //!
   TBranch        *b_ringPEs;   //!
   TBranch        *b_trueKE;   //!
   TBranch        *b_diffVtxX;   //!
   TBranch        *b_diffVtxY;   //!
   TBranch        *b_diffVtxZ;   //!
   TBranch        *b_mode;   //!
   TBranch        *b_recoToWall;   //!
   TBranch        *b_recoDWall;   //!
   TBranch        *b_recoDirX;   //!
   TBranch        *b_recoDirY;   //!
   TBranch        *b_recoDirZ;   //!

   TBranch        *b_nSubEvts2;   //!
   TBranch        *b_isHighE2;   //!
   TBranch        *b_ringPEs2;   //!
   TBranch        *b_trueKE2;   //!
   TBranch        *b_diffVtxX2;   //!
   TBranch        *b_diffVtxY2;   //!
   TBranch        *b_diffVtxZ2;   //!
   TBranch        *b_mode2;   //!
   TBranch        *b_recoToWall2;   //!
   TBranch        *b_recoDWall2;   //!
   TBranch        *b_recoDirX2;   //!
   TBranch        *b_recoDirY2;   //!
   TBranch        *b_recoDirZ2;   //!

   EnergyLookups(TTree *tree=0, TTree *tree2=0);
   virtual ~EnergyLookups();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Int_t    GetEntry2(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual Long64_t LoadTree2(Long64_t entry);
   virtual void     Init(TTree *tree, TTree *tree2);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     Show2(Long64_t entry = -1);
};

#endif

#ifdef EnergyLookups_cxx
EnergyLookups::EnergyLookups(TTree *tree,TTree *tree2) : fChain(0), fChain2(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TChain *f=new TChain("Final_Reconstruction");
      TChain *le=new TChain("Low_E");
      TChain *d=new TChain("Debug");
      TChain *f2=new TChain("Final_Reconstruction");
      TChain *le2=new TChain("Low_E");
      TChain *d2=new TChain("Debug");
      TString flavs[4] = {"numu","nue","antinumu","antinue"};
      TString hcs[4] = {"nu","antinu"};
      for(int k=0; k<2; k++){
         const char * h = hcs[k].Data();
         for (int j = 0; j < 4; j++) {
            const char *s = flavs[j].Data();
            for (int i = 1000; i < 1100; i++) {
               //Missing vector files for antinu mode:
               if (k == 1 && j == 0 && (i == 1096 || i==1037)) continue;
               if (k == 1 && j == 1 && (i == 1017 || i == 1053)) continue;
               if (k == 1 && j == 2 && i == 1078) continue;
               if (k == 1 && j == 3 && (i == 1018 || i == 1019 || i == 1053)) continue;
               char *file = Form("/data/hyperk/wchsandbox_reco/flav_%s/%s_%s_%i/%s_%s_%i_lookupreco_12in.root", s, h, s, i, h, s,i);
               if(j%2==1){
                  f->AddFile(file);
                  le->AddFile(file);
                  d->AddFile(file);
               }
               else{
                  f2->AddFile(file);
                  le2->AddFile(file);
                  d2->AddFile(file);
               }
            }
         }
      }
      f->AddFriend(le);
      f->AddFriend(f);
      f->AddFriend(d);
      f2->AddFriend(le2);
      f2->AddFriend(f2);
      f2->AddFriend(d2);
      tree = f;
      tree2 = f2;
   }
   Init(tree,tree2);
}

EnergyLookups::~EnergyLookups()
{
   if (fChain)
      delete fChain->GetCurrentFile();
   if (fChain2)
      delete fChain2->GetCurrentFile();
}

Int_t EnergyLookups::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Int_t EnergyLookups::GetEntry2(Long64_t entry)
{
// Read contents of entry.
   if (!fChain2) return 0;
   return fChain2->GetEntry(entry);
}

Long64_t EnergyLookups::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

Long64_t EnergyLookups::LoadTree2(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain2) return -5;
   Long64_t centry = fChain2->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain2->GetTreeNumber() != fCurrent2) {
      fCurrent2 = fChain2->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EnergyLookups::Init(TTree *tree, TTree *tree2)
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
   if (!tree2) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   fChain2 = tree2;
   fCurrent2 = -1;
   fChain2->SetMakeClass(1);

   fChain->SetBranchAddress("nSubevents", &nSubEvts, &b_nSubEvts);
   fChain->SetBranchAddress("isHighE", isHighE, &b_isHighE);
   fChain->SetBranchAddress("ringPEs", ringPEs, &b_ringPEs);
   fChain->SetBranchAddress("trueKE", &trueKE, &b_trueKE);
   fChain->SetBranchAddress("diffVtxX", &diffVtxX, &b_diffVtxX);
   fChain->SetBranchAddress("diffVtxY", &diffVtxY, &b_diffVtxY);
   fChain->SetBranchAddress("diffVtxZ", &diffVtxZ, &b_diffVtxZ);
   fChain->SetBranchAddress("mode", &mode, &b_mode);
   fChain->SetBranchAddress("recoToWall", recoToWall, &b_recoToWall);
   fChain->SetBranchAddress("recoDWall", recoDWall, &b_recoDWall);
   fChain->SetBranchAddress("recoDirX", recoDirX, &b_recoDirX);
   fChain->SetBranchAddress("recoDirY", recoDirY, &b_recoDirY);
   fChain->SetBranchAddress("recoDirZ", recoDirZ, &b_recoDirZ);

   fChain2->SetBranchAddress("nSubevents", &nSubEvts2, &b_nSubEvts2);
   fChain2->SetBranchAddress("isHighE", isHighE2, &b_isHighE2);
   fChain2->SetBranchAddress("ringPEs", ringPEs2, &b_ringPEs2);
   fChain2->SetBranchAddress("trueKE", &trueKE2, &b_trueKE2);
   fChain2->SetBranchAddress("diffVtxX", &diffVtxX2, &b_diffVtxX2);
   fChain2->SetBranchAddress("diffVtxY", &diffVtxY2, &b_diffVtxY2);
   fChain2->SetBranchAddress("diffVtxZ", &diffVtxZ2, &b_diffVtxZ2);
   fChain2->SetBranchAddress("mode", &mode2, &b_mode2);
   fChain2->SetBranchAddress("recoToWall", recoToWall2, &b_recoToWall2);
   fChain2->SetBranchAddress("recoDWall", recoDWall2, &b_recoDWall2);
   fChain2->SetBranchAddress("recoDirX", recoDirX2, &b_recoDirX2);
   fChain2->SetBranchAddress("recoDirY", recoDirY2, &b_recoDirY2);
   fChain2->SetBranchAddress("recoDirZ", recoDirZ2, &b_recoDirZ2);

   Notify();
}

Bool_t EnergyLookups::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EnergyLookups::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

void EnergyLookups::Show2(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain2) return;
   fChain2->Show(entry);
}

Int_t EnergyLookups::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EnergyLookups_cxx
