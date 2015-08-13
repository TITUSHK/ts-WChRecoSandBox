#define EnergyLookups_cxx
#include "EnergyLookups.h"
#include <TH2D.h>
#include <TMath.h>
#include <iostream>

void EnergyLookups::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L EnergyLookups.C
//      Root > EnergyLookups t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the fChain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   TFile * out = new TFile("energyLookups_new2.root","RECREATE");
   TH2D * nHitToWallKELookup_e     = new TH2D("nHitToWallKELookup_e"    ,"nHitToWallKELookup_e"    ,30,0,20000,30,0,2000);
   TH2D * nHitToWallVtxTrackBias_e = new TH2D("nHitToWallVtxTrackBias_e","nHitToWallVtxTrackBias_e",30,0,20000,30,0,2000);
   TH2D * nHitToWallCounts_e       = new TH2D("nHitToWallCounts_e"      ,"nHitToWallCounts_e"      ,30,0,20000,30,0,2000);
   TH2D * nHitDWallKELookup_e      = new TH2D("nHitDWallKELookup_e"     ,"nHitDWallKELookup_e"     ,30,0,20000,30,0,550);
   TH2D * nHitDWallVtxTrackBias_e  = new TH2D("nHitDWallVtxTrackBias_e" ,"nHitDWallVtxTrackBias_e" ,30,0,20000,30,0,550);
   TH2D * nHitDWallCounts_e        = new TH2D("nHitDWallCounts_e"       ,"nHitDWallCounts_e"       ,30,0,20000,30,0,550);
   TH2D * nHitToWallKELookup_mu     = new TH2D("nHitToWallKELookup_mu"    ,"nHitToWallKELookup_mu"    ,30,0,20000,30,0,2000);
   TH2D * nHitToWallVtxTrackBias_mu = new TH2D("nHitToWallVtxTrackBias_mu","nHitToWallVtxTrackBias_mu",30,0,20000,30,0,2000);
   TH2D * nHitToWallCounts_mu       = new TH2D("nHitToWallCounts_mu"      ,"nHitToWallCounts_mu"      ,30,0,20000,30,0,2000);
   TH2D * nHitDWallKELookup_mu      = new TH2D("nHitDWallKELookup_mu"     ,"nHitDWallKELookup_mu"     ,30,0,20000,30,0,550);
   TH2D * nHitDWallVtxTrackBias_mu  = new TH2D("nHitDWallVtxTrackBias_mu" ,"nHitDWallVtxTrackBias_mu" ,30,0,20000,30,0,550);
   TH2D * nHitDWallCounts_mu        = new TH2D("nHitDWallCounts_mu"       ,"nHitDWallCounts_mu"       ,30,0,20000,30,0,550);
   Long64_t nentries = fChain->GetEntries();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(mode != 1 && mode !=0) continue;
      if(nSubEvts<1) continue;
      if(!isHighE[0]) continue;
      // if (Cut(ientry) < 0) continue;
      double vtxTrackBias = diffVtxX*recoDirX[0]+diffVtxY*recoDirY[0]+diffVtxZ*recoDirZ[0];
      nHitToWallKELookup_e->Fill(ringPEs[0], recoToWall[0], trueKE);
      nHitToWallVtxTrackBias_e->Fill(ringPEs[0], recoToWall[0], vtxTrackBias);
      nHitToWallCounts_e->Fill(ringPEs[0], recoToWall[0]);
      nHitDWallKELookup_e->Fill(ringPEs[0], recoDWall[0], trueKE);
      nHitDWallVtxTrackBias_e->Fill(ringPEs[0], recoDWall[0], vtxTrackBias);
      nHitDWallCounts_e->Fill(ringPEs[0], recoDWall[0]);
   }
   Long64_t nentries2 = fChain2->GetEntries();
   Long64_t nbytes2 = 0, nb2 = 0;
   for (Long64_t jentry=0; jentry<nentries2;jentry++) {
      Long64_t ientry = LoadTree2(jentry);
      if (ientry < 0) break;
      nb2 = fChain2->GetEntry(jentry);   nbytes2 += nb2;
      if(mode2 != 1 && mode2 !=0) continue;
      if(nSubEvts2<1) continue;
      if(!isHighE2[0]) continue;
      // if (Cut(ientry) < 0) continue;
      double vtxTrackBias = diffVtxX2*recoDirX2[0]+diffVtxY2*recoDirY2[0]+diffVtxZ2*recoDirZ2[0];
      nHitToWallKELookup_mu->Fill(ringPEs2[0], recoToWall2[0], trueKE2);
      nHitToWallVtxTrackBias_mu->Fill(ringPEs2[0], recoToWall2[0], vtxTrackBias);
      nHitToWallCounts_mu->Fill(ringPEs2[0], recoToWall2[0]);
      nHitDWallKELookup_mu->Fill(ringPEs2[0], recoDWall2[0], trueKE2);
      nHitDWallVtxTrackBias_mu->Fill(ringPEs2[0], recoDWall2[0], vtxTrackBias);
      nHitDWallCounts_mu->Fill(ringPEs2[0], recoDWall2[0]);
   }
   nHitToWallKELookup_e->Divide(nHitToWallCounts_e);
   nHitToWallVtxTrackBias_e->Divide(nHitToWallCounts_e);
   nHitDWallKELookup_e->Divide(nHitDWallCounts_e);
   nHitDWallVtxTrackBias_e->Divide(nHitDWallCounts_e);
   nHitToWallKELookup_e->Write();
   nHitToWallVtxTrackBias_e->Write();
   nHitDWallKELookup_e->Write();
   nHitDWallVtxTrackBias_e->Write();
   nHitToWallKELookup_mu->Divide(nHitToWallCounts_mu);
   nHitToWallVtxTrackBias_mu->Divide(nHitToWallCounts_mu);
   nHitDWallKELookup_mu->Divide(nHitDWallCounts_mu);
   nHitDWallVtxTrackBias_mu->Divide(nHitDWallCounts_mu);
   nHitToWallKELookup_mu->Write();
   nHitToWallVtxTrackBias_mu->Write();
   nHitDWallKELookup_mu->Write();
   nHitDWallVtxTrackBias_mu->Write();
   out->Close();
}
