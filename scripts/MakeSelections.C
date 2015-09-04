#define MakeSelections_cxx
#include "MakeSelections.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <iostream>
#include <TRandom3.h>

void MakeSelections::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L MakeSelections.C
//      Root > MakeSelections t
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

   float Abspvc[100];
   //int mode;
   int Numatom=1;
   double Erec, Etrue;

   Float_t         pnu[100];   //[numnu]  GeV/c
   Float_t         dirnu[100][3];   //[numnu]
   Int_t Numbndn;
   Int_t Numbndp;
   Int_t Numfrep;
   Int_t Ibound;

   const char * s = IsAntiNu ? "RHC" : "FHC";
   TFile *out = new TFile(Form("selections_tagged_%s.root",s),"RECREATE");
   TTree *numu_mulike_nontag = MakeNtuple("numu_mulike_nontag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[0]);
   TTree *nue_mulike_nontag = MakeNtuple("nue_mulike_nontag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[1]);
   TTree *numubar_mulike_nontag = MakeNtuple("numubar_mulike_nontag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[2]);
   TTree *nuebar_mulike_nontag = MakeNtuple("nuebar_mulike_nontag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[3]);
   TTree *numu_elike_nontag = MakeNtuple("numu_elike_nontag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[0]);
   TTree *nue_elike_nontag = MakeNtuple("nue_elike_nontag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[1]);
   TTree *numubar_elike_nontag = MakeNtuple("numubar_elike_nontag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[2]);
   TTree *nuebar_elike_nontag = MakeNtuple("nuebar_elike_nontag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[3]);
   TTree *numu_mulike_ntag = MakeNtuple("numu_mulike_ntag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[0]);
   TTree *nue_mulike_ntag = MakeNtuple("nue_mulike_ntag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[1]);
   TTree *numubar_mulike_ntag = MakeNtuple("numubar_mulike_ntag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[2]);
   TTree *nuebar_mulike_ntag = MakeNtuple("nuebar_mulike_ntag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[3]);
   TTree *numu_elike_ntag = MakeNtuple("numu_elike_ntag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[0]);
   TTree *nue_elike_ntag = MakeNtuple("nue_elike_ntag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[1]);
   TTree *numubar_elike_ntag = MakeNtuple("numubar_elike_ntag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[2]);
   TTree *nuebar_elike_ntag = MakeNtuple("nuebar_elike_ntag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[3]);

   Long64_t nentries = fChain->GetEntriesFast();

   TRandom3 rand(31415);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(!isHighE[0]) continue;
      if(recoNRings[0] != 1) continue;
      Int_t pid = (recoLnLHighEMuon[0] - recoLnLHighEElectron[0] > -200) ? 13 : 11;
      Double_t recoKE = pid==13 ? recoEnergyHighEMuon[0] : recoEnergyHighEElectron[0];
      if(recoKE < 100) continue;
      if(recoKE > 2500) continue;
      Double_t vtxZ = pid==13 ? recoVtxZHighEMuon[0] : recoVtxZHighEElectron[0];
      double dWallZ = 1100-TMath::Abs(vtxZ);
      if(dWallZ < 100) continue;
      Double_t vtxX = pid==13 ? recoVtxXHighEMuon[0] : recoVtxXHighEElectron[0];
      Double_t vtxY = pid==13 ? recoVtxYHighEMuon[0] : recoVtxYHighEElectron[0];
      double recoVtxR2 = vtxX * vtxX + vtxY * vtxY;
      double dWallR = 550-TMath::Sqrt(recoVtxR2);
      if(dWallR < 100) continue;
      
      double mn = 939;
      double ml = pid ==13 ? 105 : 0.511;
      double recoE = recoKE + ml;
      double recoP = TMath::Sqrt(recoE*recoE- ml * ml);
      Double_t dirZ = pid==13 ? recoDirZHighEMuon[0] : recoDirZHighEElectron[0];
      Erec = (mn*recoE- ml * ml /2.)/(mn-recoE+recoP* dirZ);
      if(Erec > 2500 || Erec < 0) continue;
      Etrue = TMath::Sqrt(NEpvc[0][0]*NEpvc[0][0]+NEpvc[0][1]*NEpvc[0][1]+NEpvc[0][2]*NEpvc[0][2]);
      for(int i=0; i<NEnvc; i++) {
         Abspvc[i] = TMath::Sqrt(NEpvc[i][0]*NEpvc[i][0] + NEpvc[i][1]*NEpvc[i][1] + NEpvc[i][2]*NEpvc[i][2]);
         pnu[i] = Abspvc[i] / 1000.;
         for (int j = 0; j < 3; j++) {
            dirnu[i][j] = NEpvc[i][j] / Abspvc[i];
         }
      }
      if(StdHepPdg[1]==1000080160){ //Oxygen target
         Numatom = 16;
         Numbndn = 8;
         Numbndp = 8;
         Numfrep = 0;
         Ibound = 1;
      }
      else if(StdHepPdg[1]==1000010010){ //Hydrogen target
         Numatom = 1;
         Numbndn = 0;
         Numbndp = 0;
         Numfrep = 1;
         Ibound = 0;
      }
      else{
         Numatom = 0;
         Numbndn = 0;
         Numbndp = 0;
         Numfrep = 0;
         Ibound = 0;
      }

      /*
      Npvc = ntrks+2;
      //First element is incoming neutrino
      Ichvc[0]=0;
      Ipvc[0]=neutrino_id;
      double neutrino_p = neutrino_px*neutrino_px+neutrino_py*neutrino_py+neutrino_pz*neutrino_pz;
      // neutrino_px etc. seem to be just direction not magnitude of momentum, but sometimes not with modulus exactly 1
      Pvc[0][0]=neutrino_E*neutrino_px/neutrino_p;
      Pvc[0][1]=neutrino_E*neutrino_py/neutrino_p;
      Pvc[0][2]=neutrino_E*neutrino_pz/neutrino_p;
      Abspvc[0] = neutrino_E;
      //Second element is incoming nucleon
      //Need to get truth info somehow
      Ichvc[1]=0;
      Ipvc[1]=2112;
      Pvc[1][0]=0;
      Pvc[1][1]=0;
      Pvc[1][2]=0;
      Abspvc[1]=0;
      //Third element is outgoing lepton
      int i_lepton=0; while(mpid[i_lepton]<11||mpid[i_lepton]>14) i_lepton++;
      Ichvc[2]=1;
      Ipvc[2]=mpid[i_lepton];
      double m = abs(mpid[i_lepton])==11 ? 0.511 :
                 abs(mpid[i_lepton]==13) ? 105 : 0;
      double E = KE[i_lepton]+m;
      double p = TMath::Sqrt(E*E-m*m);
      double pp = TMath::Sqrt(px[i_lepton]*px[i_lepton]+py[i_lepton]*py[i_lepton]+pz[i_lepton]*pz[i_lepton]);
      Pvc[2][0]=p*px[i_lepton]/pp;
      Pvc[2][1]=p*py[i_lepton]/pp;
      Pvc[2][2]=p*pz[i_lepton]/pp;
      Abspvc[2]=p;
      //Fourth element is outgoing nucleon
      int i_nucleon =0; while(mpid[i_nucleon]!=2112&&mpid[i_nucleon]!=2212&& i_nucleon <ntrks) i_nucleon++;
      if(i_nucleon==ntrks){
         Ichvc[3] = 0;
         Ipvc[3] = 0;
         Pvc[3][0] = 0;
         Pvc[3][1] = 0;
         Pvc[3][2] = 0;
      }else {
         Ichvc[3] = 1;
         Ipvc[3] = mpid[i_nucleon];
         m = abs(mpid[i_nucleon]==2212) ? 938 :
             abs(mpid[i_nucleon]==2112) ? 940 : 0;
         E = KE[i_nucleon]+m;
         p = TMath::Sqrt(E*E-m*m);
         pp = TMath::Sqrt(px[i_lepton]*px[i_lepton]+py[i_lepton]*py[i_lepton]+pz[i_lepton]*pz[i_lepton]);
         Pvc[3][0] = p*px[i_nucleon]/pp;
         Pvc[3][1] = p*py[i_nucleon]/pp;
         Pvc[3][2] = p*pz[i_nucleon]/pp;
         Abspvc[3] = p;
      }
      //Remaining of particles
      int i=4;
      for(int i_track =0; i_track <Npvc; i_track++){
         if(i_track==i_nucleon||i_track==i_lepton) continue;
         Ichvc[i] = 1;
         Ipvc[i] = mpid[i_track];
         m = abs(mpid[i_track])==11 ? 0.511 :
                    abs(mpid[i_track]==13) ? 105 :
                    abs(mpid[i_track]==111) ? 135 :
                    abs(mpid[i_track]==211) ? 140 :
                    abs(mpid[i_track]==2212) ? 938 :
                    abs(mpid[i_track]==2112) ? 940 :
                    abs(mpid[i_track]==321) ? 494 :
                    abs(mpid[i_track]==311) ? 498 :
                    abs(mpid[i_track]==310) ? 498 :
                    abs(mpid[i_track]==130) ? 498 :
                    abs(mpid[i_track]==221) ? 548 :
                    abs(mpid[i_track]==3122) ? 1116 : 0;
         E = KE[i_track]+m;
         p = TMath::Sqrt(E*E-m*m);
         pp = TMath::Sqrt(px[i_lepton]*px[i_lepton]+py[i_lepton]*py[i_lepton]+pz[i_lepton]*pz[i_lepton]);
         Pvc[i][0] = p*px[i_track]/pp;
         Pvc[i][1] = p*py[i_track]/pp;
         Pvc[i][2] = p*pz[i_track]/pp;
         Abspvc[i] = p;
         i++;
      }
       */
      //std::cout << jentry << " " << NEipvc[0] << std::endl;
      double a = 1- dirZ * dirZ;
      Double_t dirX = pid==13 ? recoDirXHighEMuon[0] : recoDirXHighEElectron[0];
      Double_t dirY = pid==13 ? recoDirYHighEMuon[0] : recoDirYHighEElectron[0];
      double b = vtxX * dirX + vtxY * dirY;
      double c = recoVtxR2-550*550;
      double recoToWallR = (TMath::Sqrt(b*b-a*c)-b)/a;
      double recoToWallZ = 1100 - vtxZ *TMath::Abs(dirZ);
      double recoToWall = recoToWallR<recoToWallZ ? recoToWallR : recoToWallZ;
      bool neutron_tag = false;
/*      for(int iSubevent = 0; iSubevent< nSubevents && !neutron_tag; iSubevent++){
         if(ring[iSubevent]==0 && recoEnergy[iSubevent]>2 && recoEnergy[iSubevent]<10
            && recoTime[iSubevent]>1000 && recoTime[iSubevent]<100000)
            neutron_tag = true;
      }
*/
      for(int iSubevent = 0; iSubevent < nSubevents; iSubevent ++) {
         if (recoEnergy[iSubevent] > 2 && recoEnergy[iSubevent] < 10 && recoTime[iSubevent] > 1000 && recoTime[iSubevent] < 100000) {
            neutron_tag = true;
            break;
         }
      }
      bool muLike = pid==13 && recoKE > 200 && recoKE < 2000 && Erec < 1250;
      bool eLike = pid==11;
      if(neutron_tag) {
         if (muLike) {
            if (NEipvc[0] == 14)
               numu_mulike_ntag->Fill();
            else if (NEipvc[0] == -14)
               numubar_mulike_ntag->Fill();
            else if (NEipvc[0] == 12)
               nue_mulike_ntag->Fill();
            else if (NEipvc[0] == -12)
               nuebar_mulike_ntag->Fill();
         }
         else if (eLike) {
            if (NEipvc[0] == 14)
               numu_elike_ntag->Fill();
            else if (NEipvc[0] == -14)
               numubar_elike_ntag->Fill();
            else if (NEipvc[0] == 12)
               nue_elike_ntag->Fill();
            else if (NEipvc[0] == -12)
               nuebar_elike_ntag->Fill();
         }
      }
      else {
         if (muLike) {
            if (NEipvc[0] == 14)
               numu_mulike_nontag->Fill();
            else if (NEipvc[0] == -14)
               numubar_mulike_nontag->Fill();
            else if (NEipvc[0] == 12)
               nue_mulike_nontag->Fill();
            else if (NEipvc[0] == -12)
               nuebar_mulike_nontag->Fill();
         }
         else if (eLike) {
            if (NEipvc[0] == 14)
               numu_elike_nontag->Fill();
            else if (NEipvc[0] == -14)
               numubar_elike_nontag->Fill();
            else if (NEipvc[0] == 12)
               nue_elike_nontag->Fill();
            else if (NEipvc[0] == -12)
               nuebar_elike_nontag->Fill();
         }
      }
   }
   numu_mulike_ntag->Write();
   numubar_mulike_ntag->Write();
   nue_mulike_ntag->Write();
   nuebar_mulike_ntag->Write();
   numu_elike_ntag->Write();
   numubar_elike_ntag->Write();
   nue_elike_ntag->Write();
   nuebar_elike_ntag->Write();
   numu_mulike_nontag->Write();
   numubar_mulike_nontag->Write();
   nue_mulike_nontag->Write();
   nuebar_mulike_nontag->Write();
   numu_elike_nontag->Write();
   numubar_elike_nontag->Write();
   nue_elike_nontag->Write();
   nuebar_elike_nontag->Write();
   out->Close();
}

TTree *MakeSelections::MakeNtuple(TString name, float Abspvc[100], int &Numatom, double &Erec, double &Etrue, Float_t pnu[100],
                                  Float_t dirnu[100][3], Int_t &Numbndn, Int_t &Numbndp, Int_t &Numfrep, Int_t &Ibound, Int_t &flavEvents) {
   TTree * ntuple = new TTree(name, name);
   ntuple->Branch("evt",&evt,"evt/I");

   ntuple->Branch("Npvc",&NEnvc, "Npvc/I");
   ntuple->Branch("Ichvc", NEicrnvc, "Ichvc[Npvc]/I");
   ntuple->Branch("Ipvc", NEipvc, "Ipvc[Npvc]/I");
   ntuple->Branch("Abspvc",Abspvc, "Abspvc[Npvc]/F");
   ntuple->Branch("Pvc", NEpvc, "Pvc[Npvc][3]/F");
   ntuple->Branch("mode",&NEneutmode, "mode/I");

   ntuple->Branch("Crsx",&NEcrsx, "Crsx/F");
   ntuple->Branch("Crsy",&NEcrsy, "Crsy/F");
   ntuple->Branch("Crsz",&NEcrsz, "Crsz/F");
   ntuple->Branch("Crsphi",&NEcrsphi, "Crsphi/F");
   ntuple->Branch("Nvert",&NEnvert, "Nvert/I");
   ntuple->Branch("Iflgvert", NEiflgvert, "Iflgvert[Nvert]/I");
   ntuple->Branch("posvert", NEposvert, "posvert[Nvert][3]/F");
   ntuple->Branch("Nvcvert",&NEnvcvert, "Nvcvert/I");
   ntuple->Branch("Abspvert", NEabspvert, "Abspvert[Nvcvert]/F");
   ntuple->Branch("Abstpvert", NEabstpvert, "Abstpvert[Nvcvert]/F");
   ntuple->Branch("Ipvert", NEipvert, "Ipvert[Nvcvert]/I");
   ntuple->Branch("Iverti", NEiverti, "Iverti[Nvcvert]/I");
   ntuple->Branch("Ivertf", NEivertf, "Ivertf[Nvcvert]/I");
   ntuple->Branch("Dirvert", NEdirvert, "Dirvert[Nvcvert][3]/F");
   ntuple->Branch("numnu",&NEnvc, "numnu/I");
   ntuple->Branch("ipnu", NEipvc, "ipnu[numnu]/I");
   ntuple->Branch("pnu",pnu, "pnu[numnu][3]/F");
   ntuple->Branch("dirnu",dirnu, "dirnu[numnu][3]/F");
   ntuple->Branch("Numbndn",&Numbndn, "Numbndn/I");
   ntuple->Branch("Numbndp",&Numbndp, "Numbndp/I");
   ntuple->Branch("Numfrep",&Numfrep, "Numfrep/I");
   ntuple->Branch("Ibound",&Ibound, "Ibound/I");


   ntuple->Branch("Numatom",&Numatom, "Numatom/I");
   ntuple->Branch("Erec",&Erec,"Erec/D");
   ntuple->Branch("Etrue",&Etrue,"Etrue/D");

   ntuple->Branch("totFlavEvts",&flavEvents,"totFlavEvts/I");
   return ntuple;
}
