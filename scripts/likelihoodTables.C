#include "include/LikelihoodGenerator.hh"
#include "TH2D.h"
#include "TMath.h"

using namespace std;

void likelihoodTablesAll(const char *name) {
//  const int minKE = 50;
//  const int maxKE = 1050;
//  const int stepKE = 50;
  const int nKE = 21;
  int binCentresKE[nKE] = {50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050/*,3000*/};
  double binEdgesKE[nKE+1] = {25,75,125,175,225,275,325,375,425,475,525,575,625,675,725,775,825,875,925,975,1025,1075/*,4925*/};
  TFile * likelihoodTables = TFile::Open(Form("likelihood_tables_%s.root",name), "RECREATE");
  double min_pmt_dist = 0, max_pmt_dist =2460; //max depends on detector size
  double min_pmt_angle = 0, max_pmt_angle = TMath::Pi();
  const int pmt_dist_bins = 250;
  const int pmt_angle_bins = 500;
  double pmt_dist_binwidth = (max_pmt_dist - min_pmt_dist)/ pmt_dist_bins;
  double pmt_angle_binwidth = (max_pmt_angle - min_pmt_angle)/ pmt_angle_bins;
  double binEdgesPMTDist[pmt_dist_bins +1];
  double binEdgesPMTAngle[pmt_angle_bins +1];
  for(int i=0; i<= pmt_dist_bins; i++){
    binEdgesPMTDist[i]= min_pmt_dist + pmt_dist_binwidth *i;
  }
  for(int i=0; i<= pmt_angle_bins; i++){
    binEdgesPMTAngle[i]= min_pmt_angle + pmt_angle_binwidth *i;
  }
  TH3D * photons_e = new TH3D("photons_e","photons_e", pmt_dist_bins, binEdgesPMTDist, pmt_angle_bins, binEdgesPMTAngle,nKE,binEdgesKE);
  TH3D * photons_mu = new TH3D("photons_mu","photons_mu", pmt_dist_bins, binEdgesPMTDist, pmt_angle_bins,
                               binEdgesPMTAngle,nKE,binEdgesKE);
  //int pdfBins = 500;
  //double max_track_length=500;
  //TH3D * pdfs_e = new TH3D("pdfs_e","pdfs_e", pdfBins,0,max_track_length, pdfBins,0,TMath::Pi(),nKE,minKE-stepKE/2.,maxKE+stepKE/2.);
  //TH3D * pdfs_mu = new TH3D("pdfs_mu","pdfs_mu", pdfBins,0,max_track_length, pdfBins,0,TMath::Pi(),nKE,minKE-stepKE/2.,maxKE+stepKE/2.);
  for(int iKE = 0; iKE< nKE; iKE++){
    int KE = binCentresKE[iKE];
    cout << KE << "MeV:" << endl;
    TChain * ch_e = new TChain("EventTree");
    TChain * ch_mu = new TChain("EventTree");
      int fileCount = /*iKE==nKE ? 1 :*/ 5;
      for(int i=1; i<= fileCount; i++){
        ch_e->AddFile(TString::Format("/data/hyperk/wchsandbox_reco/fix/e_%iMeV/e_%iMeV_%i_FullEvent.root",KE,KE,i));
        ch_mu->AddFile(TString::Format("/data/hyperk/wchsandbox_reco/fix/mu_%iMeV/mu_%iMeV_%i_FullEvent.root",KE,KE,i));
    }
    LikelihoodGenerator * tr_e = new LikelihoodGenerator(ch_e);
    LikelihoodGenerator * tr_mu = new LikelihoodGenerator(ch_mu);
    //tr_e->Loop(true);
    //tr_mu->Loop(false);
    tr_e->Loop2();
    tr_mu->Loop2();
    //TH2D * pdf_e = (TH2D*) tr_e->h->Clone(TString::Format("pdf_e_%iMeV",KE));
    //TH2D * pdf_mu = (TH2D*) tr_mu->h->Clone(TString::Format("pdf_mu_%iMeV",KE));
    TH2D * photon_e = (TH2D*) tr_e->hPhotonCount->Clone(TString::Format("photons_e_%iMeV",KE));
    TH2D * photon_mu = (TH2D*) tr_mu->hPhotonCount->Clone(TString::Format("photons_mu_%iMeV",KE));
    //for(int ix = 1; ix<=pdfBins; ix++) for(int iy=1; iy<=pdfBins; iy++){
    //    pdfs_e->SetBinContent(ix,iy,iKE+1,pdf_e->GetBinContent(ix,iy));
    //    pdfs_mu->SetBinContent(ix,iy,iKE+1,pdf_mu->GetBinContent(ix,iy));
    //  }
    for(int ix = 1; ix<= pmt_dist_bins +1; ix++) for(int iy=1; iy<= pmt_angle_bins +1; iy++){
        double photons = photon_e->GetBinContent(ix,iy);
        if(photons < 0.00001) photons = 0.00001;
        photons_e->SetBinContent(ix,iy,iKE+1,photons);
        photons = photon_mu->GetBinContent(ix,iy);
        if(photons < 0.00001) photons = 0.00001;
        photons_mu->SetBinContent(ix,iy,iKE+1,photons);
      }
    likelihoodTables->cd();
    //pdf_e->Write();
    //pdf_mu->Write();
    photon_e->Write();
    photon_mu->Write();
  }
  likelihoodTables->cd();
  //pdfs_e->Write();
  //pdfs_mu->Write();
  photons_e->Write();
  photons_mu->Write();
  likelihoodTables->Close();
}

void likelihoodTables(int KE, const char *emu, int i, const char *name) {
  TFile * likelihoodTable = TFile::Open(Form("/data/hyperk/wchsandbox_reco/likelihoods/likelihood_table_%s_%s_%i_%i.root",name,emu,KE,i), "RECREATE");
  cout << KE << "MeV:" << endl;
  TChain * ch = new TChain("EventTree");
  ch->AddFile(TString::Format("/data/hyperk/wchsandbox_reco/fix/%s_%iMeV/%s_%iMeV_%i_FullEvent.root",emu,KE,emu,KE,i));
  LikelihoodGenerator * lg = new LikelihoodGenerator(ch);
  //lg->Loop(true);
  lg->Loop2();
  //TH2D * pdf = (TH2D*) lg->h->Clone(TString::Format("pdf_%s_%iMeV",emu,KE));
  TH2D * photon = (TH2D*) lg->hPhotonCount->Clone(TString::Format("photons_%s_%iMeV_%i",emu,KE,i));
  TH2D * time = (TH2D*) lg->hTime->Clone(TString::Format("time_%s_%iMeV_%i",emu,KE,i));
  likelihoodTable->cd();
  //pdf->Write();
  photon->Write();
  time->Write();
  likelihoodTable->Close();
}

void combineTables(const char *name){
  const int nKE = 21;
  int binCentresKE[nKE] = {50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050/*,3000*/};
  double binEdgesKE[nKE+1] = {25,75,125,175,225,275,325,375,425,475,525,575,625,675,725,775,825,875,925,975,1025,1075/*,4925*/};
  TFile * likelihoodTables = TFile::Open(Form("likelihood_tables_%s.root",name), "RECREATE");
  double min_pmt_dist = 0, max_pmt_dist =2460; //max depends on detector size
  double min_pmt_angle = 0, max_pmt_angle = TMath::Pi();
  const int pmt_dist_bins = 250;
  const int pmt_angle_bins = 500;
  double pmt_dist_binwidth = (max_pmt_dist - min_pmt_dist)/ pmt_dist_bins;
  double pmt_angle_binwidth = (max_pmt_angle - min_pmt_angle)/ pmt_angle_bins;
  double binEdgesPMTDist[pmt_dist_bins +1];
  double binEdgesPMTAngle[pmt_angle_bins +1];
  for(int i=0; i<= pmt_dist_bins; i++){
    binEdgesPMTDist[i]= min_pmt_dist + pmt_dist_binwidth *i;
  }
  for(int i=0; i<= pmt_angle_bins; i++){
    binEdgesPMTAngle[i]= min_pmt_angle + pmt_angle_binwidth *i;
  }
  TH3D * photons_e = new TH3D("photons_e","photons_e", pmt_dist_bins, binEdgesPMTDist, pmt_angle_bins, binEdgesPMTAngle,nKE,binEdgesKE);
  TH3D * photons_mu = new TH3D("photons_mu","photons_mu", pmt_dist_bins, binEdgesPMTDist, pmt_angle_bins,
                               binEdgesPMTAngle,nKE,binEdgesKE);
  TH3D * times_e = new TH3D("times_e","times_e", pmt_dist_bins, binEdgesPMTDist, pmt_angle_bins, binEdgesPMTAngle,nKE,binEdgesKE);
  TH3D * times_mu = new TH3D("times_mu","times_mu", pmt_dist_bins, binEdgesPMTDist, pmt_angle_bins, binEdgesPMTAngle,nKE,binEdgesKE);
  //int pdfBins = 500;
  //double max_track_length=500;
  //TH3D * pdfs_e = new TH3D("pdfs_e","pdfs_e", pdfBins,0,max_track_length, pdfBins,0,TMath::Pi(),nKE,minKE-stepKE/2.,maxKE+stepKE/2.);
  //TH3D * pdfs_mu = new TH3D("pdfs_mu","pdfs_mu", pdfBins,0,max_track_length, pdfBins,0,TMath::Pi(),nKE,minKE-stepKE/2.,maxKE+stepKE/2.);
  for(int iKE = 0; iKE< nKE; iKE++) {
    int KE = binCentresKE[iKE];
    cout << KE << "MeV:" << endl;
    int fileCount = /*iKE==nKE ? 1 :*/ 5;
    TH2D * hePhotonCount = new TH2D("hep", "hep", pmt_dist_bins, binEdgesPMTDist, pmt_angle_bins, binEdgesPMTAngle);
    TH2D * hmuPhotonCount = new TH2D("hmup", "hmup", pmt_dist_bins, binEdgesPMTDist, pmt_angle_bins, binEdgesPMTAngle);
    TH2D * heTime = new TH2D("het", "het", pmt_dist_bins, binEdgesPMTDist, pmt_angle_bins, binEdgesPMTAngle);
    TH2D * hmuTime = new TH2D("hmut", "hmut", pmt_dist_bins, binEdgesPMTDist, pmt_angle_bins, binEdgesPMTAngle);
    for (int i = 1; i <= fileCount; i++) {
      cout << i << endl;
      TFile *e_file = TFile::Open(Form("/data/hyperk/wchsandbox_reco/likelihoods/likelihood_table_%s_e_%i_%i.root", name, KE, i), "READ");
      TFile *mu_file = TFile::Open(Form("/data/hyperk/wchsandbox_reco/likelihoods/likelihood_table_%s_mu_%i_%i.root", name, KE, i), "READ");
      TH2D *photon_e = (TH2D *) e_file->Get(TString::Format("photons_e_%iMeV_%i", KE, i));
      TH2D *photon_mu = (TH2D *) mu_file->Get(TString::Format("photons_mu_%iMeV_%i", KE, i));
      TH2D *time_e = (TH2D *) e_file->Get(TString::Format("time_e_%iMeV_%i", KE, i));
      TH2D *time_mu = (TH2D *) mu_file->Get(TString::Format("time_mu_%iMeV_%i", KE, i));
      hePhotonCount->Add(photon_e);
      hmuPhotonCount->Add(photon_mu);
      time_e->Multiply(photon_e);
      time_mu->Multiply(photon_mu);
      heTime->Add(time_e);
      hmuTime->Add(time_mu);
      e_file->Close();
      mu_file->Close();
    }
    heTime->Divide(hePhotonCount);
    hmuTime->Divide(hmuPhotonCount);
    hePhotonCount->Scale(1./fileCount);
    hmuPhotonCount->Scale(1./fileCount);
    for(int ix = 1; ix<=pmt_dist_bins; ix++) for(int iy=1; iy<=pmt_angle_bins; iy++){
        photons_e->SetBinContent(ix,iy,iKE+1,hePhotonCount->GetBinContent(ix,iy));
        photons_mu->SetBinContent(ix,iy,iKE+1,hmuPhotonCount->GetBinContent(ix,iy));
        times_e->SetBinContent(ix,iy,iKE+1,heTime->GetBinContent(ix,iy));
        times_mu->SetBinContent(ix,iy,iKE+1,hmuTime->GetBinContent(ix,iy));
      }
    delete hePhotonCount;
    delete hmuPhotonCount;
    delete heTime;
    delete hmuTime;
//    for(int ix = 1; ix<=pmt_dist_bins; ix++) for(int iy=1; iy<=pmt_angle_bins; iy++){
//        if(photons_e->GetBinContent(ix,iy,KE+1)<0.0001)
//          photons_e->SetBinContent(ix,iy,iKE+1,0.0001);
//        if(photons_mu->GetBinContent(ix,iy,KE+1)<0.0001)
//          photons_mu->SetBinContent(ix,iy,iKE+1,0.0001);
//      }
  }
  likelihoodTables->cd();
  //pdfs_e->Write();
  //pdfs_mu->Write();
  photons_e->Write();
  photons_mu->Write();
  times_e->Write();
  times_mu->Write();
  likelihoodTables->Close();
}