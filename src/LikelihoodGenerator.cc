#include "LikelihoodGenerator.hh"
#include <TH2.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TObject.h>

ClassImp(LikelihoodGenerator)


LikelihoodGenerator::LikelihoodGenerator(TString file) : fChain(0)
{
  TTree *tree;
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(file);
  if (!f || !f->IsOpen()) {
    f = new TFile(file);
  }
  f->GetObject("EventTree",tree);

  Init(tree);
}

LikelihoodGenerator::LikelihoodGenerator(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("mu_FullEvent.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("mu_FullEvent.root");
    }
    f->GetObject("EventTree",tree);

  }
  Init(tree);
}

LikelihoodGenerator::~LikelihoodGenerator()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t LikelihoodGenerator::GetEntry(Long64_t entry)
{
// Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t LikelihoodGenerator::LoadTree(Long64_t entry)
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

void LikelihoodGenerator::Init(TTree *tree)
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
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("evt", &evt, &b_evt);
  fChain->SetBranchAddress("nphot", &nphot, &b_nphot);
  fChain->SetBranchAddress("npart", &npart, &b_npart);
  fChain->SetBranchAddress("ncapturecount", &ncapturecount, &b_ncapturecount);
  fChain->SetBranchAddress("neutroncount", &neutroncount, &b_neutroncount);
  fChain->SetBranchAddress("phot_xStart", phot_xStart, &b_phot_xStart);
  fChain->SetBranchAddress("phot_yStart", phot_yStart, &b_phot_yStart);
  fChain->SetBranchAddress("phot_zStart", phot_zStart, &b_phot_zStart);
  fChain->SetBranchAddress("phot_tStart", phot_tStart, &b_phot_tStart);
  fChain->SetBranchAddress("phot_xEnd", phot_xEnd, &b_phot_xEnd);
  fChain->SetBranchAddress("phot_yEnd", phot_yEnd, &b_phot_yEnd);
  fChain->SetBranchAddress("phot_zEnd", phot_zEnd, &b_phot_zEnd);
  fChain->SetBranchAddress("phot_tEnd", phot_tEnd, &b_phot_tEnd);
  fChain->SetBranchAddress("phot_wavelength", phot_wavelength, &b_phot_wavelength);
  fChain->SetBranchAddress("phot_processStart", phot_processStart, &b_phot_processStart);
  fChain->SetBranchAddress("phot_isScat", phot_isScat, &b_phot_isScat);
  fChain->SetBranchAddress("phot_parentid", phot_parentid, &b_phot_parentid);
  fChain->SetBranchAddress("phot_trackid", phot_trackid, &b_phot_trackid);
  fChain->SetBranchAddress("phot_hit", phot_hit, &b_phot_hit);
  fChain->SetBranchAddress("phot_capnum", phot_capnum, &b_phot_capnum);
  fChain->SetBranchAddress("part_xStart", part_xStart, &b_part_xStart);
  fChain->SetBranchAddress("part_yStart", part_yStart, &b_part_yStart);
  fChain->SetBranchAddress("part_zStart", part_zStart, &b_part_zStart);
  fChain->SetBranchAddress("part_tStart", part_tStart, &b_part_tStart);
  fChain->SetBranchAddress("part_xEnd", part_xEnd, &b_part_xEnd);
  fChain->SetBranchAddress("part_yEnd", part_yEnd, &b_part_yEnd);
  fChain->SetBranchAddress("part_zEnd", part_zEnd, &b_part_zEnd);
  fChain->SetBranchAddress("part_tEnd", part_tEnd, &b_part_tEnd);
  fChain->SetBranchAddress("part_pxStart", part_pxStart, &b_part_pxStart);
  fChain->SetBranchAddress("part_pyStart", part_pyStart, &b_part_pyStart);
  fChain->SetBranchAddress("part_pzStart", part_pzStart, &b_part_pzStart);
  fChain->SetBranchAddress("part_pxEnd", part_pxEnd, &b_part_pxEnd);
  fChain->SetBranchAddress("part_pyEnd", part_pyEnd, &b_part_pyEnd);
  fChain->SetBranchAddress("part_pzEnd", part_pzEnd, &b_part_pzEnd);
  fChain->SetBranchAddress("part_KEstart", part_KEstart, &b_part_KEstart);
  fChain->SetBranchAddress("part_KEend", part_KEend, &b_part_KEend);
  fChain->SetBranchAddress("part_processStart", part_processStart, &b_part_processStart);
  fChain->SetBranchAddress("part_processEnd", part_processEnd, &b_part_processEnd);
  fChain->SetBranchAddress("part_parentid", part_parentid, &b_part_parentid);
  fChain->SetBranchAddress("part_trackid", part_trackid, &b_part_trackid);
  fChain->SetBranchAddress("part_pid", part_pid, &b_part_pid);
  fChain->SetBranchAddress("capt_x", capt_x, &b_capt_x);
  fChain->SetBranchAddress("capt_y", capt_y, &b_capt_y);
  fChain->SetBranchAddress("capt_z", capt_z, &b_capt_z);
  fChain->SetBranchAddress("capt_t0", capt_t0, &b_capt_t0);
  fChain->SetBranchAddress("capt_E", capt_E, &b_capt_E);
  fChain->SetBranchAddress("capt_num", capt_num, &b_capt_num);
  fChain->SetBranchAddress("capt_pid", capt_pid, &b_capt_pid);
  fChain->SetBranchAddress("capt_nucleus", capt_nucleus, &b_capt_nucleus);
  fChain->SetBranchAddress("capt_nphot", capt_nphot, &b_capt_nphot);
  fChain->SetBranchAddress("capt_ngamma", capt_ngamma, &b_capt_ngamma);
  Notify();
}

Bool_t LikelihoodGenerator::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void LikelihoodGenerator::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t LikelihoodGenerator::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  return 1;
}


// First build pdfs:
//  For each photon determine position along track and emission angle and add to histogram
//  Use only photons produced from muon itself for
void LikelihoodGenerator::Loop(bool isElectron, int start, int end)
{
  std::cout << "Generating pdf... " << std::endl;
  if (fChain == 0) return;
  int nentries = fChain->GetEntries();
  if(end <= start) end = nentries;
  //h = new TH2D("h","h",100,0,500,100,angleBins);
  hPDF = new TH2D("h","h",1000,0,500,1000,0,TMath::Pi());
  Long64_t nbytes = 0, nb = 0;
  int count = end-start;
  for (Long64_t jentry=start; jentry<end;jentry++) {
    if(jentry%100==0) std::cout << jentry/100 << " " << std::endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    int ipart=0;
    while(part_parentid[ipart]!=0) ipart++;
    for(int iphot=0; iphot<nphot; iphot++){
      // Ignore photons not seen by PMT due to wavelength acceptance
      //if(phot_wavelength[iphot]>600 || phot_wavelength[iphot]<300) continue;
      double QE_Wavelengths[20] = {280., 300., 320., 340., 360., 380., 400., 420., 440., 460., 480., 500., 520., 540., 560., 580., 600., 620., 640., 660.};
      double QE_Factor[20] = {0.00, .066, .405, .801, .962, .976, 1., .957, .899, .791, .664, .550, .382, .205, .126, .069, .036, .024, .007, 0.00};
      double qe=0;
      if(phot_wavelength[iphot]<QE_Wavelengths[0]) continue;
      for(int i=0; i<19; i++) {
        if (phot_wavelength[iphot] < QE_Wavelengths[i + 1]) {
          double delta = (phot_wavelength[iphot] - QE_Wavelengths[i])/(QE_Wavelengths[i + 1] - QE_Wavelengths[i]);
          qe = (QE_Factor[i] * (1 - delta) + QE_Factor[i + 1] * delta) / 100.;
          break;
        }
      }
      if(qe<=0) continue;
      // For muons, only take photons from muon
      if(phot_parentid[iphot] != 1 && !isElectron) continue;
      // For electrons, allow all photons but only take those that reach the wall
      if(TMath::Abs(TMath::Abs(phot_zEnd[iphot])-11000)>1
          && TMath::Abs(TMath::Sqrt(phot_xEnd[iphot]*phot_xEnd[iphot]+phot_yEnd[iphot]*phot_yEnd[iphot])-5500)>1 && isElectron) continue;
      double phot_xDir = phot_xEnd[iphot]-phot_xStart[iphot];
      double phot_yDir = phot_yEnd[iphot]-phot_yStart[iphot];
      double phot_zDir = phot_zEnd[iphot]-phot_zStart[iphot];
      double phot_dist = TMath::Sqrt(phot_xDir*phot_xDir+phot_yDir*phot_yDir+phot_zDir*phot_zDir);
      phot_xDir/=phot_dist;
      phot_yDir/=phot_dist;
      phot_zDir/=phot_dist;
      double theta = TMath::ACos(phot_xDir*part_pxStart[ipart]+phot_yDir*part_pyStart[ipart]+phot_zDir*part_pzStart[ipart]);
      double track_x = phot_xStart[iphot]-part_xStart[ipart];
      double track_y = phot_yStart[iphot]-part_yStart[ipart];
      double track_z = phot_zStart[iphot]-part_zStart[ipart];
      double track_dist = TMath::Sqrt(track_x*track_x+track_y*track_y+track_z*track_z)/10.;
      hPDF->Fill(track_dist, theta, qe);
    }
  }
  
  hPDF->Scale(1./count);

  double min_track_dist = 0, max_track_dist=2460; //max depends on detector size
  double min_track_angle = 0, max_track_angle = TMath::Pi();
  const int track_dist_bins = 250;
  const int track_angle_bins = 500;
  double track_dist_binwidth = (max_track_dist-min_track_dist)/track_dist_bins;
  double track_angle_binwidth = (max_track_angle-min_track_angle)/track_angle_bins;
  min_track_dist -= 0.5*track_dist_binwidth;
  max_track_dist += 0.5*track_dist_binwidth;
  min_track_angle -= 0.5*track_angle_binwidth;
  max_track_angle += 0.5*track_angle_binwidth;
  double binEdgesTrackDist[track_dist_bins+1];
  double binEdgesTrackAngle[track_angle_bins+1];
  for(int i=0; i<=track_dist_bins; i++){
    binEdgesTrackDist[i]=min_track_dist+track_dist_binwidth*i;
  }
  for(int i=0; i<=track_angle_bins; i++){
    binEdgesTrackAngle[i]=min_track_angle+track_angle_binwidth*i;
  }
  hPhotonCount = new TH2D("h2", "h2", track_dist_bins, binEdgesTrackDist, track_angle_bins, binEdgesTrackAngle);

  std::cout << "Generating likelihood tables... " << std::endl;
  /*
  // Loop over bins in pmt_dist and pmt_track_angle
  int nbinsx = h2->GetNbinsX();
  int nbinsy = h2->GetNbinsY();
  for(int ix=1; ix<=nbinsx; ix++){
    double pmt_dist = h2->GetXaxis()->GetBinCenter(ix);
    if(ix%100==0) std::cout << ix/100 << " " << std::endl;
    for(int iy=1; iy<=nbinsy; iy++){
      double pmt_track_angle = h2->GetYaxis()->GetBinCenter(iy);
      double likelihood = 0;

      // Loop over bins in track distance, s
      int nbinss = h->GetNbinsX();
      for(int is = 1; is <=nbinss; is++){
        double track_dist = h->GetXaxis()->GetBinCenter(is);
        double theta = TMath::ATan2(pmt_dist*TMath::Sin(pmt_track_angle),pmt_dist*TMath::Cos(pmt_track_angle)-track_dist);
        int thetaBin = h->GetYaxis()->FindFixBin(theta);
        double prob = h->GetBinContent(is,thetaBin)*h->GetXaxis()->GetBinWidth(is);
        double phot_dist = pmt_dist*TMath::Sin(pmt_track_angle)/TMath::Sin(theta);
        double attenuationProb = TMath::Exp(-phot_dist/9000.);
        likelihood += prob*attenuationProb;
//        double time = (pmt_dist + N_REF*pmt_dist*TMath::Sin(pmt_track_angle)/TMath::Sin(theta))/C_VAC;
//        h3->Fill(pmt_dist,pmt_track_angle,time,prob);
      }
      h2->SetBinContent(ix,iy,likelihood);
    }
  }
  */
  //Loop over bins in pdf
  int nbins_track = hPDF->GetNbinsX();
  int nbins_theta = hPDF->GetNbinsY();
  //double theta_binwidth = h->GetYaxis()->GetBinWidth(1);
  for(int itheta=1; itheta<=nbins_theta; itheta++){
    if(itheta%100==0) std::cout << itheta/100 << " " << std::endl;
    double theta = hPDF->GetYaxis()->GetBinCenter(itheta);
    double sin_theta = TMath::Sin(theta);
    double cos_theta = TMath::Cos(theta);
    for(int itrack=1; itrack<=nbins_track; itrack++){
      double track_dist = hPDF->GetXaxis()->GetBinCenter(itrack);
      double expected_photons = hPDF->GetBinContent(itrack,itheta);

/*

      // Loop over photon dist
      double track_dist_sqr = track_dist*track_dist;
      double b = track_dist*TMath::Cos(theta);
      double max_phot_dist = -b+TMath::Sqrt(b*b+2460*2460-track_dist_sqr);
      int loop_max = max_phot_dist;
      int old_bin = -999;
      for(int phot_dist=0; phot_dist<=loop_max; phot_dist++){
        double phot_dist_sqr = phot_dist*phot_dist;
        double pmt_dist = TMath::Sqrt(track_dist_sqr+phot_dist_sqr+2.*b*phot_dist);
//        if(pmt_dist>2460) break;
        double angle = TMath::ACos((track_dist_sqr+pmt_dist*pmt_dist-phot_dist_sqr)/(2.*track_dist*pmt_dist));
        int new_pmt_dist_bin = h2->GetXaxis()->FindFixBin(pmt_dist);
        int new_angle_bin = h2->GetYaxis()->FindFixBin(angle);
        int new_bin = h2->GetBin(new_pmt_dist_bin, new_angle_bin);
        if(new_bin == old_bin) continue;
        old_bin = new_bin;
        double attenuationProb = TMath::Exp(-phot_dist/9000.);
        h2->AddBinContent(new_bin,expected_photons*attenuationProb);
      }
*/

/*
      // Loop over angle to PMT
      double min_pmt_angle = 0;
      double max_pmt_angle = theta-TMath::ASin(track_dist*sin_theta/2460.);
      int max_pmt_angle_bin = h2->GetYaxis()->FindFixBin(max_pmt_angle);
      for(int i_pmt_angle=1; i_pmt_angle <= max_pmt_angle_bin; i_pmt_angle++){
        //Loop over dist to PMT
        double min_pmt_dist = track_dist*sin_theta/TMath::Sin(theta-min_pmt_angle);
        double max_pmt_dist = track_dist*sin_theta/TMath::Sin(theta-min_pmt_angle+theta_binwidth);
        int min_pmt_dist_bin, max_pmt_dist_bin;
        if(max_pmt_dist > min_pmt_dist){
          min_pmt_dist_bin = h2->GetXaxis()->FindFixBin(min_pmt_dist);
          max_pmt_dist_bin = h2->GetXaxis()->FindFixBin(max_pmt_dist);
        }
        else{
          min_pmt_dist_bin = h2->GetXaxis()->FindFixBin(max_pmt_dist);
          max_pmt_dist_bin = h2->GetXaxis()->FindFixBin(min_pmt_dist);
        }
        for(int i_pmt_dist = min_pmt_dist_bin; i_pmt_dist <= max_pmt_dist_bin; i_pmt_dist++){
          double phot_dist = h2->GetXaxis()->GetBinCenter(i_pmt_dist)*TMath::Sin(h2->GetYaxis()->GetBinCenter(i_pmt_angle))/sin_theta;
          double attenuationProb = TMath::Exp(-phot_dist/9000.);
          h2->AddBinContent(h2->GetBin(i_pmt_dist, i_pmt_angle), expected_photons*attenuationProb);
        }
        min_pmt_angle += theta_binwidth;
      }*/



      // Loop over distance to PMT
      double track_dist_sqr = track_dist*track_dist;
      double min_pmt_dist = track_dist;
      if(cos_theta<0) min_pmt_dist *= sin_theta;
      int min_pmt_dist_bin = hPhotonCount->GetXaxis()->FindFixBin(min_pmt_dist);
      if(hPhotonCount->GetXaxis()->GetBinCenter(min_pmt_dist_bin)<min_pmt_dist) min_pmt_dist_bin++;
      int max_pmt_dist_bin = hPhotonCount->GetXaxis()->FindFixBin(2460.);
      for(int i_pmt_dist=min_pmt_dist_bin; i_pmt_dist <= max_pmt_dist_bin; i_pmt_dist++){
        double pmt_dist = hPhotonCount->GetXaxis()->GetBinCenter(i_pmt_dist);
        double pmt_dist_sqr = pmt_dist*pmt_dist;
        //Determine photon track dist
        double b = track_dist*cos_theta;
        double b2 = b*b;
        double s = TMath::Sqrt(b2+pmt_dist_sqr-track_dist_sqr);
        double phot_dist = s-b;
        double attenuationProb = TMath::Exp(-phot_dist/9000.);
        //Determine angle from vtx to pmt at this distance
        double pmt_angle = TMath::ACos((pmt_dist_sqr+track_dist_sqr-phot_dist*phot_dist)/(2*track_dist*pmt_dist));
        int i_pmt_angle = hPhotonCount->GetYaxis()->FindFixBin(pmt_angle);
        hPhotonCount->AddBinContent(hPhotonCount->GetBin(i_pmt_dist, i_pmt_angle), expected_photons*attenuationProb);
        //Check for second possible solution for angle
        if(-b<s) continue;
        phot_dist = -s-b;
        attenuationProb = TMath::Exp(-phot_dist/9000.);
        pmt_angle = TMath::ACos((pmt_dist_sqr+track_dist_sqr-phot_dist*phot_dist)/(2*track_dist*pmt_dist));
        i_pmt_angle = hPhotonCount->GetYaxis()->FindFixBin(pmt_angle);
        hPhotonCount->AddBinContent(hPhotonCount->GetBin(i_pmt_dist, i_pmt_angle), expected_photons*attenuationProb);
      }
    }
  }
  hPhotonCount->SetEntries(1); //needed because hist is filled using AddBinContent


//  h->Draw();
}



void LikelihoodGenerator::Loop2(int start,int end)
{
  std::cout << "Generating likelihood tables... " << std::endl;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  if(end <= start) end = nentries;
  //h = new TH2D("h","h",100,0,500,100,angleBins);
  hPDF = new TH2D("h","h",1000,0,500,1000,0,TMath::Pi());

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
  hPhotonCount = new TH2D("h2", "h2", pmt_dist_bins, binEdgesPMTDist, pmt_angle_bins, binEdgesPMTAngle);
  hTime = new TH2D("h3", "h3", pmt_dist_bins, binEdgesPMTDist, pmt_angle_bins, binEdgesPMTAngle);

  const int n=2;
  Long64_t nbytes = 0, nb = 0;
  int count = end-start;
  for (Long64_t jentry=start; jentry<end;jentry++) {
    if(jentry%10==0) std::cout << jentry/10 << " " << std::endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    int ipart=0;
    while(part_parentid[ipart]!=0) ipart++;
    double QE_Wavelengths[20] = {280., 300., 320., 340., 360., 380., 400., 420., 440., 460., 480., 500., 520., 540., 560., 580., 600., 620., 640., 660.};
    double QE_Factor[20] = {0.00, .066, .405, .801, .962, .976, 1., .957, .899, .791, .664, .550, .382, .205, .126, .069, .036, .024, .007, 0.00};
    for(int iphot=0; iphot<nphot; iphot+=n){ // +=n so fewer photons are used to reduce time
      if(phot_tEnd[iphot]>150) continue;
//      if(phot_parentid[iphot] != 1) continue;      
      if(phot_wavelength[iphot]<QE_Wavelengths[0]) continue;
      if(phot_wavelength[iphot]>QE_Wavelengths[19]) continue;
      int i = (int) (phot_wavelength[iphot]/20.-14);
      double qe=0;
      //for(int i=0; i<19; i++) {
      //  if (phot_wavelength[iphot] < QE_Wavelengths[i + 1]) {
          double delta = (phot_wavelength[iphot] - QE_Wavelengths[i])/(QE_Wavelengths[i + 1] - QE_Wavelengths[i]);
          qe = (QE_Factor[i] * (1 - delta) + QE_Factor[i + 1] * delta) / 100.;
      //    break;
      //  }
      //}
      //if(qe<=0) continue;
      double phot_xDir = (phot_xEnd[iphot]-phot_xStart[iphot])/10.;
      double phot_yDir = (phot_yEnd[iphot]-phot_yStart[iphot])/10.;
      double phot_zDir = (phot_zEnd[iphot]-phot_zStart[iphot])/10.;
      double phot_dist_max = TMath::Sqrt(phot_xDir*phot_xDir+phot_yDir*phot_yDir+phot_zDir*phot_zDir);
 //     cout << phot_dist_max << " " << endl;
      phot_xDir/=phot_dist_max;
      phot_yDir/=phot_dist_max;
      phot_zDir/=phot_dist_max;
      double pmt_xStart = (phot_xStart[iphot]-part_xStart[ipart])/10.;
      double pmt_yStart = (phot_yStart[iphot]-part_yStart[ipart])/10.;
      double pmt_zStart = (phot_zStart[iphot]-part_zStart[ipart])/10.;
      double pmt_xEnd = (phot_xEnd[iphot]-part_xStart[ipart])/10.;
      double pmt_yEnd = (phot_yEnd[iphot]-part_yStart[ipart])/10.;
      double pmt_zEnd = (phot_zEnd[iphot]-part_zStart[ipart])/10.;
      double a = pmt_xStart*pmt_xStart + pmt_yStart*pmt_yStart + pmt_zStart*pmt_zStart;
      double b = pmt_xStart*phot_xDir + pmt_yStart*phot_yDir + pmt_zStart*phot_zDir;
      double b2 = b*b;
      double pmt_dist_start = TMath::Sqrt(a);
      
      double a0 = pmt_xStart*part_pxStart[ipart]+pmt_yStart*part_pyStart[ipart]+pmt_zStart*part_pzStart[ipart];
      double al = phot_xDir*part_pxStart[ipart]+phot_yDir*part_pyStart[ipart]+phot_zDir*part_pzStart[ipart];

      double tStart = phot_tStart[iphot];
      double invSpeed = (phot_tEnd[iphot]-tStart)/phot_dist_max;
      //Possibly photon gets closer to vertex initially
      if(b<0){
        double pmt_dist_min = TMath::Sqrt(a- b2);
        int min_pmt_dist_bin = TMath::CeilNint(pmt_dist_min/pmt_dist_binwidth)+1;
        for(int pmt_dist_bin = TMath::FloorNint(pmt_dist_start/pmt_dist_binwidth)+1; pmt_dist_bin>=min_pmt_dist_bin; pmt_dist_bin--){
          double d = (pmt_dist_bin-1)*pmt_dist_binwidth;
          double l = -b-TMath::Sqrt(b2 -a+d*d);
          double angle = TMath::ACos((a0+l*al)/d);
          int angle_bin = hPhotonCount->GetYaxis()->FindFixBin(angle);
          hPhotonCount->AddBinContent(hPhotonCount->GetBin(pmt_dist_bin,angle_bin),qe);
          double time = tStart+l*invSpeed;
          hTime->AddBinContent(hTime->GetBin(pmt_dist_bin,angle_bin),time*qe);
        }
        pmt_dist_start = pmt_dist_min;
      }
      //Photon moves away from vertex
      int max_pmt_dist_bin;
      double pmt_dist_max = TMath::Sqrt(pmt_xEnd * pmt_xEnd + pmt_yEnd * pmt_yEnd + pmt_zEnd * pmt_zEnd);
      if(TMath::Abs(TMath::Abs(phot_zEnd[iphot])-11000)<1)
        max_pmt_dist_bin=pmt_dist_bins;
      else if(TMath::Abs(TMath::Sqrt(phot_xEnd[iphot]*phot_xEnd[iphot]+phot_yEnd[iphot]*phot_yEnd[iphot])-5500)<1)
        max_pmt_dist_bin=pmt_dist_bins;
      else {
        max_pmt_dist_bin = TMath::FloorNint(pmt_dist_max / pmt_dist_binwidth);
      }
      for(int pmt_dist_bin = TMath::CeilNint(pmt_dist_start/pmt_dist_binwidth)+1; pmt_dist_bin<=max_pmt_dist_bin; pmt_dist_bin++){
        double d = (pmt_dist_bin-1)*pmt_dist_binwidth;
        double l = -b+TMath::Sqrt(b2 -a+d*d);
        double angle = TMath::ACos((a0+l*al)/d);
        int angle_bin = hPhotonCount->GetYaxis()->FindFixBin(angle);
        double attenuationProb = (l<phot_dist_max) ? 1 : TMath::Exp(-(l-phot_dist_max)/9000.);
        hPhotonCount->AddBinContent(hPhotonCount->GetBin(pmt_dist_bin,angle_bin),qe*attenuationProb);
        double time = tStart+l*invSpeed;
        hTime->AddBinContent(hTime->GetBin(pmt_dist_bin,angle_bin),time*qe*attenuationProb);
      }
    }
  }
  hTime->SetEntries(1); //needed because hist is filled using AddBinContent
  hTime->Divide(hPhotonCount); //Average arrival time
  hPhotonCount->SetEntries(1); //needed because hist is filled using AddBinContent
  hPhotonCount->Scale(n/(double)count); //Number of photons per event, times by n because only used 1/n of the photons to save time
}

void LikelihoodGenerator::Loop3(int start,int end)
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  if(end <= start) end = nentries;
  //h = new TH2D("h","h",100,0,500,100,angleBins);
//  h = new TH2D("h","h",1000,0,500,1000,0,TMath::Pi());
  hPhotonCount = new TH2D("h2","h2",500,0,2460,500,0,TMath::Pi());
//  h3 = new TH3D("h3","h3",500,0,2460,500,0,TMath::Pi(),1000,0,100);
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=start; jentry<end;jentry++) {
    if(jentry%100==0) std::cout << jentry/100 << " " << std::endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    int ipart=0;
    while(part_parentid[ipart]!=0) ipart++;
    for(int iphot=0; iphot<nphot; iphot++){
      if(phot_wavelength[iphot]>600 || phot_wavelength[iphot]<300) continue;
      if(TMath::Abs(TMath::Abs(phot_zEnd[iphot])-11000)>1
          && TMath::Abs(TMath::Sqrt(phot_xEnd[iphot]*phot_xEnd[iphot]+phot_yEnd[iphot]*phot_yEnd[iphot])-5500)>1) continue;
      double pmt_xDir = phot_xEnd[iphot]-part_xStart[ipart];
      double pmt_yDir = phot_yEnd[iphot]-part_yStart[ipart];
      double pmt_zDir = phot_zEnd[iphot]-part_zStart[ipart];
      double pmt_dist = TMath::Sqrt(pmt_xDir*pmt_xDir+pmt_yDir*pmt_yDir+pmt_zDir*pmt_zDir);
      pmt_xDir/=pmt_dist;
      pmt_yDir/=pmt_dist;
      pmt_zDir/=pmt_dist;
      double pmt_track_angle = TMath::ACos(pmt_xDir*part_pxStart[ipart]+pmt_yDir*part_pyStart[ipart]+pmt_zDir*part_pzStart[ipart]);
      hPhotonCount->Fill(pmt_dist/10.,pmt_track_angle);
    }
  }
}
