
#include "LowEReco.hh"
#include <iostream>
#include <WChRecoLite.hh>
#include <TRandom3.h>

using namespace std;


ClassImp(LowEReco)

LowEReco::LowEReco()
{

}

LowEReco::~LowEReco()
{

}

double LowEReco::ReconstructEnergy(int nphot, const Double_t *phot_xEnd, const Double_t *phot_yEnd, const Double_t *phot_zEnd) {

    const double constant = -2.6913;//used maximum bin value
    const double gradient = 7.66944;

    std::vector<double> posX;
    std::vector<double> posY;
    std::vector<double> posZ;

    posX.clear();
    posY.clear();
    posZ.clear();


    int nphot2=1;
    // std::cout<<"test"<<std::endl;

    posX.push_back(phot_xEnd[0]);
    posY.push_back(phot_yEnd[0]);
    posZ.push_back(phot_zEnd[0]);

    for(int k=1; k<nphot; k++){
        bool notsame = true;
        double xEnd = phot_xEnd[k]; double yEnd = phot_yEnd[k]; double zEnd = phot_zEnd[k];
        // std::cout << posX.size() << std::endl;
        for(int j=0; j<posX.size(); j++){
            if(xEnd==phot_xEnd[j]&&yEnd==phot_yEnd[j]&&zEnd==phot_zEnd[j]){
                notsame=false;
                break;
            }
        }
        if(notsame){
            posX.push_back(xEnd);
            posY.push_back(yEnd);
            posZ.push_back(zEnd);
            // std::cout << posX.size() << std::endl;
            nphot2++;
        }
    }

    double reconstructedE = (nphot2*1.25 - constant)/gradient;//reconstructed E
    std::cout<<"Energy: "<< reconstructedE <<std::endl;
    return reconstructedE;
}

void LowEReco::DoLowEReco(int evt, int nhits, double * hitx, double * hity, double * hitz, double * hitt,
              double allrecoang,
              double & recoVtxX, double & recoVtxY, double & recoVtxZ, double & recoVtxTime,
              double &reconstructedDirX, double &reconstructedDirY, double &reconstructedDirZ){

    // TRandom3 numberp(42484);
    std::vector<double> fDigitX;
    std::vector<double> fDigitY;
    std::vector<double> fDigitZ;
    std::vector<double> fDigitT;
    std::vector<double> fDigitQ;
    std::vector<double> fDigitPE;
    std::vector<double> fDigitW;
    std::vector<double> fDigitV;
    std::vector<double> fDelta; // time residual

    std::vector<double> vSeedVtxX;
    std::vector<double> vSeedVtxY;
    std::vector<double> vSeedVtxZ;
    std::vector<double> vSeedVtxTime;
    std::vector<int> vSeedDigitList;


    const double N_REF=1.34; //average index of refraction
    const double C_VAC=29.9792458; //speed of light in vacuum [cm/ns]


//*********************************************************************

/*

  // Path to WCSim ROOT file
  // =======================

  //  Input should be processed first with photon mask
  // TString filename("/home/hep/wym109/volshome/hksoftware/hk-wchsandbox/analysis/anglerecon/recooutput/reco_e_20MeVvis.root");
  // TString filename("/home/hep/wym109/volshome/hksoftware/hk-wchsandbox/analysis/finalver1/recooutput/reco8inchPMT_0.root");
  TString filename("mu_out.root");
  // TString genfilename("/home/hep/wym109/volshome/hksoftware/hk-wchsandbox/analysis/anglerecon/recooutput/generatorcardfile_e_20MeVvis.root");
  TString genfilename("gun");
  // TString genfilename("generatorcardfile.root");
  TString outfilename("testout3.root");
  // Load Data
  // =========

  // WCLTREEREADER
  //**********************************************************
  // first parameter in the construct: 0=direct output of WCLite
  //                                   1=a "smeared" file made from WCL files
  // second parameter in the constructor: 1 means "include gen level tree"
  //                                      0 means no gen level tree

  WCLTreeReader *mTR = new WCLTreeReader(2,1);

  // Initialize WChRecoLite
*/  WChRecoLite* mRec = WChRecoLite::Instance();/*

  // LoadData is an overloaded function. If you specify only one name,
  // only the main tree is loaded. If you specify two files, it loads
  // both the geant output file and the file containing the gen tree
  // (in that order)

  bool isParticleGun = genfilename == "gun";
  if(isParticleGun)
    mTR->LoadData(filename);
  else
    mTR->LoadData(filename,genfilename);


    // Get positions, timing res of all PMTs
  // cout << "getting PMT info" << endl;
  // int totalPMTs = mTR->GetNpmts();
  // cout << totalPMTs << " PMTs" << endl;
  // double * pmtXall = new double[totalPMTs];
  // double * pmtYall = new double[totalPMTs];
  // double * pmtZall = new double[totalPMTs];
  // hitPMTDirX = new double[totalPMTs];
  // hitPMTDirY = new double[totalPMTs];
  // hitPMTDirZ = new double[totalPMTs];
  // double * pmtTimeResAll = new double[totalPMTs];
  // int * pmtIDall = new int[totalPMTs];
  // for(int i=0; i<totalPMTs; i++){
  //   mTR->LoadPMT(i);
  //   pmtXall[i] = mTR->get_PMTx()/10.; //pmt positions in mm, convert to cm
  //   pmtYall[i] = mTR->get_PMTy()/10.;
  //   pmtZall[i] = mTR->get_PMTz()/10.;
  //   pmtTimeResAll[i] = mTR->get_PMTtimeRes();
  //   pmtIDall[i] = mTR->get_PMTid();
  //   //cout << i << ": " <<  mTR->get_PMTid() << " (" << mTR->get_PMTx() << ", " << mTR->get_PMTy() << ", " << mTR->get_PMTz() <<")" << endl;
  //   if(TMath::Abs(pmtZall[i]-1100)<0.1){
  //     hitPMTDirX[i] = 0;
  //     hitPMTDirY[i] = 0;
  //     hitPMTDirZ[i] = -1;
  //   }
  //   else if(TMath::Abs(pmtZall[i]+1100)<0.1){
  //     hitPMTDirX[i] = 0;
  //     hitPMTDirY[i] = 0;
  //     hitPMTDirZ[i] = 1;
  //   }
  //   else{
  //     hitPMTDirX[i] = -pmtXall[i]/550.;
  //     hitPMTDirY[i] = -pmtYall[i]/550.;
  //     hitPMTDirZ[i] = 0;
  //   }
  // }

  std::cout<<"initializing treewriter"<<std::endl;
  std::cout<<"getting nentries"<<std::endl;

  int nentries = mTR->GetEntries();
  std::cout<<nentries<<std::endl;

  //define reconstrution output here
  // TString fOutputName;
  // fOutputName+="outtry.root";

  TFile f_out(outfilename,"recreate");


  Int_t evt;
  Int_t nphot;
  double recoVtxX;
  double recoVtxY;
  double recoVtxZ;
  double recoVtxTime;

  Double_t trueVtxX;
  Double_t trueVtxY;
  Double_t trueVtxZ;
  Double_t trueVtxTime;
  Double_t trueDirX;
  Double_t trueDirY;
  Double_t trueDirZ;
  Double_t trueKE;

  double DrecoVtxX;
  double DrecoVtxY;
  double DrecoVtxZ;
  double DrecoVtxTime;
  double Drecoabs;
*/
    //double allrecoang;

    // double leptonMom;
    // double leptonMomX;
    // double leptonMomY;
    // double leptonMomZ;

    double totdiffxreco=0;
    double totdiffyreco=0;
    double totdiffzreco=0;
/*
//  double reconstructedmag;
//  double reconstructeddirX;
//  double reconstructeddirY;
//  double reconstructeddirZ;

  TTree* reco_out_ntuple = new TTree("ntuple","ntuple");

  //Output variables
  reco_out_ntuple->Branch("evt",&evt,"evt/I");
  reco_out_ntuple->Branch("nphot",&nphot,"nphot/I");

  reco_out_ntuple->Branch("trueVtxX",&trueVtxX,"trueVtxX/D");
  reco_out_ntuple->Branch("trueVtxY",&trueVtxY,"trueVtxY/D");
  reco_out_ntuple->Branch("trueVtxZ",&trueVtxZ,"trueVtxZ/D");
  reco_out_ntuple->Branch("trueVtxTime",&trueVtxTime,"trueVtxTime/D");

  reco_out_ntuple->Branch("DrecoVtxX",&DrecoVtxX,"DrecoVtxX/D");
  reco_out_ntuple->Branch("DrecoVtxY",&DrecoVtxY,"DrecoVtxY/D");
  reco_out_ntuple->Branch("DrecoVtxZ",&DrecoVtxZ,"DrecoVtxZ/D");
  reco_out_ntuple->Branch("DrecoVtxTime",&DrecoVtxTime,"DrecoVtxTime/D");

  reco_out_ntuple->Branch("Drecoabs",&Drecoabs,"Drecoabs/D");
  reco_out_ntuple->Branch("reconstrucedE",&reconstructedE,"reconstrucedE/D");

  reco_out_ntuple->Branch("reconstruceddirX",&reconstructedDirX,"reconstruceddirX/D");
  reco_out_ntuple->Branch("reconstruceddirY",&reconstructedDirY,"reconstruceddirY/D");
  reco_out_ntuple->Branch("reconstruceddirZ",&reconstructedDirZ,"reconstruceddirZ/D");

  reco_out_ntuple->Branch("allrecoang",&allrecoang,"allrecoang/D");

  for(int i=0; i<nentries; i++){

    evt=i;
*/
    int fThisDigit=0;
    fDigitX.clear();
    fDigitY.clear();
    fDigitZ.clear();
    fDigitT.clear();
    fDigitW.clear();
    fDigitV.clear();
    fDigitQ.clear();
    fDigitPE.clear();
    vSeedDigitList.clear();

    int LEntry[3]={0,0,0};
    // int thisdigit[3]={0};


//    if(i%10==0) cout<<"event: "<<i<<endl;

//    mTR->LoadEvent(i);

    // ebranch->GetEntry(i);

    int nphot = nhits;


    // int npart = mTR->get_npart();

    if(!nphot){
        return;
    }
    // Int_t* phot_hit = mTR->get_phot_hit();
    // Int_t* phot_isScat = mTR->get_phot_isScat();
    // Int_t* phot_processStart = mTR->get_phot_processStart();
    // Int_t* phot_trackid = mTR->get_phot_trackid();
    // Int_t* phot_parentid = mTR->get_phot_parentid();
    // Int_t* phot_capnum = mTR->get_phot_capnum();
    // Int_t* part_processEnd = mTR->get_part_processEnd();
    // Int_t ncapturecount = mTR->get_ncapturecount();
    // Int_t intmode = mTR->get_genmode();
//    Double_t* phot_xStart = mTR->get_phot_xStart();
//    Double_t* phot_yStart = mTR->get_phot_yStart();
//    Double_t* phot_zStart = mTR->get_phot_zStart();
//    Double_t* phot_tStart = mTR->get_phot_tStart();

    Double_t* phot_xEnd = hitx;
    Double_t* phot_yEnd = hity;
    Double_t* phot_zEnd = hitz;
    // Double_t* true_time_corrected_v = mTR->get_phot_tEnd();
    // Double_t* true_time_v = mTR->get_phot_tEnd();
    Double_t* PE_time_v = hitt;
    //Double_t* photon_wavelength_v = mTR->get_phot_wavelength();

    // TRandom3 truer(0);
//    Int_t* parent = mTR->get_part_parentid();
//    int ipart=0;

//    while(parent[ipart]!=0) ipart++;

//    trueVtxX = mTR->get_part_xStart()[ipart]/10.; //convert to cm
//    trueVtxY = mTR->get_part_yStart()[ipart]/10.;
//    trueVtxZ = mTR->get_part_zStart()[ipart]/10.;
    // trueVtxTime = mTR->get_phot_tStart();

//    trueDirX = mTR->get_part_pxStart()[ipart];
//    trueDirY = mTR->get_part_pyStart()[ipart];
//    trueDirZ = mTR->get_part_pzStart()[ipart];
//    trueKE = mTR->get_part_KEstart()[ipart];
//    Int_t* part_pid = mTR->get_part_pid();


//    trueVtxTime = 0;


    //******************************************
    //Reconstruct vertex

    for(int iphot=0;iphot!=nphot;iphot++) //selection cut on tdiff
    {
        // std::cout << phot_parentid[iphot] << std::endl;

//        double t_flight = PE_time_v[iphot];
//        double xE = phot_xEnd[iphot]; double yE = phot_yEnd[iphot]; double zE = phot_zEnd[iphot];
//        double distSquare = xE*xE+yE*yE+zE*zE;

//        if( t_flight < (sqrt(distSquare)/(300./N_REF))+50.
//            && t_flight > (sqrt(distSquare)/(300./N_REF))-5. && t_flight > 0. )
//        {

            fDigitX.push_back(phot_xEnd[iphot]/10.); //!Sphere1
            fDigitY.push_back(phot_yEnd[iphot]/10.);//!Sphere1
            fDigitZ.push_back(phot_zEnd[iphot]/10.);//!Sphere1
            fDigitT.push_back(PE_time_v[iphot]);// - min_PE_time; //!Sphere1
            fDigitQ.push_back(1);
            fDigitW.push_back(/*photon_wavelength_v[iphot]*/500);
            fDigitV.push_back(C_VAC/N_REF);
            vSeedDigitList.push_back(fThisDigit);
            fThisDigit++;

//        }
    }


    std::cout<<"Photon filtering for event #"<<evt<<" has just finished. fDigits are ready."<<std::endl;

    //input digits
    mRec->SetDigits(fDigitX, fDigitY, fDigitZ, fDigitT, fDigitQ, fDigitPE, fDigitW, fDigitV, fDelta, vSeedDigitList);
    std::cout << "vSeedDigitList:" << vSeedDigitList.size() << std::endl;

    if(vSeedDigitList.size()==0) return;

    //calculate seed vertices
    mRec->CalcVertexSeeds();
    vSeedVtxX = mRec->GetSeedVtx(0);
    vSeedVtxY = mRec->GetSeedVtx(1);
    vSeedVtxZ = mRec->GetSeedVtx(2);
    vSeedVtxTime = mRec->GetSeedVtx(3);


    //now select the best vertex and save reco data
    int best_seed = mRec->SelectBestSeed();
    std::cout<<"best_seed = "<<best_seed<<"  vSeedVtxX.size() = "<<vSeedVtxX.size()<<std::endl;
    if(best_seed == -1) return;


    recoVtxX = vSeedVtxX[best_seed];
    recoVtxY = vSeedVtxY[best_seed];
    recoVtxZ = vSeedVtxZ[best_seed];
    recoVtxTime = vSeedVtxTime[best_seed];
    /*

    //Get mean as starting point for Amoeba algorithm
    int count = vSeedVtxX.size();
//    for(int iv =0; iv< count; iv++){
//        cout << vSeedVtxX[iv] << ", " << vSeedVtxY[iv] << ", " << vSeedVtxZ[iv] << ", " << vSeedVtxTime[iv] << endl;
//    }
    double meanVertex[4];
    for(int idir=0; idir<4; idir++) meanVertex[idir]=mRec->MeanVertexSeed(idir);
    cout << "Mean vtx: " << meanVertex[0] << ", " << meanVertex[1] << ", " << meanVertex[2] << ", " << meanVertex[3] << endl;
    //Get starting simplex and densities
    double p[5][4];
    mRec->GetP(p);
    double y[5];
    for(int iy=0; iy<5; iy++) y[iy]=mRec->GetSumVtx(p[iy]);
//    for(int iy=0; iy<5; iy++) {
//        cout << p[iy][0] << ", " << p[iy][1] << ", " << p[iy][2] << ", " << p[iy][3] << " : " << y[iy] << endl;
//    }
    //Run Amoeba algorithm
    int nfunc;
    mRec->Amoeba(p,y,4,0.00001,nfunc);
    recoVtxX = p[0][0];
    recoVtxY = p[0][1];
    recoVtxZ = p[0][2];
    recoVtxTime = p[0][3];
*/
/*
    DrecoVtxX = recoVtxX-trueVtxX;
    DrecoVtxY = recoVtxY-trueVtxY;
    DrecoVtxZ = recoVtxZ-trueVtxZ;
    DrecoVtxTime = recoVtxTime-trueVtxTime;
    Drecoabs = TMath::Sqrt((DrecoVtxX*DrecoVtxX)+(DrecoVtxY*DrecoVtxY)+(DrecoVtxZ*DrecoVtxZ));
*/
    //******************************************************
    //Reconstruct Direction
    for(int kphot=0;kphot!=nphot;kphot++) //selection cut on tdiff
    {
        double t_flight2 = PE_time_v[kphot];
        double xE2 = phot_xEnd[kphot]; double yE2 = phot_yEnd[kphot]; double zE2 = phot_zEnd[kphot];
        double distSquare2 = xE2*xE2+yE2*yE2+zE2*zE2;

        if( t_flight2 < (sqrt(distSquare2)/(300./N_REF))+50.
            && t_flight2 > (sqrt(distSquare2)/(300./N_REF))-5. && t_flight2 > 0. )
        {

            double trackmagreco = sqrt(((phot_xEnd[kphot]/10.)-recoVtxX)*((phot_xEnd[kphot]/10.)-recoVtxX)+((phot_yEnd[kphot]/10.)-recoVtxY)*((phot_yEnd[kphot]/10.)-recoVtxY)+((phot_zEnd[kphot]/10.)-recoVtxZ)*((phot_zEnd[kphot]/10.)-recoVtxZ));
            double diffxreco = ((phot_xEnd[kphot]/10.)-recoVtxX)/trackmagreco;
            double diffyreco = ((phot_yEnd[kphot]/10.)-recoVtxY)/trackmagreco;
            double diffzreco = ((phot_zEnd[kphot]/10.)-recoVtxZ)/trackmagreco;
            // std::cout<< "diffx: " << diffx <<std::endl;
            totdiffxreco += diffxreco;
            totdiffyreco += diffyreco;
            totdiffzreco += diffzreco;
        }
        // std::cout<< (phot_zEnd[kphot]/10.) <<std::endl;
    }

    double reconstructedMag = sqrt(totdiffxreco*totdiffxreco+totdiffyreco*totdiffyreco+totdiffzreco*totdiffzreco);
    reconstructedDirX = totdiffxreco/ reconstructedMag;
    reconstructedDirY = totdiffyreco/ reconstructedMag;
    reconstructedDirZ = totdiffzreco/ reconstructedMag;
    //******************************************************
    //Reconstruct Cherenkov Angle


    double meancrkAngle=0.;
    int countt=0;

    for(int counttriplet=0;counttriplet<1000;counttriplet++){//loop each 3 hit

        double xpos[3];
        double ypos[3];
        double zpos[3];

        TRandom3 RND((UInt_t) counttriplet);

        Double_t r0 = RND.Rndm();
        Double_t r1 = RND.Rndm();
        Double_t r2 = RND.Rndm();

        // Int_t numEntries = vSeedDigitList.size();
        LEntry[0] = TMath::FloorNint(r0*vSeedDigitList.size());
        LEntry[1] = TMath::FloorNint(r1*vSeedDigitList.size());
        LEntry[2] = TMath::FloorNint(r2*vSeedDigitList.size());
        // std::cout << LEntry0 << "  " << LEntry1 << "  " << LEntry2 << std::endl;
        if(LEntry[0]<0.||LEntry[0]>vSeedDigitList.size()) continue;//just to check
        if(LEntry[1]<0.||LEntry[1]>vSeedDigitList.size()) continue;
        if(LEntry[2]<0.||LEntry[2]>vSeedDigitList.size()) continue;

        //check not to choose same hit
        // if(thisdigit[0]==thisdigit[1]||thisdigit[0]==thisdigit[2]||thisdigit[1]==thisdigit[2]) continue;

        xpos[0] = fDigitX[LEntry[0]];
        ypos[0] = fDigitY[LEntry[0]];
        zpos[0] = fDigitZ[LEntry[0]];
        // time0 = fDigitT[LEntry0];
        xpos[1] = fDigitX[LEntry[1]];
        ypos[1] = fDigitY[LEntry[1]];
        zpos[1] = fDigitZ[LEntry[1]];
        // time1 = fDigitT[LEntry1];
        xpos[2] = fDigitX[LEntry[2]];
        ypos[2] = fDigitY[LEntry[2]];
        zpos[2] = fDigitZ[LEntry[2]];
        // time2 = fDigitT[thisdigit2];

        //check not to choose same hit
        if(xpos[0]==xpos[1]&&ypos[0]==ypos[1]&&zpos[0]==zpos[1]) continue;
        if(xpos[0]==xpos[2]&&ypos[0]==ypos[2]&&zpos[0]==zpos[2]) continue;
        if(xpos[1]==xpos[2]&&ypos[1]==ypos[2]&&zpos[1]==zpos[2]) continue;
        // std::cout << xpos[0] << "  " << ypos[0] << "  " << zpos[0] << std::endl;
        // std::cout << xpos[1] << "  " << ypos[1] << "  " << zpos[1] << std::endl;
        // std::cout << xpos[2] << "  " << ypos[2] << "  " << zpos[2] << std::endl;

        double trackmagreco[3];
        trackmagreco[0] = sqrt((xpos[0]-recoVtxX)*(xpos[0]-recoVtxX)+(ypos[0]-recoVtxY)*(ypos[0]-recoVtxY)+(zpos[0]-recoVtxZ)*(zpos[0]-recoVtxZ));
        trackmagreco[1] = sqrt((xpos[1]-recoVtxX)*(xpos[1]-recoVtxX)+(ypos[1]-recoVtxY)*(ypos[1]-recoVtxY)+(zpos[1]-recoVtxZ)*(zpos[1]-recoVtxZ));
        trackmagreco[2] = sqrt((xpos[2]-recoVtxX)*(xpos[2]-recoVtxX)+(ypos[2]-recoVtxY)*(ypos[2]-recoVtxY)+(zpos[2]-recoVtxZ)*(zpos[2]-recoVtxZ));
        // std::cout << "MAG: " << trackmagreco0 << "  " << trackmagreco1 << "  " << trackmagreco2 << std::endl;

        double unitvecX[3];
        double unitvecY[3];
        double unitvecZ[3];
        unitvecX[0] = (xpos[0]-recoVtxX)/trackmagreco[0];
        unitvecY[0] = (ypos[0]-recoVtxY)/trackmagreco[0];
        unitvecZ[0] = (zpos[0]-recoVtxZ)/trackmagreco[0];
        unitvecX[1] = (xpos[1]-recoVtxX)/trackmagreco[1];
        unitvecY[1] = (ypos[1]-recoVtxY)/trackmagreco[1];
        unitvecZ[1] = (zpos[1]-recoVtxZ)/trackmagreco[1];
        unitvecX[2] = (xpos[2]-recoVtxX)/trackmagreco[2];
        unitvecY[2] = (ypos[2]-recoVtxY)/trackmagreco[2];
        unitvecZ[2] = (zpos[2]-recoVtxZ)/trackmagreco[2];
        // std::cout << sqrt(unitvecX0*unitvecX0+unitvecY0*unitvecY0+unitvecZ0*unitvecZ0) << std::endl;
        // std::cout << sqrt(unitvecX1*unitvecX1+unitvecY1*unitvecY1+unitvecZ1*unitvecZ1) << std::endl;
        // std::cout << sqrt(unitvecX2*unitvecX2+unitvecY2*unitvecY2+unitvecZ2*unitvecZ2) << std::endl;
        double unitvecXcom = (unitvecX[0]+unitvecX[1]+unitvecX[2])/3.;
        double unitvecYcom = (unitvecY[0]+unitvecY[1]+unitvecY[2])/3.;
        double unitvecZcom = (unitvecZ[0]+unitvecZ[1]+unitvecZ[2])/3.;

        double xd = unitvecXcom-unitvecX[0];
        double yd = unitvecYcom-unitvecY[0];
        double zd = unitvecZcom-unitvecZ[0];
        // std::cout << xd << "  " << yd << "  " << zd << std::endl;

        double radiussq = sqrt(xd*xd+yd*yd+zd*zd);
        if(radiussq>1.) continue;//check
        // std::cout << "R: " << radiussq << std::endl;

        double crkAngle = TMath::ASin(radiussq);
        // recoangle->Fill(crkAngle/3.14159*180.);
        meancrkAngle+=crkAngle;

        // std::cout << "angle: " << crkAngle/3.14159*180. << std::endl;
        // std::cout << "---------------------------" << std::endl;
        countt++;
        // counttriplet++;

    }
    if(countt==0) return;//check
    allrecoang = (meancrkAngle/3.14159*180.)/countt;
    // recoangleall->Fill(allrecoang);

    // std::cout << recoangle->GetMean() << std::endl;
    // TF1 *myfit = new TF1("myfit","x*gaus(0)", 0, 90);
    // myfit->SetParameters(90,recoangle->GetMean(),recoangle->GetRMS());
    // myfit->SetParNames("Constant","Mean_value","Sigma");
    // recoangle->Fit("myfit");
    // recoangleall->Fill(recoangle->GetMean());

//  std::cout<<"TV!!: "<<trueVtxX<<" "<<trueVtxY<<" "<<trueVtxZ<<std::endl;
    std::cout<<"RECO: "<<recoVtxX<<" "<<recoVtxY<<" "<<recoVtxZ<<std::endl;
    std::cout<<"Dir: "<< reconstructedDirX <<" "<< reconstructedDirY <<" "<< reconstructedDirZ <<std::endl;
    std::cout<<"Angle: "<<allrecoang<<std::endl;

    std::cout<<"-------------------------"<<std::endl;
//  reco_out_ntuple->Fill();
//  } //end i-loop over Hits_Tree entries
//  f_out.cd();
//  reco_out_ntuple->Write();

}

