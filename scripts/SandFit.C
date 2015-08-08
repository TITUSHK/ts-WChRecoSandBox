#include "TRandom3.h"
#include <iomanip>
#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TH2D.h"
#include "TGraph2D.h"
#include "TH3D.h"
#include "TSpectrum2.h"
#include "../include/WCLTreeReader.hh"
#include "../include/HighEReco.hh"
#include "../include/LowEReco.hh"

using namespace std;

// Performs reconstruction assuming single ring high-E event
//
// filename: root output file from $WCHSANDBOX/analysis/advanced_scripts/SandboxDetSim.C
// genfilename: generatorcardfile.root file from WChSandbox
// outfilename: file to write output
// lnLfilename: file with likelihood tables
//

int* FindClusters(int nHits, double * hitTimes, int * hitCluster, int & nClusters){
    cout << "Starting cluster search" << endl;
    nClusters = 0;
    int *hitOrder = new int[nHits];
    TMath::Sort(nHits, hitTimes, hitOrder, kFALSE);
    cout << "Sorted " << nHits << " hits" << endl;
    //Loop over hits
    for(int clusterStart=0; clusterStart+9<nHits; clusterStart++){
        int clusterEnd = clusterStart+9;
        double clusterStartTime = hitTimes[hitOrder[clusterStart]];
        //Check at least 10 hits in 200ns window
        if(hitTimes[hitOrder[clusterEnd]]-clusterStartTime>200)
            continue;
        //Find cluster end time, with hit spacing of less than 10ns, allowing for at most two spacings of 10-20ns
        int nSpacings = 0;
        double clusterEndTime = clusterStartTime;
        for(clusterEnd=clusterStart; clusterEnd<nHits-1; clusterEnd++)
        {
            clusterEndTime = hitTimes[hitOrder[clusterEnd]];
            double nextHitTime = hitTimes[hitOrder[clusterEnd+1]];
            if(nextHitTime - clusterEndTime < 10)
                continue;
            if (nextHitTime - clusterEndTime < 20 && nSpacings < 2) {
                nSpacings++;
                continue;
            }
            break;
        }
        //Check at least 10 hits in cluster
        if(clusterEnd-clusterStart<9)
            continue;
        //Found a cluster
        nClusters++;
        cout << "Cluster " << nClusters << " found with " << clusterEnd-clusterStart+1 << "hits: "
             << clusterStartTime << "ns to " << clusterEndTime << "ns" << endl;
        //Set cluster number for hits in cluster
        for(int iHit = clusterStart; iHit<=clusterEnd; iHit++){
            hitCluster[hitOrder[iHit]] = nClusters;
        }
        //Look for next cluster
        clusterStart = clusterEnd;
    }
    //Count hits in each cluster
    int * clusterHitCounts = new int[nClusters+1](); for(int i=0; i<nClusters+1; i++) clusterHitCounts[i]=0;
    for(int iHit=0; iHit<nHits; iHit++) {
        clusterHitCounts[hitCluster[iHit]]++;
    }
    delete[] hitOrder;
    return clusterHitCounts;
}

void SandFit(TString filename = "out.root",
             TString genfilename = "generatorcardfile.root",
             TString outfilename = "testout2.root",
             TString lnLfilename = "likelihood_tables.root",
             TString lookupfilename = "energyLookups.root",
             int startEvent=0, int maxEvents = 0)
{

    // Path to WCSim ROOT file
    // =======================

    // Load Data
    // =========

    // WCLTREEREADER
    //**********************************************************
    // first parameter in the construct: 0=direct output of WCLite
    //                                   1=a "smeared" file made from WCL files
    //                                   2=output from SandboxDetSim
    // second parameter in the constructor: 1 means "include gen level tree"
    //                                      0 means no gen level tree

    WCLTreeReader *mTR = new WCLTreeReader(2,1);

    // LoadData is an overloaded function. If you specify only one name,
    // only the main tree is loaded. If you specify two files, it loads
    // both the geant output file and the file containing the gen tree
    // (in that order)
    cout<<"Loading data" << endl;
    bool isParticleGun = false;
    if(filename=="prod") {
        TString flavs[4] = {"numu", "nue", "antinumu", "antinue"};
        TString horns[2] = {"nu", "antinu"};
        for (int iflav = 0; iflav < 2; iflav++) {
            const char *flav = flavs[iflav].Data();
            for (int ihorn = 0; ihorn < 1; ihorn++) {
                const char *horn = horns[ihorn].Data();
                for (int i = 1000; i < 1001; i++) {
                    filename = Form("/data/hyperk/wchsandbox_reco/flav_%s/%s_%s_%i/%s_%s_%i_out_12in.root",flav,horn,flav,i,horn,flav,i);
                    genfilename = Form("/data/hyperk/wchsandbox_reco/flav_%s/%s_%s_%i/generatorcardfile.root",flav,horn,flav,i);
                    mTR->LoadData(filename,genfilename);
                }
            }
        }
    }
    else{
        if (genfilename == "gun") {
            isParticleGun = true;
            mTR->LoadData(filename);
        }
        else
            mTR->LoadData(filename, genfilename);
    }

    TFile *likelihoodTables = new TFile(lnLfilename, "READ");
    TH3D *electronPhotons = (TH3D *) likelihoodTables->Get("photons_e");
    TH3D *muonPhotons = (TH3D *) likelihoodTables->Get("photons_mu");
    TH3D *electronTimes = (TH3D *) likelihoodTables->Get("times_e");
    TH3D *muonTimes = (TH3D *) likelihoodTables->Get("times_mu");

    TFile *energyTables = new TFile(lookupfilename, "READ");
    TH2D *electronEnergyLookup = (TH2D *) energyTables->Get("nHitDWallKELookup_e");
    TH2D *muonEnergyLookup = (TH2D *) energyTables->Get("nHitDWallKELookup_mu");
    TH2D *electronVtxBiasLookup = (TH2D *) energyTables->Get("nHitDWallVtxTrackBias_e");
    TH2D *muonVtxBiasLookup = (TH2D *) energyTables->Get("nHitDWallVtxTrackBias_mu");

    // Get positions, timing res of all PMTs
    cout << "getting PMT info" << endl;
    int totalPMTs = mTR->GetNpmts();
    cout << totalPMTs << " PMTs" << endl;
    double * pmtXall = new double[totalPMTs];
    double * pmtYall = new double[totalPMTs];
    double * pmtZall = new double[totalPMTs];
    double * pmtDirXall = new double[totalPMTs];
    double * pmtDirYall = new double[totalPMTs];
    double * pmtDirZall = new double[totalPMTs];
    double * pmtTimeResAll = new double[totalPMTs];
    int * pmtIDall = new int[totalPMTs];
    for(int i=0; i<totalPMTs; i++){
        mTR->LoadPMT(i);
        pmtXall[i] = mTR->get_PMTx()/10.; //pmt positions in mm, convert to cm
        pmtYall[i] = mTR->get_PMTy()/10.;
        pmtZall[i] = mTR->get_PMTz()/10.;
        pmtTimeResAll[i] = mTR->get_PMTtimeRes();
        pmtIDall[i] = mTR->get_PMTid();
        if(TMath::Abs(pmtZall[i]-1100)<0.1){
            pmtDirXall[i] = 0;
            pmtDirYall[i] = 0;
            pmtDirZall[i] = -1;
        }
        else if(TMath::Abs(pmtZall[i]+1100)<0.1){
            pmtDirXall[i] = 0;
            pmtDirYall[i] = 0;
            pmtDirZall[i] = 1;
        }
        else{
            pmtDirXall[i] = -pmtXall[i]/550.;
            pmtDirYall[i] = -pmtYall[i]/550.;
            pmtDirZall[i] = 0;
        }
    }
    cout << "getting nEntries" << endl;

    int nEntries = mTR->GetEntries();
    cout<<nEntries<<endl;

    //define reconstruction output here

    TFile f_out(outfilename,"recreate");

    int evt;

    TTree* DebugTree = new TTree("Debug","Debug");
    TTree* LowETree = new TTree("Low_E","Low_E");
    TTree* HighEElectronTree = new TTree("High_E_Electron","High_E_Electron");
    TTree* HighEMuonTree = new TTree("High_E_Muon","High_E_Muon");
    TTree* FinalTree = new TTree("Final_Reconstruction","Final_Reconstruction");

    const int maxSubEvts = 20;
    int cluster[maxSubEvts], ring[maxSubEvts];
    int nSubevents, nClusters;

    double recoVtxXLowE[maxSubEvts];
    double recoVtxYLowE[maxSubEvts];
    double recoVtxZLowE[maxSubEvts];
    double recoTimeLowE[maxSubEvts];
    double recoDirXLowE[maxSubEvts];
    double recoDirYLowE[maxSubEvts];
    double recoDirZLowE[maxSubEvts];
    double recoChkvAngleLowE[maxSubEvts];
    double recoEnergyLowE[maxSubEvts];
    LowETree->Branch("evt",&evt,"evt/I");
    LowETree->Branch("nClusters",&nClusters,"nClusters/I");
    LowETree->Branch("nSubevents",&nSubevents,"nSubevents/I");
    LowETree->Branch("cluster", cluster,"cluster[nClusters]/I");
    LowETree->Branch("recoVtxXLowE",recoVtxXLowE,"recoVtxXLowE[nClusters]/D");
    LowETree->Branch("recoVtxYLowE",recoVtxYLowE,"recoVtxYLowE[nClusters]/D");
    LowETree->Branch("recoVtxZLowE",recoVtxZLowE,"recoVtxZLowE[nClusters]/D");
    LowETree->Branch("recoTimeLowE",recoTimeLowE,"recoTimeLowE[nClusters]/D");
    LowETree->Branch("recoDirXLowE",recoDirXLowE,"recoDirXLowE[nClusters]/D");
    LowETree->Branch("recoDirYLowE",recoDirYLowE,"recoDirYLowE[nClusters]/D");
    LowETree->Branch("recoDirZLowE",recoDirZLowE,"recoDirZLowE[nClusters]/D");
    LowETree->Branch("recoChkvAngleLowE",recoChkvAngleLowE,"recoChkvAngleLowE[nClusters]/D");
    LowETree->Branch("recoEnergyLowE",recoEnergyLowE,"recoEnergyLowE[nClusters]/D");

    double recoVtxXHighEElectron[maxSubEvts];
    double recoVtxYHighEElectron[maxSubEvts];
    double recoVtxZHighEElectron[maxSubEvts];
    double recoTimeHighEElectron[maxSubEvts];
    double recoDirXHighEElectron[maxSubEvts];
    double recoDirYHighEElectron[maxSubEvts];
    double recoDirZHighEElectron[maxSubEvts];
    double recoChkvAngleHighEElectron[maxSubEvts];
    double recoEnergyHighEElectron[maxSubEvts];
    double recoLnLHighEElectron[maxSubEvts];
    HighEElectronTree->Branch("evt",&evt,"evt/I");
    HighEElectronTree->Branch("nClusters",&nClusters,"nClusters/I");
    HighEElectronTree->Branch("nSubevents",&nSubevents,"nSubevents/I");
    HighEElectronTree->Branch("cluster", cluster,"cluster[nSubevents]/I");
    HighEElectronTree->Branch("ring",ring,"ring[nSubevents]/I");
    HighEElectronTree->Branch("recoVtxXHighEElectron",recoVtxXHighEElectron,"recoVtxXHighEElectron[nSubevents]/D");
    HighEElectronTree->Branch("recoVtxYHighEElectron",recoVtxYHighEElectron,"recoVtxYHighEElectron[nSubevents]/D");
    HighEElectronTree->Branch("recoVtxZHighEElectron",recoVtxZHighEElectron,"recoVtxZHighEElectron[nSubevents]/D");
    HighEElectronTree->Branch("recoTimeHighEElectron",recoTimeHighEElectron,"recoTimeHighEElectron[nSubevents]/D");
    HighEElectronTree->Branch("recoDirXHighEElectron",recoDirXHighEElectron,"recoDirXHighEElectron[nSubevents]/D");
    HighEElectronTree->Branch("recoDirYHighEElectron",recoDirYHighEElectron,"recoDirYHighEElectron[nSubevents]/D");
    HighEElectronTree->Branch("recoDirZHighEElectron",recoDirZHighEElectron,"recoDirZHighEElectron[nSubevents]/D");
    HighEElectronTree->Branch("recoChkvAngleHighEElectron",recoChkvAngleHighEElectron,"recoChkvAngleHighEElectron[nSubevents]/D");
    HighEElectronTree->Branch("recoEnergyHighEElectron",recoEnergyHighEElectron,"recoEnergyHighEElectron[nSubevents]/D");
    HighEElectronTree->Branch("recoLnLHighEElectron",recoLnLHighEElectron,"recoLnLHighEElectron[nSubevents]/D");

    double recoVtxXHighEMuon[maxSubEvts];
    double recoVtxYHighEMuon[maxSubEvts];
    double recoVtxZHighEMuon[maxSubEvts];
    double recoTimeHighEMuon[maxSubEvts];
    double recoDirXHighEMuon[maxSubEvts];
    double recoDirYHighEMuon[maxSubEvts];
    double recoDirZHighEMuon[maxSubEvts];
    double recoChkvAngleHighEMuon[maxSubEvts];
    double recoEnergyHighEMuon[maxSubEvts];
    double recoLnLHighEMuon[maxSubEvts];
    HighEMuonTree->Branch("evt",&evt,"evt/I");
    HighEMuonTree->Branch("nClusters",&nClusters,"nClusters/I");
    HighEMuonTree->Branch("nSubevents",&nSubevents,"nSubevents/I");
    HighEMuonTree->Branch("cluster", cluster,"cluster[nSubevents]/I");
    HighEMuonTree->Branch("ring",ring,"ring[nSubevents]/I");
    HighEMuonTree->Branch("recoVtxXHighEMuon",recoVtxXHighEMuon,"recoVtxXHighEMuon[nSubevents]/D");
    HighEMuonTree->Branch("recoVtxYHighEMuon",recoVtxYHighEMuon,"recoVtxYHighEMuon[nSubevents]/D");
    HighEMuonTree->Branch("recoVtxZHighEMuon",recoVtxZHighEMuon,"recoVtxZHighEMuon[nSubevents]/D");
    HighEMuonTree->Branch("recoTimeHighEMuon",recoTimeHighEMuon,"recoTimeHighEMuon[nSubevents]/D");
    HighEMuonTree->Branch("recoDirXHighEMuon",recoDirXHighEMuon,"recoDirXHighEMuon[nSubevents]/D");
    HighEMuonTree->Branch("recoDirYHighEMuon",recoDirYHighEMuon,"recoDirYHighEMuon[nSubevents]/D");
    HighEMuonTree->Branch("recoDirZHighEMuon",recoDirZHighEMuon,"recoDirZHighEMuon[nSubevents]/D");
    HighEMuonTree->Branch("recoChkvAngleHighEMuon",recoChkvAngleHighEMuon,"recoChkvAngleHighEMuon[nSubevents]/D");
    HighEMuonTree->Branch("recoEnergyHighEMuon",recoEnergyHighEMuon,"recoEnergyHighEMuon[nSubevents]/D");
    HighEMuonTree->Branch("recoLnLHighEMuon",recoLnLHighEMuon,"recoLnLHighEMuon[nSubevents]/D");

    int recoNRings[maxSubEvts];
    double recoVtxX[maxSubEvts];
    double recoVtxY[maxSubEvts];
    double recoVtxZ[maxSubEvts];
    double recoTime[maxSubEvts];
    double recoDirX[maxSubEvts];
    double recoDirY[maxSubEvts];
    double recoDirZ[maxSubEvts];
    double recoChkvAngle[maxSubEvts];
    double recoEnergy[maxSubEvts];
    int recoPID[maxSubEvts];
    FinalTree->Branch("evt",&evt,"evt/I");
    FinalTree->Branch("nClusters",&nClusters,"nClusters/I");
    FinalTree->Branch("nSubevents",&nSubevents,"nSubevents/I");
    FinalTree->Branch("recoNRings",recoNRings,"recoNRings[nClusters]/I");
    FinalTree->Branch("cluster", cluster,"cluster[nSubevents]/I");
    FinalTree->Branch("ring",ring,"ring[nSubevents]/I");
    FinalTree->Branch("recoVtxX",recoVtxX,"recoVtxX[nSubevents]/D");
    FinalTree->Branch("recoVtxY",recoVtxY,"recoVtxY[nSubevents]/D");
    FinalTree->Branch("recoVtxZ",recoVtxZ,"recoVtxZ[nSubevents]/D");
    FinalTree->Branch("recoTime",recoTime,"recoTime[nSubevents]/D");
    FinalTree->Branch("recoDirX",recoDirX,"recoDirX[nSubevents]/D");
    FinalTree->Branch("recoDirY",recoDirY,"recoDirY[nSubevents]/D");
    FinalTree->Branch("recoDirZ",recoDirZ,"recoDirZ[nSubevents]/D");
    FinalTree->Branch("recoChkvAngle",recoChkvAngle,"recoChkvAngle[nSubevents]/D");
    FinalTree->Branch("recoEnergy",recoEnergy,"recoEnergy[nSubevents]/D");
    FinalTree->Branch("recoPID",recoPID,"recoPID[nSubevents]/I");

    bool isHighE[maxSubEvts];
    int ringPEs[maxSubEvts];
    double trueVtxX;
    double trueVtxY;
    double trueVtxZ;
    double trueTime;
    double trueDirX;
    double trueDirY;
    double trueDirZ;
    double trueKE;
    double diffKE;
    double diffVtxX;
    double diffVtxY;
    double diffVtxZ;
    double diffTime;
    double diffDirX;
    double diffDirY;
    double diffDirZ;
    double diffVtxAbs;
    double diffDirAbs;
    double pointVtxX[maxSubEvts];
    double pointVtxY[maxSubEvts];
    double pointVtxZ[maxSubEvts];
    double pointTime[maxSubEvts];
    double trackVtxX[maxSubEvts];
    double trackVtxY[maxSubEvts];
    double trackVtxZ[maxSubEvts];
    double trackTime[maxSubEvts];
    double trackDirX[maxSubEvts];
    double trackDirY[maxSubEvts];
    double trackDirZ[maxSubEvts];
    double electronTrackCorrection[maxSubEvts];
    double muonTrackCorrection[maxSubEvts];
    int mode;
    double neutrinoE;
    double neutrinoDirX;
    double neutrinoDirY;
    double neutrinoDirZ;
    int neutrinoPID;
    int nNeutrons;
    int neutronCount;
    int nCaptures;
    const int maxPMT = 20000;
    //double trueExpectedPEHighEMuon[maxSubEvts][maxPMT];
    //double trueExpectedPEHighEElectron[maxSubEvts][maxPMT];
    //double recoExpectedPEHighEMuon[maxSubEvts][maxPMT];
    //double recoExpectedPEHighEElectron[maxSubEvts][maxPMT];
    int observedPE[maxSubEvts][maxPMT];
    int allPE[maxSubEvts][maxPMT];
    //int ring1PE[maxSubEvts][maxPMT];
    //int ring2PE[maxSubEvts][maxPMT];
    double trueToWall;
    double trueDWall;
    double recoToWall;
    double recoDWall;
    int recoCaptures;
    DebugTree->Branch("evt",&evt,"evt/I");
    DebugTree->Branch("nClusters",&nClusters,"nClusters/I");
    DebugTree->Branch("nSubevents",&nSubevents,"nSubevents/I");
    DebugTree->Branch("cluster", cluster,"cluster[nSubevents]/I");
    DebugTree->Branch("ring",ring,"ring[nSubevents]/I");
    DebugTree->Branch("isHighE",isHighE,"isHighE[nSubevents]/O");
    DebugTree->Branch("ringPEs", ringPEs,"ringPEs[nSubevents]/I");
    DebugTree->Branch("trueVtxX",&trueVtxX,"trueVtxX/D");
    DebugTree->Branch("trueVtxY",&trueVtxY,"trueVtxY/D");
    DebugTree->Branch("trueVtxZ",&trueVtxZ,"trueVtxZ/D");
    DebugTree->Branch("trueTime",&trueTime,"trueTime/D");
    DebugTree->Branch("trueDirX",&trueDirX,"trueDirX/D");
    DebugTree->Branch("trueDirY",&trueDirY,"trueDirY/D");
    DebugTree->Branch("trueDirZ",&trueDirZ,"trueDirZ/D");
    DebugTree->Branch("trueKE",&trueKE,"trueKE/D");
    DebugTree->Branch("diffKE",&diffKE,"diffKE/D");
    DebugTree->Branch("diffVtxX",&diffVtxX,"diffVtxX/D");
    DebugTree->Branch("diffVtxY",&diffVtxY,"diffVtxY/D");
    DebugTree->Branch("diffVtxZ",&diffVtxZ,"diffVtxZ/D");
    DebugTree->Branch("diffTime",&diffTime,"diffTime/D");
    DebugTree->Branch("diffDirX",&diffDirX,"diffDirX/D");
    DebugTree->Branch("diffDirY",&diffDirY,"diffDirY/D");
    DebugTree->Branch("diffDirZ",&diffDirZ,"diffDirZ/D");
    DebugTree->Branch("diffVtxAbs",&diffVtxAbs,"diffVtxAbs/D");
    DebugTree->Branch("diffDirAbs",&diffDirAbs,"diffDirAbs/D");
    DebugTree->Branch("pointVtxX",pointVtxX,"pointVtxX[nClusters]/D");
    DebugTree->Branch("pointVtxY",pointVtxY,"pointVtxY[nClusters]/D");
    DebugTree->Branch("pointVtxZ",pointVtxZ,"pointVtxZ[nClusters]/D");
    DebugTree->Branch("pointTime",pointTime,"pointTime[nClusters]/D");
    DebugTree->Branch("trackVtxX",trackVtxX,"trackVtxX[nSubevents]/D");
    DebugTree->Branch("trackVtxY",trackVtxY,"trackVtxY[nSubevents]/D");
    DebugTree->Branch("trackVtxZ",trackVtxZ,"trackVtxZ[nSubevents]/D");
    DebugTree->Branch("trackTime",trackTime,"trackTime[nSubevents]/D");
    DebugTree->Branch("trackDirX",trackDirX,"trackDirX[nSubevents]/D");
    DebugTree->Branch("trackDirY",trackDirY,"trackDirY[nSubevents]/D");
    DebugTree->Branch("trackDirZ",trackDirZ,"trackDirZ[nSubevents]/D");
    DebugTree->Branch("electronTrackCorrection",electronTrackCorrection,"electronTrackCorrection[nSubevents]/D");
    DebugTree->Branch("muonTrackCorrection",muonTrackCorrection,"muonTrackCorrection[nSubevents]/D");
    DebugTree->Branch("mode",&mode,"mode/I");
    DebugTree->Branch("neutrinoE",&neutrinoE,"neutrinoE/D");
    DebugTree->Branch("neutrinoDirX",&neutrinoDirX,"neutrinoDirX/D");
    DebugTree->Branch("neutrinoDirY",&neutrinoDirY,"neutrinoDirY/D");
    DebugTree->Branch("neutrinoDirZ",&neutrinoDirZ,"neutrinoDirZ/D");
    DebugTree->Branch("neutrinoPID",&neutrinoPID,"neutrinoPID/I");
    DebugTree->Branch("nNeutrons",&nNeutrons,"nNeutrons/I");
    DebugTree->Branch("neutronCount",&neutronCount,"neutronCount/I");
    DebugTree->Branch("nCaptures",&nCaptures,"nCaptures/I");
    DebugTree->Branch("nPMTs",&totalPMTs,"nPMTs/I");
    DebugTree->Branch("pmtID",pmtIDall,Form("pmtID[%i]/I",totalPMTs));
    DebugTree->Branch("pmtXtt",pmtXall,Form("pmtX[%i]/D",totalPMTs));
    DebugTree->Branch("pmtY",pmtYall,Form("pmtY[%i]/D",totalPMTs));
    DebugTree->Branch("pmtZ",pmtZall,Form("pmtZ[%i]/D",totalPMTs));
    DebugTree->Branch("pmtDirX",pmtDirXall,Form("pmtDirX[%i]/D",totalPMTs));
    DebugTree->Branch("pmtDirY",pmtDirYall,Form("pmtDirY[%i]/D",totalPMTs));
    DebugTree->Branch("pmtDirZ",pmtDirZall,Form("pmtDirZ[%i]/D",totalPMTs));
    DebugTree->Branch("pmtTimeRes",pmtTimeResAll,Form("pmtTimeRes[%i]/D",totalPMTs));
//    DebugTree->Branch("trueExpectedPEHighEMuon", trueExpectedPEHighEMuon,Form("trueExpectedPEHighEMuon[nSubevents][%i]/D",totalPMTs));
//    DebugTree->Branch("trueExpectedPEHighEElectron", trueExpectedPEHighEElectron,Form("trueExpectedPEHighEElectron[nSubevents][%i]/D",totalPMTs));
//    DebugTree->Branch("recoExpectedPEHighEMuon", recoExpectedPEHighEMuon,Form("recoExpectedPEHighEMuon[nSubevents][%i]/D",totalPMTs));
//    DebugTree->Branch("recoExpectedPEHighEElectron", recoExpectedPEHighEElectron,Form("recoExpectedPEHighEElectron[nSubevents][%i]/D",totalPMTs));
    DebugTree->Branch("observedPE",observedPE,Form("observedPE[nSubevents][%i]/I",totalPMTs));
    DebugTree->Branch("allPE",allPE,Form("allPE[nClusters][%i]/I",totalPMTs));
//    DebugTree->Branch("ring1PE",ring1PE,Form("ring1PE[nSubevents][%i]/I",totalPMTs));
//    DebugTree->Branch("ring2PE",ring2PE,Form("ring2PE[nSubevents][%i]/I",totalPMTs));
    DebugTree->Branch("trueToWall",&trueToWall,"trueToWall/D");
    DebugTree->Branch("trueDWall",&trueDWall,"trueDWall/D");
    DebugTree->Branch("recoToWall",&recoToWall,"recoToWall/D");
    DebugTree->Branch("recoDWall",&recoDWall,"recoDWall/D");
    DebugTree->Branch("recoCaptures",&recoCaptures,"recoCaptures/I");

    // Loop over events for reconstruction
    if(maxEvents==0 || maxEvents>nEntries) maxEvents=nEntries;
    for(int i=startEvent; i<maxEvents; i++){
        cout << "---------------------------------------" << endl;
        cout<<"Loadng event: "<<i<<endl;

        mTR->LoadEvent(i);
        evt=mTR->get_genevt();

        // Get truth info
        if(isParticleGun){ //Particle gun was used
            int* parent = mTR->get_part_parentid();
            //Find particle with no parent for vertex truth info
            int ipart=0;
            while(parent[ipart]!=0)
                ipart++;
            trueVtxX = mTR->get_part_xStart()[ipart]/10.; //convert to cm
            trueVtxY = mTR->get_part_yStart()[ipart]/10.;
            trueVtxZ = mTR->get_part_zStart()[ipart]/10.;
            trueDirX = mTR->get_part_pxStart()[ipart];
            trueDirY = mTR->get_part_pyStart()[ipart];
            trueDirZ = mTR->get_part_pzStart()[ipart];
            trueKE = mTR->get_part_KEstart()[ipart];/*
            trueEndX = mTR->get_part_xEnd()[ipart]/10.;
            trueEndY = mTR->get_part_yEnd()[ipart]/10.;
            trueEndZ = mTR->get_part_zEnd()[ipart]/10.;*/
            mode = 0;
            neutrinoE = -9999;
            neutrinoDirX = -9999;
            neutrinoDirY = -9999;
            neutrinoDirZ = -9999;
            neutrinoPID = -9999;
            nNeutrons = -9999;
            neutronCount = mTR->get_neutroncount();
            nCaptures = mTR->get_ncapturecount();
        }
        else{
            trueVtxX = mTR->get_genvtxx();
            trueVtxY = mTR->get_genvtxy();
            trueVtxZ = mTR->get_genvtxz();
            int trackNumber=0;
            int *pids = mTR->get_genpid();
            while(TMath::Abs(pids[trackNumber])<11 || TMath::Abs(pids[trackNumber])>14)
                trackNumber++;
            trueDirX = mTR->get_genpx()[trackNumber];
            trueDirY = mTR->get_genpy()[trackNumber];
            trueDirZ = mTR->get_genpz()[trackNumber];
            trueKE = mTR->get_genKE()[trackNumber];
            mode = mTR->get_genmode();
            neutrinoE = mTR->get_genE();
            neutrinoDirX = mTR->get_genbeam_px();
            neutrinoDirY = mTR->get_genbeam_py();
            neutrinoDirZ = mTR->get_genbeam_pz();
            neutrinoPID = mTR->get_genbeam_id();
            nNeutrons = mTR->get_gennneutrons();
            neutronCount = mTR->get_neutroncount();
            nCaptures = mTR->get_ncapturecount();
        }
        double trueVtxR2 = trueVtxX*trueVtxX+trueVtxY*trueVtxY;
        double trueDWallR = 550-TMath::Sqrt(trueVtxR2);
        double trueDWallZ = 1100-TMath::Abs(trueVtxZ);
        trueDWall = trueDWallR<trueDWallZ ? trueDWallR : trueDWallZ;
        double a = 1-trueDirZ*trueDirZ;
        double b = trueVtxX*trueDirX+trueVtxY*trueDirY;
        double c = trueVtxR2-550*550;
        double trueToWallR = (TMath::Sqrt(b*b-a*c)-b)/a;
        double trueToWallZ = 1100 - trueVtxZ*TMath::Abs(trueDirZ);
        trueToWall = trueToWallR<trueToWallZ ? trueToWallR : trueToWallZ;
        trueTime = 0;
        recoCaptures = 0;
//    trueDirCosTheta = trueDirZ/TMath::Sqrt(trueDirX*trueDirX+trueDirY*trueDirY+trueDirZ*trueDirZ);
        double trueDirTheta = TMath::ACos(trueDirZ /TMath::Sqrt(trueDirX * trueDirX + trueDirY * trueDirY + trueDirZ * trueDirZ));
        double trueDirPhi = TMath::ATan2(trueDirY, trueDirX);

        int totalPEs = mTR->get_nhits();
        if(totalPEs<10
       //     || trueDWall < 100 || abs(mode) != 1 || trueKE < 400 || trueKE > 1000
                ) { //ignore events with too few PEs
            cout << "----- Skipping event! -----" << endl << endl;
            nSubevents =0;
            nClusters=0;
            diffVtxX=-9999;
            diffVtxY=-9999;
            diffVtxZ=-9999;
            diffTime=-9999;
            diffDirX=-9999;
            diffDirY=-9999;
            diffDirZ=-9999;
            diffVtxAbs=-9999;
            diffDirAbs=-9999;
            diffKE=-9999;
            LowETree->Fill();
            HighEElectronTree->Fill();
            HighEMuonTree->Fill();
            FinalTree->Fill();
            DebugTree->Fill();
            continue;
        }

        int * hitCluster = new int[totalPEs](); for(int iPE =0; iPE <totalPEs; iPE++) hitCluster[iPE]=0;
        double *PEhitTimes = mTR->get_hitTime();
        int * clusterHitCounts = FindClusters(totalPEs, PEhitTimes, hitCluster, nClusters);
        cout << "Clusters found: " << nClusters << endl;

        if(nClusters <1) { //ignore events with no clusters
            cout << "----- Skipping event with no clusters! -----" << endl << endl;
            nSubevents =0;
            nClusters=0;
            diffVtxX=-9999;
            diffVtxY=-9999;
            diffVtxZ=-9999;
            diffTime=-9999;
            diffDirX=-9999;
            diffDirY=-9999;
            diffDirZ=-9999;
            diffVtxAbs=-9999;
            diffDirAbs=-9999;
            diffKE=-9999;
            LowETree->Fill();
            HighEElectronTree->Fill();
            HighEMuonTree->Fill();
            FinalTree->Fill();
            DebugTree->Fill();
            continue;
        }

        //Loop over sub-events. Cluster 0 is for hits not in any sub-event, so start at 1
        if(nClusters>maxSubEvts) nClusters =maxSubEvts;
        int subevent = 0;
        nSubevents = nClusters;
        for(int iCluster = 0; subevent< nSubevents; iCluster++, subevent++) {
            cluster[subevent] = iCluster;
            ring[subevent] = 0;
            for(int iPMT=0; iPMT<totalPMTs; iPMT++){
                observedPE[subevent][iPMT] = 0;
                allPE[iCluster][iPMT] = 0;
//                ring1PE[subevent][iPMT] = 0;
//                ring2PE[subevent][iPMT] = 0;
//                trueExpectedPEHighEElectron[subevent][iPMT] = 0;
//                trueExpectedPEHighEMuon[subevent][iPMT] = 0;
            }
            int nHits = clusterHitCounts[iCluster +1];
            cout << "Hits in sub event " << iCluster +1 << ": " << nHits << endl;
            double *hitx = new double[nHits];
            double *hity = new double[nHits];
            double *hitz = new double[nHits];
            double *hitt = new double[nHits];
            int iClusterHit = 0;
            double clusterStartTime = 99999;
            for (int iHit = 0; iHit < totalPEs; iHit++) {
                if (hitCluster[iHit] != iCluster +1) continue;
                hitx[iClusterHit] = mTR->get_hitX()[iHit];
                hity[iClusterHit] = mTR->get_hitY()[iHit];
                hitz[iClusterHit] = mTR->get_hitZ()[iHit];
                hitt[iClusterHit] = mTR->get_hitTime()[iHit];
                if(hitt[iClusterHit]<clusterStartTime) clusterStartTime=hitt[iClusterHit];
                iClusterHit++;
            }
//            for(int iHit=0; iHit<iClusterHit; iHit++){
//                hitt[iHit] -= clusterStartTime;
//            }
            // nHits energy estimate.
            LowEReco lowEReco;
            recoEnergyLowE[iCluster] = lowEReco.ReconstructEnergy(nHits, hitx, hity, hitz);
            // Do Low-E reco
            lowEReco.DoLowEReco(i, nHits, hitx, hity, hitz, hitt, recoChkvAngleLowE[iCluster],
                     recoVtxXLowE[iCluster], recoVtxYLowE[iCluster], recoVtxZLowE[iCluster], recoTimeLowE[iCluster],
                     recoDirXLowE[iCluster], recoDirYLowE[iCluster], recoDirZLowE[iCluster]);
            delete[] hitx;
            delete[] hity;
            delete[] hitz;
            delete[] hitt;
            recoVtxX[subevent] = recoVtxXLowE[iCluster];
            recoVtxY[subevent] = recoVtxYLowE[iCluster];
            recoVtxZ[subevent] = recoVtxZLowE[iCluster];
            recoTime[subevent] = recoTimeLowE[iCluster];
            recoDirX[subevent] = recoDirXLowE[iCluster];
            recoDirY[subevent] = recoDirYLowE[iCluster];
            recoDirZ[subevent] = recoDirZLowE[iCluster];
            recoChkvAngle[subevent] = recoChkvAngleLowE[iCluster];
            recoEnergy[subevent] = recoEnergyLowE[iCluster];

            cout << "--- nhits Energy estimate: " << recoEnergy[subevent];
            if(subevent == 0) cout << "  (true: " << trueKE << ")";
            cout << endl;

            // If < threshold, don't do high E reco
            const int highEThreshold = 60;
            if (recoEnergyLowE[iCluster] < highEThreshold){
                //Check for neutron capture
                if(recoEnergy[subevent] > 2 && recoEnergy[subevent] < 10 && recoTime[subevent] > 2000 && recoTime[subevent] < 100000)
                    recoCaptures++;
                recoPID[subevent] = 0;
                recoVtxXHighEElectron[subevent] = -9999;
                recoVtxYHighEElectron[subevent] = -9999;
                recoVtxZHighEElectron[subevent] = -9999;
                recoTimeHighEElectron[subevent] = -9999;
                recoDirXHighEElectron[subevent] = -9999;
                recoDirYHighEElectron[subevent] = -9999;
                recoDirZHighEElectron[subevent] = -9999;
                recoChkvAngleHighEElectron[subevent] = -9999;
                recoEnergyHighEElectron[subevent] = -9999;
                recoLnLHighEElectron[subevent] = -9999;
                recoVtxXHighEMuon[subevent] = -9999;
                recoVtxYHighEMuon[subevent] = -9999;
                recoVtxZHighEMuon[subevent] = -9999;
                recoTimeHighEMuon[subevent] = -9999;
                recoDirXHighEMuon[subevent] = -9999;
                recoDirYHighEMuon[subevent] = -9999;
                recoDirZHighEMuon[subevent] = -9999;
                recoChkvAngleHighEMuon[subevent] = -9999;
                recoEnergyHighEMuon[subevent] = -9999;
                recoLnLHighEMuon[subevent] = -9999;
                recoNRings[iCluster] = 0;
                isHighE[subevent] = false;
                continue;
            }

            //High-E reco
            HighEReco highEReco;
            highEReco.electronPhotons = electronPhotons;
            highEReco.muonPhotons = muonPhotons;
            highEReco.electronTimes = electronTimes;
            highEReco.muonTimes = muonTimes;

            int *hitPMTids = mTR->get_hitPMTid();

            //Find number of PEs from current subevent on each PMT
            double maxTime = clusterStartTime+90; //ignore all PEs after 90ns (~22m for photon in water)
            cout << "Cluster start time: " << clusterStartTime << endl;
            highEReco.nPEs = 0;
            for (int iPE = 0; iPE < totalPEs; iPE++) {
                if (PEhitTimes[iPE] > maxTime || hitCluster[iPE] != iCluster +1) continue;
                highEReco.nPEs++;
                allPE[iCluster][hitPMTids[iPE]]++;
            }
            //Now create new arrays for PMTs that were actually hit
            highEReco.nHitPMT = 0;
            for (int iPMT = 0; iPMT < totalPMTs; iPMT++) {
                highEReco.nHitPMT += (allPE[iCluster][iPMT] > 0);
            }
            cout << highEReco.nHitPMT << " PMTs with" << highEReco.nPEs << " hits" << endl;
            highEReco.hitPMTx = new double[highEReco.nHitPMT];
            highEReco.hitPMTy = new double[highEReco.nHitPMT];
            highEReco.hitPMTz = new double[highEReco.nHitPMT];
            highEReco.hitPMTDirX = new double[highEReco.nHitPMT];
            highEReco.hitPMTDirY = new double[highEReco.nHitPMT];
            highEReco.hitPMTDirZ = new double[highEReco.nHitPMT];
            highEReco.hitPMTTimeRes = new double[highEReco.nHitPMT];
            int *pmtID = new int[highEReco.nHitPMT];
            highEReco.hitPMTPEs = new int[highEReco.nHitPMT](); for(int iPMT =0; iPMT <highEReco.nHitPMT; iPMT++) highEReco.hitPMTPEs[iPMT]=0;
            int *newPMTids = new int[totalPMTs];
            for (int iPMT = 0, iPMT2 = 0; iPMT < totalPMTs; iPMT++) {
                if (allPE[iCluster][iPMT] == 0) continue;
                highEReco.hitPMTPEs[iPMT2] = allPE[iCluster][iPMT];
                highEReco.hitPMTx[iPMT2] = pmtXall[iPMT];
                highEReco.hitPMTy[iPMT2] = pmtYall[iPMT];
                highEReco.hitPMTz[iPMT2] = pmtZall[iPMT];
                highEReco.hitPMTDirX[iPMT2] = pmtDirXall[iPMT];
                highEReco.hitPMTDirY[iPMT2] = pmtDirYall[iPMT];
                highEReco.hitPMTDirZ[iPMT2] = pmtDirZall[iPMT];
                highEReco.hitPMTTimeRes[iPMT2] = pmtTimeResAll[iPMT];
                pmtID[iPMT2] = pmtIDall[iPMT];
                newPMTids[iPMT] = iPMT2;
                iPMT2++;
            }

            highEReco.hitPMT = new int[highEReco.nPEs];
            highEReco.hitT = new double[highEReco.nPEs];
            highEReco.hitRing = new int[highEReco.nPEs];
            // Loop over PEs create arrays of hit PMT, hit time, etc, for use later
            for (int iPE = 0, iPE2 = 0; iPE < totalPEs; iPE++) {
                if (PEhitTimes[iPE] > maxTime || hitCluster[iPE] != iCluster +1) continue;
                int iPMT = newPMTids[hitPMTids[iPE]];
                highEReco.hitPMT[iPE2] = iPMT;
                highEReco.hitT[iPE2] = PEhitTimes[iPE];
                highEReco.hitRing[iPE2] = 999;
                iPE2++;
            }
            delete[] newPMTids;

            //Initial point fit for starting vertex
            double recoPar[7];
            highEReco.PointFit(recoVtxX[subevent], recoVtxY[subevent], recoVtxZ[subevent], recoTime[subevent], recoChkvAngle[subevent]);
            pointVtxX[iCluster] = recoVtxX[subevent];
            pointVtxY[iCluster] = recoVtxY[subevent];
            pointVtxZ[iCluster] = recoVtxZ[subevent];
            pointTime[iCluster] = recoTime[subevent];

            int maxRings = maxSubEvts-subevent;
            double *peakTheta = new double[maxRings];
            double *peakPhi = new double[maxRings];
            int * ringPE = new int[maxRings];
            recoNRings[iCluster] = highEReco.FindRings(recoVtxX[subevent],recoVtxY[subevent],recoVtxZ[subevent],recoTime[subevent],
                                                       peakTheta,peakPhi,ringPE,maxRings);

            nSubevents += recoNRings[iCluster]-1;
            if(nSubevents > maxSubEvts) nSubevents=maxSubEvts;

            cout << "Found " << recoNRings[iCluster] << " rings." << endl;
            for(int iRing = 0; iRing < recoNRings[iCluster]; iRing++)
                cout << "Ring " << iRing+1 << ": " << ringPE[iRing] << " PEs    theta=" << peakTheta[iRing] << " phi=" << peakPhi[iRing] << endl;

            //Save all cluster hit info
            int * clusterHitPMTPEs = highEReco.hitPMTPEs;
            int * clusterHitPMT = highEReco.hitPMT;
            double * clusterHitT = highEReco.hitT;
            int clusterPEs = highEReco.nPEs;
            double * clusterHitPMTx = highEReco.hitPMTx;
            double * clusterHitPMTy = highEReco.hitPMTy;
            double * clusterHitPMTz = highEReco.hitPMTz;
            double * clusterHitPMTDirX = highEReco.hitPMTDirX;
            double * clusterHitPMTDirY = highEReco.hitPMTDirY;
            double * clusterHitPMTDirZ = highEReco.hitPMTDirZ;
            double * clusterHitPMTTimeRes = highEReco.hitPMTTimeRes;
            int * clusterPMTid = pmtID;
            int clusterHitPMTs = highEReco.nHitPMT;

            // Copy point fit seed to each ring
            for(int iRing = 1; iRing < recoNRings[iCluster]; iRing++){
                recoVtxX[subevent+iRing] = recoVtxX[subevent];
                recoVtxY[subevent+iRing] = recoVtxY[subevent];
                recoVtxZ[subevent+iRing] = recoVtxZ[subevent];
                recoTime[subevent+iRing] = recoTime[subevent];
                recoChkvAngle[subevent+iRing] = recoChkvAngle[subevent];
                for(int iPMT=0; iPMT<totalPMTs; iPMT++) observedPE[subevent+iRing][iPMT] = 0;
            }

            for(int iRing = 0; iRing < recoNRings[iCluster]; iRing++, subevent++){
                ringPEs[subevent] = ringPE[iRing];
                double recoDirTheta = peakTheta[iRing];
                double recoDirCosTheta = TMath::Cos(recoDirTheta);
                double recoDirPhi = peakPhi[iRing];

                cout << "Ring " << iRing+1 << endl;

                //Only use PEs from this ring from now on
                ring[subevent] = iRing+1;
                double *newHitT = new double[ringPEs[subevent]];
                int *newHitPMT = new int[ringPEs[subevent]];
                highEReco.hitPMTPEs = new int[clusterHitPMTs](); for(int iPMT =0; iPMT < clusterHitPMTs; iPMT++) highEReco.hitPMTPEs[iPMT]=0;
//                cout << iRing << " " << ringPEs[subevent] << ":" << endl;
                int nPEnew = 0;
                for (int iPE = 0; iPE < clusterPEs; iPE++) {
                    if (highEReco.hitRing[iPE] != iRing+1) continue;
//                    cout << iPE << " " << nPEnew << endl;
                    newHitT[nPEnew] = clusterHitT[iPE];
                    newHitPMT[nPEnew] = clusterHitPMT[iPE];
                    highEReco.hitPMTPEs[clusterHitPMT[iPE]]++;
                    nPEnew++;
                }
                highEReco.hitT = newHitT;
                highEReco.hitPMT = newHitPMT;
                highEReco.nPEs = nPEnew;

                //Now create new arrays for PMTs that were actually hit
                int nPMTnew = 0;
                for (int iPMT = 0; iPMT < clusterHitPMTs; iPMT++) {
                    nPMTnew += (highEReco.hitPMTPEs[iPMT] > 0);
                }
                cout << nPMTnew << " PMTs with hits" << endl;
                highEReco.hitPMTx = new double[nPMTnew];
                highEReco.hitPMTy = new double[nPMTnew];
                highEReco.hitPMTz = new double[nPMTnew];
                highEReco.hitPMTDirX = new double[nPMTnew];
                highEReco.hitPMTDirY = new double[nPMTnew];
                highEReco.hitPMTDirZ = new double[nPMTnew];
                highEReco.hitPMTTimeRes = new double[nPMTnew];
                pmtID = new int[nPMTnew];
                int *pmtPEsNew = new int[nPMTnew](); for(int iPMT =0; iPMT <nPMTnew; iPMT++) pmtPEsNew[iPMT]=0;
                newPMTids = new int[clusterHitPMTs];
                for (int iPMT = 0, iPMT2 = 0; iPMT < clusterHitPMTs; iPMT++) {
                    if (highEReco.hitPMTPEs[iPMT] < 1) continue;
                    pmtPEsNew[iPMT2] = highEReco.hitPMTPEs[iPMT];
                    highEReco.hitPMTx[iPMT2] = clusterHitPMTx[iPMT];
                    highEReco.hitPMTy[iPMT2] = clusterHitPMTy[iPMT];
                    highEReco.hitPMTz[iPMT2] = clusterHitPMTz[iPMT];
                    highEReco.hitPMTDirX[iPMT2] = clusterHitPMTDirX[iPMT];
                    highEReco.hitPMTDirY[iPMT2] = clusterHitPMTDirY[iPMT];
                    highEReco.hitPMTDirZ[iPMT2] = clusterHitPMTDirZ[iPMT];
                    highEReco.hitPMTTimeRes[iPMT2] = clusterHitPMTTimeRes[iPMT];
                    pmtID[iPMT2] = clusterPMTid[iPMT];
                    newPMTids[iPMT] = iPMT2;
                    iPMT2++;
                }
                delete[] highEReco.hitPMTPEs;
                highEReco.hitPMTPEs = pmtPEsNew;
                highEReco.nHitPMT = nPMTnew;
                for (int iPE = 0; iPE < highEReco.nPEs; iPE++){
                    highEReco.hitPMT[iPE] = newPMTids[highEReco.hitPMT[iPE]];
                }
                delete[] newPMTids;

                cout << "Use approx direction of ring as initial track direction" << endl;
                cout << "True vtx: (" << trueVtxX << ", " << trueVtxY << ", " << trueVtxZ << ", " << trueTime
                << ") dir: (" << TMath::Sin(trueDirTheta) * TMath::Cos(trueDirPhi) << " " <<
                TMath::Sin(trueDirTheta) * TMath::Sin(trueDirPhi) << " " << TMath::Cos(trueDirTheta) << ")" << endl;
                cout << "Reco vtx: (" << recoVtxX[subevent] << ", " << recoVtxY[subevent] << ", " << recoVtxZ[subevent] <<
                ", " << recoTime[subevent]
                << ") dir: (" << TMath::Sin(recoDirTheta) * TMath::Cos(recoDirPhi) << " " <<
                TMath::Sin(recoDirTheta) * TMath::Sin(recoDirPhi) << " " << TMath::Cos(recoDirTheta) << ")" << endl;
                recoPar[4] = recoDirCosTheta;
                recoPar[5] = recoDirPhi;

                //Fit including track:
                highEReco.TrackFit(recoVtxX[subevent], recoVtxY[subevent], recoVtxZ[subevent], recoTime[subevent],
                                   recoDirPhi, recoDirTheta,
                                   recoDirCosTheta, recoChkvAngle[subevent]);
                recoDirX[subevent] = TMath::Sin(recoDirTheta) * TMath::Cos(recoDirPhi);
                recoDirY[subevent] = TMath::Sin(recoDirTheta) * TMath::Sin(recoDirPhi);
                recoDirZ[subevent] = TMath::Cos(recoDirTheta);
                trackVtxX[subevent] = recoVtxX[subevent];
                trackVtxY[subevent] = recoVtxY[subevent];
                trackVtxZ[subevent] = recoVtxZ[subevent];
                trackTime[subevent] = recoTime[subevent];
                trackDirX[subevent] = recoDirX[subevent];
                trackDirY[subevent] = recoDirY[subevent];
                trackDirZ[subevent] = recoDirZ[subevent];

                //Energy reconstruction and correction to vertex in track direction

                double dWallZ = 1100 - TMath::Abs(recoVtxZ[subevent]);
                double dWallR = 550 - TMath::Sqrt(
                        recoVtxX[subevent] * recoVtxX[subevent] + recoVtxY[subevent] * recoVtxY[subevent]);
                double dWall = dWallR < dWallZ ? dWallR : dWallZ;
                double lookupPEs = ringPEs[subevent];
                double dWallMin = muonEnergyLookup->GetYaxis()->GetBinCenter(muonEnergyLookup->GetYaxis()->GetFirst());
                double dWallMax = muonEnergyLookup->GetYaxis()->GetBinCenter(muonEnergyLookup->GetYaxis()->GetLast());
                double PEmin = muonEnergyLookup->GetXaxis()->GetBinCenter(muonEnergyLookup->GetXaxis()->GetFirst());
                double PEmax = muonEnergyLookup->GetXaxis()->GetBinCenter(muonEnergyLookup->GetXaxis()->GetLast());
                if (dWall <= dWallMin) dWall = dWallMin + 0.01;
                else if (dWall >= dWallMax) dWall = dWallMax - 0.01;
                double correctionFactor = 1.;
                if (lookupPEs <= PEmin) {
                    correctionFactor = lookupPEs / PEmin;
                    lookupPEs = PEmin + 1;
                }
                else if (lookupPEs >= PEmax) {
                    correctionFactor = lookupPEs / PEmax;
                    lookupPEs = PEmax - 1;
                }
                recoEnergyHighEMuon[subevent] = muonEnergyLookup->Interpolate(lookupPEs, dWall)*correctionFactor;
                recoEnergyHighEElectron[subevent] = electronEnergyLookup->Interpolate(lookupPEs, dWall)*correctionFactor;
                muonTrackCorrection[subevent] = muonVtxBiasLookup->Interpolate(lookupPEs, dWall);
                electronTrackCorrection[subevent] = electronVtxBiasLookup->Interpolate(lookupPEs, dWall);

                cout << "Energy lookup:    Muon: " << recoEnergyHighEMuon[subevent] << "    Electron: " <<
                recoEnergyHighEElectron[subevent] << endl;
/*
                if (recoEnergyHighEElectron[subevent] < highEThreshold ||
                    recoEnergyHighEMuon[subevent] < highEThreshold) {
                    recoPID[subevent] = 0;
                    recoVtxXHighEElectron[subevent] = -9999;
                    recoVtxYHighEElectron[subevent] = -9999;
                    recoVtxZHighEElectron[subevent] = -9999;
                    recoTimeHighEElectron[subevent] = -9999;
                    recoDirXHighEElectron[subevent] = -9999;
                    recoDirYHighEElectron[subevent] = -9999;
                    recoDirZHighEElectron[subevent] = -9999;
                    recoChkvAngleHighEElectron[subevent] = -9999;
                    recoLnLHighEElectron[subevent] = -9999;
                    recoVtxXHighEMuon[subevent] = -9999;
                    recoVtxYHighEMuon[subevent] = -9999;
                    recoVtxZHighEMuon[subevent] = -9999;
                    recoTimeHighEMuon[subevent] = -9999;
                    recoDirXHighEMuon[subevent] = -9999;
                    recoDirYHighEMuon[subevent] = -9999;
                    recoDirZHighEMuon[subevent] = -9999;
                    recoChkvAngleHighEMuon[subevent] = -9999;
                    recoLnLHighEMuon[subevent] = -9999;
                    isHighE[subevent] = false;
                    continue;
                }
*/
                //Need to use unhit PMTs in likelihood
                for (int iPMT = 0; iPMT < highEReco.nHitPMT; iPMT++) {
                    observedPE[subevent][pmtID[iPMT]] = highEReco.hitPMTPEs[iPMT];
                }
                highEReco.nUnhitPMT = totalPMTs - highEReco.nHitPMT;
                highEReco.unhitPMTx = new double[highEReco.nUnhitPMT];
                highEReco.unhitPMTy = new double[highEReco.nUnhitPMT];
                highEReco.unhitPMTz = new double[highEReco.nUnhitPMT];
                highEReco.unhitPMTDirX = new double[highEReco.nUnhitPMT];
                highEReco.unhitPMTDirY = new double[highEReco.nUnhitPMT];
                highEReco.unhitPMTDirZ = new double[highEReco.nUnhitPMT];
                highEReco.unhitPMTTimeRes = new double[highEReco.nUnhitPMT];
                cout << "Total PMTs:" << totalPMTs << " Hit:" << highEReco.nHitPMT << " Unhit:" << highEReco.nUnhitPMT << endl;
                for (int iPMT = 0, iPMT2 = 0; iPMT < totalPMTs; iPMT++) {
                    if (observedPE[subevent][iPMT] > 0) continue;
                    highEReco.unhitPMTx[iPMT2] = pmtXall[iPMT];
                    highEReco.unhitPMTy[iPMT2] = pmtYall[iPMT];
                    highEReco.unhitPMTz[iPMT2] = pmtZall[iPMT];
                    highEReco.unhitPMTDirX[iPMT2] = pmtDirXall[iPMT];
                    highEReco.unhitPMTDirY[iPMT2] = pmtDirYall[iPMT];
                    highEReco.unhitPMTDirZ[iPMT2] = pmtDirZall[iPMT];
                    highEReco.unhitPMTTimeRes[iPMT2] = pmtTimeResAll[iPMT];
                    iPMT2++;
                }

                // Perform maximum likelihood for electron hypothesis
                cout << "Fit for electron hypothesis" << endl;
                recoVtxXHighEElectron[subevent] = recoVtxX[subevent];
                recoVtxYHighEElectron[subevent] = recoVtxY[subevent];
                recoVtxZHighEElectron[subevent] = recoVtxZ[subevent];
                recoTimeHighEElectron[subevent] = recoTime[subevent];
                recoChkvAngleHighEElectron[subevent] = recoChkvAngleLowE[iCluster];
                double recoDirPhiHighEElectron = recoDirPhi;
                double recoDirThetaHighEElectron = recoDirTheta;
                highEReco.LikelihoodFit(electronTrackCorrection[subevent], recoVtxXHighEElectron[subevent],
                                        recoVtxYHighEElectron[subevent], recoVtxZHighEElectron[subevent],
                                        recoTimeHighEElectron[subevent],
                                        recoDirPhiHighEElectron, recoDirThetaHighEElectron,
                                        recoEnergyHighEElectron[subevent], recoLnLHighEElectron[subevent], 11);
                recoDirXHighEElectron[subevent] =
                        TMath::Sin(recoDirThetaHighEElectron) * TMath::Cos(recoDirPhiHighEElectron);
                recoDirYHighEElectron[subevent] =
                        TMath::Sin(recoDirThetaHighEElectron) * TMath::Sin(recoDirPhiHighEElectron);
                recoDirZHighEElectron[subevent] = TMath::Cos(recoDirThetaHighEElectron);
//                for (int iPMT = 0, iHitPMT = 0, iUnhitPMT = 0; iPMT < totalPMTs; iPMT++) {
//                    bool isHit = allPE[subevent][iPMT] > 0;
//                    int thisID = isHit ? iHitPMT++ : iUnhitPMT++;
//                    trueExpectedPEHighEElectron[subevent][iPMT] = highEReco.ExpectedPMTPhotoelectrons(trueVtxX,
//                                                                                                     trueVtxY,
//                                                                                                     trueVtxZ,
//                                                                                                     trueDirX,
//                                                                                                     trueDirY,
//                                                                                                     trueDirZ,
//                                                                                                     thisID,
//                                                                                                     isHit,
//                                                                                                     trueKE,
//                                                                                                     11);
//                    recoExpectedPEHighEElectron[subevent][iPMT] = highEReco.ExpectedPMTPhotoelectrons(
//                            recoVtxXHighEElectron[subevent],
//                            recoVtxYHighEElectron[subevent],
//                            recoVtxZHighEElectron[subevent],
//                            recoDirXHighEElectron[subevent],
//                            recoDirYHighEElectron[subevent],
//                            recoDirZHighEElectron[subevent],
//                            thisID,
//                            isHit,
//                            recoEnergyHighEElectron[subevent],
//                            11);
//                }

                // Perform maximum likelihood for muon hypothesis
                cout << "Fit for muon hypothesis" << endl;
                recoVtxXHighEMuon[subevent] = recoVtxX[subevent];
                recoVtxYHighEMuon[subevent] = recoVtxY[subevent];
                recoVtxZHighEMuon[subevent] = recoVtxZ[subevent];
                recoTimeHighEMuon[subevent] = recoTime[subevent];
                recoChkvAngleHighEMuon[subevent] = recoChkvAngleLowE[iCluster];
                double recoDirPhiHighEMuon = recoDirPhi;
                double recoDirThetaHighEMuon = recoDirTheta;
                highEReco.LikelihoodFit(muonTrackCorrection[subevent], recoVtxXHighEMuon[subevent],
                                        recoVtxYHighEMuon[subevent], recoVtxZHighEMuon[subevent],
                                        recoTimeHighEMuon[subevent],
                                        recoDirPhiHighEMuon, recoDirThetaHighEMuon,
                                        recoEnergyHighEMuon[subevent], recoLnLHighEMuon[subevent], 13);
                recoDirXHighEMuon[subevent] = TMath::Sin(recoDirThetaHighEMuon) * TMath::Cos(recoDirPhiHighEMuon);
                recoDirYHighEMuon[subevent] = TMath::Sin(recoDirThetaHighEMuon) * TMath::Sin(recoDirPhiHighEMuon);
                recoDirZHighEMuon[subevent] = TMath::Cos(recoDirThetaHighEMuon);
//                for (int iPMT = 0, iHitPMT = 0, iUnhitPMT = 0; iPMT < totalPMTs; iPMT++) {
//                    bool isHit = allPE[subevent][iPMT] > 0;
//                    int thisID = isHit ? iHitPMT++ : iUnhitPMT++;
//                    trueExpectedPEHighEMuon[subevent][iPMT] = highEReco.ExpectedPMTPhotoelectrons(trueVtxX,
//                                                                                                 trueVtxY,
//                                                                                                 trueVtxZ,
//                                                                                                 trueDirX,
//                                                                                                 trueDirY,
//                                                                                                 trueDirZ,
//                                                                                                 thisID,
//                                                                                                 isHit,
//                                                                                                 trueKE,
//                                                                                                 13);
//                    recoExpectedPEHighEMuon[subevent][iPMT] = highEReco.ExpectedPMTPhotoelectrons(
//                            recoVtxXHighEMuon[subevent],
//                            recoVtxYHighEMuon[subevent],
//                            recoVtxZHighEMuon[subevent],
//                            recoDirXHighEMuon[subevent],
//                            recoDirYHighEMuon[subevent],
//                            recoDirZHighEMuon[subevent],
//                            thisID,
//                            isHit,
//                            recoEnergyHighEMuon[subevent],
//                            13);
//                }

                cout << "Reco  LnL_mu - LnL_e = " << recoLnLHighEMuon[subevent] - recoLnLHighEElectron[subevent] << endl;
                cout << endl << endl;

                if (recoLnLHighEMuon[subevent] - recoLnLHighEElectron[subevent] > 0) {
                    recoPID[subevent] = 13;
                    recoVtxX[subevent] = recoVtxXHighEMuon[subevent];
                    recoVtxY[subevent] = recoVtxYHighEMuon[subevent];
                    recoVtxZ[subevent] = recoVtxZHighEMuon[subevent];
                    recoTime[subevent] = recoTimeHighEMuon[subevent];
                    recoDirX[subevent] = recoDirXHighEMuon[subevent];
                    recoDirY[subevent] = recoDirYHighEMuon[subevent];
                    recoDirZ[subevent] = recoDirZHighEMuon[subevent];
                    recoChkvAngle[subevent] = recoChkvAngleHighEMuon[subevent];
                    recoEnergy[subevent] = recoEnergyHighEMuon[subevent];
                }
                else {
                    recoPID[subevent] = 11;
                    recoVtxX[subevent] = recoVtxXHighEElectron[subevent];
                    recoVtxY[subevent] = recoVtxYHighEElectron[subevent];
                    recoVtxZ[subevent] = recoVtxZHighEElectron[subevent];
                    recoTime[subevent] = recoTimeHighEElectron[subevent];
                    recoDirX[subevent] = recoDirXHighEElectron[subevent];
                    recoDirY[subevent] = recoDirYHighEElectron[subevent];
                    recoDirZ[subevent] = recoDirZHighEElectron[subevent];
                    recoChkvAngle[subevent] = recoChkvAngleHighEElectron[subevent];
                    recoEnergy[subevent] = recoEnergyHighEElectron[subevent];
                }
                isHighE[subevent] = true;

                delete[] highEReco.hitPMTx;
                delete[] highEReco.hitPMTy;
                delete[] highEReco.hitPMTz;
                delete[] highEReco.hitPMTDirX;
                delete[] highEReco.hitPMTDirY;
                delete[] highEReco.hitPMTDirZ;
                delete[] highEReco.hitPMTTimeRes;
                delete[] pmtID;
                delete[] highEReco.hitT;
                delete[] highEReco.hitPMT;
                delete[] highEReco.hitPMTPEs;
                delete[] highEReco.unhitPMTx;
                delete[] highEReco.unhitPMTy;
                delete[] highEReco.unhitPMTz;
                delete[] highEReco.unhitPMTDirX;
                delete[] highEReco.unhitPMTDirY;
                delete[] highEReco.unhitPMTDirZ;
                delete[] highEReco.unhitPMTTimeRes;
            }
            subevent--;
            delete[] highEReco.hitRing;
            delete[] clusterHitPMTPEs;
            delete[] clusterHitPMT;
            delete[] clusterHitT;
            delete[] clusterHitPMTx;
            delete[] clusterHitPMTy;
            delete[] clusterHitPMTz;
            delete[] clusterHitPMTDirX;
            delete[] clusterHitPMTDirY;
            delete[] clusterHitPMTDirZ;
            delete[] clusterHitPMTTimeRes;
            delete[] clusterPMTid;
            delete[] peakTheta;
            delete[] peakPhi;
            delete[] ringPE;
//    break;
        }
        trueDirX = TMath::Sin(trueDirTheta) * TMath::Cos(trueDirPhi);
        trueDirY = TMath::Sin(trueDirTheta) * TMath::Sin(trueDirPhi);
        trueDirZ = TMath::Cos(trueDirTheta);
        diffVtxX = recoVtxX[0]-trueVtxX;
        diffVtxY = recoVtxY[0]-trueVtxY;
        diffVtxZ = recoVtxZ[0]-trueVtxZ;
        diffVtxAbs = TMath::Sqrt(diffVtxX*diffVtxX + diffVtxY*diffVtxY + diffVtxZ*diffVtxZ);
        diffTime = recoTime[0]-trueTime;
        diffDirX = recoDirX[0]-trueDirX;
        diffDirY = recoDirY[0]-trueDirY;
        diffDirZ = recoDirZ[0]-trueDirZ;
        diffDirAbs = TMath::ACos(recoDirX[0]*trueDirX + recoDirY[0]*trueDirY + recoDirZ[0]*trueDirZ);
        diffKE= recoEnergy[0]-trueKE;
        double recoVtxR2 = recoVtxX[0]*recoVtxX[0]+recoVtxY[0]*recoVtxY[0];
        double recoDWallR = 550-TMath::Sqrt(recoVtxR2);
        double recoDWallZ = 1100-TMath::Abs(recoVtxZ[0]);
        recoDWall = recoDWallR<recoDWallZ ? recoDWallR : recoDWallZ;
        a = 1-recoDirZ[0]*recoDirZ[0];
        b = recoVtxX[0]*recoDirX[0]+recoVtxY[0]*recoDirY[0];
        c = recoVtxR2-550*550;
        double recoToWallR = (TMath::Sqrt(b*b-a*c)-b)/a;
        double recoToWallZ = 1100 - recoVtxZ[0]*TMath::Abs(recoDirZ[0]);
        recoToWall = recoToWallR<recoToWallZ ? recoToWallR : recoToWallZ;
        cout << "True vtx: (" << trueVtxX << "," << trueVtxY << "," << trueVtxZ << ")  time: " << trueTime
             << "  dir: (" << trueDirX << "," << trueDirY << "," << trueDirZ << ")  Energy:" << trueKE << endl;
        cout << "Reco vtx: (" << recoVtxX[0] << "," << recoVtxY[0] << "," << recoVtxZ[0] << ")  time: " << recoTime[0]
        << "  dir: (" << recoDirX[0] << "," << recoDirY[0] << "," << recoDirZ[0] << ")  Energy:" << recoEnergy[0] << endl;
        cout << "Diff vtx: (" << diffVtxX << "," << diffVtxY << "," << diffVtxZ << ")  time: " << diffTime
        << "  dir: (" << diffDirX << "," << diffDirY << "," << diffDirZ << ")  Energy:" << diffKE << endl;
        cout << "Diff vtx abs: " << diffVtxAbs << "     Diff dir abs: " << diffDirAbs*180./TMath::Pi() << endl << endl;
        LowETree->Fill();
        HighEElectronTree->Fill();
        HighEMuonTree->Fill();
        FinalTree->Fill();
        DebugTree->Fill();
        delete[] clusterHitCounts;
        delete[] hitCluster;
    } //end i-loop over Hits_Tree entries

    f_out.cd();
    LowETree->Write();
    HighEElectronTree->Write();
    HighEMuonTree->Write();
    FinalTree->Write();
    DebugTree->Write();
    delete mTR;
    delete LowETree;
    delete HighEElectronTree;
    delete HighEMuonTree;
    delete FinalTree;
    delete DebugTree;
//
    delete[] pmtXall;
    delete[] pmtYall;
    delete[] pmtZall;
    delete[] pmtDirXall;
    delete[] pmtDirYall;
    delete[] pmtDirZall;
    delete[] pmtTimeResAll;
    delete[] pmtIDall;
    likelihoodTables->Close();
    energyTables->Close();
}