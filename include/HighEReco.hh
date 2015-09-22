#ifndef HIGHERECO_HH
#define HIGHERECO_HH

#include "TObject.h"
#include "TH3D.h"

class HighEReco : public TObject {

public:

    static const double C_VAC;
    static const double N_REF;
    static const double C_WAT;
//    static const double factor;
    static const double multiRingFactor;

    double *hitPMTx;
    double *hitPMTy;
    double *hitPMTz;
    double *hitPMTDirX;
    double *hitPMTDirY;
    double *hitPMTDirZ;
    double *hitPMTTimeRes;
    int *hitPMTPEs;
    int nHitPMT;
    double *unhitPMTx;
    double *unhitPMTy;
    double *unhitPMTz;
    double *unhitPMTDirX;
    double *unhitPMTDirY;
    double *unhitPMTDirZ;
    double *unhitPMTTimeRes;
    int nUnhitPMT;
    double *hitT;
    int *hitRing;
    int *hitPMT;
    int nPEs;
    TH3D * electronPhotons;
    TH3D * muonPhotons;
    TH3D * electronTimes;
    TH3D * muonTimes;

    HighEReco();
    ~HighEReco();
    int FindRings(double vtxX, double vtxY, double vtxZ, double vtxT, double *thetaPeaks, double *phiPeaks,
                  int *ringPEs, int maxRings, bool useTrack = true, double *hough = 0);
    TH3D *FindTracks(double vtxX, double vtxY, double vtxZ);
    void PointFit(double &recoVtxX, double &recoVtxY, double &recoVtxZ, double &recoT, double &cherenkovAngle);
    void TrackFit(double &recoVtxX, double &recoVtxY, double &recoVtxZ, double &recoT, double &recoDirPhi,
                  double &recoDirTheta, double &recoDirCosTheta, double &cherenkovAngle);
    void LikelihoodFit(double &trackCorrection, double &recoVtxX, double &recoVtxY,
                       double &recoVtxZ, double &recoT, double &recoDirPhi,
                       double &recoDirTheta, double &recoKE, double &recoLnL, int ipnu);
    double PointBadness(const double *par);
    double TrackBadness(const double *par);
    double TrackGoodness(const double *par, int ring = 0, double factor = 1);
    double ElectronLnLikelihood(const double *par);
    double MuonLnLikelihood(const double *par);
    double FullTimeOfFlight(double vtxX, double vtxY, double vtxZ, double dirX, double dirY, double dirZ, int pmt, double cherenkovAngle);
    double ExpectedPMTPhotoelectrons(double vtxX, double vtxY, double vtxZ, double dirX, double dirY, double dirZ, int pmt, bool isHit, double kineticEnergy, int ipnu, double * expectedTime = 0);

private:

    double Goodness(const double *timesOfFlight, double time = -1000, double factor = 1, int ring = 0);
    double PointTimeOfFlight(double vtxX, double vtxY, double vtxZ, int pmt);
    double PointChkvAngle(double vtxX, double vtxY, double vtxZ, double dirX, double dirY, double dirZ, int pmt);
    double TrackChkvAngle(double vtxX, double vtxY, double vtxZ, double dirX, double dirY, double dirZ, int pmt, double tof);
    double PMTlnLikelihood(double vtxX, double vtxY, double vtxZ, double dirX, double dirY, double dirZ, int pmt, bool isHit, double kineticEnergy, int ipnu);
    double LnLikelihood(const double *par, int ipnu, bool total, bool print);
    double ElectronLnLikelihoodTotal(const double *par);
    double MuonLnLikelihoodTotal(const double *par);
    int FindRing(double vtxX, double vtxY, double vtxZ, double vtxT, double &thetaPeak, double &phiPeak,
                 int ringNumber = 0, bool useTrack = true, double *hough = 0);

ClassDef(HighEReco,0)

};

#endif
