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
                double FullTimeOfFlight( const double& vtxX, const double& vtxY,const double& vtxZ,const double& dirX, const double& dirY, const double& dirZ, const int& pmt,const double& cherenkovAngle);
//		double FullTimeOfFlight(double vtxX, double vtxY,double vtxZ,double dirX, double dirY, double dirZ, int pmt,double cherenkovAngle);
		double ExpectedPMTPhotoelectrons(double vtxX, double vtxY, double vtxZ, double dirX, double dirY, double dirZ, int pmt, bool isHit, double kineticEnergy, int ipnu, double * expectedTime = 0);

	private:

		double Goodness(const double *timesOfFlight, double time = -1000, double factor = 1, int ring = 0);
		//double PointTimeOfFlight( double vtxX, double vtxY, double vtxZ, int pmt);
		double PointTimeOfFlight(const double& vtxX, const double& vtxY, const double& vtxZ, const int& pmt);
		//double PointChkvAngle(double vtxX, double vtxY, double vtxZ, double dirX,  double dirY, double dirZ, int pmt);
		double PointChkvAngle(const double& vtxX,const double& vtxY,const double& vtxZ, const double& dirX, const double& dirY, const double& dirZ, const int& pmt);
		double TrackChkvAngle( const double& vtxX,const double& vtxY, const double& vtxZ, const double& dirX, const double& dirY, const double& dirZ, const int& pmt, const double& tof);
		//double TrackChkvAngle( double vtxX,double vtxY, double vtxZ, double dirX, double dirY, double dirZ, int pmt, double tof);
		//double PMTlnLikelihood( double vtxX, double vtxY,double vtxZ, double dirX, double dirY, double dirZ, int pmt, bool isHit, double kineticEnergy, int ipnu);

		double PMTlnLikelihood(const double& vtxX,const double& vtxY,const double& vtxZ,const double& dirX,const double& dirY,const double& dirZ,const int& pmt,const bool& isHit,const double& kineticEnergy,const int& ipnu);
		double LnLikelihood(const double *par, int ipnu, bool total);
		double ElectronLnLikelihoodTotal(const double *par);
		double MuonLnLikelihoodTotal(const double *par);
		int FindRing(double vtxX, double vtxY, double vtxZ, double vtxT, double &thetaPeak, double &phiPeak,
				int ringNumber = 0, bool useTrack = true, double *hough = 0);



		ClassDef(HighEReco,0)

};

//'''''''''''''''''''''''''''''''''''''''''''''''''''''//
/////////////////////////////////////////////////////////
//.....................................................//


// Calculated photon emission angle from vertex to PMT assuming point source of photon
inline double  HighEReco::PointChkvAngle(const double& vtxX,const double& vtxY,const double& vtxZ, const double& dirX, const double& dirY, const double& dirZ, const int& pmt){
double vtxToPMTx = hitPMTx[pmt]-vtxX;
double vtxToPMTy = hitPMTy[pmt]-vtxY;
double vtxToPMTz = hitPMTz[pmt]-vtxZ;
double pmtDist = TMath::Sqrt(vtxToPMTx*vtxToPMTx+vtxToPMTy*vtxToPMTy+vtxToPMTz*vtxToPMTz);
return TMath::ACos((dirX*vtxToPMTx+dirY*vtxToPMTy+dirZ*vtxToPMTz)/pmtDist);
}



// Calculated time of flight from vertex to PMT assuming point source of photons
inline double  HighEReco::PointTimeOfFlight(const double& vtxX, const double& vtxY,const double& vtxZ,const int& pmt){
double vtxToPMTx = hitPMTx[pmt]-vtxX;
double vtxToPMTy = hitPMTy[pmt]-vtxY;
double vtxToPMTz = hitPMTz[pmt]-vtxZ;
double pmtDist = TMath::Sqrt(vtxToPMTx*vtxToPMTx+vtxToPMTy*vtxToPMTy+vtxToPMTz*vtxToPMTz);
return pmtDist/C_WAT;
}



// Calculated time of flight from vertex to PMT including lepton track
inline double HighEReco::TrackChkvAngle(const double& vtxX,const double& vtxY, const double& vtxZ,const double& dirX,const double& dirY,const double& dirZ,const int& pmt,const double& tof){
double vtxToPMTx = hitPMTx[pmt]-vtxX;
double vtxToPMTy = hitPMTy[pmt]-vtxY;
double vtxToPMTz = hitPMTz[pmt]-vtxZ;
//Distance from vertex to PMT:
double pmtDist = TMath::Sqrt(vtxToPMTx*vtxToPMTx+vtxToPMTy*vtxToPMTy+vtxToPMTz*vtxToPMTz);
double cosPMTtrackAngle = (dirX*vtxToPMTx+dirY*vtxToPMTy+dirZ*vtxToPMTz)/pmtDist;
double sinPMTtrackAngle = TMath::Sqrt(1-cosPMTtrackAngle*cosPMTtrackAngle);
double A = (C_VAC*tof/pmtDist-cosPMTtrackAngle)/sinPMTtrackAngle;
if(A*A-N_REF*N_REF+1<0)
return 0;
return TMath::ATan((A-N_REF* TMath::Sqrt(A*A-N_REF*N_REF+1))/(N_REF*N_REF-A*A));
}

inline double HighEReco::FullTimeOfFlight(const double& vtxX, const double& vtxY,const double& vtxZ,const double& dirX,const double& dirY,const double& dirZ,const int& pmt,const double& cherenkovAngle){
double vtxToPMTx = hitPMTx[pmt]-vtxX;
double vtxToPMTy = hitPMTy[pmt]-vtxY;
double vtxToPMTz = hitPMTz[pmt]-vtxZ;
double sinChkvAngle = TMath::Sin(cherenkovAngle);
double pmtDist = TMath::Sqrt(vtxToPMTx*vtxToPMTx+vtxToPMTy*vtxToPMTy+vtxToPMTz*vtxToPMTz);
double pmtTrackAngle = TMath::ACos((vtxToPMTx*dirX+vtxToPMTy*dirY+vtxToPMTz*dirZ)/pmtDist);
double trackDist = pmtDist* TMath::Sin(cherenkovAngle-pmtTrackAngle)/sinChkvAngle;
return trackDist/C_VAC + pmtDist*TMath::Sin(pmtTrackAngle)/(C_WAT*sinChkvAngle);
}

inline double HighEReco::PMTlnLikelihood(const double& vtxX,const double& vtxY,const double& vtxZ,const double& dirX,const double& dirY,const double& dirZ,const int& pmt,const bool& isHit,const double& kineticEnergy,const int& ipnu){
double expectedPhotoelectrons = ExpectedPMTPhotoelectrons(vtxX, vtxY, vtxZ, dirX, dirY, dirZ, pmt, isHit, kineticEnergy, ipnu);
double observedPhotoelectrons = hitPMTPEs[pmt];
double l = observedPhotoelectrons* TMath::Log(expectedPhotoelectrons)-expectedPhotoelectrons;
return l;
}


#endif
