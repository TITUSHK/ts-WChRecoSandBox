#include <TMath.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TH3.h>
#include "HighEReco.hh"

using namespace std;

ClassImp(HighEReco)

HighEReco::HighEReco()
{

}

HighEReco::~HighEReco()
{

}

const double HighEReco::C_VAC = 29.9792458;
const double HighEReco::N_REF = 1.34;
const double HighEReco::factor = 1.;
const double HighEReco::C_WAT = C_VAC/N_REF;
const double HighEReco::multiRingFactor = 0.09;


// General goodness function based on timing residuals of each PE
// timesOfFlight : array of TOF calculated for each PMT
// time          : time of interaction (mean correctedTime is used if time < -999)

double HighEReco::Goodness(const double *timesOfFlight, double time){
    // Loop to calculate PE times corrected for TOF to PMT
    double *correctedTimes = new double[nPEs];
    for(int iPE=0; iPE<nPEs; iPE++)
    {
        int pmt = hitPMT[iPE];
        correctedTimes[iPE] = hitT[iPE]-timesOfFlight[pmt];
    }
    if(time < -999){
        for(int iPE=0; iPE<nPEs; iPE++)
            time += correctedTimes[iPE];
        time /= nPEs;
    }
    double tRes = hitPMTTimeRes[0]; //assume all PMTs the same timing resolution

    // Sum over PEs to get goodness
    double result = 0;
    for(int iPE=0; iPE<nPEs; iPE++)
    {
        result += TMath::Exp(-0.5* TMath::Power((correctedTimes[iPE]-time)/(factor*tRes),2));
    }
    delete[] correctedTimes;
    return result/nPEs;
}

// Calculated time of flight from vertex to PMT assuming point source of photons
double  HighEReco::PointTimeOfFlight(double vtxX, double vtxY, double vtxZ, int pmt){
    double vtxToPMTx = hitPMTx[pmt]-vtxX;
    double vtxToPMTy = hitPMTy[pmt]-vtxY;
    double vtxToPMTz = hitPMTz[pmt]-vtxZ;
    double pmtDist = TMath::Sqrt(vtxToPMTx*vtxToPMTx+vtxToPMTy*vtxToPMTy+vtxToPMTz*vtxToPMTz);
    return pmtDist/C_WAT;
}

// Calculated time of flight from vertex to PMT including lepton track
double HighEReco::FullTimeOfFlight(double vtxX, double vtxY, double vtxZ, double dirX, double dirY, double dirZ, int pmt, double cherenkovAngle){
    double vtxToPMTx = hitPMTx[pmt]-vtxX;
    double vtxToPMTy = hitPMTy[pmt]-vtxY;
    double vtxToPMTz = hitPMTz[pmt]-vtxZ;
    double sinChkvAngle = TMath::Sin(cherenkovAngle);
    //Distance from vertex to PMT:
    double pmtDist = TMath::Sqrt(vtxToPMTx*vtxToPMTx+vtxToPMTy*vtxToPMTy+vtxToPMTz*vtxToPMTz);
    //Angle between track direction and vertex-PMT direction:
    double pmtTrackAngle = TMath::ACos((vtxToPMTx*dirX+vtxToPMTy*dirY+vtxToPMTz*dirZ)/pmtDist);
    //Distance along track traveled before reaching point at Cherenkov angle to PMT:
    double trackDist = pmtDist* TMath::Sin(cherenkovAngle-pmtTrackAngle)/sinChkvAngle;
    //Time of flight of lepton along track plus time of flight of photon to PMT:
    return trackDist/C_VAC + pmtDist*TMath::Sin(pmtTrackAngle)/(C_WAT*sinChkvAngle);
}

// Calculated photon emission angle from vertex to PMT assuming point source of photons
double  HighEReco::PointChkvAngle(double vtxX, double vtxY, double vtxZ, double dirX, double dirY, double dirZ, int pmt){
    double vtxToPMTx = hitPMTx[pmt]-vtxX;
    double vtxToPMTy = hitPMTy[pmt]-vtxY;
    double vtxToPMTz = hitPMTz[pmt]-vtxZ;
    double pmtDist = TMath::Sqrt(vtxToPMTx*vtxToPMTx+vtxToPMTy*vtxToPMTy+vtxToPMTz*vtxToPMTz);
    return TMath::ACos((dirX*vtxToPMTx+dirY*vtxToPMTy+dirZ*vtxToPMTz)/pmtDist);
}

// Calculated time of flight from vertex to PMT including lepton track
double HighEReco::TrackChkvAngle(double vtxX, double vtxY, double vtxZ, double dirX, double dirY, double dirZ, int pmt, double tof){
    double vtxToPMTx = hitPMTx[pmt]-vtxX;
    double vtxToPMTy = hitPMTy[pmt]-vtxY;
    double vtxToPMTz = hitPMTz[pmt]-vtxZ;
    //Distance from vertex to PMT:
    double pmtDist = TMath::Sqrt(vtxToPMTx*vtxToPMTx+vtxToPMTy*vtxToPMTy+vtxToPMTz*vtxToPMTz);
    double cosPMTtrackAngle = (dirX*vtxToPMTx+dirY*vtxToPMTy+dirZ*vtxToPMTz)/pmtDist;
    double sinPMTtrackAngle = TMath::Sqrt(1-cosPMTtrackAngle*cosPMTtrackAngle);
    //To find cherenkov angle, solve equation for TOF: t=d/(c*sin(thetaC)) * (sin(thetaC - alpha) + n_ref*sin(alpha))
    //Solution is:
    // tan(thetaC)=(A-n_ref*sqrt(A^2-n_ref^2+1))/(n_ref^2-A^2)
    //where
    //     A = (c*t/d-cos(alpha))/sin(alpha)
    double A = (C_VAC*tof/pmtDist-cosPMTtrackAngle)/sinPMTtrackAngle;
    // Must have A^2-n_ref^2+1>0 for physical solution.
    if(A*A-N_REF*N_REF+1<0)
        return 0;
    return TMath::ATan((A-N_REF* TMath::Sqrt(A*A-N_REF*N_REF+1))/(N_REF*N_REF-A*A));
}

//Goodness of fit for point fitter
double HighEReco::PointGoodness(const double *par){
    double vtxZ = par[0];
    double vtxY = par[1];
    double vtxX = par[2];
    if(vtxX*vtxX+vtxY*vtxY>550*550) // No vertex outside detector radius
        return 0;
    // Loop to find expected time of flight from vertex to each PMT
    double* timesOfFlight = new double[nHitPMT];
    for(int iPMT=0; iPMT< nHitPMT; iPMT++)
        timesOfFlight[iPMT] = PointTimeOfFlight(vtxX, vtxY, vtxZ, iPMT);
    double result = Goodness(timesOfFlight,0);
    delete[] timesOfFlight;
    return -result; //negative for minimise not maximise
}

// Goodness of fit including the track length for given vertex and direction
double HighEReco::TrackGoodness(const double *par){
    double vtxZ = par[4];
    double vtxY = par[5];
    double vtxX = par[6];
    double time = par[0];
    if(vtxX*vtxX+vtxY*vtxY>550*550)
        return 0;
    double dirCosTheta = par[2];
    double dirSinTheta = TMath::Sqrt(1-dirCosTheta*dirCosTheta);
    double dirPhi = par[3];
    double dirX = dirSinTheta* TMath::Cos(dirPhi);
    double dirY = dirSinTheta* TMath::Sin(dirPhi);
    double dirZ = dirCosTheta;
    double cherenkovAngle = par[1];
    double* timesOfFlight = new double[nHitPMT];
    // Loop to find expected time of flight from vertex to PMT for given Cherenkov angle
    for(int iPMT=0; iPMT< nHitPMT; iPMT++)
        timesOfFlight[iPMT] = FullTimeOfFlight(vtxX, vtxY, vtxZ, dirX, dirY, dirZ, iPMT, cherenkovAngle);
    double result = Goodness(timesOfFlight, time);
    //Penalise cherenkov angles far from 42 deg
    delete[] timesOfFlight;
    //cout << vtxX << " " << vtxY << " " << vtxZ << " " << time << " " << dirCosTheta << " " << dirPhi << " " << cherenkovAngle << " " << nHitPMT << " " << result << endl;
    return -result; //negative for minimise not maximise
}

TH3D* HighEReco::FindTracks(double vtxX, double vtxY, double vtxZ){
    double cherenkovAngle = TMath::ACos(1.0/N_REF);
    double sinChkvAngle = TMath::Sin(cherenkovAngle);
    double cosChkvAngle = TMath::Cos(cherenkovAngle);
    const int xBins = 100;
    const int yBins = 100;
    const int zBins = 200;
    //double track[xBins][yBins][zBins];
    TH3D *trackHist = new TH3D("trackHist","trackHist",xBins,-550,550,yBins,-550,550,zBins,-1100,1100);
    // Loop over PEs to find angle between track and direction from vertex to PMT assuming fixed Cherenkov angle
    for(int iPE=0; iPE<nPEs; iPE++){
        int pmt = hitPMT[iPE];
        double vtxToPMTx = hitPMTx[pmt]-vtxX;
        double vtxToPMTy = hitPMTy[pmt]-vtxY;
        double vtxToPMTz = hitPMTz[pmt]-vtxZ;
        //Distance from vertex to PMT:
        double pmtDist = TMath::Sqrt(vtxToPMTx*vtxToPMTx+vtxToPMTy*vtxToPMTy+vtxToPMTz*vtxToPMTz);

        //To find angle alpha, solve equation for TOF: t=d/(c*sin(thetaC)) * (sin(thetaC - alpha) + n_ref*sin(alpha))
        //Solution is:
        // tan(alpha)=(A-B*sqrt(A^2-B^2+1))/(B^2-A^2)
        //where
        // A=(n_ref-cos(thetaC))/sin(thetaC)
        //  is the difference between TOF of lepton and TOF of photon divided by distance to PMT
        // B=c*t/d
        //  is the ratio of c times TOF and distance from vertex to PMT
        double A = (N_REF-cosChkvAngle)/sinChkvAngle;
        double B = C_VAC*hitT[iPE]/pmtDist;
        // Must have B>1 and A^2-B^2+1>0 for physical solution.
        if(B<1 || A*A-B*B+1<0)
            continue;
        double pmtTrackAngle= TMath::ATan((A-B* TMath::Sqrt(A*A-B*B+1))/(B*B-A*A));
        double trackDist = pmtDist* TMath::Sin(cherenkovAngle-pmtTrackAngle)/sinChkvAngle;
        //For single PE with fixed pmtTrackAngle there is circle of possible directions
        //Scale number of points since circle with large pmtTrackAngle will cover more bins
        //Radius of circle on sphere is proportional to sin of half opening angle
        int circlePoints = (int) (1000* TMath::Sin(pmtTrackAngle));
        for(int iPhi=0; iPhi<circlePoints; iPhi++){
            //Calculate point on circle pmtTrackAngle away from vertex-PMT direction
            double Phi = iPhi*2.* TMath::Pi()/circlePoints;
            double cosPhi = TMath::Cos(Phi);
            double sinPhi = TMath::Sin(Phi);
            double cosAlpha = TMath::Cos(pmtTrackAngle);
            double sinAlpha = TMath::Sin(pmtTrackAngle);
            double phi0 = TMath::ATan2(vtxToPMTy,vtxToPMTx);
            double theta0 = TMath::ACos(vtxToPMTz/pmtDist);
            double cosphi0 = TMath::Cos(phi0);
            double sinphi0 = TMath::Sin(phi0);
            double costheta0 = TMath::Cos(theta0);
            double sintheta0 = TMath::Sin(theta0);
            double x = costheta0*cosphi0*sinAlpha*cosPhi - sinphi0*sinAlpha*sinPhi + sintheta0*cosphi0*cosAlpha;
            double y = sinphi0*costheta0*sinAlpha*cosPhi + cosphi0*sinAlpha*sinPhi + sintheta0*sinphi0*cosAlpha;
            double z = costheta0*cosAlpha-sintheta0*sinAlpha*cosPhi;
            //double trackDist = pmtDist*TMath::Sin(cherenkovAngle-pmtTrackAngle)/sinChkvAngle;
            x *= trackDist;
            y *= trackDist;
            z *= trackDist;
            trackHist->Fill(vtxX+x,vtxY+y,vtxZ+z);
        }
    }
    return trackHist;
}

// Find rings by assuming photons emitted at cherenkov angle from vertex
int HighEReco::FindRing(double vtxX, double vtxY, double vtxZ, double vtxT, double &thetaPeak, double &phiPeak,
                        int ringNumber, bool useTrack, int *hough) {
    double cherenkovAngle = TMath::ACos(1.0/N_REF);
    double sinChkvAngle = TMath::Sin(cherenkovAngle);
    double cosChkvAngle = TMath::Cos(cherenkovAngle);
    double A = (N_REF-cosChkvAngle)/sinChkvAngle; //See below for explanation of use
    const int nBins = 2500;
    double theta[nBins];
    double phi[nBins];
    //Distribute bins evenly around unit sphere using spiral method
    //Calculate shift to correct for end bins covering larger area
    double shift = 0.5;
    //Use Newton's method with initial guess 0.5 to solve equation:
    //  a*sin(a) - pi/N = 0
    //where
    //  a=sqrt(pi/(N-x))
    //  x is correct shift
    //only need one iteration for good enough estimate
    double a= TMath::Sqrt(TMath::Pi()/(nBins-shift));
    shift -= (a* TMath::Sin(a)- TMath::Pi()/nBins)*2.0*(nBins-shift)/(a*(TMath::Sin(a)+a* TMath::Cos(a)));
    double dz = 2.0/(nBins-shift); //Change in z=cos(theta) for consecutive bins
    double s = TMath::Sqrt(4* TMath::Pi()/(nBins-shift)); //Approx distance between bins
    double sOverdz = s/dz;
    bool isSaved = hough !=0;
    if(!isSaved) hough = new int[2500];
    //Calculate theta and phi position of each bin according to position on spiral
    for(int iBin=0; iBin<nBins; iBin++){
        theta[iBin] = TMath::ACos(1-dz*(iBin+0.5*(1-shift)));
        phi[iBin] = fmod(theta[iBin]*sOverdz+ TMath::Pi(), TMath::TwoPi())- TMath::Pi();
        hough[iBin] = 0;
    }
    // Loop over PEs to add circle of possible directions for each PE
    for(int iPE=0; iPE<nPEs; iPE++){
        if(hitRing[iPE]<ringNumber) continue; //ignore PEs identified as part of previous ring
        int pmt = hitPMT[iPE];
        double vtxToPMTx = hitPMTx[pmt]-vtxX;
        double vtxToPMTy = hitPMTy[pmt]-vtxY;
        double vtxToPMTz = hitPMTz[pmt]-vtxZ;
        //Distance from vertex to PMT:
        double pmtDist = TMath::Sqrt(vtxToPMTx*vtxToPMTx+vtxToPMTy*vtxToPMTy+vtxToPMTz*vtxToPMTz);

        //Just use cherenkov angle if not using track
        double pmtTrackAngle = cherenkovAngle;
        //Otherwise, use track calculation
        if(useTrack){
            //To find angle, solve TOF equation for alpha:
            //  t=d/(c*sin(thetaC)) * (sin(thetaC - alpha) + n_ref*sin(alpha))
            //Solution is:
            //  tan(alpha)=(A-B*sqrt(A^2-B^2+1))/(B^2-A^2)
            //where
            //  A=(n_ref-cos(thetaC))/sin(thetaC)
            //  B=c*t/d
            double B = C_VAC*(hitT[iPE]-vtxT)/pmtDist;
            // Must have B>1 and A^2-B^2+1>0 for physical solution.
            if(B<1 || A*A-B*B+1<0)
                continue;
            pmtTrackAngle= TMath::ATan((A-B* TMath::Sqrt(A*A-B*B+1))/(B*B-A*A));
        }

        //For single PMT PE with fixed pmtTrackAngle there is circle of possible directions
        //Scale number of points since circle with large pmtTrackAngle will cover more bins
        //Radius of circle on sphere is proportional to sin of half opening angle
        int circlePoints = (int) (1000* TMath::Sin(pmtTrackAngle));

        //Pre-calculate useful quantities
        double cosAlpha = TMath::Cos(pmtTrackAngle);
        double sinAlpha = TMath::Sin(pmtTrackAngle);
        //phi0 and theta0 are for vector pointing from vertex to PMT
        double phi0 = TMath::ATan2(vtxToPMTy,vtxToPMTx);
        double theta0 = TMath::ACos(vtxToPMTz/pmtDist);
        double cosphi0 = TMath::Cos(phi0);
        double sinphi0 = TMath::Sin(phi0);
        double costheta0 = TMath::Cos(theta0);
        double sintheta0 = TMath::Sin(theta0);
        //int prevBin = -1;
        for(int iPhi=0; iPhi<circlePoints; iPhi++){
            //Calculate point on circle pmtTrackAngle away from vertex-PMT direction (theat0,phi0)
            double Phi = iPhi*2.* TMath::Pi()/circlePoints;
            double cosPhi = TMath::Cos(Phi);
            double sinPhi = TMath::Sin(Phi);
            double x = costheta0*cosphi0*sinAlpha*cosPhi - sinphi0*sinAlpha*sinPhi + sintheta0*cosphi0*cosAlpha;
            double y = sinphi0*costheta0*sinAlpha*cosPhi + cosphi0*sinAlpha*sinPhi + sintheta0*sinphi0*cosAlpha;
            double z = costheta0*cosAlpha-sintheta0*sinAlpha*cosPhi;
            double newTheta = TMath::ACos(z);
            double newPhi = TMath::ATan2(y,x);
            // Find correct bin
            int m = (int) ((newTheta*sOverdz-newPhi)/ TMath::TwoPi()+0.5);
            double binTheta = (m* TMath::TwoPi()+newPhi)/sOverdz;
            int bin = (int) ((1.0- TMath::Cos(binTheta))/dz+0.5*shift);
            //if(bin==prevBin) continue;
            // Add to bin
            hough[bin]++;
            //prevBin=bin;
        }
    }
    //Find peak (just use bin with highest value)
    int peakBin = 0;
    int peakCount = 0;
    for(int iBin=0; iBin<nBins; iBin++){
        if(hough[iBin]>peakCount){
            peakBin=iBin;
            peakCount= hough[iBin];
        }
    }
    thetaPeak = theta[peakBin];
    phiPeak = phi[peakBin];
    double dirX = TMath::Sin(thetaPeak)* TMath::Cos(phiPeak);
    double dirY = TMath::Sin(thetaPeak)* TMath::Sin(phiPeak);
    double dirZ = TMath::Cos(thetaPeak);
    //Identify PEs from this ring
    int nRingPEs=0;
    for(int iPE=0; iPE<nPEs; iPE++){
        int pmt = hitPMT[iPE];
        if(hitRing[iPE]<ringNumber) continue;
        //calculate expected time of flight
        double expectedTime = FullTimeOfFlight(vtxX, vtxY, vtxZ, dirX, dirY, dirZ, pmt, cherenkovAngle);
        if(TMath::Abs(hitT[iPE]-vtxT-expectedTime)<2.0*hitPMTTimeRes[pmt]){
            hitRing[iPE]=ringNumber;
            nRingPEs++;
        }
    }
    cout << "Ring " << ringNumber << " found with direction: (" << dirX << ", " << dirY << ", " << dirZ << ") (" << thetaPeak << ", " << phiPeak << ") PEs: " << nRingPEs << endl;

//  TGraph2D *gr = new TGraph2D();
//  for(int iBin=0; iBin<nBins; iBin++)
//    gr->SetPoint(iBin,phi[iBin],theta[iBin],hough[iBin]);
//  return gr;
    if(!isSaved) delete[] hough;
    return nRingPEs;
}

int HighEReco::FindRings(double vtxX, double vtxY, double vtxZ, double vtxT, double *thetaPeaks, double *phiPeaks,
                         int *ringPEs, int maxRings, bool useTrack, int*hough) {
    double * tmpThetaPeaks = new double[maxRings];
    double * tmpPhiPeaks = new double[maxRings];
    int * tmpRingPEs = new int[maxRings];
    //Keep looking for rings until few enough PEs in next ring
    int nrings;
    int *saveHough = hough;
    for(nrings=0; nrings<maxRings; nrings++){
        tmpRingPEs[nrings]= FindRing(vtxX, vtxY, vtxZ, vtxT, tmpThetaPeaks[nrings], tmpPhiPeaks[nrings], nrings + 1,
                                     useTrack, saveHough);
        saveHough = 0;
        if(tmpRingPEs[nrings]<0.09*tmpRingPEs[0]) break;
    }

    if(nrings <2){
        ringPEs[0] = tmpRingPEs[0];
        thetaPeaks[0] = tmpThetaPeaks[0];
        phiPeaks[0] = tmpPhiPeaks[0];
    }
    else {
        //Redistribute hits to rings
/*
        double cherenkovAngle = TMath::ACos(1.0 / N_REF);
        double *dirX = new double[nrings];
        double *dirY = new double[nrings];
        double *dirZ = new double[nrings];
        for (int iRing = 0; iRing < nrings; iRing++) {
            dirX[iRing] = TMath::Sin(tmpThetaPeaks[iRing]) * TMath::Cos(tmpPhiPeaks[iRing]);
            dirY[iRing] = TMath::Sin(tmpThetaPeaks[iRing]) * TMath::Sin(tmpPhiPeaks[iRing]);
            dirZ[iRing] = TMath::Cos(tmpThetaPeaks[iRing]);
        }
        for (int iPE = 0; iPE < nPEs; iPE++) {
            int pmt = hitPMT[iPE];
            int ring = hitRing[iPE]-1;
            double expectedTime = FullTimeOfFlight(vtxX, vtxY, vtxZ, dirX[ring], dirY[ring], dirZ[ring], pmt,
                                                   cherenkovAngle);
            for (int newRing = ring+1; newRing < nrings; newRing++) {
                double newExpectedTime = FullTimeOfFlight(vtxX, vtxY, vtxZ, dirX[newRing], dirY[newRing], dirZ[newRing],
                                                          pmt, cherenkovAngle);
                if (TMath::Abs(newExpectedTime - hitT[iPE] - vtxT) < TMath::Abs(expectedTime - hitT[iPE] - vtxT)) {
                    tmpRingPEs[hitRing[iPE]-1]--;
                    tmpRingPEs[newRing]++;
                    hitRing[iPE] = newRing + 1;
                    expectedTime = newExpectedTime;
                }
            }
        }
        delete[] dirX;
        delete[] dirY;
        delete[] dirZ;
*/
        //Reorder rings to be in order of number of hits
        for (int iRing = 0; iRing < nrings; iRing++) ringPEs[iRing] = tmpRingPEs[iRing];
        std::sort(ringPEs, ringPEs + nrings, std::greater<int>());
        int *newOrder = new int[nrings];
        for (int newRingIndex = 0; newRingIndex < nrings; newRingIndex++) {
            int oldRingIndex = 0;
            while (tmpRingPEs[oldRingIndex] != ringPEs[newRingIndex]) oldRingIndex++;
            thetaPeaks[newRingIndex] = tmpThetaPeaks[oldRingIndex];
            phiPeaks[newRingIndex] = tmpPhiPeaks[oldRingIndex];
            tmpRingPEs[oldRingIndex] = -1;
            newOrder[oldRingIndex] = newRingIndex;
//            cout << oldRingIndex+1 << " " << newRingIndex+1 << endl;
        }
        for (int iPE = 0; iPE < nPEs; iPE++) {
            if(hitRing[iPE]<=nrings)
            hitRing[iPE] = newOrder[hitRing[iPE]-1] + 1;
        }
        delete[] newOrder;
    }
/*
    for (int iRing = 0; iRing < nrings; iRing++){
        ringPEs[iRing] = tmpRingPEs[iRing];
        thetaPeaks[iRing] = tmpThetaPeaks[iRing];
        phiPeaks[iRing] = tmpPhiPeaks[iRing];
    }
*/
    delete[] tmpThetaPeaks;
    delete[] tmpPhiPeaks;
    delete[] tmpRingPEs;
    return nrings;
}

// Calculate expected number of photoelectrons based on Cherenkov pattern
double HighEReco::ExpectedPMTPhotoelectrons(double vtxX, double vtxY, double vtxZ, double dirX, double dirY, double dirZ, int pmt, bool isHit, double kineticEnergy, int ipnu, double *expectedTime){
    //cout << "test" << hitPMTDirX[0]<<endl;i
    TH3D * photonsLookup = (ipnu == 11) ? electronPhotons : muonPhotons;
    TH3D * timesLookup = (ipnu == 11) ? electronTimes : muonTimes;
    double vtxToPMTx, vtxToPMTy, vtxToPMTz;
    if (isHit) {
        vtxToPMTx = hitPMTx[pmt] - vtxX;
        vtxToPMTy = hitPMTy[pmt] - vtxY;
        vtxToPMTz = hitPMTz[pmt] - vtxZ;
    }
    else {
        vtxToPMTx = unhitPMTx[pmt] - vtxX;
        vtxToPMTy = unhitPMTy[pmt] - vtxY;
        vtxToPMTz = unhitPMTz[pmt] - vtxZ;
    }
    //Distance from vertex to PMT:
    double pmtDist = TMath::Sqrt(vtxToPMTx*vtxToPMTx+vtxToPMTy*vtxToPMTy+vtxToPMTz*vtxToPMTz);
    //Angle between track direction and vertex-PMT direction:
    double pmtTrackAngle = TMath::ACos((vtxToPMTx*dirX+vtxToPMTy*dirY+vtxToPMTz*dirZ)/pmtDist);
    double minDist = photonsLookup->GetXaxis()->GetBinCenter(1);
    if(pmtDist <= minDist) pmtDist = minDist+0.0001;
    double maxDist = photonsLookup->GetXaxis()->GetBinCenter(photonsLookup->GetNbinsX());
    if(pmtDist >= maxDist) pmtDist = maxDist-0.0001;
    double minAngle = photonsLookup->GetYaxis()->GetBinCenter(1);
    if(pmtTrackAngle <= minAngle) pmtTrackAngle = minAngle+0.0001;
    double maxAngle = photonsLookup->GetYaxis()->GetBinCenter(photonsLookup->GetNbinsY());
    if(pmtTrackAngle >= maxAngle) pmtTrackAngle = maxAngle-0.0001;
    double minKE = photonsLookup->GetZaxis()->GetBinCenter(1);
    if(kineticEnergy <= minKE) kineticEnergy= minKE+0.0001;
    double maxKE = photonsLookup->GetZaxis()->GetBinCenter(photonsLookup->GetNbinsZ());
    if(kineticEnergy >= maxKE) kineticEnergy = maxKE-0.0001;
//    if(pmt==0) cout << pmtDist << " " << pmtTrackAngle << " " << kineticEnergy << endl;

    double expectedPhotons = photonsLookup->Interpolate(pmtDist, pmtTrackAngle, kineticEnergy);//*100000000;//Remove this factor when likelihood tables normalisation is fixed
    if(expectedTime) *expectedTime = timesLookup->Interpolate(pmtDist, pmtTrackAngle, kineticEnergy);

    double QE = 2.*0.22; //WChSandBox has half as many photons so need so double QE

    //Account for solid angle covered by PMT at given distance
    double PMTsize = 12*2.54/2.; //12 inch diameter
    double solidAngleFactor = (1.-pmtDist/ TMath::Sqrt(pmtDist*pmtDist+PMTsize*PMTsize)); // solid angle covered by PMT
    double halfAngleBinWidth = photonsLookup->GetYaxis()->GetBinWidth(photonsLookup->GetYaxis()->FindFixBin(pmtTrackAngle))/2.;
    solidAngleFactor /= 2*TMath::Sin(pmtTrackAngle)*TMath::Sin(halfAngleBinWidth); //solid angle covered by bin in theta

    //Account for reduction of solid angle covered by PMT not directly facing towards vertex
    double cosPMTangle = isHit ?
                         -(vtxToPMTx* hitPMTDirX[pmt] + vtxToPMTy* hitPMTDirY[pmt] + vtxToPMTz* hitPMTDirZ[pmt])/pmtDist :
                         -(vtxToPMTx* unhitPMTDirX[pmt] + vtxToPMTy* unhitPMTDirY[pmt] + vtxToPMTz* unhitPMTDirZ[pmt])/pmtDist;

    double expectedPhotoelectrons = expectedPhotons*QE*solidAngleFactor*cosPMTangle;
    if(expectedPhotoelectrons<0.00001) expectedPhotoelectrons = 0.00001;
    return expectedPhotoelectrons;
}

// Calculate likelihood of number of photoelectrons based on Poisson
double HighEReco::PMTlnLikelihood(double vtxX, double vtxY, double vtxZ, double dirX, double dirY, double dirZ, int pmt, bool isHit, double kineticEnergy, int ipnu){
    double expectedPhotoelectrons = ExpectedPMTPhotoelectrons(vtxX, vtxY, vtxZ, dirX, dirY, dirZ, pmt, isHit, kineticEnergy, ipnu);
    double observedPhotoelectrons = hitPMTPEs[pmt];
    //if(expectedPhotoelectrons < 2 && observedPhotoelectrons < 2) return 0;
    double l = observedPhotoelectrons* TMath::Log(expectedPhotoelectrons)-expectedPhotoelectrons;
    return l;
}

// Calculate total likelihood for all PMTs
double HighEReco::LnLikelihood(const double *par, int ipnu, bool total, bool print) {
    double vtxZ = par[1];
    double vtxY = par[2];
    double vtxX = par[3];
    double time = par[4];
    double dirCosTheta = par[5];
    double dirSinTheta = TMath::Sqrt(1-dirCosTheta*dirCosTheta);
    double dirPhi = par[6];
    double dirX = dirSinTheta* TMath::Cos(dirPhi);
    double dirY = dirSinTheta* TMath::Sin(dirPhi);
    double dirZ = dirCosTheta;
    double kineticEnergy = par[7];
    double trackOffset = par[0];
    vtxZ -= dirZ*trackOffset;
    vtxY -= dirY*trackOffset;
    vtxX -= dirX*trackOffset;
    if(vtxX*vtxX+vtxY*vtxY>550*550 || TMath::Abs(vtxZ)>1100
        || vtxX!=vtxX || vtxY!=vtxY || vtxZ!=vtxZ || time!=time || dirCosTheta!=dirCosTheta
        || dirPhi!=dirPhi || kineticEnergy!=kineticEnergy || trackOffset!=trackOffset
         ){
        return -9999999999;
    }
    double lnLikelihood = 0;
    double *hitExpectedPE = new double[nHitPMT];
    double *hitExpectedT = new double[nHitPMT];
    double totalExpectedPE = 0;
    double totalObservedPE = 0;
    //Hit PMTs
    for(int iPMT = 0; iPMT< nHitPMT; iPMT++){
        hitExpectedPE[iPMT] = ExpectedPMTPhotoelectrons(vtxX,vtxY,vtxZ,dirX,dirY,dirZ,iPMT,true,kineticEnergy,ipnu,hitExpectedT+iPMT);
        totalExpectedPE += hitExpectedPE[iPMT];
        totalObservedPE += hitPMTPEs[iPMT];
        if(print)cout<<iPMT<<" "<<hitPMTx[iPMT]<<" "<<hitPMTy[iPMT]<<" "<<hitPMTz[iPMT]<<" "<<hitPMTDirX[iPMT]<<" "<<hitPMTDirY[iPMT]<<" "<<hitPMTDirZ[iPMT]<<" "<<hitExpectedPE[iPMT]<<" "<<hitPMTPEs[iPMT]<<endl;
    }
    //Unhit PMTs
    double *unhitExpectedPE = new double[nUnhitPMT];
    for(int iPMT = 0; iPMT< nUnhitPMT; iPMT++){
        unhitExpectedPE[iPMT] = ExpectedPMTPhotoelectrons(vtxX,vtxY,vtxZ,dirX,dirY,dirZ,iPMT,false,kineticEnergy,ipnu);
        if(print)cout<<iPMT<<" "<<unhitPMTx[iPMT]<<" "<<unhitPMTy[iPMT]<<" "<<unhitPMTz[iPMT]<<" "<<unhitPMTDirX[iPMT]<<" "<<unhitPMTDirY[iPMT]<<" "<<unhitPMTDirZ[iPMT]<<" "<<unhitExpectedPE[iPMT]<<" "<<0<<endl;
    }
    //Include only PMTs with highest expected # of hits (add 10% as many unhit PMTs as hit PMTs)
    int n;// = (nHitPMT + 1)/10;
//    cout << totalExpectedPE << " " << n << " ";
//    if(n>nUnhitPMT)
        n=nUnhitPMT;
//    else //std::nth_element sorts the list such that the nth element is in correct place and all smaller elements are before and all larger elements are after
//        nth_element(unhitExpectedPE,unhitExpectedPE+nUnhitPMT-n,unhitExpectedPE+nUnhitPMT);
    for(int iPMT = nUnhitPMT-n; iPMT<nUnhitPMT; iPMT++) //use PMTs with larger expected PE than nth element
        totalExpectedPE += unhitExpectedPE[iPMT];
    if(print) cout << totalExpectedPE << " " << totalObservedPE << endl;
    if(!total){
        double normalisation = 1./totalExpectedPE;
        for(int iPMT = 0; iPMT< nHitPMT; iPMT++){
            //Poisson:
            //lnLikelihood += hitPMTPEs[iPMT]*TMath::Log(hitExpectedPE[iPMT]);

            //Multinomial:
            lnLikelihood += hitPMTPEs[iPMT]*TMath::Log(hitExpectedPE[iPMT]*normalisation);

            //Negative Binomial, r=2:
            //const double r = 1000;
            //double p = 1/(1+r/(normalisation*hitExpectedPE[iPMT]));
            //lnLikelihood += hitPMTPEs[iPMT]*TMath::Log(p)+r*TMath::Log(1-p);
        }
        //Poisson:
        //lnLikelihood -= totalExpectedPE;
        for(int iPE = 0; iPE<nPEs; iPE++){
            //Gaussian for time:
            const double sigma = 4.;
            lnLikelihood -= 0.5*TMath::Power((hitT[iPE]-hitExpectedT[hitPMT[iPE]])/(sigma+hitPMTTimeRes[hitPMT[iPE]]),2);
        }
    } else {
        //Total PE term
        if(ipnu == 13 && totalExpectedPE>600) totalExpectedPE = 2.*totalExpectedPE - 600.;
        else if(ipnu==11 && totalExpectedPE>500/0.6) totalExpectedPE = 1.6*totalExpectedPE - 500.;
        lnLikelihood += totalObservedPE* TMath::Log(totalExpectedPE) - totalExpectedPE;
    }
    delete[] hitExpectedPE;
    delete[] hitExpectedT;
    delete[] unhitExpectedPE;
    return lnLikelihood;
}

double HighEReco::ElectronLnLikelihood(const double *par){
    return -LnLikelihood(par, 11, false, false);
}

double HighEReco::MuonLnLikelihood(const double *par){
//    for(int i =0; i<8; i++) cout << par[i] << " ";
    double result = -LnLikelihood(par, 13, false, false);
//    cout << result << endl;
    return result;
}

double HighEReco::ElectronLnLikelihoodTotal(const double *par){
    return -LnLikelihood(par, 11, true, false);
}

double HighEReco::MuonLnLikelihoodTotal(const double *par){
    return -LnLikelihood(par, 13, true, false);
}

void HighEReco::PointFit(double &recoVtxX, double &recoVtxY, double &recoVtxZ, double &recoT, double &cherenkovAngle) {
    double * timesOfFlight= new double[nHitPMT];
    cout<<"Init vtx: ("<<recoVtxX<<" "<<recoVtxY<<" "<<recoVtxZ<<") t:"  << 0. << endl;

    cout << "Point fit to find vertex" << endl;
    // Set up Minuit
    ROOT::Math::Minimizer*minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(1000000);
    minimizer->SetTolerance(0.0001);
    minimizer->SetPrintLevel(0);
    //ROOT::Math::GradFunctor f(&PointGoodness, &PointGoodnessDeriv, 4);
    ROOT::Math::Functor f(this, &HighEReco::PointGoodness, 3);
    minimizer->SetFunction(f);
    minimizer->SetLimitedVariable(0,"vtxZ",recoVtxZ,1,-1100,1100);
    minimizer->SetLimitedVariable(1,"vtxY",recoVtxY,1,-550,550);
    minimizer->SetLimitedVariable(2,"vtxX",recoVtxX,1,-550,550);
    // Perform minimization
    minimizer->Minimize();
    recoVtxZ = minimizer->X()[0];
    recoVtxY = minimizer->X()[1];
    recoVtxX = minimizer->X()[2];
    cout<<"Reco vtx: ("<<recoVtxX<<" "<<recoVtxY<<" "<<recoVtxZ<<")" << endl;
    cherenkovAngle = TMath::ACos(1.0/N_REF);

    // Find interaction time from mean PE timing residual
    cout << "Estimate interaction time" << endl;
    recoT = 0;//Loop over PMTs to find expected TOF from reco vertex
    for(int iPMT=0; iPMT< nHitPMT; iPMT++)
    {
        double vtxToPMTx = hitPMTx[iPMT]-recoVtxX;
        double vtxToPMTy = hitPMTy[iPMT]-recoVtxY;
        double vtxToPMTz = hitPMTz[iPMT]-recoVtxZ;
        double pmtDist = TMath::Sqrt(vtxToPMTx*vtxToPMTx+vtxToPMTy*vtxToPMTy+vtxToPMTz*vtxToPMTz);
        timesOfFlight[iPMT] = N_REF*pmtDist/C_VAC;
    }
    //Loop over PEs to find mean timing residual (hit time minus TOF)
    for(int iPE=0; iPE<nPEs; iPE++)
    {
        int pmt = hitPMT[iPE];
        recoT += hitT[iPE]-timesOfFlight[pmt];
    }
    recoT /= nPEs;
    cout << "Reco t: " << recoT << endl;
    delete[] timesOfFlight;
}

void HighEReco::TrackFit(double &recoVtxX, double &recoVtxY, double &recoVtxZ, double &recoT, double &recoDirPhi,
              double &recoDirTheta, double &recoDirCosTheta, double &cherenkovAngle) {
    //double recoPar2[7]={recoT,cherenkovAngle,recoDirCosTheta,recoDirPhi,recoVtxZ,recoVtxY,recoVtxX};
//    TrackGoodness(recoPar2);
//  trueG = -TrackGoodness(truePar);
//  cout<<"   Goodness True: "<<trueG<<"  Reco: "<<recoG<<endl;
    cout << "Fit vertex including track" << endl;
    //Set up minimizer:
    //ROOT::Math::Minimizer*minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    ROOT::Math::Minimizer*minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(1000000);
    minimizer->SetTolerance(0.0001);
    minimizer->SetPrintLevel(0);
    ROOT::Math::Functor f2(this,&HighEReco::TrackGoodness,7);
    minimizer->SetFunction(f2);
    cout<<"Init vtx: ("<<recoVtxX<<" "<<recoVtxY<<" "<<recoVtxZ<<") t:" << recoT
    << " dir: ("<< TMath::Sin(recoDirTheta)* TMath::Cos(recoDirPhi)<<" "<< TMath::Sin(recoDirTheta)* TMath::Sin(recoDirPhi)<<" "<<
    TMath::Cos(recoDirTheta)<<")"<< endl;
    minimizer->SetLimitedVariable(0,"T",recoT,1.,-100,100);
    minimizer->SetLimitedVariable(1,"chkvAngle",cherenkovAngle,0.1,0.35,1.1);
    minimizer->SetLimitedVariable(2,"dirCosTheta",recoDirCosTheta,0.1,-0.99999,0.99999);
    minimizer->SetLimitedVariable(3,"dirPhi",recoDirPhi,0.1,-TMath::Pi(), TMath::Pi());
    minimizer->SetLimitedVariable(4,"vtxZ",recoVtxZ,1.,-1100,1100);
    minimizer->SetLimitedVariable(5,"vtxY",recoVtxY,1.,-550,550);
    minimizer->SetLimitedVariable(6,"vtxX",recoVtxX,1.,-550,550);
    //Loop to iterate fitting direction then vertex
    for(int iteration = 0; iteration<2; iteration++){
        //Fix cherenkov angle
//        minimizer->FixVariable(1);
        //Fix vertex
        minimizer->FixVariable(4);
        minimizer->FixVariable(5);
        minimizer->FixVariable(6);
        //Perform fit
        cout << "Fit dir" << endl;
        minimizer->Minimize();
        //Release vertex
        minimizer->ReleaseVariable(4);
        minimizer->ReleaseVariable(5);
        minimizer->ReleaseVariable(6);
        //Fix direction
        minimizer->FixVariable(2);
        minimizer->FixVariable(3);
        //Perform fit
        cout << "Fit vtx" << endl;
        minimizer->Minimize();
        //Release direction
        minimizer->ReleaseVariable(2);
        minimizer->ReleaseVariable(3);
    }
    //Get fitted values
    recoVtxZ = minimizer->X()[4];
    recoVtxY = minimizer->X()[5];
    recoVtxX = minimizer->X()[6];
    recoT = minimizer->X()[0];
    recoDirCosTheta = minimizer->X()[2];
    recoDirPhi = minimizer->X()[3];
    recoDirTheta= TMath::ACos(recoDirCosTheta);
    cherenkovAngle = minimizer->X()[1];
    //Print
//  cout<<"True vtx: ("<<trueVtxX<<" "<<trueVtxY<<" "<<trueVtxZ<<") t" << trueT
//      <<" dir: ("<< TMath::Sin(trueDirTheta)* TMath::Cos(trueDirPhi)<<" "<< TMath::Sin(trueDirTheta)* TMath::Sin(trueDirPhi)<<" "<<
//  TMath::Cos(trueDirTheta)<<")"<<endl;
    cout<<"Reco vtx: ("<<recoVtxX<<" "<<recoVtxY<<" "<<recoVtxZ<<") t" << recoT
    <<" dir: ("<< TMath::Sin(recoDirTheta)* TMath::Cos(recoDirPhi)<<" "<< TMath::Sin(recoDirTheta)* TMath::Sin(recoDirPhi)<<" "<<
    TMath::Cos(recoDirTheta)<<")"
    <<" thetaC: " << cherenkovAngle*180.0/ TMath::Pi() << endl;
    //cout<<"rErr vtx: ("<<recoErrVtxX<<" "<<recoErrVtxY<<" "<<recoErrVtxZ<<")"<<" dir: (c"<<recoErrDirCosTheta<<" "<<recoErrDirPhi<<")"<<endl;
    //double recoPar4[7]={recoT,cherenkovAngle,recoDirCosTheta,recoDirPhi,recoVtxZ,recoVtxY,recoVtxX};
//    TrackGoodness(recoPar4);
//  cout<<"Goodness True: "<<trueG<<"  Reco: "<<recoG<<endl<<endl;

}

void HighEReco::LikelihoodFit(double &trackCorrection, double &recoVtxX, double &recoVtxY,
                   double &recoVtxZ, double &recoT, double &recoDirPhi,
                   double &recoDirTheta, double &recoKE, double &recoLnL, int ipnu) {
  //Set up minimizer:
    ROOT::Math::Minimizer*minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(1000000);
//    minimizer->SetTolerance(0.0001);
    minimizer->SetPrintLevel(0);

    ROOT::Math::Functor func;
    if(ipnu==11) {
        func = ROOT::Math::Functor(this,&HighEReco::ElectronLnLikelihood,8);
    }
    else{
        func = ROOT::Math::Functor(this,&HighEReco::MuonLnLikelihood,8);
    }
    minimizer->SetFunction(func);
    double recoDirCosTheta = TMath::Cos(recoDirTheta);
    cout<<"Init vtx: ("<<recoVtxX<<" "<<recoVtxY<<" "<<recoVtxZ<<") t:" << recoT
    << " dir: ("<< TMath::Sin(recoDirTheta)* TMath::Cos(recoDirPhi)<<" "<< TMath::Sin(recoDirTheta)* TMath::Sin(recoDirPhi)<<" "<<
    recoDirCosTheta <<")(" << recoDirCosTheta << "," << recoDirPhi << ") KE: " << recoKE << " Offset: " << trackCorrection << endl;
//    double initPar[8]={trackCorrection,recoVtxZ,recoVtxY,recoVtxX,recoT,recoDirCosTheta,recoDirPhi,recoKE};
//    LnLikelihood(initPar, ipnu,false,true);
    if(recoKE < 50) recoKE = 51;
    if(recoKE > 3000) recoKE = 2999;
    minimizer->SetLimitedVariable(0,"trackOffset",trackCorrection,10,-1000,1000);
    minimizer->SetLimitedVariable(1,"vtxZ",recoVtxZ,1.,-1100,1100);
    minimizer->SetLimitedVariable(2,"vtxY",recoVtxY,1.,-550,550);
    minimizer->SetLimitedVariable(3,"vtxX",recoVtxX,1.,-550,550);
    minimizer->SetLimitedVariable(4,"T",recoT,1.,-100,100);
    minimizer->SetLimitedVariable(5,"dirCosTheta", recoDirCosTheta,0.1,-0.99999,0.99999);
    minimizer->SetLimitedVariable(6,"dirPhi",recoDirPhi,0.1,-TMath::Pi(), TMath::Pi());
    minimizer->SetLimitedVariable(7,"kineticEnergy",recoKE,10,50,3000);
    //Fix direction, vertex
//    minimizer->FixVariable(0);
    minimizer->FixVariable(1);
    minimizer->FixVariable(2);
    minimizer->FixVariable(3);
    minimizer->FixVariable(4);
    minimizer->FixVariable(5);
    minimizer->FixVariable(6);
    minimizer->FixVariable(7);
    //Perform fit
//ls    minimizer->Minimize();
//    minimizer->ReleaseVariable(7);
//    minimizer->Minimize();
    minimizer->ReleaseVariable(5);
    minimizer->ReleaseVariable(6);
    minimizer->Minimize();
    minimizer->ReleaseVariable(1);
    minimizer->ReleaseVariable(2);
    minimizer->ReleaseVariable(3);
    minimizer->FixVariable(0);
    minimizer->Minimize();
/*
    ROOT::Math::Functor fElectronTotal(&ElectronLnLikelihoodTotal,8);
    minimizer->SetFunction(fElectronTotal);
    minimizer->FixVariable(0);
    minimizer->FixVariable(1);
    minimizer->FixVariable(2);
    minimizer->FixVariable(3);
    minimizer->FixVariable(4);
    minimizer->FixVariable(5);
//    minimizer->FixVariable(6);
    minimizer->FixVariable(7);
    minimizer->ReleaseVariable(6);
    minimizer->Minimize();
*/
    recoVtxZ = minimizer->X()[1];
    recoVtxY = minimizer->X()[2];
    recoVtxX = minimizer->X()[3];
    recoT = minimizer->X()[4]- minimizer->X()[0]/C_VAC;
    recoDirCosTheta = minimizer->X()[5];
    recoDirPhi = minimizer->X()[6];
    recoDirTheta= TMath::ACos(recoDirCosTheta);
    recoKE = minimizer->X()[7];
    recoVtxZ -= minimizer->X()[0]*recoDirCosTheta;
    recoVtxY -= minimizer->X()[0]* TMath::Sin(recoDirTheta)* TMath::Sin(recoDirPhi);
    recoVtxX -= minimizer->X()[0]* TMath::Sin(recoDirTheta)* TMath::Cos(recoDirPhi);
//  cout<<"True vtx: ("<<trueVtxX<<" "<<trueVtxY<<" "<<trueVtxZ<<") t" << trueT
//      <<" dir: ("<< TMath::Sin(trueDirTheta)* TMath::Cos(trueDirPhi)<<" "<< TMath::Sin(trueDirTheta)* TMath::Sin(trueDirPhi)<<" "<<
//  TMath::Cos(trueDirTheta)<<")"
//      <<" KE: " << trueKE << endl;
    cout<<"reco vtx: ("<<recoVtxX<<" "<<recoVtxY<<" "<<recoVtxZ<<") t" << recoT
    <<" dir: ("<< TMath::Sin(recoDirTheta)* TMath::Cos(recoDirPhi)<<" "<<
    TMath::Sin(recoDirTheta)* TMath::Sin(recoDirPhi)<<" "<< recoDirCosTheta <<")"
    <<" KE: " << recoKE << endl;
    double recoPar[8]={0,recoVtxZ,recoVtxY,recoVtxX,recoT,recoDirCosTheta,recoDirPhi, recoKE};
    recoLnL = LnLikelihood(recoPar, ipnu, false, false);
//  trueElectronLnL = -func(trueParL, true);
//  cout << "Log Likelihood (electron) True: "<< trueElectronLnL <<"  Reco: " << recoElectronLnL << endl <<endl;
}


