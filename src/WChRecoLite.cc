/* Andrey Elagin, September 5, 2014
 * LightReco is a light-weight standalone vertex and (in the near future) directionality
 * reconsruction code for 0vbb-decay events in liquid scintillator. The code is based on
 * quadruplet-based vertex-finding method by Michael Smy. Many lines are directly copied
 * from WCSimAnalysis package.
 *
 * See README.txt for instructions on how to use and notes on significant updates.
 */


#include "WChRecoLite.hh"
#include "TObject.h"
#include "WCSimRecoDigit.hh"
#include "WCSimRecoObjectTable.hh"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TError.h"
#include "TMinuit.h"
#include "TRandom.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>

#include <math.h>

//#include "help_func.C"

ClassImp(WChRecoLite)

static WChRecoLite* fgWChRecoLite = 0;

static void vertex_time_lnl(Int_t&, Double_t*, Double_t& f, Double_t* par, Int_t)
{
    Double_t vtx_time = par[0];
    Double_t fom = 0.0;
//  time_fit_itr();
    fgWChRecoLite->TimePropertiesLnL(vtx_time,fom);
    f = -fom; // note: need to maximize this fom
    return;
}


WChRecoLite* WChRecoLite::Instance()
{
    if( !fgWChRecoLite ){
        fgWChRecoLite = new WChRecoLite();
    }
    /*
    if( !fgWChRecoLite ){
      assert(fgWChRecoLite);
    }
    */
    if( fgWChRecoLite ){

    }

    return fgWChRecoLite;
}


WChRecoLite::WChRecoLite()
{

    R_SPHERE=550; //sphere diameter [cm]??????? was 650
    R_LENGTH=1100; //length of detector
    N_REF=1.34; //average index of refraction
    C_VAC=29.9792458; //speed of light in vacuum [cm/ns]
    NSeedsTarget=400; //number of quadruplets
    TSIGMA=0.5; //total time spread (including detector TTS chromatic dispersions)
    fNDigits=0;
    fThisDigit=0;
    fLastEntry=0;
    fCounter=0;
    fMinTime=0;

    fBaseFOM=100.0; //Figure of merit. Borrowed from WCSim: the higher it is the better
    meanTime=0.;
    seedTime=0.;

    //this is for diagnostics, not finished yet
    fFOM = new TFile("fFOM.root","recreate");
    hT = new TH1F("hT","hT",100,-10,10);
    hDT0 = new TH1F("hDT0","hDT0",100,-5,5);
    hDT = new TH1F("hDT","hDT",100,-5,5);

    hX = new TH1F("hX","hX",100,-550,550);
    hY = new TH1F("hY","hY",100,-550,550);
    hZ = new TH1F("hZ","hZ",100,-1100,1100);
    hFOM = new TH1F("hFOM","hFOM",50,100.,200.);

}

WChRecoLite::~WChRecoLite()
{
    this->Reset();

    fFOM->Close();
    delete hT;
    delete hDT0;
    delete hDT;
    delete hX;
    delete hY;
    delete hZ;
    delete hFOM;
}


void WChRecoLite::Reset()
{
    fDigitX.clear();
    fDigitY.clear();
    fDigitZ.clear();
    fDigitT.clear();
    fDigitQ.clear();
    fDigitPE.clear();
    fDigitW.clear();
    fDigitV.clear();
    fDelta.clear();

    vSeedVtxX.clear();
    vSeedVtxY.clear();
    vSeedVtxZ.clear();
    vSeedVtxTime.clear();
    vSeedDigitList.clear();
}


// This function solves system of 4 equations with 4 unknowns to find the seed vertex for each quadruple
int WChRecoLite::FindVertex(Double_t x0, Double_t y0, Double_t z0, Double_t t0, Double_t x1, Double_t y1, Double_t z1, Double_t t1, Double_t x2, Double_t y2, Double_t z2, Double_t t2, Double_t x3, Double_t y3, Double_t z3, Double_t t3, Double_t& vxm, Double_t& vym, Double_t& vzm, Double_t& vtm, Double_t& vxp, Double_t& vyp, Double_t& vzp, Double_t& vtp)
{
    vxm = -99999.9;
    vym = -99999.9;
    vzm = -99999.9;
    vtm = -99999.9;

    vxp = -99999.9;
    vyp = -99999.9;
    vzp = -99999.9;
    vtp = -99999.9;

    // speed of light in water
    // =======================
    Double_t c = C_VAC/N_REF;

    // causality checks
    // ================
    if( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0) >= c*c*(t1-t0)*(t1-t0)
        && (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) >= c*c*(t2-t1)*(t2-t1)
        && (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) + (z3-z2)*(z3-z2) >= c*c*(t3-t2)*(t3-t2)
        && (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0) >= c*c*(t2-t0)*(t2-t0)
        && (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) + (z3-z1)*(z3-z1) >= c*c*(t3-t1)*(t3-t1)
        && (x3-x0)*(x3-x0) + (y3-y0)*(y3-y0) + (z3-z0)*(z3-z0) >= c*c*(t3-t0)*(t3-t0) ){

        // [Note: for causality, require that |x_{i}-x_{j}| >= c*|t_{i}-t_{j}|
        //        for each pair of points]
        Double_t dx1 = x1-x0;  Double_t dy1 = y1-y0;  Double_t dz1 = z1-z0;  Double_t dt1 = c*(t1-t0);
        Double_t dx2 = x2-x0;  Double_t dy2 = y2-y0;  Double_t dz2 = z2-z0;  Double_t dt2 = c*(t2-t0);
        Double_t dx3 = x3-x0;  Double_t dy3 = y3-y0;  Double_t dz3 = z3-z0;  Double_t dt3 = c*(t3-t0);

        Double_t epsilon = 1.0e-7;

        // check that points don't all lie in a plane
        if( !( fabs(dx1)<epsilon && fabs(dx2)<epsilon && fabs(dx3)<epsilon )
            && !( fabs(dy1)<epsilon && fabs(dy2)<epsilon && fabs(dy3)<epsilon )
            && !( fabs(dz1)<epsilon && fabs(dz2)<epsilon && fabs(dz3)<epsilon )
            && !( fabs(dx1)<epsilon && fabs(dy1)<epsilon && fabs(dz1)<epsilon )
            && !( fabs(dx2)<epsilon && fabs(dy2)<epsilon && fabs(dz2)<epsilon )
            && !( fabs(dx3)<epsilon && fabs(dy3)<epsilon && fabs(dz3)<epsilon ) ){

            // [Note: this is a problem for detectors with flat faces!]

            Double_t Mdata[9] = { dx1+epsilon, dy1, dz1,
                                  dx2, dy2+epsilon, dz2,
                                  dx3, dy3, dz3+epsilon };

            Double_t Qdata[3] = { 0.5*( dx1*dx1 + dy1*dy1 + dz1*dz1 - dt1*dt1 ),
                                  0.5*( dx2*dx2 + dy2*dy2 + dz2*dz2 - dt2*dt2 ),
                                  0.5*( dx3*dx3 + dy3*dy3 + dz3*dz3 - dt3*dt3 ) };

            Double_t Tdata[3] = { dt1,
                                  dt2,
                                  dt3 };

            TMatrixD M(3,3,Mdata);
            TMatrixD Q(3,1,Qdata);
            TMatrixD T(3,1,Tdata);

            if( M.Determinant() != 0.0 ){

                TMatrixD A(3,1);
                TMatrixD B(3,1);

                M.Invert();
                A.Mult(M,T);
                B.Mult(M,Q);

                Double_t ax = A(0,0);
                Double_t ay = A(1,0);
                Double_t az = A(2,0);

                Double_t bx = B(0,0);
                Double_t by = B(1,0);
                Double_t bz = B(2,0);

                Double_t ab = ax*bx + ay*by + az*bz;
                Double_t a2 = ax*ax + ay*ay + az*az;
                Double_t b2 = bx*bx + by*by + bz*bz;

                Double_t qa = a2-1.0;
                Double_t qb = 2.0*ab;
                Double_t qc = b2;

                // check for solutions
                if( qb*qb-4.0*qa*qc>0.0 ){

                    // The common vertex is given by a quadratic equation, which has two solutions.
                    // Typically, one solution corresponds to photons travelling forwards in time,
                    // and the other solution corresponds to photons travelling backwards in time.
                    // However, sometimes there appear to be two valid solutions.

                    Double_t ctm = ( -qb - sqrt(qb*qb-4.0*qa*qc) ) / ( 2.0*qa );
                    Double_t ctp = ( -qb + sqrt(qb*qb-4.0*qa*qc) ) / ( 2.0*qa );

                    Double_t tm = t0 + ctm/c;
                    Double_t xm = x0 + ctm*ax + bx;
                    Double_t ym = y0 + ctm*ay + by;
                    Double_t zm = z0 + ctm*az + bz;
                    Bool_t foundVertexM = 0;

                    if( tm<t0 && tm<t1
                        && tm<t2 && tm<t3 ){
                        vxm = xm;
                        vym = ym;
                        vzm = zm;
                        vtm = tm;
                        foundVertexM = 1;
                    }

                    Double_t tp = t0 + ctp/c;
                    Double_t xp = x0 + ctp*ax + bx;
                    Double_t yp = y0 + ctp*ay + by;
                    Double_t zp = z0 + ctp*az + bz;
                    Bool_t foundVertexP = 0;

                    if( tp<t0 && tp<t1
                        && tp<t2 && tp<t3 ){
                        vxp = xp;
                        vyp = yp;
                        vzp = zp;
                        vtp = tp;
                        foundVertexP = 1;
                    }

                }
                else
                {
                    // std::cout << "qb*qb-4.0*qa*qc<0.0." << std::endl;
                }
            }
            else
            {
                // std::cout << "M.Determinant() == 0.0. " << std::endl;
            }

        }
        else
        {
            // std::cout << "The 4 points lie on a plane. " << std::endl;
        }
    }
    else
    {
        // std::cout << "Causality check for this 4-group of digits has not been passed. " << std::endl;
    }


    return 0;
}



// Randomly select photon hit (x,y,z,t) for a quadruple
Int_t WChRecoLite::ChooseNextDigit(Double_t& xpos, Double_t& ypos, Double_t& zpos, Double_t& time)
{
    xpos=0.; ypos=0.; zpos=0.; time=0.;
    // ROOT random number generator
    Double_t r = RND.Rndm();

    //  std::cout<<"I'm inside ChooseNextDitig"<<std::endl;
    // pseudo-random number generator
    Int_t numEntries = vSeedDigitList.size();
/*  cout<<"fCounter = "<<fCounter<<"   fLastEntry = "<<fLastEntry<<endl;
  fCounter++;
  if( fCounter>=fNDigits ) fCounter = 0;
  fThisDigit = vSeedDigitList.at(fLastEntry);

  Double_t t0 = 0.5 + fDigitT[fCounter] - fMinTime;
  Double_t q0 = 0.5 + fDigitQ[fCounter];

  Double_t t1 = 0.5 + fDigitT[fThisDigit] - fMinTime;
  Double_t q1 = 0.5 + fDigitQ[fThisDigit];

  Double_t tq = 100.0*(t0*q0+t1*q1);
  Double_t r = tq - TMath::Floor(tq);
*/
//  r = gRandom->Uniform(); // Christoph Aberle, August 14, 2013: use of a proper RN generator since I saw that quadruplets were duplicated with the pseudo-random number generator used in the lines above

    fLastEntry = (Int_t)(r*numEntries);

    //  std::cout<<"fLastEntry = "<<fLastEntry<<"   r = "<<r<<"   numEntries = "<<numEntries<<std::endl;

    // return the new digit
    fThisDigit = vSeedDigitList.at(fLastEntry);
    //  std::cout<<"fThisDigit = "<<fThisDigit<<std::endl;
    xpos = fDigitX[fThisDigit];
    ypos = fDigitY[fThisDigit];
    zpos = fDigitZ[fThisDigit];
    time = fDigitT[fThisDigit];
    // std::cout<<"xpos = "<<xpos<<"   ypos = "<<ypos<<"   zpos = "<<zpos<<"   time = "<<time<<std::endl;
    return fThisDigit;
}


// Make a quadruple
Int_t WChRecoLite::ChooseNextQuadruple(Double_t& x0, Double_t& y0, Double_t& z0, Double_t& t0, Double_t& x1, Double_t& y1, Double_t& z1, Double_t& t1, Double_t& x2, Double_t& y2, Double_t& z2, Double_t& t2, Double_t& x3, Double_t& y3, Double_t& z3, Double_t& t3)
{
    int code=0; // 0 -if OK, 1 -if failed to chose 4 different digits
    int digit0=0;
    int digit1=-1;
    int digit2=-2;
    int digit3=-3;
    int counter1=0;
    int counter2=0;
    int counter3=0;
    digit0 = ChooseNextDigit(x0,y0,z0,t0);
//  digit1 = ChooseNextDigit(x1,y1,z1,t1);
//  digit2 = ChooseNextDigit(x2,y2,z2,t2);
//  digit3 = ChooseNextDigit(x3,y3,z3,t3);

//check that selected digits for the quadruple are not identical
//if they are after 100 attempts then be it, the quadruple will not
//be used later
    while( (digit1<0 || digit1==digit0) && counter1<100)
    {
        digit1 = ChooseNextDigit(x1,y1,z1,t1);
        counter1++;
    }

    while( (digit2<0 || digit2==digit0 || digit2==digit1) && counter2<100)
    {
        digit2 = ChooseNextDigit(x2,y2,z2,t2);
        counter2++;
    }

    while( (digit3<0 || digit3==digit0 || digit3==digit1 || digit3==digit2) && counter3<100)
    {
        digit3 = ChooseNextDigit(x3,y3,z3,t3);
        counter3++;
    }

    if(counter1>=100 || counter2>=100 || counter3>=100) code=1;

    return code;
}


// Calculate NSeedsTarget verticies
Int_t WChRecoLite::CalcVertexSeeds()
{
    // reset list of seeds
    // ===================
    vSeedVtxX.clear();
    vSeedVtxY.clear();
    vSeedVtxZ.clear();
    vSeedVtxTime.clear();

    Double_t x0 = 0.0;
    Double_t y0 = 0.0;
    Double_t z0 = 0.0;
    Double_t t0 = 0.0;

    Double_t x1 = 0.0;
    Double_t y1 = 0.0;
    Double_t z1 = 0.0;
    Double_t t1 = 0.0;

    Double_t x2 = 0.0;
    Double_t y2 = 0.0;
    Double_t z2 = 0.0;
    Double_t t2 = 0.0;

    Double_t x3 = 0.0;
    Double_t y3 = 0.0;
    Double_t z3 = 0.0;
    Double_t t3 = 0.0;

    double fVtxX1, fVtxY1, fVtxZ1, fVtxTime1;
    double fVtxX2, fVtxY2, fVtxZ2, fVtxTime2;

    int counter=0;
    std::cout<<"I'm inside CalcVertexSeeds"<<std::endl;

    while( vSeedVtxX.size()<NSeedsTarget && counter<100*NSeedsTarget )
    {
        //    std::cout<<"counter = "<<std::endl;
        ChooseNextQuadruple(x0,y0,z0,t0,
                            x1,y1,z1,t1,
                            x2,y2,z2,t2,
                            x3,y3,z3,t3);
        /*
        std::cout<<"counter = "<<counter<<std::endl;
        std::cout << "   digit0: (x,y,z,t)=(" << x0 << "," << y0 << "," << z0 << "," << t0 << ") " << std::endl;
        std::cout << "   digit1: (x,y,z,t)=(" << x1 << "," << y1 << "," << z1 << "," << t1 << ") " << std::endl;
        std::cout << "   digit2: (x,y,z,t)=(" << x2 << "," << y2 << "," << z2 << "," << t2 << ") " << std::endl;
        std::cout << "   digit3: (x,y,z,t)=(" << x3 << "," << y3 << "," << z3 << "," << t3 << ") " << std::endl;
        */

        FindVertex(x0,y0,z0,t0,
                   x1,y1,z1,t1,
                   x2,y2,z2,t2,
                   x3,y3,z3,t3,
                   fVtxX1,fVtxY1,fVtxZ1,fVtxTime1,
                   fVtxX2,fVtxY2,fVtxZ2,fVtxTime2);
        /*
         std::cout << "   result1: (x,y,z,t)=(" << fVtxX1 << "," << fVtxY1 << "," << fVtxZ1 << "," << fVtxTime1 << ") " << std::endl
                   << "   result2: (x,y,z,t)=(" << fVtxX2 << "," << fVtxY2 << "," << fVtxZ2 << "," << fVtxTime2 << ") " << std::endl;
         std::cout<<"-------------"<<std::endl;
        */

        if(fVtxX1==-99999.9 && fVtxX2==-99999.9) //no solutions try anouther quadruple
        {
            counter++;
            continue;
        }


        bool inside_det;
        if(sqrt(fVtxX1*fVtxX1+fVtxY1*fVtxY1)<R_SPHERE && fVtxZ1<R_LENGTH) inside_det=true; else inside_det=false;
        // std::cout<<"Solution1: inside_det = "<<inside_det<<"  X = "<<fVtxX1<<"  Y = "<<fVtxY1<<"  Z = "<<fVtxZ1<<"  T = "<<fVtxTime1<<std::endl;
//    if(!inside_det)
//    {
//      std::cout<<"Solution outside detector"<<std::endl;
//      continue;
//    }

        // add first digit
        if( inside_det ){
            vSeedVtxX.push_back(fVtxX1);
            vSeedVtxY.push_back(fVtxY1);
            vSeedVtxZ.push_back(fVtxZ1);
            vSeedVtxTime.push_back(fVtxTime1);
            // std::cout << "New vertex seed 1: x= " << fVtxX1 << ", " << fVtxY1 << ", " << fVtxZ1 << ", " <<fVtxTime1 << std::endl;
        }


        if(sqrt(fVtxX2*fVtxX2+fVtxY2*fVtxY2)<R_SPHERE && fVtxZ2<R_LENGTH) inside_det=true; else inside_det=false;
        // std::cout<<"Solution2: inside_det = "<<inside_det<<"  X = "<<fVtxX2<<"  Y = "<<fVtxY2<<"  Z = "<<fVtxZ2<<"  T = "<<fVtxTime2<<std::endl;
//    if(!inside_det)
//    {
//      std::cout<<"Solution outside detector"<<std::endl;
//      continue;
//    }

        // add second digit
        if( inside_det ){
            vSeedVtxX.push_back(fVtxX2);
            vSeedVtxY.push_back(fVtxY2);
            vSeedVtxZ.push_back(fVtxZ2);
            vSeedVtxTime.push_back(fVtxTime2);
            // std::cout << "New vertex seed 2: x= " << fVtxX2 << ", " << fVtxY2 << ", " << fVtxZ2 << ", " <<fVtxTime2 << std::endl;
        }

        counter++;

    }

    std::cout << "The number of calculated seeds in vSeedVtxX, vSeedVtxY, vSeedVtxZ, vSeedVtxTime for this event is = " << vSeedVtxX.size() << std::endl;


    return 0;
}


double WChRecoLite::Amotry(double p[5][4], double y[], int ndim, double psum[], const int ihi, const double fac)
{
    int j;
    double fac1,fac2,ytry;
    double ptry[ndim];
    fac1=(1.0-fac)/ndim;
    fac2=fac1-fac;
    for (j=0;j<ndim;j++)
        ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    ytry=GetSumVtx(ptry);
    if (ytry < y[ihi]) {
        y[ihi]=ytry;
        for (j=0;j<ndim;j++) {
            psum[j] += ptry[j]-p[ihi][j];
            p[ihi][j]=ptry[j];
        }
    }
//    std::cout << ptry[0] << ", " << ptry[1] << ", " << ptry[2] << ", " << ptry[3] << " : " << ytry << std::endl;
    return ytry;
}

void get_psum(double p[5][4], double psum[], int ndim)
{
    int i,j;
    double sum;
    int mpts=ndim+1;
    for (j=0;j<ndim;j++) {
        for (sum=0.0,i=0;i<mpts;i++)
            sum += p[i][j];
        psum[j]=sum;
    }
}

void SWAP(double &a, double &b)
{
    double dum=a; a=b; b=dum;
}

void WChRecoLite::Amoeba(double p[5][4], double y[], int ndim, const double ftol, int &nfunk)
{
    const int NMAX=5000;
    const double TINY=1.0e-10;
    int i,ihi,ilo,inhi,j;
    double rtol,ysave,ytry;
    int mpts=ndim+1;
    double psum[ndim];
    nfunk=0;
    get_psum(p,psum,ndim);
    for (;;) {
        ilo=0;
        ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
        for (i=0;i<mpts;i++) {
            if (y[i] <= y[ilo]) ilo=i;
            if (y[i] > y[ihi]) {
                inhi=ihi;
                ihi=i;
            } else if (y[i] > y[inhi] && i != ihi) inhi=i;
        }
        rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
        if (rtol < ftol) {
            SWAP(y[0],y[ilo]);
            for (i=0;i<ndim;i++) SWAP(p[0][i],p[ilo][i]);
            std::cout << p[0][0] << ", " << p[0][1] << ", " << p[0][2] << ", " << p[0][3] << " : " << y[0] << std::endl;
            break;
        }
        if (nfunk >= NMAX){
            std::cout << "NMAX exceeded" << std::endl;
            break;
        }
        nfunk += 2;
        ytry=Amotry(p,y,ndim,psum,ihi,-1.0);
        if (ytry <= y[ilo])
            ytry=Amotry(p,y,ndim,psum,ihi,2.0);
        else if (ytry >= y[inhi]) {
            ysave=y[ihi];
            ytry=Amotry(p,y,ndim,psum,ihi,0.5);
            if (ytry >= ysave) {
                for (i=0;i<mpts;i++) {
                    if (i != ilo) {
                        for (j=0;j<ndim;j++)
                            p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
                        y[i]=GetSumVtx(psum);
                    }
                }
                nfunk += ndim;
                get_psum(p,psum,ndim);
            }
        } else --nfunk;
    }
}

double fomtotal;
double weightedRecoX;
double weightedRecoY;
double weightedRecoZ;
std::vector<double> fomarray;

// Calculates FOM by looking at time residuals
// make sure fDelta has been filled before calling this function
Int_t WChRecoLite::TimePropertiesLnL(double & vtx_time, double & fom)
{
    double A = 1.0 / ( 2.0*TSIGMA*TMath::Sqrt(0.5*TMath::Pi()) );
    double Preal=0.;
    double P=0.;
    fom=0.;
    double chi2=0.;
    double ndof=0.0;

    // std::cout<<"fDelta.size() = "<<fDelta.size()<<std::endl;
    for( Int_t idigit=0; idigit<fDigitX.size(); idigit++ )
    {
        double delta = fDelta[idigit] - vtx_time;
        // std::cout << delta << std::endl;
        Preal = A*exp(-(delta*delta)/(2.0*TSIGMA*TSIGMA));
        // P = (1.0-Pnoise)*Preal + Pnoise;
        // chi2 += -2.0*log(P);
        // std::cout<<"A = "<<A<<"   delta = "<<delta<<"   TSIGMA = "<<TSIGMA<<"   Preal = "<<Preal<<std::endl;
        // std::cout<<"delta = "<<delta<<"   Preal = "<<Preal<<"   fDelta[idigit] = "<<fDelta[idigit]<<"   vtx_time = "<<vtx_time<<std::endl;

        // chi2 += -2.0*log(Preal);
        // chi2 += 2.0*(-Preal);
        chi2 += Preal;
        ndof += 1.0;
    }
    // std::cout<<"chi2 = "<<chi2<<"   ndof = "<<ndof<<std::endl;
    // if( ndof>0.0 ){
    // fom = fBaseFOM - 1.0*chi2/ndof;
    // fom = fBaseFOM - 100.0*chi2/ndof;
    fom = chi2;
    fomtotal +=fom;
    // std::cout << "fom=" << fom << std::endl;
//    fom = chi2/ndof;
    // }


    return 0;
}


//Minuit manipulation are decribed here
void WChRecoLite::FitPointTimePropertiesLnL(Double_t& fit_time, Double_t& fom)
{
    Int_t err = 0;
    Int_t flag = 0;

//  Double_t seedTime = meanTime;

//  Double_t fitTime = 0.0;
    Double_t fitTimeErr = 0.0;

    TMinuit* fMinuitTimeFit = new TMinuit();
    fMinuitTimeFit->SetMaxIterations(5000);

    Double_t* arglist = new Double_t[10];
    arglist[0]=1;  // 1: standard minimization
    // 2: try to improve minimum

    fMinuitTimeFit->mncler();
    fMinuitTimeFit->SetFCN(vertex_time_lnl);
    fMinuitTimeFit->mnexcm("SET STR",arglist,1,err);
    fMinuitTimeFit->mnparm(0,"vtx_time",seedTime,0.1,-100.0,100.0,err);  //Negative times have to be possible since the true time is at zero
    flag = fMinuitTimeFit->Migrad();
    fMinuitTimeFit->GetParameter(0,fit_time,fitTimeErr);

    //  std::cout <<"fitTime = " << fit_time << ", fitTimeErr = " << fitTimeErr << std::endl;

    delete [] arglist;
    delete fMinuitTimeFit;

    return;
}

Int_t WChRecoLite::SelectBestSeed(int evt_num)
{
    std::ostringstream oss;
    oss<<evt_num;
    std::string iter_evt = oss.str();

    TH1F* ht = (TH1F*) hT->Clone(("hT_"+iter_evt).c_str());
    TH1F* hdt0 = (TH1F*) hDT0->Clone(("hDT0_"+iter_evt).c_str());
    TH1F* hdt = (TH1F*) hDT->Clone(("hDT_"+iter_evt).c_str());

    TH1F* hx = (TH1F*) hX->Clone(("hX_"+iter_evt).c_str());
    TH1F* hy = (TH1F*) hY->Clone(("hY_"+iter_evt).c_str());
    TH1F* hz = (TH1F*) hZ->Clone(("hZ_"+iter_evt).c_str());
    TH1F* hfom = (TH1F*) hFOM->Clone(("hFOM_"+iter_evt).c_str());

    Int_t bestSeed = -1;
    Double_t bestFOM = -1.0;

    fomtotal = 0.;
    weightedRecoX = 0.;
    weightedRecoY = 0.;
    weightedRecoZ = 0.;
    fomarray.clear();

    for(int i=0;i!=vSeedVtxX.size();++i)
    {
        // vSeedVtxX.clear();
        // std::cout << "SEEDvtx0=" << vSeedVtxX[i] << std::endl;
        hx->Fill(vSeedVtxX[i]);
        hy->Fill(vSeedVtxY[i]);
        hz->Fill(vSeedVtxZ[i]);
        ht->Fill(vSeedVtxTime[i]);

        // loop over digits
        // ================
        double Swx=0.;
        double Sw=0.;
//    double meanTime=0.;
        fDelta.clear();
        for( Int_t idigit=0; idigit<fDigitX.size(); idigit++ )
        {
            Double_t dx = fDigitX[idigit]-vSeedVtxX[i];
            Double_t dy = fDigitY[idigit]-vSeedVtxY[i];
            Double_t dz = fDigitZ[idigit]-vSeedVtxZ[i];
            Double_t ds = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
            //need to check what is the proper variable out of the two listed below
            Double_t time0 = fDigitT[idigit] - 0; //this is what was done in WCSim for the JINST paper
            Double_t time1 = fDigitT[idigit] - vSeedVtxTime[i];

            double fPointResidual0 = time0 - ds/(C_VAC/N_REF);
            double fPointResidual = time1 - ds/(C_VAC/N_REF);
// TEMP test: use true velocity:
            // double fPointResidual0 = time0 - ds/fDigitV[idigit];
            // double fPointResidual = time - ds/fDigitV[idigit];
//      cout<<"Lambda = "<<fDigitW[idigit]<<"   index_ref = "<<C_VAC/fDigitV[idigit]<<endl;
            hdt0->Fill(fPointResidual0);
            hdt->Fill(fPointResidual);

            fDelta.push_back(fPointResidual0);
            double weight = 1.0/(TSIGMA*TSIGMA);
            Swx += time1*weight; // here is some room for upgrade id TSIGMA is not always the same
            Sw += weight;

        }

        meanTime=Swx/Sw;
        seedTime=vSeedVtxTime[i];

        double vtx_time=0.;
        double fom=0.;
        //    FitPointTimePropertiesLnL(vtx_time,fom); //do time fit
        TimePropertiesLnL(vtx_time,fom); //calculate fom for fitted time value (possibly this can be done during the fit)
        hfom->Fill(fom);
        // std::cout << "vSeedVtxX[i]:" << vSeedVtxX[i] << "  " << fom << "  " << fomtotal << std::endl;

        fomarray.push_back(fom);


        // std::cout << "--------------------:" << fomarray[i] << std::endl;

        if( fom>bestFOM ){
            bestSeed = i;
            bestFOM = fom;
        }

        // std::cout<<"bestseed" <<bestSeed<<"  fom = "<<fom<<"   vtx_time = "<<vtx_time<<std::endl;
    }
    std::cout<<"fomtotal: " << fomtotal << std::endl;

    for(int iweight=0;iweight!=vSeedVtxX.size();++iweight){
        double seedweight= fomarray[iweight]/fomtotal;
        weightedRecoX += vSeedVtxX[iweight]*seedweight;
        weightedRecoY += vSeedVtxY[iweight]*seedweight;
        weightedRecoZ += vSeedVtxZ[iweight]*seedweight;
        // std::cout <<  seedweight << "  " << weightedRecoX << "  " << vSeedVtxX[iweight] << std::endl;

    }
    std::cout << "--------------------:" << weightedRecoX << "  " << weightedRecoY << "  " << weightedRecoZ << std::endl;

    fFOM->cd();
    ht->Write();
    hdt0->Write();
    hdt->Write();

    hx->Write();
    hy->Write();
    hz->Write();

    hfom->Write();
    // hx->Fit("gaus");
    // hy->Fit("gaus");
    // hz->Fit("gaus");

    return bestSeed;
}

double WChRecoLite::WeightedSumSeed(int xyz)//do select best seed first
{
    if(xyz==0) return weightedRecoX;
    if(xyz==1) return weightedRecoY;
    if(xyz==2) return weightedRecoZ;

    // else{
    //   continue;
    // }
}

double WChRecoLite::MeanVertexSeed(int dir)
{
    double meanVertex = 0;
    std::vector<double> vSeedVtx;
    if(dir == 0) vSeedVtx=vSeedVtxX;
    else if(dir == 1) vSeedVtx=vSeedVtxY;
    else if(dir == 2) vSeedVtx=vSeedVtxZ;
    else if(dir == 3) vSeedVtx=vSeedVtxTime;
    else return 0;
    int count = vSeedVtx.size();
    for(int i =0; i< count; i++){
        meanVertex += vSeedVtx[i];
    }
    meanVertex /= count;
}

void WChRecoLite::SetDigits(std::vector<double> iDigitX,std::vector<double> iDigitY,std::vector<double> iDigitZ,std::vector<double> iDigitT,std::vector<double> iDigitQ,std::vector<double> iDigitPE,std::vector<double> iDigitW,std::vector<double> iDigitV,std::vector<double> iDelta,std::vector<int> iSeedDigitList)
{
    this->Reset();

    fDigitX=iDigitX;
    fDigitY=iDigitY;
    fDigitZ=iDigitZ;
    fDigitT=iDigitT;
    fDigitQ=iDigitQ;
    fDigitPE=iDigitPE;
    fDigitW=iDigitW;
    fDigitV=iDigitV;
    fDelta=iDelta;
    vSeedDigitList=iSeedDigitList;
}


std::vector<double> WChRecoLite::GetSeedVtx(int wvert){

    if(wvert==0) return vSeedVtxX;
    if(wvert==1) return vSeedVtxY;
    if(wvert==2) return vSeedVtxZ;
    if(wvert==3) return vSeedVtxTime;
}

void WChRecoLite::GetP(double p[5][4]){

    //double **p = new double[5][4];
    // double sumvtx[5]={0.};
    double sphere_rad = 30.;
    int j, k;
    // double ddx[vSeedVtxX.size()], ddy[vSeedVtxX.size()], ddz[vSeedVtxX.size()], ddt[vSeedVtxX.size()];
    // std::cout
    for(j=0;j<4;j++) {
        p[0][j] = MeanVertexSeed(j);
    }
    for(j=0;j<4;j++){
        for (k=0;k<4;k++){
            // p[k][j] = mRec->MeanVertexSeed(j);//mean xyzt
            p[j+1][k] = p[0][k];//mean xyzt
        }
        p[j+1][j] = p[j+1][j] + sphere_rad;//mean xyzt+300
    }
}

double WChRecoLite::GetSumVtx(double vertex[4]){

    // double p[5][4];
    double sumvtx=0;
    double sphere_rad = 30.;
    int j, k;
    double ddx[vSeedVtxX.size()], ddy[vSeedVtxY.size()], ddz[vSeedVtxZ.size()], ddt[vSeedVtxTime.size()];

    // for(j=0;j<4;j++){
    //   for (k=0;k<5;k++){
    //     // p[k][j] = mRec->MeanVertexSeed(j);//mean xyzt
    //     p[k][j] = MeanVertexSeed(j);//mean xyzt
    //   }
    //   p[j+1][j] = MeanVertexSeed(j) + sphere_rad;//mean xyzt+300
    // }

    for(int i=0;i!=vSeedVtxX.size();++i){
        ddx[i]=(fabs(vertex[0]-vSeedVtxX[i]))/sphere_rad;
        ddy[i]=(fabs(vertex[1]-vSeedVtxY[i]))/sphere_rad;
        ddz[i]=(fabs(vertex[2]-vSeedVtxZ[i]))/sphere_rad;
        ddt[i]=(fabs(vertex[3]-vSeedVtxTime[i]))/(sphere_rad/10);
        // std::cout << ddx[i] <<" " <<ddy[i] << "  " << ddz[i] << "  " << ddt[i] << std::endl;
        if(ddx[i]*ddx[i]+ddy[i]*ddy[i]+ddz[i]*ddz[i]+ddt[i]*ddt[i]<1.){
            // std::cout << "inside sphere!"<< std::endl;
            sumvtx=sumvtx-1.;
        }
    }

    return sumvtx;
}