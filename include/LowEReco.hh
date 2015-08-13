#ifndef LOWERECO_HH
#define LOWERECO_HH

#include "TObject.h"

class LowEReco : public TObject {

public:

    LowEReco();
    ~LowEReco();

    double ReconstructEnergy(int nphot, const Double_t *phot_xEnd, const Double_t *phot_yEnd, const Double_t *phot_zEnd);
    void DoLowEReco(int evt, int nhits, double * hitx, double * hity, double * hitz, double * hitt,
                              double allrecoang,
                              double & recoVtxX, double & recoVtxY, double & recoVtxZ, double & recoVtxTime,
                              double &reconstructedDirX, double &reconstructedDirY, double &reconstructedDirZ);

        private:

ClassDef(LowEReco,0)

};

#endif
