//Header class for dark pulse generation helper classes
//Creating the dark pulses is done in the photon masks
//Methods declared here are from darkPulseHelpers.cc

#ifndef DARK_PULSE_HELPERS
#define DARK_PULSE_HELPERS

//Used to decide when in the time window to place a hit
double randomTime(double start, double stop);
//Given the cumulative probs (see below) calculates using random numbers how many dark noise hits
int numPulses(double cumulativeProbs[], int maxNumHits);
//Used to calculate the likleyhood of diffrent each possible number of dark noise hits
double * cumulativeProbs(double darkPulseRate, double eventWindowLength, int maxNumHits);

#endif

