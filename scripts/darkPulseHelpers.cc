// Contains all helper functions for dark pulse generation
// Creating the dark pulses is done in the photon masks
#include "TMath.h" //has a power class and an exponent class and a factorial class
#include "darkPulseHelpers.h"
#include "TRandom3.h"
#include <iostream>



// probability method; to be called once when window time and dark pulse rates are known
// Returns the cumulative probablility array
// The cumulative probablility of a number of hits is the probablility of having
// bettween 0 and that number of hits inclusive.
// We then generate a random number bettween 0 and 1, and identify the number of hits
// corrisponing to this random number by which cumulative probablilities it falls bettween.
// As there could theoreticaly be an infinite number of hits we chose a maximum number
// of hits beond which we will not calculate cumulative probablilities
// darkPulseRate in Hz, eventWindowLength is receved in nano seconds then converted to seconds
double * cumulativeProbs(double darkPulseRate, double eventWindowLength, int maxNumHits){
	double * cumulativeProbs = new double[maxNumHits + 1];
	double lengthInSeconds = eventWindowLength*TMath::Power(10.,-9); 
	double mu = darkPulseRate*lengthInSeconds; //the poission distribution variable
	//Its easyest to calculate the probablity of 0 hits separatly
	cumulativeProbs[0] = TMath::Exp(-mu);
	for(int i = 1; i <= maxNumHits; i++){
		//We use a poisson distribution v
		double probIHits = TMath::Exp(-mu)*TMath::Power(mu,i)/TMath::Factorial(i);
		double cumulative = cumulativeProbs[i-1] + probIHits;
		//Sometimes the really small numbers generate infinities, no idea why.
		//Just take them out. The results should tend to 1 anyway.
		if(cumulative>1) cumulative = 1; 
		cumulativeProbs[i] = cumulative;
	}
	return cumulativeProbs;

}


//Decides using a random number generator how many dark pulses the pmt sees
// this takes forever for some reason related to TRandom3
int numPulses(double cumulativeProbs[], int maxNumHits)
{
	static TRandom3 numberp(42484);
	int numHits = 0;
	//Get a random number bettween 0 and 1
	double random = numberp.Rndm();

	while(random > cumulativeProbs[numHits]){
 		numHits++;
		if (numHits > maxNumHits){
			return maxNumHits + 1; 
			std::cout << "Error: maxNumHits has been reached. (If this happens frequently then raise maxNumHits in the photon mask) <<<<<<"<< std::endl;
		}
	}
	return numHits;
	
}

//Picks a random time in the time window 
// this takes forever for some reason related to TRandom3
double randomTime(double start, double stop)
{
	static TRandom3 numberp(42484);
	//Get a random number bettween 0 and 1
	double random = numberp.Rndm();

	double length = stop - start;
	return start + length*random;
	
}

