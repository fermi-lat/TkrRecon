/**
* @class TkrNeuralNetPR
*
* @brief Class definition for the Neural Net PatRec Transient Data Object.
*
* last modified 3/02
*
* @authors b. allgood and w. atwood
*
* $Header$
*/

#ifndef __TKRNEURALNETPR_H
#define __TKRNEURALNETPR_H

#include "src/PatRec/NeuralNet/TkrNeuralNet.h"
#include "TkrRecon/PatRec/TkrPatRecon.h"

class TkrNeuralNetPR : public TkrPatRecon
{
public:
    TkrNeuralNetPR(ITkrGeometrySvc* pTkrGeo) {pGeometry = pTkrGeo;}
   ~TkrNeuralNetPR() {}

    TkrCandidates* doPatRecon(TkrClusters* pTkrClus, double energy, 
                              Point position)
	{return new TkrNeuralNet(pGeometry, pTkrClus, energy, position);}

private:
    ITkrGeometrySvc* pGeometry;
};

#endif  // __TKRNEURALNETPR_H