//
//********************************************************************
//
// Class definition for candidate pattern reconstruction track
// Contains the information necessary to run the track fit
//
// Note that TkrBase contains the base information needed to seed
// the Kalman Filter fit routines, everything else is extra
//
// Under development!
// 1/15/02 Tracy Usher
//
//
#ifndef __TkrPatCand_H
#define __TkrPatCand_H

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "TkrRecon/Track/TkrBase.h"
#include "TkrRecon/PatRec/TkrPatCandHit.h"

class TkrPatCand: public TkrBase
{    
public:
    
    TkrPatCand(int layer, int tower, const Ray& testRay);
   ~TkrPatCand();

    //Provide a method to writeout the contents of the class
    void writeOut(MsgStream& log) const; 

    //Method to store hits into vector
    void addCandHit(TkrCluster*   pCluster);
    void addCandHit(TkrPatCandHit candHit);

    //Access to some information regarding hits (if they exist)
    int              lastLayer();

    //Provide access to the vector of hit candidates
    int              numPatCandHits()    {return m_hits.size();}
    TkrPatCandHit*   getCandHit(int idx) {return &m_hits[idx];}

    //Provide access to a vector iterator (do we want this?)
    CandHitVectorPtr getCandHitPtr()     {return m_hits.begin();}
           
private:
    //For sorting the hits
    friend bool      sortHits(const TkrPatCandHit o1, const TkrPatCandHit o2) {return o1 < o2;}

    CandHitVector m_hits;
};

//Following typedefs for containing track candidate objects
typedef std::vector<TkrPatCand*>            CandTrkVector;
typedef std::vector<TkrPatCand*>::iterator  CandTrkVectorPtr;


#endif
