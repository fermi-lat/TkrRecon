#ifndef __TkrPatCand_H
#define __TkrPatCand_H
/** 
* @class TkrPatCand
*
* @brief Contains the candidate track information needed for the track fit
*        Note that everything needed to seed the fitter (in find hits mode)
*        can be obtained by the methods defined in TkrRecInfo.
*
* @author The Tracking Software Group
*/
#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "TkrRecon/Track/TkrRecInfo.h"
#include "TkrRecon/PatRec/TkrPatCandHit.h"

class TkrPatCand: public TkrRecInfo
{    
public:
    
    TkrPatCand(int layer, int tower, double energy, double quality, const Ray& testRay);
   ~TkrPatCand();

    /// Define the TkrRecInfo routines
    double        getQuality()                      const {return m_quality;};
    double        energy(TrackEnd end = Start)      const {return m_energy;}
    int           layer(TrackEnd end = Start)       const {return m_firstLayer;}
    int           tower(TrackEnd end = Start)       const {return m_itower;}
    Point         position(TrackEnd end = Start)    const {return m_position;}
    Vector        direction(TrackEnd end = Start)   const {return m_direction;}
    Ray           ray(TrackEnd end = Start)         const {return Ray(position(),direction());}
    TkrFitPar     TrackPar(TrackEnd end = Start)    const;
    double        TrackParZ(TrackEnd end = Start)   const {return m_position.z();}
    TkrFitMatrix  TrackCov(TrackEnd end = Start)    const; 
    bool          empty()                           const {return m_firstLayer >= 0;}

    //Provide a method to writeout the contents of the class
    void writeOut(MsgStream& log) const; 

    //Method to store hits into vector
    void addCandHit(TkrCluster*   pCluster);
    void addCandHit(TkrPatCandHit candHit);

    //Methods to set member variables if need be
    void setEnergy(double energy) {m_energy = energy;}

    //Access to some information regarding hits (if they exist)
    int              lastLayer();

    //Provide access to the vector of hit candidates
    int              numPatCandHits()    {return m_hits.size();}
    TkrPatCandHit*   getCandHit(int idx) {return &m_hits[idx];}
    TkrPatCandHit*   getFoLPlane(TrackEnd end = Start);

    //Provide access to a vector iterator (do we want this?)
    CandHitVectorPtr getCandHitPtr()     {return m_hits.begin();}
           
private:
    //For sorting the hits
    friend bool      sortHits(const TkrPatCandHit o1, const TkrPatCandHit o2) {return o1 < o2;}

    CandHitVector m_hits;

    //Information mandated by TkrRecInfo
    Point         m_position;
    Vector        m_direction;
    double        m_energy;
    double        m_quality;
    int           m_firstLayer;
    int           m_itower;    
};

//Following typedefs for containing track candidate objects
typedef std::vector<TkrPatCand*>            CandTrkVector;
typedef std::vector<TkrPatCand*>::iterator  CandTrkVectorPtr;


#endif
