
#ifndef __TkrFitTrack_H
#define __TkrFitTrack_H 1

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "TkrRecon/Track/TkrBase.h"
#include "TkrRecon/TrackFit/TkrFitPlane.h"
#include "TkrRecon/Cluster/TkrCluster.h"

/** 
* @class TkrFitTrack
*
* @brief Contains the data members which specify a TKR fit track, and access methods
*
* Adapted from original version of Bill Atwood
*
* @author Bill Atwood, Tracy Usher, Leon Rochester
*/

class TkrFitTrack: public TkrBase
{    
public:
    /// Constructor with arguments
    /**
    * Instatiate a TkrFitTrack class and initialize
    * @param layer   The Layer Number of the first hit
    * @param tower   The tower number of this first hit
    * @param energy  Initial energy to assign to the track (and use in fit)
    * @param testRay Initial position and direction of candidate track
    */    
    TkrFitTrack(int layer, int tower, double energy, const Ray& testRay);
   ~TkrFitTrack();
    
    /// Utilities 
    void   clear();
    void   writeOut(MsgStream& log) const; 
    void   draw(gui::DisplayRep& v);

    /// Is there a track?
    bool          empty()                  const;

    /// Access to primary quantities
    inline double getChiSquare()           const {return m_chisq;}
    inline double getChiSquareSmooth()     const {return m_chisqSmooth;}
    inline double getQuality()             const {return m_Q;};
    inline double getKalThetaMS()          const {return m_KalThetaMS;}
    inline double getKalEnergy()           const {return m_KalEnergy;}
    inline double getScatter()             const {return m_rmsResid;}

    /// Access errors at track start
    double        getErrorXPosition()      const;
    double        getErrorXSlope()         const;
    double        getErrorYPosition()      const;
    double        getErrorYSlope()         const;

    /// Access to derived information on gaps
    int           getNumGaps()             const;
    int           getNumXGaps()            const {return m_Xgaps;}
    int           getNumYGaps()            const {return m_Ygaps;}
    int           getNumXFirstGaps ()      const {return m_XistGaps;}
    int           getNumYFirstGaps ()      const {return m_YistGaps;}

    /// Access to derived information on kinks
    double        getKink(int iplane)      const;
    double        getKinkNorma(int iplane) const;

    /// Access to plane information
    TkrFitPlane   getFirstPlane()          const;
    TkrFitPlane   getLastPlane()           const;
    int           getNumHits()             const {return m_hits.size();}

    TkrFitPlaneConPtr hitIterConst()  {return m_hits.begin();}
    TkrFitPlaneConPtr hitIterEnd()    {return m_hits.end();}
      
    enum Status {EMPTY, FOUND, CRACK}; 

    /// Members below are protected - this gives access to the
    /// classes which fit the tracks (and inherit from this class)
protected:	
    /// Visualization routines (move to own algorithm?)
    void   drawChiSq(gui::DisplayRep& v, TkrFitHit::TYPE typ);
    void   drawTrack(gui::DisplayRep& v, TkrFitHit::TYPE typ);

    // Track quality information
    double m_chisq;
    double m_chisqSmooth;
    double m_rmsResid;
    double m_Q;

    // Information from the fitter
    double m_KalEnergy;
    double m_KalThetaMS;

    // Can this information be defined in methods?
    int    m_Xgaps;
    int    m_Ygaps;
    int    m_XistGaps;
    int    m_YistGaps;

    // Vector containing hits for each plane traversed
    TkrFitPlaneCol m_hits;
};

//Following typedefs for containing fit track objects
typedef std::vector<TkrFitTrack*>            TkrFitTrackCol;
typedef std::vector<TkrFitTrack*>::iterator  TkrFitTrackColPtr;

#endif
