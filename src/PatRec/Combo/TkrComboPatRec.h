/**
  * @class TkrComboPatRec
  *
  * @brief A Combinatoric Pattern recognition for GLAST
  *
  * 01-Dec-2001
  * This is meant to be close to what has been used historically 
  * to find tracks in GLAST.  
  *
  * It uses two basic methods: 
  * 1) Events seeded with an energy centroid in the Cal and an energy
  * 2) Events without cal information
  * If the incoming Cal Point is null its assumed that there is no cal 
  * information
  *
  * @author Bill Atwood, SCIPP/UCSC
  *
  * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/Combo/TkrComboPatRec.h,v 1.16 2002/08/30 18:35:40 atwood Exp $
*/

#ifndef __TKRCOMBOPATREC_H
#define __TKRCOMBOPATREC_H

#include "Event/Recon/TkrRecon/TkrPatCandCol.h"
#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "src/Track/TkrControl.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "src/TrackFit/KalFitTrack/KalFitTrack.h"
#include "geometry/Ray.h"

#include <vector>

using namespace Event;

class TkrComboPatRec : public TkrPatCandCol
{
public:
    TkrComboPatRec(ITkrGeometrySvc* pTkrGeo, TkrClusterCol* pClusters, double CalEnergy, Point CalPosition);
        ~TkrComboPatRec() {};

private:
    class Candidate
    {
    public:
        Candidate(TkrClusterCol* clusters,
                  ITkrGeometrySvc* geometry,
                  int layer, int twr, double e, 
                  Point x, Vector t, float d, float s, int g, int top);
        ~Candidate();

        /// Access
        void setDeflection(float def) {m_deflection = def;}
        void setSigma(float sig)      {m_sigma      = sig;}
        void setQuality(float Q)      {m_qual       = Q;}
        void setGap(int gs)           {m_gap        = gs;}
        void setConEnergy(double e)   {m_ConEnergy  = e;}
        int  adjustType(int incr);
        float  deflection()    const {return m_deflection;}
        float  sigma()         const {return m_sigma;}
        float  quality()       const {return m_qual;}
        int    gap()           const {return m_gap;} 
        int    type()          const {return m_type;}
        double conEnergy()     const {return m_ConEnergy;}
        KalFitTrack *track()         {return m_track;} 
        
    private:    
        float m_deflection;    // End point deflection of line of 1st two hits
        float m_sigma;         // Number of sigma deflection corresponds to
        float m_qual;          // Resulting track Quality
        int m_gap;             // Size of gap (hopefully zero!) 
        int m_type;            // Track type: 100*PR-Mode + 10*Energy-Mode + Lead_Hits
        double m_ConEnergy;    // Constraind energy results
        KalFitTrack *m_track;  // The trial track fit
    };

    /// Access methods for getting the individual candidate tracks
    typedef std::vector<Candidate*> CandidateList; 
    typedef std::vector<Candidate*>::iterator iterator;

    CandidateList& candidates() {return m_candidates;}
    iterator begin()            {return m_candidates.begin();}
    iterator end()              {return m_candidates.end();}

    /// Major Sub sections 
    void searchCandidates(double CalEnergy, Point CalPosition);
    void setEnergies(double CalEnergy); 
    void loadOutput(); 
    
    /// Internal drivers
    void findBlindCandidates();
    void findCalCandidates();

    /// Internal utilities
    float findNextHit(int, Ray&, float&);
    bool  incorporate(Candidate*);
    
    /// Data members
    CandidateList m_candidates;  // Internal list of found hypothesises

    Point m_Pcal;      // Calorimeter seed point
    Point m_nextHit;   // Space point transfer space
    double m_cut;      // Sigma cut used 
    double m_energy;   // Energy used to compute errors
    double m_arclen;   // arclength transfer space 
    int m_BestHitCount;// highest hit count on a track this event
    int m_TopLayer;    // Upper most layer in which a track was found
    int m_firstLayer;  // Find first hit layer once

    /// Pointers to clusters, geometry, and control parameters
    ITkrGeometrySvc* m_tkrGeo;
    TkrClusterCol*   m_clusters;
    TkrControl*      m_control; 

};

#endif