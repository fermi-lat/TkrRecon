/*
	This class defines the tree part of the "link and tree" pattern 
	recognition algorithm. 
	Tracy Usher Nov 29, 2000
*/
#ifndef __TKRCOMBOPATREC_H
#define __TKRCOMBOPATREC_H

#include "Event/Recon/TkrRecon/TkrPatCandCol.h"
#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "src/TrackFit/KalFitTrack/KalFitTrack.h"
#include "TkrRecon/Track/TkrBase.h"
#include "geometry/Ray.h"

#include <vector>

using namespace TkrRecon;

class TkrComboPatRec : public TkrPatCandCol
{
public:
	TkrComboPatRec(ITkrGeometrySvc* pTkrGeo, TkrClusterCol* pClusters, double CalEnergy, Point CalPosition);
   ~TkrComboPatRec();
private:

    class Candidate: public TkrBase 
    {
    public:
        Candidate(int layer, int twr, double e, 
                   Point x, Vector t, float d, float s, int g);
        ~Candidate() {};

        // access
        void setDeflection(float def) {m_deflection = def;}
        void setSigma(float sig)      {m_sigma      = sig;}
        void setQuality(float Q)      {m_qual       = Q;}
        void setGap(int gs)           {m_gap        = gs;}
        float  deflection()    const {return m_deflection;}
        float  sigma()         const {return m_sigma;}
        float  quality()       const {return m_qual;}
        int    gap()           const {return m_gap;} 
        
    private:    
        float m_deflection;    // End point deflection of line of 1st two hits
        float m_sigma;         // Number of sigma deflection corresponds to
        float m_qual;          // Resulting track Quality
        int m_gap;             // Size of gap (hopefully zero!) 
    };


    // Access methods for getting the individual candidate tracks
    typedef std::vector<Candidate> CandidateList; 
    typedef std::vector<Candidate>::const_iterator const_iterator;

    const CandidateList& candidates() const {return m_candidates;}
    const_iterator begin()            const {return m_candidates.begin();}
    const_iterator end()              const {return m_candidates.end();}

    void searchCandidates(double CalEnergy, Point CalPosition);
    
    // internal drivers
    void findBlindCandidates();
    void findCalCandidates();

    // internal utilities
    float findNextHit(int, float, Ray&, float&, int&);
    void  incorporate(Candidate&);
    
    // data members
    CandidateList  m_candidates;  // List of found hypothesises
    TkrFitCol      m_tracks;      // List of attempted fits

    Point m_Pcal;      // Calorimeter seed point
    Point m_nextHit;   // Space point transfer space
    double m_cut;      // Sigma cut used 
    double m_energy;   // Energy used to compute errors
    double m_arclen;   // arclength transfer space 

};

#endif