//---------------------------------------------------
//   TkrComboPR Implementation
//
//   Search for track candidates by combinatorics
// 
//    W. Atwood, Nov., 2001  
//---------------------------------------------------

#ifndef __TkrComboPR_H
#define __TkrComboPR_H 1

#include <vector>
#include "TkrRecon/Track/TkrBase.h"


class TkrComboPR 
{
    
public:
    
    // constructor
    TkrComboPR(double cut, double calEne = 0., Point calHit = Point(0.,0.,0.));
    ~TkrComboPR() {}
    
    // access
    int numCandidates() {return (int) m_candidates.size();}

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
    
private:
    
    // internal drivers
    void findBlindCandidates();
    void findCalCandidates();

    // internal utilities
    float findNextHit(int, float, Ray&, float&, int&);
    void  incorporate(Candidate);
    
    // data members
    CandidateList m_candidates;  // List of found hypothesises

    Point m_Pcal;      // Calorimeter seed point
    Point m_nextHit;   // Space point transfer space
    double m_cut;      // Sigma cut used 
    double m_energy;   // Energy used to compute errors
    double m_arclen;   // arclength transfer space 
    
};

#endif
