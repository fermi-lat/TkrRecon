//---------------------------------------------------
//   GFcandidates
//
//     Search for GFxxx candidates (xxx = particle or gamma)
//---------------------------------------------------

#ifndef __GFcandidates_H
#define __GFcandidates_H 1

#include <vector>
#include "src/PatRec/Combo/GFgamma.h"
#include "TkrRecon/Track/GFcontrol.h"

//##########################################################
class GFcandidates 
//##########################################################
{
public:
    
    enum type {TRACK,PARTICLE,PAIR,GAMMA};
    
public:
    
    // construction
    GFcandidates(type t, double ene, double cut, Point Pend = Point(0.,0.,0.),
                   Point Pini = Point(0.,0.,0.));
    ~GFcandidates() {}
    void clear();
    
    // access
    int numCandidates() {return (int) m_candidates.size();}
    
    std::vector<GFdata> m_candidates;
    
    // utilities
    static GFdata GFconstructor(GFcandidates::type , double ene, int ilayer, 
        const Ray testRay);
    
private:
    
    // internal drivers
    void ini();
//    bool findCandidates();
    
    // lower level itnernal drivers
    bool findCandidates(std::vector<GFdata>& candidates, 
        const GFdata& candidate,
        double ene, GFcandidates::type);
    
    bool findSeedCandidates(std::vector<GFdata>& candidates, 
        GFcandidates::type, double ene); 
    bool findSeedCandidates(std::vector<GFdata>& candidates, 
        GFcandidates::type , double ene, int iplane, int itower = 0);
    
    // internal utilities
    static void incorporate(std::vector<GFdata>& candidates, const GFdata);
    Point createPend(int ilayer, const Point& PIni);
    
private:
    
    GFcandidates::type m_type;
    GFcandidates::type m_seedtype;
    
    Point m_Pini;
    Point m_Pend;
    
    double m_eneCandidate;
    
};


#endif
