
#ifndef __SIRECOBJS_H
#define __SIRECOBJS_H 1

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataObject.h"
#include "src/PatRec/Combo/GFcandidates.h"
#include "TkrRecon/Track/TkrFitTrack.h"
#include "GlastEvent/Recon/ISiRecObjs.h"

#include "gui/DisplayRep.h"

//----------------------------------------------
//
//   SiRecObjs
//
//   Transient Storage Data
//----------------------------------------------
//   It contains the high level recostructed objects
//   of the tracker (gammas and tracks)
//----------------------------------------------
//             J.A Hernando, Santa Cruz 02/29/00
//----------------------------------------------

extern const CLID& CLID_SiRecObjs;
//const static CLID CLID_SiRecObjs = 254;

/*!
SiRecObjs container and service class of the Tracker Reconstructed objects (gamma and tracks)
*/
//##########################################################
class SiRecObjs : public ISiRecObjs
//##########################################################
{
    
public:
    
    //! constructor - ini the list
    SiRecObjs(ITkrGeometrySvc* pTkrGeo, TkrClusters* pTkrClus, double CalEnergy, Point CalPosition);
    //! destructor - delete the pointers to the data
    virtual ~SiRecObjs() {clear();}
    
    // GAUDI members to be use by the converters
    static const CLID& classID() {return CLID_SiRecObjs;}
    virtual const CLID& clID() const {return classID();}
    
    //! add a GFgamma to the list
    void addGamma(GFgamma* g)       {m_GFgammaList.push_back(g);}
    //! add a GFparticle to the list
    void addParticle(TkrFitTrack* p) {m_trackList.push_back(p);}
    
    void searchGammas(double CalEnergy, Point CalPosition);
        //! returns number of gammas
    int numGammas() const			{return m_GFgammaList.size();}
    //! returns pointer to GFgamma in i position
    GFgamma*  Gamma(int i)   const              {return m_GFgammaList[i];}
    
    void searchParticles(double CalEnergy, Point CalPosition);
    //! returns number of GFparticles
    int numParticles() const			{return m_GFparticleList.size();}
     //! returns number of GFparticles
    int numTracks() const			{return m_trackList.size();}
    //! returns pointer to GFparticle in position i
    GFparticle* Particle(int i) const		{return m_GFparticleList[i];}
    
    //! draws the SiRecObjs
    void update(gui::DisplayRep& v) {draw(v);}
    
    //! clear the list of SiRecObjs
    virtual void clear();
    //! empty method
    virtual void make() {}
    //! write out the information of the SiRecObjs
    virtual void writeOut(MsgStream& log) const;
    
    //new methods required for the interface
    //! Get the X slope of the ith GFparticle
    double getXGFparticleSlope(int i) {return m_GFparticleList[i]->getXGFtrack()->k_direction().x();}
    //! Get the Y slope of the ith GFparticle
    double getYGFparticleSlope(int i) { return m_GFparticleList[i]->getYGFtrack()->k_direction().y();}
    //! get the vextex of the ith Gamma
    Point getGammaVertex(int i) { return m_GFgammaList[i]->vertex(); }
    //! get the direction vector of the ith Gamma
    Vector getGammaDirection(int i) { return Vector(0,0,-1.); }
    //! Get the X slope of the ith GFgamma
    double getXGFgammaSlope(int i) {return m_GFgammaList[i]->getPair(TkrCluster::X)->k_direction().x(); }
    //! Get the Y slope of the ith GFgamma
    double getYGFgammaSlope(int i) { return m_GFgammaList[i]->getPair(TkrCluster::Y)->k_direction().y(); }
    //	void draw(GraphicsRep& v);
    void draw(gui::DisplayRep& v);
    
private:
    
    //! ini the lists
    virtual void ini();
    //! draws the objects

    
private:
    
    //! vector with the list of GFparticle reconstructed
    std::vector<TkrFitTrack*> m_trackList;
    //! vector with the list of GFparticle reconstructed
    std::vector<GFparticle*> m_GFparticleList;
    //! vector with the list of GFgamma reconstructed
    std::vector<GFgamma*> m_GFgammaList;
};

#endif
