//
//-----------------------------------------------------------------
//
//  TkrVertex
//
//  Class definition for a found vertex 
//  ** Test Version **
//
//  Adapted from TkrFitTrack.h
//  Tracy Usher March 1, 2002
//
//-----------------------------------------------------------------
//
#ifndef __TkrVertex_H
#define __TkrVertex_H 1

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "TkrRecon/Track/TkrRecInfo.h"
#include "TkrRecon/Track/TkrFitTrack.h"

class TkrVertex: public TkrRecInfo
{    
public:
    
    TkrVertex( int layer, int tower, double energy, double quality, const Ray& testRay);
   ~TkrVertex() {}

    /// Define the TkrBase routines
    double        getQuality()                      const {return m_quality;};
    double        energy(TrackEnd end = Start)      const {return m_energy;}
    int           layer(TrackEnd end = Start)       const {return m_firstLayer;}
    int           tower(TrackEnd end = Start)       const {return m_itower;}
    Point         position(TrackEnd end = Start)    const {return m_position;}
    Vector        direction(TrackEnd end = Start)   const {return m_direction;}
    Ray           ray(TrackEnd end = Start)         const {return Ray(position(),direction());}
    TkrFitPar     TrackPar(TrackEnd end = Start)    const {return m_vertexPar;}
    double        TrackParZ(TrackEnd end = Start)   const {return m_position.z();}
    TkrFitMatrix  TrackCov(TrackEnd end = Start)    const {return m_vertexCov;}
    bool          empty()                           const {return m_firstLayer >= 0;}

    // Add tracks to the list
    void addTrack(TkrFitTrack* pTrack) {m_tracks.push_back(pTrack);}
    
    // How many tracks in the vertex?
    int  getNumTracks()                {return m_tracks.size();}

    // Pointers to track info
    TkrFitTrackColPtr getTrackPtr()    {return m_tracks.begin();}
    TkrFitTrackColPtr getTrackEnd()    {return m_tracks.end();}

    /// Utilities 
    void writeOut(MsgStream& log) const; 
    
private:
    TkrFitPar    m_vertexPar;
    TkrFitMatrix m_vertexCov;
    Point        m_position;
    Vector       m_direction;
    double       m_energy;
    double       m_quality;
    int          m_firstLayer;
    int          m_itower;        
    
    TkrFitTrackCol m_tracks;
};

//Following typedefs for containing fit track objects
typedef std::vector<TkrVertex*>            TkrVertexVec;
typedef std::vector<TkrVertex*>::iterator  TkrVertexVecPtr;

#endif
