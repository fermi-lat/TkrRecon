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
#include "TkrRecon/Track/TkrBase.h"
#include "TkrRecon/Track/TkrFitTrack.h"

class TkrVertex: public TkrBase
{    
public:
    
    TkrVertex( int layer, int tower, double energy, const Ray& testRay);
    ~TkrVertex() {}

    // Add tracks to the list
    void addTrack(TkrFitTrack* pTrack) {m_tracks.push_back(pTrack);}
    
    // How many tracks in the vertex?
    int  getNumTracks()                {return m_tracks.size();}

    // Pointers to track info
    TkrFitTrackColPtr getTrackPtr()    {return m_tracks.begin();}
    TkrFitTrackColPtr getTrackEnd()    {return m_tracks.end();}

    /// Utilities 
    void writeOut(MsgStream& log) const; 
        
protected:	
    
private:
    TkrFitTrackCol m_tracks;
};

//Following typedefs for containing fit track objects
typedef std::vector<TkrVertex*>            TkrVertexVec;
typedef std::vector<TkrVertex*>::iterator  TkrVertexVecPtr;

#endif
