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
    
    /// Utilities 
    void writeOut(MsgStream& log) const; 
        
protected:	
    
private:
    TkrVector m_tracks;
};

//Following typedefs for containing fit track objects
typedef std::vector<TkrVertex*>            TkrVertexVec;
typedef std::vector<TkrVertex*>::iterator  TkrVertexVecPtr;

#endif
