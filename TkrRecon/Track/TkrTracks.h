
#ifndef TKRTRACKS_H
#define TKRTRACKS_H

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "TkrRecon/Track/TkrFitTrack.h"
#include "GaudiKernel/DataObject.h"

extern const CLID& CLID_TkrTracks;

//
//------------------------------------------------------------------------
//
// TkrTracks
//
// Transient Data Object for the fit, reconstructed tracks of Tracker Reconstruction
//
//------------------------------------------------------------------------
//

class TkrTracks : public DataObject
{
public:
    TkrTracks() {}
	TkrTracks(ITkrGeometrySvc* pTkrGeo, TkrClusters* pTkrClus);
   ~TkrTracks();

    //Provide ability to output some information about self...
    void writeOut(MsgStream& log) const;
    
    //! draws the individual tracks
    void update(gui::DisplayRep& v) {draw(v);}
    
    //	void draw(GraphicsRep& v);
    void draw(gui::DisplayRep& v);

	//! GAUDI members to be use by the converters
	static const CLID&  classID()           {return CLID_TkrTracks;}
	virtual const CLID& clID()        const {return classID();}

    //How many reconstructed tracks are there?
    int                 getNumTracks() const {return m_Tracks.size();}

    //Access to tracks through an iterator
    TkrVectorPtr        getTrackPtr()        {return m_Tracks.begin();}

    //Access to tracks by index
    TkrFitTrack*        getTrack(int idx)    {return m_Tracks[idx];}

    //Add tracks to the list
    void                addTrack(TkrFitTrack* track);

private:
    TkrVector m_Tracks;
};

#endif