
#ifndef PatRecTracks_h
#define PatRecTracks_h

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "TkrRecon/PatRec/TkrPatCand.h"
#include "GaudiKernel/DataObject.h"

extern const CLID& CLID_TkrCandidates;

//
//------------------------------------------------------------------------
//
// TkrCandidates
//
// Transient Data Object for the Pattern Recognition state of Tracker Reconstruction
//
//------------------------------------------------------------------------
//

class TkrCandidates : public DataObject
{
public:
	TkrCandidates();
   ~TkrCandidates();

    //Provide ability to output some information about self...
    void writeOut(MsgStream& log) const;

	//! GAUDI members to be use by the converters
	static const CLID&  classID()           {return CLID_TkrCandidates;}
	virtual const CLID& clID()        const {return classID();}

    //How many track candidates are there?
    int                 getNumCands() const {return m_candTracks.size();}

    //Access to tracks through an iterator
    CandTrkVectorPtr    getTrackPtr()       {return m_candTracks.begin();}

    //Access to tracks by index
    TkrPatCand*         getTrack(int idx)   {return m_candTracks[idx];}

    //Add tracks to the list
    void                addTrack(TkrPatCand* candTrack);

private:
    CandTrkVector m_candTracks;
};

#endif