/*
	Code to implement the PatRecTracks class
	Tracy Usher Nov 27, 2000
*/

#include "TkrRecon/PatRec/TkrCandidates.h"

TkrCandidates::TkrCandidates()
{
    m_candTracks.clear();

	return;
}

TkrCandidates::~TkrCandidates()
{
    int              numCands = getNumCands();
    CandTrkVectorPtr cands    = getTrackPtr();

    while(numCands--) delete *cands++;

	return;
}


void TkrCandidates::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG << " --- TkrCandidates::writeOut --- " << endreq;
    log << MSG::DEBUG << "     Number of candidate tracks: " << getNumCands() << endreq;

    int trackNo = 0;

    while(trackNo < getNumCands())
    {
        log << MSG::DEBUG << "     Track Number: " << trackNo << endreq;
        m_candTracks[trackNo++]->writeOut(log);
    }

    return;
}

void TkrCandidates::addTrack(TkrPatCand* candTrack)
{
    m_candTracks.push_back(candTrack);

    return;
}