/*
	Code to implement the PatRecTracks class
	Tracy Usher Nov 27, 2000
*/

#include "TkrRecon/Track/TkrTracks.h"
#include "TkrRecon/Track/GFtutor.h"

TkrTracks::TkrTracks(ITkrGeometrySvc* pTkrGeo, TkrClusters* pTkrClus)
{
    // Store cluster and geometry information for the subclasses
    GFtutor::load(pTkrClus, pTkrGeo);

    m_Tracks.clear();

	return;
}

TkrTracks::~TkrTracks()
{
    int               numTracks = getNumTracks();
    TkrFitTrackColPtr trks      = getTrackPtr();

    while(numTracks--) delete *trks++;
}


void TkrTracks::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG << " --- TkrTracks::writeOut --- " << endreq;
    log << MSG::DEBUG << "     Number of tracks: " << getNumTracks() << endreq;

    int trackNo = 0;

    while(trackNo < getNumTracks())
    {
        log << MSG::DEBUG << "     Track Number: " << trackNo << endreq;
        m_Tracks[trackNo++]->writeOut(log);
    }

    return;
}

void TkrTracks::draw(gui::DisplayRep& v)
{
    v.setColor("blue");

    int numTracks = getNumTracks();

    if (numTracks > 0) 
    {
        int trkIdx = 0;

        while(trkIdx < numTracks)
        {
            TkrFitTrack* track = getTrack(trkIdx++);

            track->draw(v);
        }
    }

    return;
}

void TkrTracks::addTrack(TkrFitTrack* Track)
{
    m_Tracks.push_back(Track);

    return;
}