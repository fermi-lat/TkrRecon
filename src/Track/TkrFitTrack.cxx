//------------------------------------------------------------------------------
//
//     TkrFitTrack
//
//     Implementation of the Tracker Fit Track output class
//				
//------------------------------------------------------------------------------


#include "TkrRecon/Track/TkrFitTrack.h"
#include "TkrRecon/Track/GFtutor.h"
#include "TkrRecon/Track/GFcontrol.h"

// Constructor takes no arguments and basically just initializes to
// a null state. In the class' current incarnation, it is expected to
// be inherited from a class which actually does the track fit. 
TkrFitTrack::TkrFitTrack()
{ 
    m_chisq       = -9999.;
    m_chisqSmooth = -9999.;
    m_rmsResid    = 0.;
    m_Q           = -1e6;
    m_KalEnergy   = 0.;
    m_KalThetaMS  = 0.;
    m_hits.clear();
}

TkrFitTrack::~TkrFitTrack()
{
    clear();
}

bool TkrFitTrack::empty() const
{
    bool empty = false;
    if (layer()        < 0)                         empty = true;
    if (getNumHits()   < GFcontrol::minSegmentHits) empty = true;
    if (getChiSquare() < 0.)                        empty = true;

    return empty;
}

void TkrFitTrack::clear()
{   
    m_hits.clear();
    
    m_Q            = -1e6;
}

Vector TkrFitTrack::direction(TrackEnd end) const
{
    TkrFitPar trkPar = getFoLPlane(end).getHit(TkrFitHit::SMOOTH).getPar();

    if   (end == TkrRecInfo::Start) return Vector(-trkPar.getXSlope(),-trkPar.getYSlope(),-1.).unit();
    else                            return Vector( trkPar.getXSlope(), trkPar.getYSlope(), 1.).unit();
}



void TkrFitTrack::writeOut(MsgStream& log) const
{
    
    log << MSG::DEBUG << " --- TkrFitTrack::writeOut --- "            << endreq;
    log << MSG::DEBUG << " Position      = " << position().x() << " " <<position().y() << " " << position().z() << endreq;
    log << MSG::DEBUG << " Direction     = " << direction().x() << " " << direction().y() << " " << direction().z() << endreq;
    log << MSG::DEBUG << " Energy        = " << energy() << endreq;
    log << MSG::DEBUG << " first Layer   = " << layer() << endreq;
    log << MSG::DEBUG << " Tower         = " << tower() << endreq;
    log << MSG::DEBUG << " quality       = " << getQuality()       << endreq;
    log << MSG::DEBUG << " num m_hits    = " << getNumHits()       << endreq;
}


TkrFitPlane TkrFitTrack::getFoLPlane(TrackEnd end) const
{
    if (m_hits.size() == 0) 
    {
        return TkrFitPlane();
    }
    else
    {
        if (end == TkrRecInfo::Start) return m_hits.front();
        else                          return m_hits.back();
    }
}
