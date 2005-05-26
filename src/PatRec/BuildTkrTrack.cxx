//
//

#include "BuildTkrTrack.h"

BuildTkrTrack::BuildTkrTrack(const ITkrGeometrySvc* tkrGeo) : m_tkrGeom(tkrGeo)
{
    // Set up control
    m_control = TkrControl::getPtr();

    return;
}

Event::TkrTrack* BuildTkrTrack::makeNewTkrTrack(Point startPos, 
                                                Vector startDir, 
                                                double energy, 
                                                std::vector<const Event::TkrCluster*>& clusterVec)
{
    Event::TkrTrack* track = new Event::TkrTrack();
            
    track->setInitialPosition(startPos);
    track->setInitialDirection(startDir);
    track->setInitialEnergy(energy);
    track->setStatusBit(Event::TkrTrack::LATENERGY);

    for(std::vector<const Event::TkrCluster*>::iterator clusIter = clusterVec.begin(); 
        clusIter != clusterVec.end(); 
        clusIter++)
    {
        const Event::TkrCluster* cluster = *clusIter;

        Event::TkrTrackHit* trackHit = makeTkrTrackHit(cluster);

        track->push_back(trackHit);
    }

    setFirstHitParams(track);

    return track;
}

Event::TkrTrackHit* BuildTkrTrack::makeTkrTrackHit(const Event::TkrCluster* cluster)
{
    Event::TkrTrackHit* trackHit = new Event::TkrTrackHit(const_cast<Event::TkrCluster*>(cluster), 
                                                          cluster->getTkrId(),
                                                          cluster->position().z(),   
                                                          0., 0., 0., 0., 0.);

    // Retrieve a reference to the measured parameters (for setting)
    Event::TkrTrackParams& params = trackHit->getTrackParams(Event::TkrTrackHit::MEASURED);

    // Set measured track parameters
    params(Event::TkrTrackParams::xPosIdx) = cluster->position().x();
    params(Event::TkrTrackParams::xSlpIdx) = 0.;
    params(Event::TkrTrackParams::yPosIdx) = cluster->position().y();
    params(Event::TkrTrackParams::ySlpIdx) = 0.;

    int measIdx = trackHit->getParamIndex(Event::TkrTrackHit::SSDMEASURED,    Event::TkrTrackParams::Position);
    int nonmIdx = trackHit->getParamIndex(Event::TkrTrackHit::SSDNONMEASURED, Event::TkrTrackParams::Position);

    double sigma     = m_tkrGeom->siResolution();
    double sigma_alt = m_tkrGeom->trayWidth()/sqrt(12.);

    params(measIdx,measIdx) = sigma * sigma;
    params(nonmIdx,nonmIdx) = sigma_alt * sigma_alt;

    // Last: set the hit status bits
    unsigned int status_bits = Event::TkrTrackHit::HITONFIT | Event::TkrTrackHit::HASMEASURED |
                               Event::TkrTrackHit::HITISSSD | Event::TkrTrackHit::HASVALIDTKR;

    if (cluster->getTkrId().getView() == idents::TkrId::eMeasureX) status_bits |= Event::TkrTrackHit::MEASURESX;
    else                                                           status_bits |= Event::TkrTrackHit::MEASURESY;

    trackHit->setStatusBit((Event::TkrTrackHit::StatusBits)status_bits);

    return trackHit;
}

bool BuildTkrTrack::setFirstHitParams(Event::TkrTrack* track)
{
    bool successful = false;

    // We can only do this if the track contains hits already
    if (track->size() > 0)
    {
        Event::TkrTrackHit* trackHit = track->front();

        trackHit->setEnergy(track->getInitialEnergy());
        
        int    measIdx   = trackHit->getParamIndex(Event::TkrTrackHit::SSDMEASURED,    Event::TkrTrackParams::Position);
        int    nonmIdx   = trackHit->getParamIndex(Event::TkrTrackHit::SSDNONMEASURED, Event::TkrTrackParams::Position);
        int    xSlpIdx   = Event::TkrTrackParams::xSlpIdx;
        int    ySlpIdx   = Event::TkrTrackParams::ySlpIdx;
        double sigma     = m_tkrGeom->siResolution();
        double sigma_alt = m_tkrGeom->trayWidth()/sqrt(12.);

        Vector startDir  = track->getInitialDirection();
        Point  startPos  = track->getInitialPosition();

        double x_slope   = startDir.x()/startDir.z();
        double y_slope   = startDir.y()/startDir.z();
        Event::TkrTrackParams firstParams(startPos.x(), x_slope, startPos.y(), y_slope,
                                          5., 0., 0., 0., 0., 0., 0., 5., 0., 0.);

        firstParams(measIdx)          = trackHit->getMeasuredPosition(Event::TkrTrackHit::MEASURED);
        firstParams(measIdx, measIdx) = sigma * sigma;
        firstParams(nonmIdx, nonmIdx) = sigma_alt * sigma_alt;
        firstParams(xSlpIdx, xSlpIdx) = m_control->getIniErrSlope() * m_control->getIniErrSlope();
        firstParams(ySlpIdx, ySlpIdx) = m_control->getIniErrSlope() * m_control->getIniErrSlope();

        // Now do the same for the FILTERED params
        Event::TkrTrackParams& filtPar = trackHit->getTrackParams(Event::TkrTrackHit::FILTERED);
        filtPar = firstParams;

        // And now do the same for the PREDICTED params
        Event::TkrTrackParams& predPar = trackHit->getTrackParams(Event::TkrTrackHit::PREDICTED);
        predPar = firstParams;

        // Last: set the hit status bits
        unsigned int status_bits = trackHit->getStatusBits();

        status_bits |= Event::TkrTrackHit::HASPREDICTED 
                    |  Event::TkrTrackHit::HASFILTERED 
                    |  Event::TkrTrackHit::HASVALIDTKR;

        // Update the TkrTrackHit status bits
        trackHit->setStatusBit((Event::TkrTrackHit::StatusBits)status_bits);

        successful = true;
    }

    return successful;
}
