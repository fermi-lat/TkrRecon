//------------------------------------------------------------------------------
//
//     TrackFitUtils
//
//     Provides utility routines for the generic Kalman Filter track fit
//
//     Adapted from existing TkrRecon code, mostly from KalFitter originally
//     authored by Bill Atwood. 
//              
//------------------------------------------------------------------------------


#include "TrackFitUtils.h" 
#include "src/Track/TkrControl.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrFailureModeSvc.h"
#include "src/Track/TkrControl.h"
#include "idents/TowerId.h"

#include <cmath>

//using namespace Event;

//-----------------------------------------------------
// 
//   TrackFitUtils
//
//-----------------------------------------------------

TrackFitUtils::TrackFitUtils(ITkrGeometrySvc* tkrGeom, IFitHitEnergy* hitEnergy) :
m_tkrGeom(tkrGeom),
m_tkrFail(m_tkrGeom->getTkrFailureModeSvc()),
m_hitEnergy(hitEnergy)
{
    // Purpose and Method: Constructor - Initialization for TrackFitUtils

    // Inputs:  Pointer to Detector Geometery, pointer to TrkClusters
    // Outputs: No output, object creation
    // Dependencies: None
    // Restrictions and Caveats:  None

    // Set up control
    m_control = TkrControl::getPtr();
    // vector for track compare
    m_hitVec.resize(m_tkrGeom->numPlanes());
}


void TrackFitUtils::flagAllHits(Event::TkrTrack& track, int iflag)
{
    // Purpose and Method: Flag all clusters on a given track as having been used
    // Inputs: Track, flag that is passed on the TkrCluster (!= 0 means flagged) 
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    Event::TkrTrackHitVecItr hitPtr = track.begin();

    while(hitPtr != track.end()) {
        Event::TkrTrackHit* hitplane = *hitPtr++; 

        if (!(hitplane->getStatusBits() & Event::TkrTrackHit::HITONFIT)) continue;

        Event::TkrClusterPtr cluster = hitplane->getClusterPtr();

        cluster->flag(iflag);
    }
}

void TrackFitUtils::unFlagAllHits(Event::TkrTrack& track)
{
    // Purpose and Method: Unflag all clusters on a given track as having been used
    // Inputs: Track
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    Event::TkrTrackHitVecItr hitPtr = track.begin();

    while(hitPtr != track.end()) {
        Event::TkrTrackHit* hitplane = *hitPtr++; 

        if (!(hitplane->getStatusBits() & Event::TkrTrackHit::HITONFIT)) continue;

        Event::TkrClusterPtr cluster = hitplane->getClusterPtr();

        cluster->unflag();
    }
}  

void TrackFitUtils::unFlagHit(Event::TkrTrack& track, int num)
{
    // Purpose and Method: Unflag a specfic cluster on a given track
    // Inputs: Track, Index of cluster to unflag
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    Event::TkrTrackHitVecItr hitPtr = track.begin();
    hitPtr += num;
    Event::TkrTrackHit* hitplane     = *hitPtr;

    if (hitplane->getStatusBits() & Event::TkrTrackHit::HITONFIT)
    {
        Event::TkrClusterPtr cluster = hitplane->getClusterPtr();

        cluster->unflag();
    }
}  

void TrackFitUtils::finish(Event::TkrTrack& track)
{
    // Purpose and Method: Kalman clean-up.  Translates fit parameters
    //         back to direction cosines and a point, computes the 
    //         chisquareds, computes the Kalman Energy (MS determined energy)
    //         ** Original code taken from KalFitter **
    // Inputs: None
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    // Compute the fit variables  
    if (track.getChiSquareSmooth() >= 0) 
    {
        int                 nplanes    = track.getNumHits();
        Event::TkrTrackHit* firstPlane = track[0];
        Point               x0         = firstPlane->getPoint(Event::TkrTrackHit::SMOOTHED);
        Vector              dir        = firstPlane->getDirection(Event::TkrTrackHit::SMOOTHED);

        track.setInitialPosition(x0);
        track.setInitialDirection(dir);
        track.setScatter(0.);

        int    last_XLayer      = -1; 
        int    last_YLayer     = -1;
        int    num_xPlanes      =  0;
        int    num_yPlanes      =  0;
        int    Xgaps            =  0;
        int    Ygaps            =  0;
        double rmsResid         =  0.;
        double start_energy     = track.getInitialEnergy(); 
        double cos_inv          = 1./fabs(track.getInitialDirection().z()); 
        double z0               = x0.z(); 
        double rad_len          = 0.; 
        int    plane_count      = 1; 
        int    numSegmentPoints = 0; 
        bool    quit_first      = false; 

        // Loop over the hits on the track. Note planes are layers of SSDs
        for(Event::TkrTrackHitVecItr hitPtr = track.begin(); hitPtr != track.end(); hitPtr++)
        {
            Event::TkrTrackHit* hit = *hitPtr;
            rad_len += hit->getRadLen(); 

            if(plane_count > 4 && !quit_first) {
                double arc_len  = (z0- hit->getZPlane())*cos_inv; 
                double theta_ms = 13.6/start_energy * sqrt(rad_len) *
                    (1. + .038*log(rad_len));
                double plane_err = cos_inv*arc_len*theta_ms/1.7321; 
                quit_first  = plane_err > 2.*m_tkrGeom->siStripPitch();
            }

            if ((hit->getStatusBits()& Event::TkrTrackHit::HITONFIT) && hit->validCluster()) { 

                if(!quit_first) numSegmentPoints++;

                int this_plane = m_tkrGeom->trayToPlane(hit->getTkrId().getTray(),hit->getTkrId().getBotTop());
                bool xPlane = hit->getTkrId().getView() == idents::TkrId::eMeasureX; 

                if(hit->validSmoothedHit()) {
                    double x  = hit->getMeasuredPosition(Event::TkrTrackHit::SMOOTHED);
                    double xm = hit->getMeasuredPosition(Event::TkrTrackHit::MEASURED);

                    int this_layer, view;
                    m_tkrGeom->planeToLayer(this_plane, this_layer, view);  

                    if (xPlane) {
                        num_xPlanes++;
                        if(last_XLayer >= 0) {
                            Xgaps += last_XLayer - this_layer - 1; 
                            if(num_xPlanes < 3 || !quit_first) track.setNumXFirstGaps(Xgaps);
                        }
                        last_XLayer = this_layer; 
                    }
                    else {
                        num_yPlanes++; 
                        if(last_YLayer>= 0) {
                            Ygaps += last_YLayer- this_layer - 1;
                            if(num_yPlanes < 3 || !quit_first) track.setNumYFirstGaps(Ygaps);
                        } 
                        last_YLayer= this_layer;
                    }
                    rmsResid+= (x-xm)*(x-xm);
                    plane_count++; 
                }
            }
        }
        rmsResid=sqrt(rmsResid/(1.*plane_count));

        track.setScatter(rmsResid);
        track.setNumSegmentPoints(numSegmentPoints);
        track.setNumXHits(num_xPlanes);
        track.setNumYHits(num_yPlanes); 
        track.setNumXGaps(Xgaps);
        track.setNumYGaps(Ygaps);

        // Energy calculations
        eneDetermination(track);

        // Segment Calculation
        if (track.getChiSquareFilter() >= 0) {
            track.setChiSqSegment(computeChiSqSegment(track, track.getNumSegmentPoints()));
            track.setQuality(computeQuality(track));
        }   

        // Compute the radiation lengths to the calorimeter front face
        double arc_min = (track.getInitialPosition().z() - m_tkrGeom->calZTop())/fabs(track.getInitialDirection().z()); 
        IKalmanParticle* TkrFitPart = m_tkrGeom->getPropagator();
        TkrFitPart->setStepStart(x0, dir, arc_min);
        track.setTkrCalRadLen(TkrFitPart->radLength()); 
    }
}

double TrackFitUtils::computeQuality(const Event::TkrTrack& track) const
{ 
    // Purpose and Method: Ascribes a "quality" for each track.
    //         Somewhat arbitrary - main components are number
    //         of hits and Smoothed Chisq.  Long tracks are 
    //         high quality - tracks are penalized by the size 
    //         of chisquared. 
    //         ** Original code taken from KalFitter **
    // Inputs: None
    // Outputs: Quality parameter
    // Dependencies: None
    // Restrictions and Caveats:  None

    // Determine plane number 
    idents::TkrId hitId = track[0]->getTkrId();
    int           layer = m_tkrGeom->trayToBiLayer(hitId.getTray(),hitId.getBotTop());

    // Calc. How many hits are possible?
    //int num_max = 2*(m_tkrGeom->numPlanes() - layer);
    int num_max = 2*(layer+1); //Layers start at zero...
    if(num_max > 16) num_max = 16;

    // Don't allow more then 8 of each projection
    int num_Hits  = (track.getNumXHits() <= 8) ? track.getNumXHits() : 8;
    num_Hits     += (track.getNumYHits() <= 8) ? track.getNumYHits() : 8;

    // Scale to max. allowed 
    float hit_count_factor = (1.*num_Hits)/(1.*num_max);

    // Overall factors are to make this ~ match older def's 
    double quality = 64.*hit_count_factor - 2.*sqrt(track.getChiSquareSmooth()); 


    //    double quality = 4*(m_nxHits+m_nyHits-4. - (m_Xgaps+m_Ygaps)) 
    //                    -2.*sqrt(m_chisqSmooth); 
    return quality;
}

void TrackFitUtils::eneDetermination(Event::TkrTrack& track)
{
    // Purpose and Method:Computes the track energy from the amount
    //     of multiple scattering alongthe track. (refered to as 
    //     the KalEnergy). 
    //         ** Original code taken from KalFitter and written by Bill Atwood **
    // Inputs: None
    // Outputs: sets kalEnergy and its error for this track
    // Dependencies: None
    // Restrictions and Caveats:  None

    double totalRad  = 0.;
    double eneSum    = 0.;
    double thetaSum  = 0.;
    double sX        = 0.;
    double sY        = 0.; 
    double radLen    = 0.; 
    double count     = 0.; 
    double eSumCount = 0.;
    double tSumCount = 0.;
    double x_cls_size = 1.;
    double y_cls_size = 1.;

    Vector t0(0.,0.,0.); 

    // Keep track of the current bilayer
    idents::TkrId hitId = track[0]->getTkrId();
    int old_BiLayer_Id = m_tkrGeom->trayToBiLayer(hitId.getTray(),hitId.getBotTop());

    for (Event::TkrTrackHitVecItr planeIter = track.begin(); planeIter != track.end(); planeIter++)
    { 
        // Get the last cluster size for range estimation
        const Event::TkrTrackHit& plane   = **planeIter;

        if (!((plane.getStatusBits() & Event::TkrTrackHit::HITONFIT)&& plane.validSmoothedHit())) continue;

        // Valid hit, get cluster and id info
        Event::TkrClusterPtr cluster = plane.getClusterPtr();
        idents::TkrId        hitId   = plane.getTkrId();
        int                  biLayer = m_tkrGeom->trayToBiLayer(hitId.getTray(),hitId.getBotTop());

        if (hitId.getView() == idents::TkrId::eMeasureX) x_cls_size = cluster->size();
        else                                             y_cls_size = cluster->size();

        if( biLayer == old_BiLayer_Id) 
        {
            sX += plane.getTrackParams(Event::TkrTrackHit::SMOOTHED).getxSlope();
            sY += plane.getTrackParams(Event::TkrTrackHit::SMOOTHED).getySlope();
            radLen += plane.getRadLen(); 
            count += 1.; 
            if(planeIter != track.end()) continue;
        }
        totalRad += radLen;    
        Vector t1 = Vector(-sX/count, -sY/count, -1.).unit();

        if(t0.magnitude() < .1) { // Trap for first plane
            t0 = t1; 
            sX = plane.getTrackParams(Event::TkrTrackHit::SMOOTHED).getxSlope();
            sY = plane.getTrackParams(Event::TkrTrackHit::SMOOTHED).getySlope();
            radLen = plane.getRadLen();
            count =1.;
            old_BiLayer_Id = biLayer; 
            continue; 
        }

        double e_factor = exp(-totalRad);        
        double t0t1 = t0*t1;
        double theta = acos(t0t1);
        double rl_factor = radLen;

        eSumCount += 1.; 
        tSumCount += 1./rl_factor; 
        eneSum += (theta * e_factor)*(theta * e_factor)/rl_factor; 
        thetaSum += theta * theta /rl_factor; 

        // Reset for next X,Y measuring plane
        old_BiLayer_Id = biLayer; 
        t0 = t1; 
        sX = plane.getTrackParams(Event::TkrTrackHit::SMOOTHED).getxSlope();
        sY = plane.getTrackParams(Event::TkrTrackHit::SMOOTHED).getySlope();
        radLen = plane.getRadLen();
        count =1.;
    }

    // Set a max. energy based on range 
    //        - use cluster size as indicator of range-out
    double prj_size_x = m_tkrGeom->siThickness()*fabs(sX)/
        m_tkrGeom->siStripPitch()             + 1.;
    double prj_size_y = m_tkrGeom->siThickness()*fabs(sY)/
        m_tkrGeom->siStripPitch()             + 1.;
    double range_limit = 100000.;  // 100 GeV max... 
    if((x_cls_size - prj_size_x) > 2 || (y_cls_size - prj_size_y) > 2) {
        range_limit = totalRad * 50.; // 10 MeV = 15% rad. len 
    }

    double kalEnergy = track.getInitialEnergy();
    double kalEneErr = 10000.;
    double thetaMS   = 0.;

    if (eSumCount > 0)
    {
        double e_inv = sqrt(eneSum  /2./eSumCount);

        kalEnergy = 13.6 / e_inv;
        if(kalEnergy < m_control->getMinEnergy()/3.) kalEnergy = m_control->getMinEnergy()/3.; 
        if(kalEnergy > range_limit) kalEnergy = range_limit;
        kalEneErr = kalEnergy/sqrt(eSumCount);
        thetaMS   = sqrt(thetaSum/2./tSumCount);
    }
    else
    {
        // just a place to set a breakpoint
        thetaMS = 0.0;
    }

    track.setKalThetaMS(thetaMS);

    track.setKalEnergy(kalEnergy); //Units MeV
    track.setKalEnergyError(kalEneErr);
}

double TrackFitUtils::computeChiSqSegment(const Event::TkrTrack& track, int nhits, Event::TkrTrackHit::ParamType typ)
{
    // Purpose and Method: Computes the chisquared for the first
    //            portion of the track
    //         ** Original code taken from KalFitter **
    // Inputs: no. of hits to include and the chisquared type (FILTER or SMOOTH)
    // Outputs: chisquared/D.F.
    // Dependencies: None
    // Restrictions and Caveats:  None

    double chi2 = 0;

    if (typ == Event::TkrTrackHit::FILTERED)
    {
        for (int ihit=0; ihit<nhits; ihit++) chi2 += track[ihit]->getChiSquareFilter();
    }
    else
    {
        for (int ihit=0; ihit<nhits; ihit++) chi2 += track[ihit]->getChiSquareSmooth();
    }

    chi2 /= (1.*nhits);
    return chi2;
}

void TrackFitUtils::setSharedHitsStatus(Event::TkrTrack& track, int maxShare)
{
    // Purpose and Method: Sets the shared hit status for a given track
    //         ** Original code taken from KalFitter **
    // Inputs: The track t set
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    // Hits are shared depending on cluster size 
    // and track direction
    Event::TkrTrackHitVecItr plane = track.begin();

    int i_Hit = 0; 
    int i_share = 0;

    while(plane != track.end()) 
    {
        // First hits (x & y) are shared: gamma conversion vertex
        if(i_Hit < 2 && (*plane)->validCluster()) { 
            unFlagHit(track, i_Hit);
            i_Hit++;
            i_share++;
            plane++;
            continue;
        }
        if(!(*plane)->validCluster()) {
            plane++;
            i_Hit++;
            continue;
        }

        // For the rest - unflag according to Cluster size and Trajectory
        Event::TkrClusterPtr cluster = (*plane)->getClusterPtr();

        double slope = (*plane)->getMeasuredSlope(Event::TkrTrackHit::FILTERED);

        double cls_size = cluster->size();        
        double prj_size = m_tkrGeom->siThickness()*fabs(slope)
            /m_tkrGeom->siStripPitch() + 1.;
        if(cls_size> prj_size) 
        {
            unFlagHit(track, i_Hit);
            i_share++;
        }
        if(i_share >= maxShare) break; 
        i_Hit++;
        plane++;
    }

    return;
}

int TrackFitUtils::numUniqueHits(Event::TkrTrack& track1,
                                 Event::TkrTrack& track2,
                                 int minUnique) const
{
    // Alternative to compareTracks(), above
    // This one counts unique hits on the shorter track, 
    // by looking at each layer (once!)

    int size = m_hitVec.size();
    m_hitVec.assign(size, 0);

    int nHit1 = track1.getNumFitHits();
    int nHit2 = track2.getNumFitHits();
    Event::TkrTrack& trackToCompare = (nHit1<nHit2 ? track1 : track2);
    Event::TkrTrack& trackToLoad    = (nHit1<nHit2 ? track2 : track1);
    Event::TkrTrackHitVecItr hitPtr;
    Event::TkrTrackHit* hit;
    int plane;

    /* debug stuff... keep this for a while
    hitPtr = track1.begin();
    hit = *hitPtr;
    plane = m_tkrGeom->getPlane(hit->getTkrId());
    int nxHits = track1.getNumXHits();
    int nyHits = track1.getNumYHits();
    std::cout << "Track 1 1stPlane " << plane << 
        ", lastPlane " << plane-track1.getNumHits()+1<< std::endl;
    std::cout << "nXhits " << nxHits << " nYhits " << nyHits 
        << " nFitHits " << track1.getNumFitHits() << std::endl;
    // now look at each hit:
    for (; hitPtr!= track1.end(); ++hitPtr,--plane) {
        hit = *hitPtr;
        int status = hit->getStatusBits();
        std::cout <<"Plane " << plane << " validCl " << hit->validCluster();
        if (hit->validCluster()) {
            bool measuresX = ((status & Event::TkrTrackHit::MEASURESX)!=0);
            bool measuresY = ((status & Event::TkrTrackHit::MEASURESY)!=0);
            std::cout << "measX/Y " << measuresX << " " << measuresY ;
        }
        std::cout << std::endl;
    }

    hitPtr = track2.begin();
    hit = *hitPtr;
    plane = m_tkrGeom->getPlane(hit->getTkrId());
    std::cout << "Track 2 1stPlane " << plane << 
        ", lastPlane " << plane-track2.getNumHits()+1<< std::endl;
    */

    int nUnique = 0;
    hitPtr = trackToLoad.begin();
    plane = m_tkrGeom->getPlane((*hitPtr)->getTkrId());

    for (; hitPtr!= trackToLoad.end(); ++hitPtr,--plane) {
        hit = *hitPtr;
        /*
        if(plane<0 || plane>m_tkrGeom->numPlanes()-1) {
            idents::TkrId tkrId = hit->getTkrId();
            std::cout << "tray/bottop" << tkrId.getTray() << " " << tkrId.getBotTop() << std::endl;
            std::cout << "plane " << plane <<" bad!" << std::endl;
        }
        */
        m_hitVec[plane] = hit;
    }

    hitPtr = trackToCompare.begin();
    Event::TkrTrackHit* hit1 = *hitPtr;
    plane = m_tkrGeom->getPlane(hit1->getTkrId());
    for (; hitPtr!=trackToCompare.end(); ++hitPtr,--plane) {
        hit1 = *hitPtr;
        /*
        if(plane<0 || plane>m_tkrGeom->numPlanes()-1) {
            std::cout << "plane " << plane <<" bad!" << std::endl;
        }
        */
        hit = m_hitVec[plane];

        // if the other hit is not there, hit is unique
        nUnique++;
        if (hit!=0) {
            bool valid  = hit->validCluster();
            bool valid1 = hit1->validCluster();
            // if both clusters are valid and equal, hit is shared
            if (valid && valid1) {
                Event::TkrClusterPtr cl = hit->getClusterPtr();
                Event::TkrClusterPtr cl1 = hit1->getClusterPtr();
                if(cl==cl1) nUnique--;
            } 
            // if both clusters are invalid, hit is possibly "shared"
            // could check distance between to be more sure...
            else if (!valid && !valid1) nUnique--;
        }
        // bail early if 3rd arg is supplied, and there are enuf hits
        if (nUnique>=minUnique) return nUnique;
    }
    return nUnique;
}


int TrackFitUtils::compareTracks(Event::TkrTrack& track1, 
                                 Event::TkrTrack& track2,
                                 int stopAfter) const
{
    // stopAfter defaults to a large number, so loops go to completion
    //  if method is called with only 2 arguments

    int num_sharedHits = 0;
    /*
    int num_compares = 0;
    int minHits = std::min(track1.getNumHits(), track2.getNumHits());
    //std::cout << "Begin compareTrack for " << minHits << " hits" << std::endl;
    */
    // Loop over the hits on the track1
    Event::TkrTrackHitVecItr hitPtr1, hitPtr2;  
    hitPtr1 = track1.begin();
    for (; hitPtr1 != track1.end(); ++hitPtr1) {
        Event::TkrTrackHit* hit1 = *hitPtr1;
        if(!hit1->validCluster()) continue;
        Event::TkrClusterPtr pClus1 = hit1->getClusterPtr();

        // Loop over the hits on track2
        hitPtr2 = track2.begin();
        for (;hitPtr2 != track2.end(); ++hitPtr2) {
            Event::TkrTrackHit* hit2 = *hitPtr2;
            //num_compares++;
            if(!hit2->validCluster()) continue;
            if(pClus1 == hit2->getClusterPtr()) {
                //std::cout << "CT: plane " 
                //    << m_tkrGeom->getPlane(hit2->getTkrId())<<std::endl;
                num_sharedHits++;
                if (num_sharedHits>stopAfter) {
                  //std::cout << num_compares << " comparisons, left early" 
                  //  << std::endl;
                    return num_sharedHits;
                }
                // found a match for this hit, on to the next one
                break;
            }
        }
    }
    //std::cout << num_sharedHits << " overlaps, " << num_compares 
    //    << " comparisons" << std::endl;
    return num_sharedHits;
}

double TrackFitUtils::firstKinkNorm(Event::TkrTrack& track)
{
    Event::TkrTrackHitVecItr hitPtr = track.begin();
    Event::TkrTrackHit* hit1 = *hitPtr;
    Vector t0 = hit1->getDirection(Event::TkrTrackHit::SMOOTHED);

    int layer0 = m_tkrGeom->getLayer(hit1->getTkrId());
    double energy = hit1->getEnergy(); 
    double rad_len = 0.; 
    Vector t1;
    hitPtr++; 
    while(hitPtr != track.end()) {
        Event::TkrTrackHit* hit = *hitPtr;
        int next_layer = m_tkrGeom->getLayer(hit->getTkrId()); 
        if(abs(next_layer - layer0) > 1) break;
        t1 = hit->getDirection(Event::TkrTrackHit::SMOOTHED);
        rad_len += hit->getRadLen();
        hitPtr++;
    }
    double kink_angle = acos(t0*t1);
    double mscat_angle = 13.6 * sqrt(rad_len) * ( 1 + .038*log(rad_len))*
        1.414 /energy; 
    return kink_angle/mscat_angle;
}


