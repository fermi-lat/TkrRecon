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

using namespace Event;

//-----------------------------------------------------
// 
//   TrackFitUtils
//
//-----------------------------------------------------

TrackFitUtils::TrackFitUtils(Event::TkrClusterCol* clusters, ITkrGeometrySvc* geo, IFitHitEnergy* hitEnergy) :
                             m_clusters(clusters),
                             m_tkrGeo(geo),
                             m_tkrFail(m_tkrGeo->getTkrFailureModeSvc()),
                             m_hitEnergy(hitEnergy)
{
    // Purpose and Method: Constructor - Initialization for TrackFitUtils
  
    // Inputs:  Pointer to Detector Geometery, pointer to TrkClusters
    // Outputs: No output, object creation
    // Dependencies: None
    // Restrictions and Caveats:  None

    // Set up control
    m_control = TkrControl::getPtr();
}


void TrackFitUtils::flagAllHits(const TkrKalFitTrack& track, int iflag)
{
   // Purpose and Method: Flag all clusters on a given track as having been used
   // Inputs: Track, flag that is passed on the TkrCluster (!= 0 means flagged) 
   // Outputs: None
   // Dependencies: None
   // Restrictions and Caveats:  None

    TkrFitPlaneConPtr hitPtr = track.begin();

    while(hitPtr != track.end()) {
        const TkrFitPlane& hitplane = *hitPtr++; 
        m_clusters->flagHit(hitplane.getIDHit(), iflag);
    }
}

void TrackFitUtils::unFlagAllHits(const TkrKalFitTrack& track)
{
   // Purpose and Method: Unflag all clusters on a given track as having been used
   // Inputs: Track
   // Outputs: None
   // Dependencies: None
   // Restrictions and Caveats:  None

    TkrFitPlaneConPtr hitPtr = track.begin();

    while(hitPtr != track.end()) {
        const TkrFitPlane& hitplane = *hitPtr++; 
        m_clusters->unflagHit(hitplane.getIDHit());
    }
}  

void TrackFitUtils::unFlagHit(const TkrKalFitTrack& track, int num)
{
    // Purpose and Method: Unflag a specfic cluster on a given track
    // Inputs: Track, Index of cluster to unflag
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    TkrFitPlaneConPtr hitPtr = track.begin();
    hitPtr += num;
    const TkrFitPlane& hitplane     = *hitPtr;
    m_clusters->unflagHit(hitplane.getIDHit());
}  

TkrKalFitTrack* TrackFitUtils::newFitTrack(TkrPatCand& patCand)
{
    // Purpose and Method: create a new TkrKalFitTrack object from a pattern recognition track
    // Inputs: Pointer to the pattern recognition track to reproduce
    // Outputs: A pointer to the new TkrKalFitTrack object
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    //Go through each candidate and pass to the fitter
    //int    iniLayer = patCand.getLayer();
    //int    iniTower = patCand.getTower();
    Ray    testRay  = patCand.getRay();
    double energy   = patCand.getEnergy();
    int    type     = (int)(patCand.getQuality()); 
    int    numHits  = patCand.numPatCandHits();
        
    TkrKalFitTrack* track  = new Event::TkrKalFitTrack();

    track->setType(type); 
    track->setInitialPosition(testRay.position());
    track->setInitialDirection(testRay.direction());
    track->setStartEnergy(energy);
        
    //Now fill the hits from the pattern track
    Event::CandHitVectorPtr candPtr = patCand.getHitIterBegin();
    while(numHits--)
    {
        const Event::TkrPatCandHit* candHit = *candPtr++;
        addMeasHit(*track, patCand, *candHit);
    }

    return track;
}


void TrackFitUtils::addMeasHit(TkrKalFitTrack& track, const TkrPatCand& patCand, const TkrPatCandHit& candHit)
{
    // Purpose and Method: Add a measured hit (TkrCluster) from pattern recognition
    //                     to the given track 
    // Inputs: track, hit from Pat. Rec. Track
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    track.push_back(newMeasPlane(candHit, m_hitEnergy->initialHitEnergy(patCand, candHit, track.getStartEnergy())));

    if (candHit.View() == TkrCluster::X) track.setNumXHits(track.getNumXHits()+1);
    else                                 track.setNumYHits(track.getNumYHits()+1);
    
    return;
}

TkrFitPlane TrackFitUtils::newMeasPlane(const TkrPatCandHit& candHit, const double energy)
{
    // Purpose and Method: Add a measured hit (TkrCluster) from pattern recognition
    //                     to the given track 
    // Inputs: track, hit from Pat. Rec. Track
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    Point            planePos  = candHit.Position();
    TkrCluster::view planeView = candHit.View();
    int              clusIdx   = candHit.HitIndex();
    int              towerIdx  = m_clusters->getHit(clusIdx)->tower();
    TkrFitHit        measHit   = makeMeasHit(m_clusters->position(clusIdx), planeView);

    TkrFitPlane      newPlane(clusIdx, towerIdx, candHit.PlaneIndex(), energy, 
                              planePos.z(), measHit, planeView);
    
    return newPlane;
}

void TrackFitUtils::addNewHit(TkrFitPlane& plane, TkrFitHit::TYPE type, TkrFitPar& statePar, TkrFitMatrix& stateCovMat)
{
    // Purpose and Method: Add a new set of track parameters and their covariance matrix to a track hit
    // Inputs: The plane to add info to, the type of parameters, the parameters and their covariance matrix
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    TkrFitHit    stateHit(type,statePar,stateCovMat);

    plane.setHit(stateHit);

    return;
}

void TrackFitUtils::updateMaterials(TkrFitPlane& plane, TkrFitMatrix& Qmat, double radLen, double actDist, double energy)
{
    // Purpose and Method: Updates the materials parameters in a given plane
    // Inputs: The plane to update, the multiple scattering error matrix, radiation lengths, active 
    //         distance and hit energy
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    plane.setRadLen(radLen);
    plane.setActiveDist(actDist);

    // Change energy on the condition it is not negative...
    if (double newEnergy = m_hitEnergy->updateHitEnergy(energy,radLen) >= 0.) plane.setEnergy(newEnergy);

    plane.setQmaterial(Qmat);

    return;
}

TkrFitHit TrackFitUtils::makeMeasHit(const Point& x0, const TkrCluster::view& planeView)
{
    // Purpose and Method: Creates a new measured TkrFitHit which can be added to a track
    // Inputs: The hit coordinates and the hit view (orientation)
    // Outputs: A TkrFitHit object for a measured point
    // Dependencies: None
    // Restrictions and Caveats:  None
    const double oneOverSqrt12 = 1./sqrt(12.);

    double sigma     = m_tkrGeo->siResolution();
    double sigma_alt = m_tkrGeo->trayWidth() * oneOverSqrt12;
    double cx,cy;

    if (planeView == TkrCluster::X)
    {
        cx = sigma * sigma;
        cy = sigma_alt * sigma_alt;
    }
    else
    {
        cy = sigma * sigma;
        cx = sigma_alt * sigma_alt;
    }

    TkrFitPar measPar(x0.x(), 0., x0.y(), 0.);
    TkrFitMatrix measCov(1);

    measCov(1,1) = cx;
    measCov(3,3) = cy;

    TkrFitHit measHit(TkrFitHit::MEAS, measPar, measCov);

    return measHit;
}

void TrackFitUtils::finish(TkrKalFitTrack& track)
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
    if (track.getChiSquare() >= 0) 
    {
        int          nplanes    = track.getNumHits();
        TkrFitPlane& firstPlane = track[0];
        double       x          = firstPlane.getHit(TkrFitHit::SMOOTH).getPar().getXPosition();
        double       y          = firstPlane.getHit(TkrFitHit::SMOOTH).getPar().getYPosition();
        double       z          = firstPlane.getZPlane(); 
        Point        x0         = Point(x,y,z);
        
        double       x_slope    = firstPlane.getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
        double       y_slope    = firstPlane.getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
        Vector       dir        = Vector(-1.*x_slope,-1.*y_slope,-1.).unit();
        
        track.setInitialPosition(x0);
        track.setInitialDirection(dir);
        //track.setChiSquare(track.getChiSquare() / (nplanes-4.)); // 1 measurement per plane - 4 parameters in 3D fit
        //track.setChiSquareSmooth(track.getChiSquareSmooth() / (nplanes-4.));  
        //track.setChiSquare(track.getChiSquare() / (2.*nplanes-4.)); // 1 measurement per plane - 4 parameters in 3D fit
        //track.setChiSquareSmooth(track.getChiSquareSmooth() / (2.*nplanes-4.));  
        track.setScatter(0.);
        
        TkrFitPlaneColPtr hitPtr = track.begin();
        
        int last_Xplane      = -1; 
        int last_Yplane      = -1;
        int num_xPlanes      =  0;
        int num_yPlanes      =  0;
        int Xgaps            =  0;
        int Ygaps            =  0;
        double rmsResid      =  0.;
        double start_energy  = track.getEnergy(); 
        double cos_inv       = 1./fabs(track.getDirection().z()); 
        double z0            = 0; 
        double rad_len       = 0.; 
        int plane_count      = 0; 
        int numSegmentPoints = 0; 
        bool quit_first      = false; 
        while(hitPtr != track.end()) {
            plane_count++; 
            const TkrFitPlane& hit = *hitPtr++;
            
            if(plane_count > 2 && !quit_first) {
                if(plane_count == 3) z0 = hit.getZPlane();
                rad_len += hit.getRadLen(); 
                double arc_len  = (z0- hit.getZPlane())*cos_inv; 
                double theta_ms = 13.6/start_energy * sqrt(rad_len) *
                    (1. + .038*log(rad_len));
                double plane_err = cos_inv*arc_len*theta_ms/1.7321; 
                quit_first  = plane_err > 2.*m_tkrGeo->siStripPitch();
            }
            if(!quit_first) numSegmentPoints++;
            
            int this_plane = hit.getIDPlane();
            bool xPlane = (hit.getProjection() == TkrCluster::X); 
            double xm; 
            if (xPlane) {
                num_xPlanes++;
                x  = hit.getHit(TkrFitHit::SMOOTH).getPar().getXPosition();
                xm = hit.getHit(TkrFitHit::MEAS).getPar().getXPosition();
                if(last_Xplane > 0) {
                    Xgaps += this_plane-last_Xplane-1; 
                    if(num_xPlanes < 3 || !quit_first) {
                        track.setNumXFirstGaps(Xgaps);
                    }
                }
                last_Xplane = this_plane; 
            }
            else {
                num_yPlanes++; 
                x  = hit.getHit(TkrFitHit::SMOOTH).getPar().getYPosition();
                xm = hit.getHit(TkrFitHit::MEAS).getPar().getYPosition();
                if(last_Yplane >= 0) {
                    Ygaps += this_plane-last_Yplane-1;
                    if(num_yPlanes < 3 || !quit_first) {
                        track.setNumYFirstGaps(Ygaps);
                    }
                } 
                last_Yplane = this_plane;
            }
            rmsResid+= (x-xm)*(x-xm);
        }
        rmsResid=sqrt(rmsResid/(1.*nplanes));
        
        track.setScatter(rmsResid);
        track.setNumSegmentPoints(numSegmentPoints);
        track.setNumXGaps(Xgaps);
        track.setNumYGaps(Ygaps);
        
        // Energy calculations
        eneDetermination(track);
        
        // Segment Calculation
        if (track.getChiSquare() >= 0) 
        {
            track.setChiSqSegment(computeChiSqSegment(track, track.getNumSegmentPoints()));
            track.setQuality(computeQuality(track));
        }   
        
        // Compute the radiation lengths to the calorimeter front face
        double arc_min = (track.getInitialPosition().z() - m_tkrGeo->calZTop())/fabs(track.getInitialDirection().z()); 
        IKalmanParticle* TkrFitPart = m_tkrGeo->getPropagator();
        TkrFitPart->setStepStart(x0, dir, arc_min);
        track.setTkrCalRadLen(TkrFitPart->radLength()); 
    }
    

}

double TrackFitUtils::computeQuality(const TkrKalFitTrack& track) const
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
    
    // Calc. How many hits are possible?
    int num_max = m_tkrGeo->numPlanes() - 2*track.getLayer();
    if(num_max > 16) num_max = 16;
    
    // Don't allow more than 8 of each projection
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

void TrackFitUtils::eneDetermination(TkrKalFitTrack& track)
{
    // Purpose and Method:Computes the track energy from the amount
    //     of multiple scattering alongthe track. (refered to as 
    //     the KalEnergy). 
    //         ** Original code taken from KalFitter and written by Bill Atwood **
    // Inputs: None
    // Outputs: sets kalEnergy and its error for this track
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    int nplanes = track.size()-2; // Assume last 2 hits are x,y pair
    // No new info. here - using SMOOTHED
    
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
    int old_Plane_Id = track[0].getIDPlane(); 
    
    for (int iplane = 0; iplane < nplanes; iplane++) { 
        // Get the last cluster size for range estimation
        TkrCluster::view hit_proj = track[iplane].getProjection();
        int hit_Id = track[iplane].getIDHit();
        if(hit_proj == TkrCluster::X) {
            x_cls_size = m_clusters->size(hit_Id);
        }
        else {
            y_cls_size = m_clusters->size(hit_Id);
        }
        if( track[iplane].getIDPlane() == old_Plane_Id) {
            sX += track[iplane].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
            sY += track[iplane].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
            radLen += track[iplane].getRadLen(); 
            count += 1.; 
            if(iplane != nplanes-1) continue;
        }
        totalRad += radLen;    
        Vector t1 = Vector(-sX/count, -sY/count, -1.).unit();
        
        if(t0.magnitude() < .1) { // Trap for first plane
            t0 = t1; 
            sX = track[iplane].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
            sY = track[iplane].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
            radLen = track[iplane].getRadLen();
            count =1.;
            old_Plane_Id = track[iplane].getIDPlane(); 
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
        old_Plane_Id = track[iplane].getIDPlane(); 
        t0 = t1; 
        sX = track[iplane].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
        sY = track[iplane].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
        radLen = track[iplane].getRadLen();
        count =1.;
    }
    
    // Set a max. energy based on range 
    //        - use cluster size as indicator of range-out
    double prj_size_x = m_tkrGeo->siThickness()*fabs(sX)/
        m_tkrGeo->siStripPitch()             + 1.;
    double prj_size_y = m_tkrGeo->siThickness()*fabs(sY)/
        m_tkrGeo->siStripPitch()             + 1.;
    double range_limit = 100000.;  // 100 GeV max... 
    if((x_cls_size - prj_size_x) > 2 || (y_cls_size - prj_size_y) > 2) {
        range_limit = totalRad * 50.; // 10 MeV = 15% rad. len 
    }

    double kalEnergy = track.getStartEnergy();
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

double TrackFitUtils::computeChiSqSegment(const TkrKalFitTrack& track, int nhits, TkrFitHit::TYPE typ)
{
    // Purpose and Method: Computes the chisquared for the first
    //            portion of the track
    //         ** Original code taken from KalFitter **
    // Inputs: no. of hits to include and the chisquared type (FILTER or SMOOTH)
    // Outputs: chisquared/D.F.
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    double chi2 = 0;
    int ihit =0;
    for (ihit =0; ihit < nhits; ihit++) {
        chi2 += track[ihit].getDeltaChiSq(typ);
    }
    chi2 /= (1.*nhits);
    return chi2;
}

void TrackFitUtils::setSharedHitsStatus(TkrKalFitTrack& track)
{
    // Purpose and Method: Sets the shared hit status for a given track
    //         ** Original code taken from KalFitter **
    // Inputs: The track t set
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    // Hits are shared depending on cluster size 
    // and track direction
    TkrFitPlaneConPtr pln_pointer = track.begin();
                
    int i_Hit = 0; 
    int i_share = 0;
    while(pln_pointer != track.end() && i_Hit < 6) 
    {
        // First 2 hits (x & y) are shared
        if(i_Hit < 2) 
        { 
            unFlagHit(track, i_Hit);
            i_Hit++;
            i_share++;
            pln_pointer++;
            continue;
        }
        // For the rest - unflag according to Cluster size and Trajectory
        Event::TkrFitPlane plane = *pln_pointer;
        Event::TkrCluster::view hit_proj = plane.getProjection();
        Event::TkrFitPar tkr_par = plane.getHit(Event::TkrFitHit::FIT).getPar();
        double slope = tkr_par.getYSlope();
        if(hit_proj == Event::TkrCluster::X) 
        {
            slope = tkr_par.getXSlope();
        }        
        int hit_Id = plane.getIDHit();;
        double cls_size = m_clusters->size(hit_Id);        
        double prj_size = m_tkrGeo->siThickness()*fabs(slope)
                         /m_tkrGeo->siStripPitch() + 1.;
        if(cls_size> prj_size) 
        {
            unFlagHit(track, i_Hit);
            i_share++;
        }
        if(i_share >= 5) break; 
        i_Hit++;
        pln_pointer++;
    }

    return;
}


/*

TkrFitHit TrackFitUtils::initialFitHit(const TkrFitPar& initialPar, const TkrFitMatrix& baseCovMat)
{  
    // Purpose and Method: Set parameters for the first hit
    //          Errors are somewhat arbitrary. 
    // Inputs: Input Parameters from Pat.Rec.
    // Outputs: a TkrFitHit
    // Dependencies: None
    // Restrictions and Caveats:  None

    //  The first error is arbitrary to a degree
    TkrFitMatrix first_errors; 
    first_errors(2,2) = m_control->getIniErrSlope() * m_control->getIniErrSlope();
    first_errors(4,4) = m_control->getIniErrSlope() * m_control->getIniErrSlope();
    
    TkrFitHit hitf(TkrFitHit::FIT, initialPar, baseCovMat + first_errors);
    
    return hitf;
}

TkrFitMatrix TrackFitUtils::computeMeasCov(const TkrFitPar& newPars, const TkrFitMatrix& oldCovMat, 
                                           const TkrCluster& cluster)
{
    // Purpose and Method: Set parameters for the first hit
    //          Errors are somewhat arbitrary. 
    // Inputs: Input Parameters from Pat.Rec.
    // Outputs: a TkrFitHit
    // Dependencies: None
    // Restrictions and Caveats:  None
   
    // Compute the Measurement covariance taking into account the 
    // Local track slope
    const double oneOverSqrt12 = 1. / sqrt(12.);

    // The following sets the error to the slop between the track 
    // and the cluster over sqrt(12). It protects against getting
    // too small.

    TkrFitMatrix newCov(1);
    ////double min_err = 0.707 * m_tkrGeo->siResolution(); 
    double min_err = m_tkrGeo->siResolution(); 
    double clusWid = const_cast<TkrCluster&>(cluster).size();

    if(cluster.v() == TkrCluster::X) 
    {
        double x_slope  = newPars.getXSlope();
        double wid_proj = fabs(x_slope * m_tkrGeo->siThickness());
        double wid_cls  = clusWid * m_tkrGeo->siStripPitch();
        double error    = (wid_cls - wid_proj) * oneOverSqrt12;
        ////double error    = errorFactor(clusWid, x_slope) * m_tkrGeo->siStripPitch() * oneOverSqrt12;
        ////if (error == 0.)
        ////{
        ////    error = (wid_cls - wid_proj) * oneOverSqrt12;
        ////    error = (error > min_err) ? error : min_err; 
        ////}
        ////error = m_tkrGeo->siResolution();

        error       = (error > min_err) ? error : min_err; 
        newCov(1,1) = error * error;
        newCov(3,3) = oldCovMat(3,3);
    }
    else 
    {
        double y_slope  = newPars.getYSlope();
        double wid_proj = fabs(y_slope * m_tkrGeo->siThickness());
        double wid_cls  = clusWid * m_tkrGeo->siStripPitch();
        double error    = (wid_cls - wid_proj) * oneOverSqrt12;
        ////double error    = errorFactor(clusWid, y_slope) * m_tkrGeo->siStripPitch() * oneOverSqrt12;
        ////if (error == 0.)
        ////{
        ////    error = (wid_cls - wid_proj) * oneOverSqrt12;
        ////    error = (error > min_err) ? error : min_err; 
        ////}
        ////error = m_tkrGeo->siResolution();
    
        error = (error > min_err) ? error : min_err; 
        newCov(1,1) = oldCovMat(1,1);
        newCov(3,3) = error * error;
    }

    return newCov;
}

double TrackFitUtils::errorFactor(double strips, double slope) 
{
    // Kludgey code that returns the factor by which to multiply width/sqrt(12)

    // strips is the number of strips in the cluster
    // slope is the slope of the track in the measuring view


    double stripAspect = 0.57;  // 228/400
    double absSlope = fabs(slope/stripAspect);

    // calculation below is done in units of strips
    // absSlope = 1 is the slope that crosses one strip exactly

    // This is an empircal fit the the "measured" errors as a function
    // of number of strips in cluster and slope of the track.
    // The "1.05" accounts for the threshold.
    // Not sure about the 0.6 yet.

    // For now, return zero if outside the range of applicability.
    // For clusters wider than expected, maybe max(0.707, fabs(meas - projected)) is a good guess.
    // For clusters narrower than expected, there must be missing strips,
    // so the error should also be larger, perhaps again max(0.707, fabs(meas-projected))

    // actually, we can do better... most of the 2nd case are tracks going through the edge
    // a wafer, so we can "fix" them post facto.

    bool oldErrors = false;

    double factor = 0.0;
    double eps0 = 0.0; // use this to extent or restrict the valid range for 1-strip clusters
    double eps1 = 0.0; // ditto for the rest of the clusters
    double loSlope, hiSlope, peakSlope;
    double loPar1, hiPar1, peakDev;

    if(oldErrors) { return 1.0;}

    int nStrips = floor(strips+.01);  // just to be safe
    if (nStrips==1) {
        if (absSlope<1.5+eps0) factor = 1 - 0.52*absSlope;
    } else if (nStrips<11) {
        if (nStrips==2) {
            loSlope = .5 ; hiSlope = 2.5; peakSlope = 1.61;
            peakDev = 0.97; loPar1 = .613; hiPar1 = .697;
        } else if (nStrips==3) {
            loSlope = 1.8 ; hiSlope = 3.5; peakSlope = 2.78;
            peakDev = 0.97; loPar1 = .600; hiPar1 = .759;
        } else if (nStrips==4) {
            loSlope = 3.0 ; hiSlope = 4.6; peakSlope = 3.80;
            peakDev = 0.90; loPar1 = .691; hiPar1 = .755;
        } else if (nStrips==5) {
            loSlope = 4.2 ; hiSlope = 5.6; peakSlope = 4.83;
            peakDev = 0.94; loPar1 = .769; hiPar1 = .819;
        } else if (nStrips>=6) {
            double nm6 = 1.03*(nStrips - 6);
            loSlope = 5.0 + nm6 ; hiSlope = 6.6 + nm6; peakSlope = 5.88 + nm6;
            peakDev = 0.96; loPar1 = .714; hiPar1 = .851;
        }
        if (absSlope>loSlope-eps1 && absSlope < peakSlope) {
            factor = peakDev - loPar1*(peakSlope - absSlope);
        } else if (absSlope>peakSlope && absSlope < hiSlope + eps1 ) {
            factor = peakDev - hiPar1*(absSlope - peakSlope);
        }
    }

    //double factor = 0.0;
    //
    //int nStrips = floor(strips+.01);  // just to be safe
    //if (nStrips==1) {
    //    if (absSlope<1.4) factor = 1 - 0.5*absSlope;
    //} else if (nStrips==2) {
    //    if (absSlope>.4 && absSlope<2.6) factor = 0.9 - 0.6*fabs(absSlope-1.6);
    //} else if (nStrips<11) {
    //    if (fabs(absSlope-(2.7+1.05*(nStrips-3)))<1.) 
    //        factor = 0.9 - 0.6*fabs(absSlope - (2.7 + 1.05*(nStrips-3)));
    //}

    return factor;
}
*/
