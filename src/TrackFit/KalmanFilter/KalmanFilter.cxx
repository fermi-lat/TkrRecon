
//----------------------------------------------------------------------
//    
//    Implementation of the Kalman Filter Functions 
//               TkrFitPlane
//
//      Original due to Jose Hernando-Angel circa 1997-1999
//      Re-written to combine both X and Y projections (2001)
// 
//      Bill Atwood, SCIPP/UCSC, Nov. 2001
//-----------------------------------------------------------------------

#include "KalmanFilter.h"
#include "Utilities/TkrException.h"

using namespace Event;

namespace {
    enum paramIndex {XPOS=1, XSLOPE=2, YPOS=3, YSLOPE=4};
}



KalmanFilter::KalmanFilter(TkrClusterCol* clusters, ITkrGeometrySvc* geo)
{
    m_clusters   = clusters;
    m_tkrGeo     = geo;

    m_radLength  = 0.;
    m_activeDist = 0;
    m_control    = TkrControl::getPtr();
}

//-------------------------------------
//  Kalman Functions
//-------------------------------------

TkrFitHit KalmanFilter::predicted(TkrFitPlane& start, TkrFitHit::TYPE typ, 
                                  int &nlayers, int klayer, double &zend,
                                  double &arc_min)
{
    // Extrapolate the hit to the next SSD layer - maybe an X or a Y
    //     Note: distance will be > arc_min;   arc_len used returned as arc_min
    // Returns (nlayers)  the number of GLAST Planes crossed 
    // Note: nlayers = 0 for crossing between X and Y in same plane

    TkrFitHit    hit     = start.getHit(typ);
    
    TkrFitPar    pp      = hit.getPar();
    TkrFitMatrix Ck      = hit.getCov();
    
    double       ene     = start.getEnergy();
        if( ene <=0 || ene > 1e6 ){
        throw( TkrException("bad energy in KalmanFilter") );
    }

    double       x_slope = pp.getXSlope();   
    double       y_slope = pp.getYSlope(); 
    Vector       dir_ini = Vector(-x_slope, -y_slope, -1.).unit();

    double       x0      = pp.getXPosition();
    double       y0      = pp.getYPosition();
    double       z0      = start.getZPlane();
    Point        x_ini(x0,y0,z0);

    int          nsteps  = klayer - start.getIDPlane();

    double       down    = nsteps < 0 ? +1. : -1.;

    IKalmanParticle* TkrFitPart = m_tkrGeo->getPropagator();
    TkrFitPart->setStepStart(x_ini, dir_ini, arc_min);
    if(TkrFitPart->trackToNextPlane()) 
    {
        AXIS planeProjection = TkrCluster::Y;
        if(TkrFitPart->isXPlane()) planeProjection = TkrCluster::X; 

        start.setNextProj(planeProjection);

        arc_min = TkrFitPart->arcLength();
        nlayers = static_cast<int> (arc_min*fabs(dir_ini.z())/29.0); //mm
    }
    else 
    {
        nlayers = -1;
        return TkrFitHit(); 
   }
    
    zend    = TkrFitPart->position().z(); 
    double relDeltaZ = down*fabs(arc_min*dir_ini.z());

    TkrFitMatrix F(relDeltaZ);
                          
    m_Qmaterial = TkrFitPart->mScat_Covr(ene, arc_min); 
    pp = F*pp;
    if (down == -1.)     Ck=(F*(Ck*F.T()))+m_Qmaterial;
    else if (down == +1) Ck=(F*(Ck+m_Qmaterial)*F.T());
    
    TkrFitHit hitpred(TkrFitHit::PRED, pp, Ck);

    m_radLength  = TkrFitPart->radLength();
    m_activeDist = TkrFitPart->insideActArea();
 
    
    return hitpred;
} 

TkrFitHit KalmanFilter::predicted(TkrFitPlane& start, TkrFitHit::TYPE typ, 
                                  int /*klayer*/, double &zend,
                                  double &arc_min)
{
    // Extrapolate by arc_min
    // Returns (nlayers)  the number of GLAST Planes crossed 
    // Note: nlayers = 0 for crossing between X and Y in same plane

    TkrFitHit    hit     = start.getHit(typ);
    
    TkrFitPar    pp      = hit.getPar();
    TkrFitMatrix Ck      = hit.getCov();
    
    double       ene     = start.getEnergy();
    if( ene <=0 || ene > 1e6 ){
        throw( TkrException("bad energy in KalmanFilter") );
    }
    double       x_slope = pp.getXSlope();   
    double       y_slope = pp.getYSlope(); 
    Vector       dir_ini = Vector(-x_slope, -y_slope, -1.).unit();

    double       x0      = pp.getXPosition();
    double       y0      = pp.getYPosition();
    double       z0      = start.getZPlane();
    Point        x_ini(x0,y0,z0);

    double       down    = -1.;

    IKalmanParticle* TkrFitPart = m_tkrGeo->getPropagator();
    TkrFitPart->setStepStart(x_ini, dir_ini, arc_min);
    if(arc_min >= 0) {
        AXIS planeProjection = TkrCluster::Y;
        if(TkrFitPart->isXPlane()) planeProjection = TkrCluster::X; 
        start.setNextProj(planeProjection);
    }
    else 
    {
        return TkrFitHit(); 
    }
    
    zend    = TkrFitPart->position().z(); 
    double relDeltaZ = down*fabs(arc_min*dir_ini.z());

    TkrFitMatrix F(relDeltaZ);
                          
    m_Qmaterial = TkrFitPart->mScat_Covr(ene, arc_min); 
    pp = F*pp;
    if (down == -1.)     Ck=(F*(Ck*F.T()))+m_Qmaterial;
    else if (down == +1) Ck=(F*(Ck+m_Qmaterial)*F.T());
    
    TkrFitHit hitpred(TkrFitHit::PRED, pp, Ck);

    m_radLength  = TkrFitPart->radLength();
    m_activeDist = TkrFitPart->insideActArea(); 
    
    return hitpred;
} 

TkrFitHit KalmanFilter::predicted(TkrFitPlane& start, TkrFitPlane& kplaneNext) 
{  
    // Extrapolate the hit to the next Kplane; 
    // Note that GLAST trajectory is evolving in negative z coordinates.
    // The energy is necessary to compute the MS error.
    
    TkrFitHit    hit     = start.getHit(TkrFitHit::FIT);
    
    TkrFitPar    pp      = hit.getPar();
    TkrFitMatrix Ck      = hit.getCov();
    
//    double       ene     = start.getEnergy();
//        if( ene <=0 || ene > 1e6 ){
//        throw( TkrException("bad energy in KalmanFilter") );
//    }
    double ene = 0.;
    try
    {
        ene = start.getEnergy();
    }
    catch(...)
    {
        throw( TkrException("bad energy in KalmanFilter") );
    }

    double       x_slope = pp.getXSlope(); 
    double       y_slope = pp.getYSlope(); 
    Vector       dir_ini = Vector(-x_slope, -y_slope, -1.).unit();

    double       x0      = pp.getXPosition();
    double       y0      = pp.getYPosition();
    double       z0      = start.getZPlane();
    Point        x_ini(x0,y0,z0); 

    int          iplane  = start.getIDPlane();
    int          nsteps  = kplaneNext.getIDPlane() - iplane;
    double       down    = -1.;
    if (nsteps <0 ) 
    {
        down = +1.; // going up;
        dir_ini = -dir_ini;
    }
    

    double       deltaZ  = kplaneNext.getZPlane() - start.getZPlane();
    double       arc_len = fabs(deltaZ/dir_ini.z()); 

    IKalmanParticle* TkrFitPart = m_tkrGeo->getPropagator();
    TkrFitPart->setStepStart(x_ini, dir_ini, arc_len);

    double relDeltaZ = down * fabs(deltaZ);

    TkrFitMatrix F(relDeltaZ);
                          
    m_Qmaterial = TkrFitPart->mScat_Covr(ene, arc_len); 
    pp = F*pp;
    if (down == -1.)     Ck=(F*(Ck*F.T()))+m_Qmaterial;
    else if (down == +1) Ck=(F*(Ck+m_Qmaterial)*F.T());
    
    // store the matrix with the material
    kplaneNext.setQmaterial(m_Qmaterial);

    TkrFitHit hitpred(TkrFitHit::PRED, pp, Ck);

    m_radLength  = TkrFitPart->radLength();
    m_activeDist = TkrFitPart->insideActArea(); 

    return hitpred;
    
} 

TkrFitHit KalmanFilter::filter(TkrFitPlane& filterPlane)
{
    // Filter = Weighting the hits & prediction with their cov matrices

   // Re-compute the meas. cov. matrix taking into account cls-size and
    // track slopes
    TkrFitHit    pred_hit   = filterPlane.getHit(TkrFitHit::PRED);
    computeMeasCov(filterPlane, pred_hit.getPar());  
    
    TkrFitHit    meas_hit   = filterPlane.getHit(TkrFitHit::MEAS);
    TkrFitPar    p_meas     = meas_hit.getPar();

    // The following incorporates the meas.projection matrix H
    TkrFitMatrix G(1); // The meas. weight matrix
    if(filterPlane.getProjection()==TkrCluster::X) {
        G(XPOS,XPOS) = 1./meas_hit.getCov().getcovX0X0();
        G(YPOS,YPOS) = 0.;
    }
    else {
        G(XPOS,XPOS) = 0.;
        G(YPOS,YPOS) = 1./meas_hit.getCov().getcovY0Y0();
    }

    TkrFitMatrix Ckpred    = pred_hit.getCov();
    int          i_error;
    TkrFitMatrix Ckpredinv = Ckpred;
    Ckpredinv.invert(i_error); 
    TkrFitPar p_pred=pred_hit.getPar();
    
 // TkrFitMatrix Ck        = (Ckpredinv+H*(G*H));  
    TkrFitMatrix Ck        = (Ckpredinv+G);
    Ck.invert(i_error);
 // TkrFitPar pk=((Ck*Ckpredinv)*p_pred)+((Ck*(H*G))*p_meas);
    TkrFitPar pk=((Ck*Ckpredinv)*p_pred)+((Ck*G)*p_meas);
    TkrFitHit hitfit(TkrFitHit::FIT, pk, Ck);
    
    return hitfit;   
}

TkrFitHit KalmanFilter::smoother(TkrFitPlane& start, const TkrFitPlane& kplaneLast)
{
    // Apply the smoother to a Fitted hit!
    TkrFitHit    hitf0    = start.getHit(TkrFitHit::FIT);
    
    TkrFitHit    hitpred  = kplaneLast.getHit(TkrFitHit::PRED);
    TkrFitHit    hitsm    = kplaneLast.getHit(TkrFitHit::SMOOTH);
    
    double       distance = kplaneLast.getZPlane() - start.getZPlane();
    TkrFitMatrix F(distance);

    TkrFitMatrix Ck       = hitf0.getCov();

    TkrFitMatrix Ck1pred  = hitpred.getCov();
    TkrFitMatrix Ck1sm    = hitsm.getCov();
    
    TkrFitPar    pk       = hitf0.getPar();
    TkrFitPar    pk1pred  = hitpred.getPar();
    TkrFitPar    pk1sm    = hitsm.getPar();
 
    int i_error; 
    TkrFitMatrix Ck1predinv = Ck1pred;
    Ck1predinv.invert(i_error); 
    TkrFitMatrix A=Ck*((F.T())*Ck1predinv);

    TkrFitPar    psm=pk+A*(pk1sm-pk1pred);
    TkrFitMatrix Csm=Ck+A*((Ck1sm-Ck1pred)*A.T());
    
    TkrFitHit new_hit_sm(TkrFitHit::SMOOTH,psm,Csm);
    
    return new_hit_sm; 
}

void KalmanFilter::computeMeasCov(TkrFitPlane& plane, TkrFitPar pred_pars)
{
    // Compute the Measurement covariance taking into account the 
    // Local track slope

    // Get the measure hit, the prediction and the cluster
    TkrFitHit    meas_hit      = plane.getHit(TkrFitHit::MEAS);
    int id_Cls = plane.getIDHit();
    double size_Cls = (*m_clusters)[id_Cls]->size();

    TkrFitMatrix newCov(1);

    int measured, other;
    double covOther;
    double slope;

    if(plane.getProjection()==TkrCluster::X) {
        slope = pred_pars.getXSlope();
        measured = XPOS;
        other    = YPOS;
        covOther = meas_hit.getCov().getcovY0Y0();
    } else {
        slope = pred_pars.getYSlope();
        measured = YPOS;
        other    = XPOS;
        covOther = meas_hit.getCov().getcovX0X0();
    }

    double error = getError(size_Cls, slope);
    
    newCov(measured, measured) = error*error;
    newCov(other, other)       = covOther;

    TkrFitHit newMeas(TkrFitHit::MEAS, meas_hit.getPar(), newCov);
    plane.setHit(newMeas);   
    return;
}

double KalmanFilter::getError(double strips, double slope) const 
{
    // Kludgey code that returns the error, depending on errorType
    // 0 -> same as before
    // 1 -> first attempt at slope-dependent errors
    // 2 -> second attempt

    // strips is the number of strips in the cluster
    // slope is the slope of the track in the measuring view

    double stripAspect = m_tkrGeo->siThickness()/m_tkrGeo->siStripPitch();
    double absSlope = fabs(slope*stripAspect);

    // calculation below is done in units of strips
    // absSlope = 1 is the slope that crosses one strip exactly

    // For clusters narrower than expected, there must be missing strips,
    // so the error should also be larger, perhaps again max(sqrt(.5), fabs(meas-projected-1))

    // actually, we could do better... most of these are tracks going through the edge
    // a wafer, so we could "fix" them post facto.

    double error;
    double factor = 0.0;
    double minErr = m_tkrGeo->siResolution(); 
    double oneOverSqrtTwelve = 1./sqrt(12.);
    double clusterWidth  = strips*m_tkrGeo->siStripPitch();
    double projectedWidth = fabs(slope)*m_tkrGeo->siThickness();
    int    errorType = m_control->getErrorType();
    int    nStrips = (int) strips+.01;  // just to be safe

    if(errorType==0) {

        error = (clusterWidth - projectedWidth)*oneOverSqrtTwelve;
        error = std::max(error, minErr);

    }  else if (errorType==1) {

        // for the error types below, the error depends on the number of strips in the cluster.
        // For each number of strips, the error is triangular in slope, peaking roughly at
        // nStrips - 1. The base of the triangle is about 2 strips wide.
        // Outside of these limits, there is a minimum error of stripWidth/sqrt(24),
        // except for wide strips, where the error is increased by the excess width.

        if (nStrips==1) {
            if (absSlope<1.4) factor = 1 - 0.5*absSlope;
        } else if (nStrips==2) {
            if (absSlope>.4 && absSlope<2.6) factor = 0.9 - 0.6*fabs(absSlope-1.6);
        } else if (nStrips<11) {
            if (fabs(absSlope-(2.7+1.05*(nStrips-3)))<1.) 
                factor = 0.9 - 0.6*fabs(absSlope - (2.7 + 1.05*(nStrips-3)));
        }
        if (factor==0) {
            double delta = clusterWidth - projectedWidth - 0.5;
            factor = std::max(fabs(delta), 1.);
        } 
        error = factor*minErr;
    
    } else {

        double eps0 = -0.1; // use this to extent or restrict the valid range for 1-strip clusters
        double eps1 = -0.1; // ditto for the rest of the clusters
        double loSlope, hiSlope, peakSlope;
        double loPar1, hiPar1, peakDev;

        if (nStrips==1) {
            if (absSlope<1.5+eps0) factor = 1 - 0.52*absSlope;
        } else if (nStrips<11) {
            if (nStrips==2) {
                loSlope = .4 ; hiSlope = 2.5; peakSlope = 1.61;
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
        if (factor==0) {
            double delta = clusterWidth - projectedWidth - 0.5;
            factor = std::max(fabs(delta), 1.);
        } 

        error = factor*minErr;
    }

    return error;
}
