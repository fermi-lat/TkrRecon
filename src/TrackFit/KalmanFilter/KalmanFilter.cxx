
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

KalmanFilter::KalmanFilter(TkrClusterCol* clusters, ITkrGeometrySvc* geo)
{
    m_clusters   = clusters;
    m_tkrGeo     = geo;

    m_radLength  = 0.;
    m_activeDist = 0;
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
        G(1,1) = 1./meas_hit.getCov().getcovX0X0();
        G(3,3) = 0.;
    }
    else {
        G(1,1) = 0.;
        G(3,3) = 1./meas_hit.getCov().getcovY0Y0();
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

    // The following sets the error to the slop between the track 
    // and the cluster over sqrt(12). It protects against getting
    // too small.

    TkrFitMatrix newCov(1);
    double min_err = m_tkrGeo->siResolution();   

    if(plane.getProjection()==TkrCluster::X) {
        double size_Cls = m_clusters->size(id_Cls);
        double x_slope = pred_pars.getXSlope();
        double wid_proj = fabs(x_slope*m_tkrGeo->siThickness());
        double wid_cls  = size_Cls*m_tkrGeo->siStripPitch();
        double error    = (wid_cls - wid_proj)/3.4641;
        error = (error > min_err) ? error : min_err; 
        newCov(1,1) = error*error;
        newCov(3,3) = meas_hit.getCov().getcovY0Y0();
    }
    else {
        double size_Cls = m_clusters->size(id_Cls);
        double y_slope = pred_pars.getYSlope();
        double wid_proj = fabs(y_slope*m_tkrGeo->siThickness());
        double wid_cls  = size_Cls*m_tkrGeo->siStripPitch();
        double error    = (wid_cls - wid_proj)/3.4641;
        error = (error > min_err) ? error : min_err; 
        newCov(1,1) = meas_hit.getCov().getcovX0X0();
        newCov(3,3) = error*error;
    }
    TkrFitHit newMeas(TkrFitHit::MEAS, meas_hit.getPar(), newCov);
    plane.setHit(newMeas);   
    return;
    
}
