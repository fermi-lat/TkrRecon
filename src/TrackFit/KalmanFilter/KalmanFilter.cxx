
//----------------------------------------------------------------------
//    
//    Implementation of the Kalman Filter Functions 
//               TkrFitPlane
//
//      Original due to Jose Hernando-Angel circa 1997-1999
//      Re-written to combine both X and Y projections (2001) 
//      
//-----------------------------------------------------------------------

#include "KalmanFilter.h"
#include "TkrRecon/GaudiAlg/TkrReconAlg.h"
#include "GlastSvc/Reco/IKalmanParticle.h"


//-------------------------------------
//  Kalman Functions
//-------------------------------------

TkrFitHit KalmanFilter::predicted(TkrFitPlane& start, TkrFitHit::TYPE typ, int &nlayers, int klayer, double &zend,
                           double &arc_min)
{
    // Extrapolate the hit to the next SSD layer - maybe an X or a Y
    // Returns (nlayers)  the number of GLAST Planes crossed 
    // Note: nlayers = 0 for crossing between X and Y in same plane

    TkrFitHit    hit     = start.getHit(typ);
    
    TkrFitPar    pp      = hit.getPar();
    TkrFitMatrix Ck      = hit.getCov();
    
    double       ene     = start.getEnergy();
    double       x_slope = pp.getXSlope();   
    double       y_slope = pp.getYSlope(); 
    Vector       dir_ini = Vector(-x_slope, -y_slope, -1.).unit();

    double       x0      = pp.getXPosition();
    double       y0      = pp.getYPosition();
    double       z0      = start.getZPlane();
    Point        x_ini(x0,y0,z0);

    int          nsteps  = klayer - start.getIDPlane();

    double       down    = nsteps < 0 ? +1. : -1.;

    double       arc_len = nlayers * 32.6/fabs(dir_ini.z()); //mm

    IKalmanParticle* TkrFitPart = TkrReconAlg::m_KalParticle;
    TkrFitPart->setStepStart(x_ini, dir_ini, arc_min);
    if(TkrFitPart->trackToNextPlane()) 
    {
        AXIS planeProjection = TkrCluster::Y;
        if(TkrFitPart->isXPlane()) planeProjection = TkrCluster::X; 

        start.setNextProj(planeProjection);

        arc_len = TkrFitPart->arcLength();
        nlayers = arc_len*fabs(dir_ini.z())/29.0; //mm
    }
    else 
    {
        nlayers = -1;
        return TkrFitHit(); 
    }
    
    zend    = TkrFitPart->position().z(); 
    arc_min = arc_len;

    double relDeltaZ = down*fabs(arc_len*dir_ini.z());

    TkrFitMatrix F(relDeltaZ);
	                      
    TkrFitMatrix Q = TkrFitPart->mScat_Covr(ene, arc_len); 
    pp = F*pp;
    if (down == -1.)     Ck=(F*(Ck*F.T()))+Q;
    else if (down == +1) Ck=(F*(Ck+Q)*F.T());
    
    TkrFitHit hitpred(TkrFitHit::PRED, pp, Ck);
    
    return hitpred;
} 

TkrFitHit KalmanFilter::predicted(TkrFitPlane& start, TkrFitHit::TYPE typ, int nsteps)
{
    // Extrapolate the hit to the plane (delta)
    // Note that GLAST trajectory evolutes along negative z values.
    
    TkrFitHit    hit     = start.getHit(typ);
    
    TkrFitPar    pp      = hit.getPar();
    TkrFitMatrix Ck      = hit.getCov();
    
    double       ene     = start.getEnergy();
    double       x_slope = pp.getXSlope(); 
    double       y_slope = pp.getYSlope(); 
    Vector       dir_ini = Vector(-x_slope, -y_slope, -1.).unit();

    double       x0      = pp.getXPosition();
    double       y0      = pp.getYPosition();
    double       z0      = start.getZPlane();
    Point        x_ini(x0,y0,z0); 

    int          iplane  = start.getIDPlane();
    double       down    = -1.;
    if (nsteps <0 ) 
    {
        down = +1.; // going up;
        dir_ini = -dir_ini;
    }
    
    double       deltaZ  = down*fabs(nsteps)*32.575; //mm and also Very Bad ... need to re-engineer
    double       arc_len = fabs(deltaZ/dir_ini.z()); 

    IKalmanParticle* TkrFitPart = TkrReconAlg::m_KalParticle;
    TkrFitPart->setStepStart(x_ini, dir_ini, arc_len);

    double relDeltaZ = down * deltaZ;

    TkrFitMatrix F(relDeltaZ);
	                      
    TkrFitMatrix Q = TkrFitPart->mScat_Covr(ene, arc_len); 
    pp = F*pp;
    if (down == -1.)     Ck=(F*(Ck*F.T()))+Q;
    else if (down == +1) Ck=(F*(Ck+Q)*F.T());
    
    TkrFitHit hitpred(TkrFitHit::PRED, pp, Ck);
    
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

    IKalmanParticle* TkrFitPart = TkrReconAlg::m_KalParticle;
    TkrFitPart->setStepStart(x_ini, dir_ini, arc_len);

    double relDeltaZ = down * fabs(deltaZ);

    TkrFitMatrix F(relDeltaZ);
	                      
    TkrFitMatrix Q = TkrFitPart->mScat_Covr(ene, arc_len); 
    pp = F*pp;
    if (down == -1.)     Ck=(F*(Ck*F.T()))+Q;
    else if (down == +1) Ck=(F*(Ck+Q)*F.T());
    
    // store the matrix with the material
    kplaneNext.setQmaterial(Q);

    TkrFitHit hitpred(TkrFitHit::PRED, pp, Ck);
    
    return hitpred;
    
} 

TkrFitHit KalmanFilter::filter(TkrFitPlane& filterPlane)
{
    // Weight the hits with their cov matricex
    
    TkrFitHit    hit1      = filterPlane.getHit(TkrFitHit::MEAS);
    TkrFitPar    pmeas     = hit1.getPar();
    double       sigXX     = hit1.getCov().getcovX0X0();
    double       sigYY     = hit1.getCov().getcovY0Y0();

    TkrFitMatrix H(1);
    TkrFitMatrix G(1);
    G(1,1) = 1./sigXX;
    G(3,3) = 1./sigYY;
    
    TkrFitHit    hit2      = filterPlane.getHit(TkrFitHit::PRED);
    TkrFitMatrix Ckpred    = hit2.getCov();
    int          i_error;
    TkrFitMatrix Ckpredinv = Ckpred;
    Ckpredinv.invert(i_error); 
    TkrFitPar ppred=hit2.getPar();
    

    TkrFitMatrix Ck        = (Ckpredinv+H*(G*H));  
    Ck.invert(i_error);
//    TkrFitPar pk=((Ck*Ckpredinv)*ppred)+((Ck*G)*pmeas);
    TkrFitPar pk=((Ck*Ckpredinv)*ppred)+((Ck*(H*G))*pmeas);

    TkrFitHit hitfit(TkrFitHit::FIT, pk, Ck);
    
    return hitfit;
    
}

TkrFitHit KalmanFilter::smoother(TkrFitPlane& start, const TkrFitPlane& kplaneLast)
{
    // Apply the smoother to a Fitted hit!
    TkrFitHit    hitf0    = start.getHit(TkrFitHit::FIT);
    
    TkrFitHit    hitpred  = kplaneLast.getHit(TkrFitHit::PRED);
    TkrFitHit    hitsm    = kplaneLast.getHit(TkrFitHit::SMOOTH);
    
    // double distance=abs(kplaneLast.getZPlane()-getZPlane());

    // double distance=isteps*GFtutor::traySpacing();
    // TkrFitMatrix F(1.,-1.*abs(distance),0.,1.);
 
    double       distance = +kplaneLast.getZPlane()-start.getZPlane();
    // double distance=isteps*GFtutor::traySpacing();
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
    
    TkrFitHit newhitsm(TkrFitHit::SMOOTH,psm,Csm);
    
    return newhitsm;
}

//-------------------------------------------
