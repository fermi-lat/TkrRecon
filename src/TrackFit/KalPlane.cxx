
// $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalPlane.cxx,v 1.6 2002/03/29 02:09:43 lsrea Exp $

//----------------------------------------------------------------------
//    
//    Implementation of the Kalman Filter Functions 
//               KalPlane
//
//      Original due to Jose Hernando-Angel circa 1997-1999
//      Re-written to combine both X and Y projections (2001) 
//      
//-----------------------------------------------------------------------

#include "TkrRecon/TrackFit/KalPlane.h"
#include "TkrRecon/GaudiAlg/TkrReconAlg.h"
#include "GlastSvc/Reco/IKalmanParticle.h"

void KalPlane::removeHit()
{
	m_IDHit = 0;
	m_IDTower = 0;
	KalPar pnull(0.,0., 0.,0.);
	KalMatrix covnull(HepMatrix(4,4,0));
        KalHit temp(KalHit::MEAS,pnull,covnull);
	setHit(temp);
	clean();
}

void KalPlane::clean()
{
	KalPar pnull(0.,0.,0.,0.);
	KalMatrix covnull(HepMatrix(4,4,0));
	setHit(KalHit(KalHit::PRED,pnull,covnull));
	setHit(KalHit(KalHit::FIT,pnull,covnull));
	setHit(KalHit(KalHit::SMOOTH,pnull,covnull));
}

void KalPlane::clear()
{
	KalPar pnull(0.,0.,0.,0.);
	KalMatrix covnull(HepMatrix(4,4,0));

	m_eneplane = 0.;
	m_IDHit = 0xffffffff;
	m_IDPlane  = -1;
	m_IDTower = -1;	
	m_zplane = 0.;

	setHit(KalHit(KalHit::MEAS,pnull,covnull));
	setHit(KalHit(KalHit::PRED,pnull,covnull));
	setHit(KalHit(KalHit::FIT,pnull,covnull));
	setHit(KalHit(KalHit::SMOOTH,pnull,covnull));
}

void KalPlane::setHit(const KalHit& hit)
{
    KalHit::TYPE type;
    switch (type=hit.getType()){
    case KalHit::PRED:
	m_hitpred=hit;
	break;
    case KalHit::MEAS:
	m_hitmeas=hit;
	break;
    case KalHit::FIT:
	m_hitfit=hit;
	break;
    case KalHit::SMOOTH:
	m_hitsmooth=hit;
	break;
    case KalHit::UNKNOWN:
        break;
    }   
}

KalHit KalPlane::getHit(KalHit::TYPE type) const
{  
    switch (type){
    case KalHit::PRED:
	return KalHit(m_hitpred);
    case KalHit::MEAS:
	return KalHit(m_hitmeas);
    case KalHit::FIT:
	return KalHit(m_hitfit);
    case KalHit::SMOOTH:
	return KalHit(m_hitsmooth);
    case KalHit::UNKNOWN:
        break;
    } 
    return KalHit();
}

Point KalPlane::getPoint(KalHit::TYPE type)const
{
    KalHit hit=getHit(type);
    return Point(hit.getPar().getXPosition(),
                 hit.getPar().getYPosition(),getZPlane());
}

double KalPlane::getDeltaChiEne(KalHit::TYPE type)const
{
    KalHit hit=getHit(type);
    double delparX=m_hitmeas.getPar().getXPosition()-hit.getPar().getXPosition();
    double delparY=m_hitmeas.getPar().getYPosition()-hit.getPar().getYPosition();
    double sigma2X=m_hitmeas.getCov().getcovX0X0();
    double sigma2Y=m_hitmeas.getCov().getcovY0Y0();
    
    double variance=(delparX*delparX)/sigma2X + (delparY*delparY)/sigma2Y;
    return variance;
}

void KalPlane::setDeltaEne(double ene)

{       
    double radlen = getRadLen();
    double factor = exp(-1.*radlen);

    setEnergy(ene*factor);
}

double KalPlane::getSigma(KalHit::TYPE type) const
{
    double sigma = 1e6;
    KalHit hit=getHit(type);
    KalHit hitmeas = getHit(KalHit::MEAS);
    double delX=hit.getPar().getXPosition()-hitmeas.getPar().getXPosition();
    double delY=hit.getPar().getYPosition()-hitmeas.getPar().getYPosition();
    double sigma2X=hit.getCov().getcovX0X0();
    double sigma2Y=hit.getCov().getcovY0Y0();
    
    sigma=(delX*delX)/sigma2X + (delY*delY)/sigma2Y;
    return sigma;
}

double KalPlane::getDeltaChiSq(KalHit::TYPE type) const
{  
    KalHit hit=getHit(type);
    double delparX=m_hitmeas.getPar().getXPosition()-hit.getPar().getXPosition();
    double delparY=m_hitmeas.getPar().getYPosition()-hit.getPar().getYPosition();
    double sigma2X=m_hitmeas.getCov().getcovX0X0();
    double sigma2Y=m_hitmeas.getCov().getcovY0Y0();
    
    double chi2=(delparX*delparX)/sigma2X + (delparY*delparY)/sigma2Y;
    return chi2;
}

//-------------------------------------
//  Kalman Functions
//-------------------------------------

KalHit KalPlane::predicted(KalHit::TYPE typ, int &nlayers, int klayer, double &zend,
                           double &arc_min)
{
    // Extrapolate the hit to the next SSD layer - maybe an X or a Y
    // Returns (nlayers)  the number of GLAST Planes crossed 
    // Note: nlayers = 0 for crossing between X and Y in same plane

    KalHit hit=getHit(typ);
    
    KalPar pp=hit.getPar();
    KalMatrix Ck=hit.getCov();
    
    double ene =getEnergy();
    double x_slope = pp.getXSlope();   
    double y_slope = pp.getYSlope(); 
    Vector dir_ini = Vector(-x_slope, -y_slope, -1.).unit();

    double x0 = pp.getXPosition();
    double y0 = pp.getYPosition();
    double z0 = getZPlane();
    Point   x_ini(x0,y0,z0);

    int nsteps=klayer-getIDPlane();
    double down = -1.;
    if (nsteps <0 ) down = +1.; // going up;

    double arc_len = nlayers * 32.6/fabs(dir_ini.z()); //mm
    //std::auto_ptr<IKalmanParticle> 
    //        kalPart(TkrReconAlg::m_gismoSvc->kalmanParticle(x_ini, dir_ini, arc_min));
    IKalmanParticle* kalPart = TkrReconAlg::m_KalParticle;
    kalPart->setStepStart(x_ini, dir_ini, arc_min);
    if(kalPart->trackToNextPlane()) {
        if(kalPart->isXPlane()) m_projPlus = TkrCluster::X; 
        else                    m_projPlus = TkrCluster::Y;
        arc_len = kalPart->arcLength();
        nlayers = arc_len*fabs(dir_ini.z())/29.0; //mm
    }
    else {
        nlayers = -1;
        return KalHit(); 
    }   
    zend = kalPart->position().z(); 
    arc_min = arc_len;
    double relDeltaZ = down*fabs(arc_len*dir_ini.z());

    KalMatrix F(relDeltaZ);
	                      
    KalMatrix Q = kalPart->mScat_Covr(ene, arc_len); 
    pp = F*pp;
    if (down == -1.)     Ck=(F*(Ck*F.T()))+Q;
    else if (down == +1) Ck=(F*(Ck+Q)*F.T());
    
    KalHit hitpred(KalHit::PRED, pp, Ck);
    
    return hitpred;
    
} 

KalHit KalPlane::predicted(KalHit::TYPE typ, int nsteps)
{
    // Extrapolate the hit to the plane (delta)
    // Note that GLAST trajectory evolutes along negative z values.
    
    KalHit hit=getHit(typ);
    
    KalPar pp=hit.getPar();
    KalMatrix Ck=hit.getCov();
    
    double ene =getEnergy();
    double x_slope = pp.getXSlope(); 
    double y_slope = pp.getYSlope(); 
    Vector dir_ini = Vector(-x_slope, -y_slope, -1.).unit();
    //Vector dir_ini(-x_slope, -y_slope, -1.);
    //dir_ini.unit();
    double x0 = pp.getXPosition();
    double y0 = pp.getYPosition();
    double z0 = getZPlane();
    Point   x_ini(x0,y0,z0); 

    int iplane=getIDPlane();
    double down = -1.;
    if (nsteps <0 ) {
        down = +1.; // going up;
        dir_ini = -dir_ini;
    }
    
    double deltaZ=down*fabs(nsteps)*32.575; //mm and also Very Bad ... need to re-engineer
    
    double arc_len = fabs(deltaZ/dir_ini.z()); 
    //std::auto_ptr<IKalmanParticle> 
    //        kalPart(TkrReconAlg::m_gismoSvc->kalmanParticle(x_ini, dir_ini, arc_len));
    IKalmanParticle* kalPart = TkrReconAlg::m_KalParticle;
    kalPart->setStepStart(x_ini, dir_ini, arc_len);
    double relDeltaZ = down * deltaZ;

    KalMatrix F(relDeltaZ);
	                      
    KalMatrix Q = kalPart->mScat_Covr(ene, arc_len); 
    pp = F*pp;
    if (down == -1.)     Ck=(F*(Ck*F.T()))+Q;
    else if (down == +1) Ck=(F*(Ck+Q)*F.T());
    
    KalHit hitpred(KalHit::PRED, pp, Ck);
    
    return hitpred;
   
} 

KalHit KalPlane::predicted(KalPlane& kplaneNext) 
{  
    // Extrapolate the hit to the next Kplane; 
    // Note that GLAST trajectory is evolving in negative z coordinates.
    // The energy is necessary to compute the MS error.
    
    KalHit hit=getHit(KalHit::FIT);
    
    KalPar pp=hit.getPar();
    KalMatrix Ck=hit.getCov();
    
    double ene =getEnergy();
    double x_slope = pp.getXSlope(); 
    double y_slope = pp.getYSlope(); 
    Vector dir_ini = Vector(-x_slope, -y_slope, -1.).unit();
    //Vector dir_ini(-x_slope, -y_slope, -1.);
    //dir_ini.unit();
    double x0 = pp.getXPosition();
    double y0 = pp.getYPosition();
    double z0 = getZPlane();
    Point   x_ini(x0,y0,z0); 

    int iplane=getIDPlane();
    int nsteps=kplaneNext.getIDPlane()-getIDPlane();;
    double down = -1.;
    if (nsteps <0 ) {
        down = +1.; // going up;
        dir_ini = -dir_ini;
    }
    

    double deltaZ=kplaneNext.getZPlane()-getZPlane();
    
    double arc_len = fabs(deltaZ/dir_ini.z()); 
    //std::auto_ptr<IKalmanParticle> 
    //        kalPart(TkrReconAlg::m_gismoSvc->kalmanParticle(x_ini, dir_ini, arc_len));
    IKalmanParticle* kalPart = TkrReconAlg::m_KalParticle;
    kalPart->setStepStart(x_ini, dir_ini, arc_len);
    double relDeltaZ = down * fabs(deltaZ);

    KalMatrix F(relDeltaZ);
	                      
    KalMatrix Q = kalPart->mScat_Covr(ene, arc_len); 
    pp = F*pp;
    if (down == -1.)     Ck=(F*(Ck*F.T()))+Q;
    else if (down == +1) Ck=(F*(Ck+Q)*F.T());
    
    // store the matrix with the material
    kplaneNext.m_Qmaterial = Q;

    KalHit hitpred(KalHit::PRED, pp, Ck);
    
    return hitpred;
    
} 

KalHit KalPlane::filter()
{
    // Weight the hits with their cov matricex
    
    KalHit hit1=getHit(KalHit::MEAS);
    KalMatrix H(1);
    double sigXX = hit1.getCov().getcovX0X0();
    double sigYY = hit1.getCov().getcovY0Y0();

    KalMatrix G(1);
    G(1,1) = 1./sigXX;
    G(3,3) = 1./sigYY;
    KalPar pmeas=hit1.getPar();
    
    KalHit hit2=getHit(KalHit::PRED);
    KalMatrix Ckpred=hit2.getCov();
    int i_error;
    KalMatrix Ckpredinv = Ckpred;
    Ckpredinv.invert(i_error); 
    KalPar ppred=hit2.getPar();
    

    KalMatrix Ck=(Ckpredinv+H*(G*H));  
    Ck.invert(i_error);
//    KalPar pk=((Ck*Ckpredinv)*ppred)+((Ck*G)*pmeas);
    KalPar pk=((Ck*Ckpredinv)*ppred)+((Ck*(H*G))*pmeas);

    KalHit hitfit(KalHit::FIT, pk, Ck);
    
    return hitfit;
    
}

KalHit KalPlane::smoother(const KalPlane& kplaneLast)
{
    // Apply the smoother to a Fited hit!
    KalHit hitf0=getHit(KalHit::FIT);
    
    KalHit hitpred=kplaneLast.getHit(KalHit::PRED);
    KalHit hitsm=kplaneLast.getHit(KalHit::SMOOTH);
    
    // double distance=abs(kplaneLast.getZPlane()-getZPlane());

    // double distance=isteps*GFtutor::traySpacing();
    // KalMatrix F(1.,-1.*abs(distance),0.,1.);
 
    double distance=+kplaneLast.getZPlane()-getZPlane();
    // double distance=isteps*GFtutor::traySpacing();
    KalMatrix F(distance);

    KalMatrix Ck=hitf0.getCov();

    KalMatrix Ck1pred=hitpred.getCov();
    KalMatrix Ck1sm=hitsm.getCov();
    
    KalPar pk=hitf0.getPar();
    KalPar pk1pred=hitpred.getPar();
    KalPar pk1sm=hitsm.getPar();
 
    int i_error; 
    KalMatrix Ck1predinv = Ck1pred;
    Ck1predinv.invert(i_error); 
    KalMatrix A=Ck*((F.T())*Ck1predinv);

    KalPar psm=pk+A*(pk1sm-pk1pred);
    KalMatrix Csm=Ck+A*((Ck1sm-Ck1pred)*A.T());
    
    KalHit newhitsm(KalHit::SMOOTH,psm,Csm);
    
    return newhitsm;
}

//-------------------------------------------
