//----------------------------------------------------------------------
//    
//    Implementation of the Kalman Filter Functions 
//               TkrFitPlane
//
//      Original due to Jose Hernando-Angel circa 1997-1999
//      Re-written to combine both X and Y projections (2001) 
//      
//-----------------------------------------------------------------------

#include "TkrRecon/TrackFit/TkrFitPlane.h"
#include "GlastSvc/Reco/IKalmanParticle.h"

void TkrFitPlane::removeHit()
{
	m_IDHit   = 0;
	m_IDTower = 0;
	TkrFitPar pnull(0.,0., 0.,0.);
	TkrFitMatrix covnull(HepMatrix(4,4,0));
        TkrFitHit temp(TkrFitHit::MEAS,pnull,covnull);
	setHit(temp);
	clean();
}

void TkrFitPlane::clean()
{
	TkrFitPar pnull(0.,0.,0.,0.);
	TkrFitMatrix covnull(HepMatrix(4,4,0));
	setHit(TkrFitHit(TkrFitHit::PRED,pnull,covnull));
	setHit(TkrFitHit(TkrFitHit::FIT,pnull,covnull));
	setHit(TkrFitHit(TkrFitHit::SMOOTH,pnull,covnull));
}

void TkrFitPlane::clear()
{
	TkrFitPar pnull(0.,0.,0.,0.);
	TkrFitMatrix covnull(HepMatrix(4,4,0));

	m_eneplane = 0.;
	m_IDHit = 0xffffffff;
	m_IDPlane  = -1;
	m_IDTower = -1;	
	m_zplane = 0.;

	setHit(TkrFitHit(TkrFitHit::MEAS,pnull,covnull));
	setHit(TkrFitHit(TkrFitHit::PRED,pnull,covnull));
	setHit(TkrFitHit(TkrFitHit::FIT,pnull,covnull));
	setHit(TkrFitHit(TkrFitHit::SMOOTH,pnull,covnull));
}

void TkrFitPlane::setHit(const TkrFitHit& hit)
{
    TkrFitHit::TYPE type;
    switch (type=hit.getType()){
    case TkrFitHit::PRED:
	m_hitpred=hit;
	break;
    case TkrFitHit::MEAS:
	m_hitmeas=hit;
	break;
    case TkrFitHit::FIT:
	m_hitfit=hit;
	break;
    case TkrFitHit::SMOOTH:
	m_hitsmooth=hit;
	break;
    case TkrFitHit::UNKNOWN:
        break;
    }   
}

TkrFitHit TkrFitPlane::getHit(TkrFitHit::TYPE type) const
{  
    switch (type){
    case TkrFitHit::PRED:
	return TkrFitHit(m_hitpred);
    case TkrFitHit::MEAS:
	return TkrFitHit(m_hitmeas);
    case TkrFitHit::FIT:
	return TkrFitHit(m_hitfit);
    case TkrFitHit::SMOOTH:
	return TkrFitHit(m_hitsmooth);
    case TkrFitHit::UNKNOWN:
        break;
    } 
    return TkrFitHit();
}

Point TkrFitPlane::getPoint(TkrFitHit::TYPE type)const
{
    TkrFitHit hit=getHit(type);
    return Point(hit.getPar().getXPosition(),
                 hit.getPar().getYPosition(),getZPlane());
}

double TkrFitPlane::getDeltaChiEne(TkrFitHit::TYPE type)const
{
    TkrFitHit hit=getHit(type);
    double delparX=m_hitmeas.getPar().getXPosition()-hit.getPar().getXPosition();
    double delparY=m_hitmeas.getPar().getYPosition()-hit.getPar().getYPosition();
    double sigma2X=m_hitmeas.getCov().getcovX0X0();
    double sigma2Y=m_hitmeas.getCov().getcovY0Y0();
    
    double variance=(delparX*delparX)/sigma2X + (delparY*delparY)/sigma2Y;
    return variance;
}

void TkrFitPlane::setDeltaEne(double ene)

{       
    double radlen = getRadLen();
    double factor = exp(-1.*radlen);

    setEnergy(ene*factor);
}

double TkrFitPlane::getSigma(TkrFitHit::TYPE type) const
{
    double sigma = 1e6;
    TkrFitHit hit=getHit(type);
    TkrFitHit hitmeas = getHit(TkrFitHit::MEAS);
    double delX=hit.getPar().getXPosition()-hitmeas.getPar().getXPosition();
    double delY=hit.getPar().getYPosition()-hitmeas.getPar().getYPosition();
    double sigma2X=hit.getCov().getcovX0X0();
    double sigma2Y=hit.getCov().getcovY0Y0();
    
    sigma=(delX*delX)/sigma2X + (delY*delY)/sigma2Y;
    return sigma;
}

double TkrFitPlane::getDeltaChiSq(TkrFitHit::TYPE type) const
{  
    TkrFitHit hit=getHit(type);
    double delparX=m_hitmeas.getPar().getXPosition()-hit.getPar().getXPosition();
    double delparY=m_hitmeas.getPar().getYPosition()-hit.getPar().getYPosition();
    double sigma2X=m_hitmeas.getCov().getcovX0X0();
    double sigma2Y=m_hitmeas.getCov().getcovY0Y0();
    
    double chi2=(delparX*delparX)/sigma2X + (delparY*delparY)/sigma2Y;
    return chi2;
}
