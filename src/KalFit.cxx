
// $Header: /nfs/slac/g/glast/ground/cvs/tb_recon/src/reconstruction/KalFit.cxx,v 1.9 2000/10/26 17:24:43 burnett Exp $

//----------------------------------------------------------------------
//    
//    Implementation of the Kalman Filter Functions
//               KalTrack - doFit() 
//               KalPlane, KalHit, KalMatrix, KalPar  - related Objects 
//
//               J.A. Hernando, B. Atwood.  Santa Cruz,  7/8/98,
//-----------------------------------------------------------------------

#include "TkrRecon/KalFit.h"
#include "geometry/Ray.h"
//#include "Event/compConf.h"
#include "TkrRecon/GFcontrol.h"
#include <cmath>
// JAH, HA, error compiling with MVC++, 05/22/00
// using std::abs;

bool CONTROL_setDeltaEne = false;
double CONTROL_MaxEne = 2.;
int CONTROL_maxPlanesEne = 8;

//-------------------------------------------------------------------
//  Access Data Functions (keep compatible with LSQFit)
//-------------------------------------------------------------------

//KalTrack::~KalTrack() {}

//#####################################
void KalTrack::clear()
//#####################################
{
	kplanelist.clear();
	ini();
}
/*
//#####################################
void KalTrack::writeOut(std::ostream& out) const
//#####################################
{
	if (numDataPoints() <3 || chiSquare() < 0) return;

	out << " --- KalTrack Information --- " << "\n";
	out << " IniEne " << iniEnergy() << "\n";
	out << " nhits  " << numDataPoints() << "\n";
	out << " ngaps  " << numGaps()       << "\n";
	out << " chiSQFit " << chiSquare() << "\n";
	out << " chiSQSmooth " << chiSquareSmooth()       << "\n";
	out << " chiSQsegment " << chiSquareSegment()       << "\n";

	out << " x      " << position(0.)    << "\n";
	out << " Slope  " << slope()         << "\n";
	out << " errX   " << errorPosition() << "\n";
	out << " errSlp " << errorSlope()    << "\n";
	out << " errSlpV" << errorSlopeAtVertex() << "\n";
	out << " ECene  " << KalEnergy()     << "\n";

	std::cout << " --> KalPlanes:" << "\n";
	for (int i = 0; i < 2; i++) kplanelist[i].writeOut(out);
}
*/
//#########################################
void KalTrack::drawChiSq(gui::DisplayRep& v, SiCluster::view axis, KalHit::TYPE typ)
//#########################################
{
 	
    // Two dimensional Projection Track Drawing
    int nplanes = kplanelist.size();
    for (int iplane=0; iplane<nplanes; iplane++) {
		double x = kplanelist[iplane].getHit(typ).getPar().getPosition();
		double Ox = kplanelist[iplane].getOrthPar().getPosition(); 
		double z0  = kplanelist[iplane].getZPlane()+0.01;
		double x0,y0;
		if (axis == SiCluster::X) {
			x0 = x;
			y0 = Ox;
		} else {
			x0 = Ox;
			y0 = x;
		}

		double delta= kplanelist[iplane].getDeltaChiSq(typ);
		double xl,xr,yl,yr;
		if (axis == SiCluster::X) {
			xl = x0-0.5*delta;
			xr = x0+0.5*delta;
			yl = y0;
			yr = y0;
		} else {
			xl = x0;
			xr = x0;
			yl = y0-0.5*delta;
			yr = y0+0.5*delta;
		}		
    	v.moveTo(Point(xl,yl,z0));
        v.lineTo(Point(xr,yr,z0));
    }
}

//#########################################
void KalTrack::drawTrack(gui::DisplayRep& v, SiCluster::view axis, KalHit::TYPE typ)
//#########################################
{
    // Two dimensional Projection Track Drawing
    //unused:	double XLIMIT = -95.;
    int nplanes = kplanelist.size();
    for (int iplane=0; iplane<nplanes-1; iplane++) {
	double x = kplanelist[iplane].getHit(typ).getPar().getPosition();
	double Ox = kplanelist[iplane].getOrthPar().getPosition(); 
//	Ox = XLIMIT;
	double z0  = kplanelist[iplane].getZPlane()+0.01;
	double x0,y0;
	if (axis == SiCluster::X) {
	    x0 = x;
	    y0 = Ox;
	} else {
	    x0 = Ox;
            y0 = x;
	}

	double tanx, tany;
	double tanX  = kplanelist[iplane].getHit(typ).getPar().getSlope();
		double tanOX = kplanelist[iplane].getOrthPar().getSlope();
//		tanOX = 0.;
		if (axis == SiCluster::X){
			tanx = tanX;
			tany = tanOX;
		} else {
			tanx = tanOX;
			tany = tanX;
		}

		Point origin(x0,y0,z0);
		Vector dir = Vector(-1.*tanx,-1.*tany,-1.);
		double mag = dir.mag();
		dir = dir*(1./mag);
	
		Ray segment(origin,dir);
		double zstep=kplanelist[iplane+1].getZPlane()-z0;
		double cosz=dir.z();
		v.moveTo(segment.position(0.));
        v.lineTo(segment.position(zstep/cosz));
    }
}
//#####################################
int KalTrack::numGaps() const
//#####################################
{
    int numGaps =0;
    if (numDataPoints() == 0) return numGaps;

    numGaps =1+ kplanelist.back().getIDPlane() -
	kplanelist.front().getIDPlane() - numDataPoints();

    return numGaps;
}

//#####################################
int KalTrack::compareFits(KalTrack& ktrack)
//#####################################
{
    int numComData=0;
    if (kplanelist.size()==0||ktrack.kplanelist.size()==0) 
	return numComData;
    
    for (unsigned int i=0;i<kplanelist.size();i++){
	for (unsigned int j=0; j<ktrack.kplanelist.size();j++){
	    if (kplanelist[i].getIDHit()==
		ktrack.kplanelist[j].getIDHit()) numComData++;
	}
    }
    return numComData;
}

//#####################################
Point KalTrack::getHit(unsigned ipos)const
//#####################################
{
    if (ipos<kplanelist.size()){
	return kplanelist[ipos].getPoint(KalHit::MEAS);
    }
    //  else return Point(FLT_MAX,FLT_MAX,FLT_MAX);
    return Point(-999999.,-999999.,-999999.);
}
//#####################################
unsigned KalTrack::getHitIndex(unsigned ipos)const
//#####################################
{
    unsigned index= 0xFFFFFFFF;
    if (ipos<kplanelist.size()) 
	index=kplanelist[ipos].getIDHit();
    return index;
}
/*
//#####################################
void KalTrack::printOn(std::ostream &) const
//#####################################
{}
*/ 
//################################################
double KalTrack::positionAtZ(double const z) const
//################################################ 
{
    if( kplanelist.empty() ) return 99;
    const KalPlane& plane = *kplanelist.begin();
    double  x = plane.getHit(KalHit::SMOOTH).getPar().getPosition(),
		z0  = plane.getZPlane()+0.01,
		tanX  = plane.getHit(KalHit::SMOOTH).getPar().getSlope();
    return x+(z-z0)*tanX;
}

//#####################################
void KalTrack::setIniEnergy(double e) 
//#####################################
{
	m_energy0 = e;
	for (unsigned int iplane = 0; iplane < kplanelist.size(); iplane++)
		kplanelist[iplane].setEnergy(m_energy0);
} 
//-------------------------------------
//   Kalman Track Private
//-------------------------------------

//#####################################
KalTrack::KalTrack()
: m_energy0(0), m_x0(0), m_slopeX(0)
, m_chisq(1e6), m_chisqSmooth(1e6) 
, m_KalEnergy(0), m_KalThetaMS(0), m_rmsResid(0)
, m_numSegmentPoints(0), m_chisqSegment(0)
//#####################################
{
}
//#####################################
double KalTrack::doFit()
//#####################################
{
    ini();
    
    int nplanes=kplanelist.size();
    if (nplanes<=2) return m_chisq;
    
    m_chisq       = 0.;
	m_chisqSmooth = 0.;
    
    // Generate the initial hit to start the Kalman Filter
    //----------------------------------------------------
    KalHit hitf=generateFirstFitHit();
    kplanelist[0].setHit(hitf); 
    
    //  Filter 
    //------------
    int iplane = 0;  // to be compatible with new scoping rules for (MSC_VER)
    for (iplane = 0 ; iplane<nplanes-1;iplane++){
		filterStep(iplane);
        if(iplane > 0) m_chisq+=kplanelist[iplane+1].getDeltaChiSq(KalHit::FIT);
    }
    
    // Smoother
    //---------
    KalHit hitsm=(kplanelist[nplanes-1].getHit(KalHit::FIT)).changeType(KalHit::SMOOTH);
    kplanelist[nplanes-1].setHit(hitsm);
    m_chisqSmooth=kplanelist[nplanes-1].getDeltaChiSq(KalHit::SMOOTH);
    
    for (iplane=nplanes-2; iplane>=0;iplane--){
	KalHit hitsm=kplanelist[iplane].smoother(kplanelist[iplane+1]);
        kplanelist[iplane].setHit(hitsm);
	m_chisqSmooth+=kplanelist[iplane].getDeltaChiSq(KalHit::SMOOTH);                
    }
   
    // End the Calculations
    //---------------------
    finish();

    return m_chisq;
}

//#####################################
void KalTrack::ini()
//#####################################
{
    m_energy0  = 0.;
    m_x0       = 0.;
    m_slopeX   = 0.;
    m_rmsResid = 0.;
    m_KalEnergy = 0.;
    m_chisq    = 1e6;
    m_chisqSmooth = 1e6;
    m_KalThetaMS = 0.;
    m_rmsResid = 0.;
    m_numSegmentPoints = 0;
    m_chisqSegment = 1e6;
    
    unsigned int iplane =0;
    for (;iplane < kplanelist.size();iplane++) {
        kplanelist[iplane].clean();
    }
    
}
//#####################################
void KalTrack::finish()
//#####################################
{
    // Compute the fit variables  
    if (m_chisq>=0){
		int nplanes = kplanelist.size();
        m_x0=kplanelist[0].getHit(KalHit::SMOOTH).getPar().getPosition();
        m_slopeX=kplanelist[0].getHit(KalHit::SMOOTH).getPar().getSlope();
        m_rmsResid=-1.;
        m_chisq=m_chisq/(1.*nplanes-2.);
		m_chisqSmooth/=(1.*nplanes-2.);
        m_rmsResid=0.;
		int iplane = 0;
        for (iplane=0;iplane<nplanes;iplane++){
            double res=kplanelist[iplane].
                getHit(KalHit::SMOOTH).getPar().getPosition()-
                kplanelist[iplane].
                getHit(KalHit::MEAS).getPar().getPosition();
            m_rmsResid+=res*res;
        }
        m_rmsResid=sqrt(m_rmsResid/(1.*nplanes));
    }

    //   Energy calculations
    eneDetermination();

    // Segment Calculation
    if (m_chisq>=0){
	m_numSegmentPoints = computeNumSegmentPoints();
	m_chisqSegment = computeChiSqSegment(m_numSegmentPoints);
    }	
}

//#######################################################
void KalTrack::filterStep(int iplane) 
//#######################################################
{
    KalHit hitp = kplanelist[iplane].predicted(kplanelist[iplane+1]);
    kplanelist[iplane+1].setHit(hitp);
    KalHit hitf1 = kplanelist[iplane+1].filter();
    kplanelist[iplane+1].setHit(hitf1);

    if (CONTROL_setDeltaEne) 
        kplanelist[iplane+1].setDeltaEne(kplanelist[iplane].getEnergy());
        
}

//##########################################
KalHit KalTrack::generateFirstFitHit()
//##########################################
{
    
    int nplanes=kplanelist.size();
    if (nplanes<2) 
	std::cout << "ERROR - Kaltrack::generateFirstFitHit" << '\n';
    
    double slopeguess=(kplanelist[1].getHit(KalHit::MEAS).getPar().getPosition()-
	kplanelist[0].getHit(KalHit::MEAS).getPar().getPosition())/
	(kplanelist[1].getZPlane()-kplanelist[0].getZPlane());
    double x0guess=kplanelist[0].getHit(KalHit::MEAS).getPar().getPosition();
    
    
    double orth_slope = kplanelist[1].getOrthPar().getSlope();
    double cosZ = sqrt(1./(1. + slopeguess*slopeguess + orth_slope*orth_slope));
    double energy = kplanelist[1].getEnergy();
    if (energy == 0.) energy = m_energy0;
    // if (energy < 0.01) energy =0.01;
    double theta0=KalPlane::theta0ms(energy,
	cosZ,KalPlane::radLen(kplanelist[0].getIDPlane()));
    //  The first error is arbitrary to a degree: choose 10*(ms + 1-plane hit precision)
    
    double thetaerror=10.*(theta0*(1+slopeguess*slopeguess) + .0025);
    //	if (thetaerror>3.1416) thetaerror=3.1416;
    thetaerror*=thetaerror;
    KalMatrix m(0.,0.,0.,thetaerror); 
    //KalMatrix m(0.,0.,0.,0.01);
    
    KalPar parguess(x0guess,slopeguess);
    KalHit hitf(KalHit::FIT,parguess,
	(kplanelist[0].getHit(KalHit::MEAS)).getCov()+m);
    
    return hitf;
}

//##########################################
void KalTrack::eneDetermination()
//##########################################
{
    int nplanes = numDataPoints();
    int iplane = 2;
    double totalRad = 0.;
    double eneSum = 0.;
    double thetaSum = 0.;
    for (iplane = 2; iplane < nplanes; iplane++) {
        int idplane = kplanelist[iplane].getIDPlane();
        double chie = kplanelist[iplane].getDeltaChiEne(KalHit::PRED);
        double eta = ((chie-1.)*(chie-1.)*(chie-1.))/((chie+2.)*(chie+2.));
        eta = sqrt(fabs(eta));

        double slope = kplanelist[iplane].getHit(KalHit::SMOOTH).getPar().getSlope();
        double slope_orth = kplanelist[iplane].getOrthPar().getSlope();
        double cosZ= sqrt(1/(1.+ slope*slope+slope_orth*slope_orth));
        double cosX= sqrt(1/(1.+ slope*slope));
        
        double radlen = KalPlane::radLen(idplane)/cosZ;
        totalRad += radlen;
        double factor = 1./(2.-exp(-1.*totalRad));
        // double factor = 1./(1.+totalRad);

        double sigma = sqrt(kplanelist[iplane].getHit(KalHit::MEAS).getCov().getsiga());
        double distance = abs(kplanelist[iplane].getZPlane()-kplanelist[iplane-1].getZPlane());

        // factor - energy loss factor
        // sigma/distance - dimensions
        // eta - undimesnionless parameters
        // cosX, etx - geometrical factors
        double theta = factor*(sigma/distance)*eta*(cosX*cosZ*sqrt(cosZ));
        theta /= (sqrt(radlen)*(1+0.038*log(radlen)));
        thetaSum += theta;

    }
    m_KalThetaMS = thetaSum/(nplanes-2.);
    m_KalEnergy = 0.0136/m_KalThetaMS;
    double radlen = KalPlane::radLen(kplanelist[0].getIDPlane());
    m_KalThetaMS *= sqrt(radlen)*(1+0.038*log(radlen));
}

//##########################################
double KalTrack::errorPosition() const
//##########################################
{
	double error = 0.;
	if (m_chisqSmooth == 0) return error;
	error = sqrt(kplanelist[0].getHit(KalHit::SMOOTH).getCov().getsiga());
	return error;
}
//##########################################
double KalTrack::errorSlope() const
//##########################################
{
	double error = 0.;
	if (m_chisqSmooth == 0) return error;
	error = sqrt(kplanelist[0].getHit(KalHit::SMOOTH).getCov().getsigb());
	return error;
}
//##########################################
double KalTrack::errorSlopeAtVertex() const
//##########################################
{
	double error = 0.;
	if (m_chisqSmooth == 0) return error;

	double slope     = m_slopeX;
	double ene       = kplanelist[0].getEnergy();
	KalMatrix Cov    = kplanelist[0].getHit(KalHit::SMOOTH).getCov();
	double orthSlope = kplanelist[0].getOrthPar().getSlope();
	int iplane       = kplanelist[0].getIDPlane();
	double radlen0   = KalPlane::radLen(iplane);
	KalMatrix Q=KalPlane::Q(ene,slope,orthSlope,0.5*radlen0);

	Cov = Cov + Q;

	error = sqrt(Cov.getsigb());
	return error;
}

//##########################################
int KalTrack::computeNumSegmentPoints(KalHit::TYPE typ)
//##########################################
{
	
    unsigned int num_ist=0;
    // OSF alpha detected m_energy0 == -NaN here. 
    // Paul_Kunz@slac.stanford.edu
#ifndef _MSC_VER
    if ( !isnan(m_energy0) ) { 			
#endif
        // potential square root of negative number here.
        // tburnett@u.washington.edu
        num_ist = m_energy0>0? static_cast<unsigned int>( 4.*sqrt(m_energy0)): 1000;
#ifndef _MSC_VER
    } else {
      num_ist = 1000;
    }
#endif
    if (num_ist <= 2) num_ist = 3;
    if (num_ist > kplanelist.size()) num_ist = kplanelist.size(); 
			
    return num_ist; 
	/*int npoints = 1;
	int iplane = 1;
	for (iplane = 1; iplane< kplanelist.size(); iplane++) {
		double sigb0 = sqrt(kplanelist[iplane-1].getHit(typ).getCov().getsigb());
		double sigb1 = sqrt(kplanelist[iplane].getHit(typ).getCov().getsigb());
		double delta = (sigb0-sigb1)/sigb0;
		if (abs(delta)>0.15) npoints++;
		else break;
	}
	if (npoints<3) npoints = 3;
	return npoints; */

}
//##########################################
double KalTrack::computeChiSqSegment(int nhits, KalHit::TYPE typ)
//##########################################
{
	double chi2 = 0;
    int ihit =0;
	for (ihit =0; ihit < nhits; ihit++) {
		chi2 += kplanelist[ihit].getDeltaChiSq(typ);
	}
	chi2 /= (nhits-2.);
	return chi2;
}
//##########################################
double KalTrack::kink(int iplane) const
//##########################################
{
	double kink = 0.;
	if (iplane<0 || iplane == numDataPoints()) return kink;

	double slope0 = kplanelist[iplane].getHit(KalHit::SMOOTH).getPar().getSlope();
	double slope1 = kplanelist[iplane+1].getHit(KalHit::SMOOTH).getPar().getSlope();

	kink = slope1-slope0;

	return kink;
}
//##########################################
double KalTrack::kinkNorma(int iplane) const
//##########################################
{
	double k = 0.;
	k = kink(iplane);

	if (iplane <0 || iplane >= numDataPoints()) return k;

	double errorSQ = kplanelist[iplane+1].getHit(KalHit::SMOOTH).getCov().getsigb();
	errorSQ += kplanelist[iplane+1].getQmaterial().getsigb();

	double error = sqrt(errorSQ);

	double kinkNorma =  0.;
	if (error != 0.) kinkNorma = fabs(k)/error;

	return kinkNorma;

}
//-------------------------------------
//   Kalman Plane
//-------------------------------------

//#####################################
void KalPlane::writeOut(std::ostream& out) const
//#####################################
{
	out << " --- KalPlane Information --- " << "\n";
	out << " IDhit    " << getIDHit() << "\n";
	out << " IDplane  " << getIDPlane() << "\n";
	out << " IDtower  " << getIDTower() << "\n";
	out << " z        " << getZPlane() << "\n";
	out << " xmeas    " << getHit(KalHit::MEAS).getPar().getPosition() << "\n";
	out << " xpred    " << getHit(KalHit::PRED).getPar().getSlope() << "\n";
	out << " errXpred " << sqrt(getHit(KalHit::PRED).getCov().getsiga()) << "\n";
	
	out << " xSm     " << getHit(KalHit::SMOOTH).getPar().getPosition() << "\n";
	out << " SlopeSm " << getHit(KalHit::SMOOTH).getPar().getSlope() << "\n";
	out << " errxSm  " << sqrt(getHit(KalHit::SMOOTH).getCov().getsiga()) << "\n";
	out << " errSlSm " << sqrt(getHit(KalHit::SMOOTH).getCov().getsigb()) << "\n";
	out << " chiSqFit " << getDeltaChiSq(KalHit::FIT) << "\n";
	out << " chiSqSm  " << getDeltaChiSq(KalHit::SMOOTH) << "\n";
	out << " chiSqEne " << getDeltaChiEne(KalHit::PRED)<< "\n";

}
//#####################################
void KalPlane::removeHit()
//#####################################
{
	m_IDHit = 0;
	m_IDTower = 0;
	KalPar pnull(0.,0.);
	KalMatrix covnull(0.,0.,0.);
	setHit(KalHit(KalHit::MEAS,pnull,covnull));
	clean();
}
//#####################################
void KalPlane::clean()
//#####################################
{
	KalPar pnull(0.,0.);
	KalMatrix covnull(0.,0.,0.);
	setHit(KalHit(KalHit::PRED,pnull,covnull));
	setHit(KalHit(KalHit::FIT,pnull,covnull));
	setHit(KalHit(KalHit::SMOOTH,pnull,covnull));
}
//#####################################
void KalPlane::clear()
//#####################################
{
	KalPar pnull(0.,0.);
	KalMatrix covnull(0.,0.,0.);

	m_eneplane = 0.;
	m_IDHit = 0xffffffff;
	m_IDPlane  = -1;
	m_IDTower = -1;
	m_OrthPar = pnull;	
	m_zplane = 0.;

	setHit(KalHit(KalHit::MEAS,pnull,covnull));
	setHit(KalHit(KalHit::PRED,pnull,covnull));
	setHit(KalHit(KalHit::FIT,pnull,covnull));
	setHit(KalHit(KalHit::SMOOTH,pnull,covnull));
}
//#####################################
void KalPlane::setHit(const KalHit& hit)
//#####################################
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

//#####################################
KalHit KalPlane::getHit(KalHit::TYPE type) const
//#####################################
{
    KalHit hit(KalHit::UNKNOWN);
    
    switch (type){
    case KalHit::PRED:
	hit=m_hitpred;
	break;
    case KalHit::MEAS:
	hit=m_hitmeas;
	break;
    case KalHit::FIT:
	hit=m_hitfit;
	break;
    case KalHit::SMOOTH:
	hit=m_hitsmooth;
        break;
    case KalHit::UNKNOWN:
        break;
    } 
    
    return hit;
}

//#####################################
Point KalPlane::getPoint(KalHit::TYPE type)const
//#####################################
{
    KalHit hit=getHit(type);
    return Point(hit.getPar().getPosition(),0.,getZPlane());
}

//#####################################
double KalPlane::getDeltaChiEne(KalHit::TYPE type)const
//#####################################
{
    KalHit hit=getHit(type);
    double delpar=m_hitmeas.getPar().getPosition()-hit.getPar().getPosition();
    double sigma2=m_hitmeas.getCov().getsiga();
    
    double variance=(delpar*delpar)/sigma2;
    return variance;
}

//#####################################
void KalPlane::setDeltaEne(double ene)
//#####################################
{
    double slope = getHit(KalHit::PRED).getPar().getSlope();
    double slope_orth = getOrthPar().getSlope();
    double cosZ= sqrt(1/(1.+ slope*slope+slope_orth*slope_orth));
    double cosX= sqrt(1/(1.+ slope*slope));
        
    double radlen = KalPlane::radLen(getIDPlane())/cosZ;
    double factor = exp(-1.*radlen);

    setEnergy(ene*factor);
}

//#####################################
double KalPlane::getSigma(KalHit::TYPE type) const
//#####################################
{

	double sigma = 1e6;

    KalHit hit = getHit(type);
	KalHit hitmeas = getHit(KalHit::MEAS);
	
	double deltap = hit.getPar().getPosition()-hitmeas.getPar().getPosition();
	double error2 = hit.getCov().getsiga();
	sigma = fabs(deltap)/sqrt(error2);

	return sigma;
}
//#####################################
double KalPlane::getDeltaChiSq(KalHit::TYPE type) const
//#####################################
{
    
    KalHit hit(KalHit::UNKNOWN);

    switch (type){
    case KalHit::FIT:
	hit=m_hitfit;
	break;
    case KalHit::SMOOTH:
        hit=m_hitsmooth;
        break;
    case KalHit::MEAS:
    case KalHit::PRED:
    case KalHit::UNKNOWN:
        break;
    }
    
	/*
    KalMatrix V=hitmeas.getCov();
    KalMatrix H(1.,0.,0.,0.);
    KalMatrix Ck=hit.getCov();
    KalMatrix R=V-(H*(Ck*H));
    double xval=1./R.getsiga();
    
    double pmeas=hitmeas.getPar().getPosition();
    double pk=hit.getPar().getPosition();
    double delpar=pmeas-pk;
    
    double deltachi2=delpar*xval*delpar;
    
    return deltachi2; */

    double delpar =   m_hitmeas.getPar().getPosition() -hit.getPar().getPosition();
    double variance = m_hitmeas.getCov().getsiga() -    hit.getCov().getsiga();

    // prevent division by zero: what should it do then?
	double chi2 = 0.;
	if (variance>0.) chi2 = (delpar*delpar)/variance;
	else if (variance<0.) chi2 = 1e6;

    return chi2;
}

//-------------------------------------
//  Kalman Functions
//-------------------------------------

//#######################################################
KalHit KalPlane::predicted(KalHit::TYPE typ, double zEnd, int klayer)
//#######################################################
{
    // Extrapolate the hit to the plane klayer, with a position zEnd;
    // This functions handels UP and DOWN predictions
    // The normal GLAST trajectories go down (negative z)

    KalHit hit=getHit(typ);
    
    KalPar pp=hit.getPar();
    KalMatrix Ck=hit.getCov();
    
    double ene =getEnergy();
    double slope = pp.getSlope(); 
    double orth_slope = getOrthPar().getSlope();
    
    int nsteps=klayer-getIDPlane();
    double down = -1.;
    if (nsteps <0 ) down = +1.; // going up;
    nsteps = abs(nsteps);
    double deltaZ=zEnd-getZPlane();
    int iplane=getIDPlane();
    double movedDeltaZ=0.;
    
    for (int istep=0; istep< abs(nsteps); istep++){	
	// double relDeltaZ= GFtutor::traySpacing();
	// if (istep==nsteps-1) relDeltaZ=abs(deltaZ)-movedDeltaZ;
	// movedDeltaZ+=relDeltaZ;
	
	// KalMatrix F(1.,-1.*relDeltaZ,0.,1.);

		double relDeltaZ = down * GFtutor::trayGap();
		if (istep == abs(nsteps) -1) relDeltaZ = deltaZ - movedDeltaZ;
		movedDeltaZ+=relDeltaZ;

		KalMatrix F(1.,relDeltaZ,0.,1.);
		KalMatrix Q=KalPlane::Q(ene, slope, orth_slope, KalPlane::radLen(iplane+1));
		pp=F*pp;
		if (down == -1.) Ck=(F*(Ck*F.transpose()))+Q;
		else if (down == +1) Ck=(F*(Ck+Q)*F.transpose());
	}
    
//    if (movedDeltaZ != deltaZ) std::cout << " error KalPlane predicted ";
    
    KalHit hitpred(KalHit::PRED, pp, Ck);
    
    return hitpred;
    
} 
//#######################################################
KalHit KalPlane::predicted(KalHit::TYPE typ, int nsteps)
//#######################################################
{
    // Extrapolate the hit to the plane (delta)
    // Note that GLAST trajectory evolutes along negative z values.
    
    KalHit hit=getHit(typ);
    
    KalPar pp=hit.getPar();
    KalMatrix Ck=hit.getCov();
    
    double ene =getEnergy();
    double slope = pp.getSlope(); 
    double orth_slope = getOrthPar().getSlope();
    
    double down = -1;
    if (nsteps < 0) down = +1;
    double deltaZ=down*fabs(nsteps)*GFtutor::trayGap();
    
    int iplane=getIDPlane();
    double movedDeltaZ=0.;
    
    
    for (int istep=0; istep<abs(nsteps); istep++){

	double relDeltaZ = down*GFtutor::trayGap();
	if (istep == abs(nsteps) -1) relDeltaZ = deltaZ - movedDeltaZ;
	movedDeltaZ+=relDeltaZ;

	KalMatrix F(1.,relDeltaZ,0.,1.);
	KalMatrix Q=KalPlane::Q(ene, slope, orth_slope, KalPlane::radLen(iplane+1));
	pp=F*pp;
	if (down == -1.) Ck=(F*(Ck*F.transpose()))+Q;
	else if (down == +1) Ck=(F*(Ck+Q)*F.transpose());
    }

//    if (movedDeltaZ != deltaZ) std::cout << " error KalPlane predicted ";
    
    KalHit hitpred(KalHit::PRED, pp, Ck);
    
    return hitpred;
   
} 
//#######################################################
KalHit KalPlane::predicted(KalPlane& kplaneNext) 
//#######################################################
{
    
    // Extrapolate the hit to the next Kplane; 
    // Note that GLAST trajectory is evolving in negative z coordinates.
    // The energy is necessary to compute the MS error.
    
    KalHit hit=getHit(KalHit::FIT);
    
    KalPar pp=hit.getPar();
    KalMatrix Ck=hit.getCov();
    
    double ene= getEnergy();
    double slope = pp.getSlope(); 
    double orth_slope = getOrthPar().getSlope();
    
    int nsteps=kplaneNext.getIDPlane()-getIDPlane();

    double down = -1.;
    if (nsteps <0) down = +1.;

    double deltaZ=kplaneNext.getZPlane()-getZPlane();
    double movedDeltaZ =0.;
    int iplane=getIDPlane();

    KalMatrix Qaccumul(0.,0.,0.);

    for (int istep=0; istep<abs(nsteps); istep++){

	double relDeltaZ = down*GFtutor::trayGap();
	if (istep == abs(nsteps) -1) relDeltaZ = deltaZ - movedDeltaZ;
	movedDeltaZ+=relDeltaZ;

	KalMatrix F(1.,relDeltaZ,0.,1.);
	KalMatrix Q=KalPlane::Q(ene, slope, orth_slope, KalPlane::radLen(iplane+1));

	pp=F*pp;
	if (down == -1.) Ck=(F*(Ck*F.transpose()))+Q;
	else if (down == +1) Ck=(F*(Ck+Q)*F.transpose());
    
	if (down == -1.) {
	    Qaccumul = Q + F*(Qaccumul*F.transpose());
	}
            if (down == +1.) Qaccumul = F*((Qaccumul+Q)*F.transpose());  
    }
    // store the matrix with the material
    kplaneNext.m_Qmaterial = Qaccumul;

    KalHit hitpred(KalHit::PRED, pp, Ck);
    
    return hitpred;
    
} 

//#######################################################
KalHit KalPlane::filter()
//#######################################################
{
    // Weight the hits with their cov matricex
    
    KalHit hit1=getHit(KalHit::MEAS);
    KalMatrix H(1.,0.,0.,0.);
    KalMatrix G(1./hit1.getCov().getsiga(),0.,0.);
    KalPar pmeas=hit1.getPar();
    
    KalHit hit2=getHit(KalHit::PRED);
    KalMatrix Ckpred=hit2.getCov();
    KalMatrix Ckpredinv=Ckpred.invert();
    KalPar ppred=hit2.getPar();
    
    KalMatrix Ck=(Ckpredinv+H*(G*H)).invert();   
    KalPar pk=((Ck*Ckpredinv)*ppred)+((Ck*G)*pmeas);
    
    KalHit hitfit(KalHit::FIT, pk, Ck);
    
    return hitfit;
    
}

//#######################################################
KalHit KalPlane::smoother(const KalPlane& kplaneLast)
//#######################################################
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
    KalMatrix F(1.,distance,0.,1.);

    KalMatrix Ck=hitf0.getCov();

    KalMatrix Ck1pred=hitpred.getCov();
    KalMatrix Ck1sm=hitsm.getCov();
    
    KalPar pk=hitf0.getPar();
    KalPar pk1pred=hitpred.getPar();
    KalPar pk1sm=hitsm.getPar();
    
    KalMatrix A=Ck*((F.transpose())*Ck1pred.invert());
    
    KalPar psm=pk+A*(pk1sm-pk1pred);
    KalMatrix Csm=Ck+A*((Ck1sm-Ck1pred)*A.transpose());
    
    KalHit newhitsm(KalHit::SMOOTH,psm,Csm);
    
    return newhitsm;
}

//-------------------------------------
// KalPlane Static Functions
//-------------------------------------
//#####################################
double KalPlane::theta0ms(double ene, double cosZ, double radlen0)
//#####################################
{
    double radlen=radlen0/cosZ;
    double theta0=(0.0136/ene)*sqrt(radlen)*(1+0.038*log(radlen));
    //	double theta0 = radlen0;  // dimensionless case
    return theta0;
}

//#####################################
double KalPlane::radLen(int kplane)
//#####################################
{
    double foamx0=625., siliconx0=9.36;
    
    int ipln = GFtutor::numPlanes() - 1 - kplane; 
    double radlen=(GFtutor::convRadLen(ipln)+GFtutor::trayGap()/foamx0+
	2.*0.02/18.+2.*GFtutor::siThickness()/siliconx0);
    
    //	radlen=GFtutor::convRadLen(); //dimensionless case
    
    return radlen;
}

//#####################################
KalMatrix KalPlane::Q(double ene, double slope, double orth_slope, double radlen0)
//#####################################
{
    // WARNING - this functions uses the angle theta instead of the
    // slope==tan(theta) BUT it returns in Q the error for the slope
    double cosZ = sqrt(1./(1.+ slope*slope + orth_slope*orth_slope));
    double cosX = sqrt(1./(1.+ slope*slope));
    double theta0=KalPlane::theta0ms(ene,cosZ,radlen0);
    // double ocos=(1./(cosX*CosX));
    double ocos=(1./(cosX*cosZ));
    double jacovian=ocos*ocos;
    KalMatrix Q(0.,theta0*theta0*jacovian,0.);
    //  dimensionless case
    //	double theta0=KalPlane::theta0ms(ene,slope,radlen);
    //	KalMatrix Q(0.,theta0*theta0,0.);
    return Q;
}
//-------------------------------------------
//   Kalman Hit 
//-------------------------------------------

//#######################################################
KalHit KalHit::changeType(TYPE typ)
//#######################################################
{
    
    KalHit hit;
    
    hit.type=typ;
    hit.par=par;
    hit.cov=cov;
    
    return hit;
}
//-------------------------------------------
//   Kalman Parameter and Matrix 
//-------------------------------------------

int KalMatrix::kfdim=2;

//#####################################
double operator*(const KalPar& pa, const KalPar& pb)
//#####################################
{
    double po=pa.position*pb.position+pa.slope*pb.slope;	
    return po;
}

//#####################################
KalPar operator+(const KalPar& pa, const KalPar& pb)
//#####################################
{
    KalPar po(pa.position+pb.position,
	pa.slope+pb.slope);	
    return po;
}

//#####################################
KalPar operator-(const KalPar& pa, const KalPar& pb)
//#####################################
{
    KalPar po(pa.position-pb.position,
	pa.slope-pb.slope);	
    return po;
}

//#####################################
KalPar operator*(const KalMatrix& m, const KalPar& p)
//#####################################
{
    KalPar po(m.a[0][0]*p.position+m.a[0][1]*p.slope,
	m.a[1][0]*p.position+m.a[1][1]*p.slope);	
    return po;
}

//#####################################
KalMatrix operator*(const KalMatrix& ma, const KalMatrix& mb)
//#####################################
{
    KalMatrix mc;
    
    for (int i=0;i<KalMatrix::kfdim;i++){
	for (int j=0;j<KalMatrix::kfdim;j++){	
	    for (int k=0;k<KalMatrix::kfdim;k++){
		mc.a[i][j]+=ma.a[i][k]*mb.a[k][j];
	    }
	}
    }
    
    return mc;
}

//#####################################
KalMatrix operator+(const KalMatrix& ma, const KalMatrix& mb)
//#####################################
{
    KalMatrix mc;
    
    for (int i=0;i<KalMatrix::kfdim;i++){
	for (int j=0;j<KalMatrix::kfdim;j++){	
	    mc.a[i][j]=ma.a[i][j]+mb.a[i][j];	
	}
    }
    
    
    return mc;
}

//#####################################
KalMatrix operator-(const KalMatrix& ma, const KalMatrix& mb)
//#####################################
{
    KalMatrix mc;
    
    for (int i=0;i<KalMatrix::kfdim;i++){
	for (int j=0;j<KalMatrix::kfdim;j++){	
	    mc.a[i][j]=ma.a[i][j]-mb.a[i][j];	
	}
    }
    
    return mc;
}

//#####################################
KalMatrix KalMatrix::invert()
//#####################################
{
    KalMatrix mc;
    
    double det=a[0][0]*a[1][1]-a[0][1]*a[1][0];
    if (det!=0.){
	mc.a[0][0]=a[1][1]/det;
	mc.a[1][0]=-1.*a[1][0]/det;
	mc.a[0][1]=-1.*a[0][1]/det;
	mc.a[1][1]=a[0][0]/det;	
    }
    
    return mc;
}

//#####################################
KalMatrix KalMatrix::transpose()
//#####################################
{
    KalMatrix mc;
    
    for (int i=0;i<KalMatrix::kfdim;i++){	
	for (int j=0;j<KalMatrix::kfdim;j++){	
	    mc.a[i][j]=a[j][i];	
	}
    }
    
    return mc;
}

