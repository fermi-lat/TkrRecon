
// $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalTrack.cxx,v 1.7 2002/03/30 20:40:54 lsrea Exp $

//----------------------------------------------------------------------
//    
//    Implementation of the Kalman Filter Functions
//               KalTrack 
//
//      Original due to Jose Hernando-Angel circa 1997-1999
//      Re-written to combine both X and Y projections (2001) 
//      
//-----------------------------------------------------------------------

#include "TkrRecon/TrackFit/KalTrack.h"
#include "geometry/Ray.h"
#include "TkrRecon/Track/GFcontrol.h"
#include "TkrRecon/GaudiAlg/TkrReconAlg.h"
#include "GlastSvc/Reco/IKalmanParticle.h"
#include <cmath>
#ifdef __CYGWIN__
# include <ieeefp.h>
#endif
bool CONTROL_setDeltaEne = false;
double CONTROL_MaxEne = 2.;
int CONTROL_maxPlanesEne = 8;


void KalTrack::clear()
{
    kplanelist.clear();
    ini();
}

void KalTrack::drawChiSq(gui::DisplayRep& v,  KalHit::TYPE typ)
{
    KalHit::TYPE fit = KalHit::SMOOTH;
    int nplanes = kplanelist.size();
    for (int iplane=0; iplane<nplanes; iplane++) {
        KalPlane::AXIS prj = kplanelist[iplane].getProjection();
        double x0, y0, z0, xl, xr, yl, yr;
        double delta= kplanelist[iplane].getDeltaChiSq(typ);
        if(prj == TkrCluster::X){
            x0 = kplanelist[iplane].getHit(typ).getPar().getXPosition();
            y0 = kplanelist[iplane].getHit(fit).getPar().getYPosition(); 
            z0 = kplanelist[iplane].getZPlane()+0.1;
            xl = x0-0.5*delta;
            xr = x0+0.5*delta;
            yl = y0;
            yr = y0;
        } 
        else {
            x0 = kplanelist[iplane].getHit(fit).getPar().getXPosition();
            y0 = kplanelist[iplane].getHit(typ).getPar().getYPosition(); 
            z0 = kplanelist[iplane].getZPlane()+0.1;
            xl = x0;
            xr = x0;
            yl = y0-0.5*delta;
            yr = y0+0.5*delta;
        }		
        v.moveTo(Point(xl,yl,z0));
        v.lineTo(Point(xr,yr,z0));
    }
}

void KalTrack::drawTrack(gui::DisplayRep& v, KalHit::TYPE typ)
{
    KalHit::TYPE fit = KalHit::SMOOTH;
    int nplanes = kplanelist.size();
    for (int iplane=0; iplane<nplanes-1; iplane++) {
        KalPlane::AXIS prj = kplanelist[iplane].getProjection();
        double x0, y0, z0;

		KalHit::TYPE xtyp, ytyp;
		xtyp = (prj == TkrCluster::X ? typ : fit);
		ytyp = (prj == TkrCluster::X ? fit : typ);

		// this sets up the track segment to the next plane
    
        x0 = kplanelist[iplane].getHit(xtyp).getPar().getXPosition();
        y0 = kplanelist[iplane].getHit(ytyp).getPar().getYPosition(); 
        z0 = kplanelist[iplane].getZPlane();
 
		double tanx = kplanelist[iplane].getHit(typ).getPar().getXSlope();
        double tany = kplanelist[iplane].getHit(typ).getPar().getYSlope();
        
        Point origin(x0,y0,z0);
        Vector dir = Vector(-1.*tanx,-1.*tany,-1.).unit();
        
        Ray segment(origin,dir);
        double zstep=kplanelist[iplane+1].getZPlane()-z0;
        double cosz=dir.z();

        // this sets up the dotted line from the lower part of the extrapolated track
		//  to the next hit.


		prj = kplanelist[iplane].getNextProj();

		xtyp = (prj == TkrCluster::X ? typ : fit);
		ytyp = (prj == TkrCluster::X ? fit : typ);

        x0 = kplanelist[iplane+1].getHit(xtyp).getPar().getXPosition();
        y0 = kplanelist[iplane+1].getHit(ytyp).getPar().getYPosition(); 
        z0 = kplanelist[iplane+1].getZPlane();

        Point p(x0, y0, z0);

		// do them in this order, so that the connection doesn't cover the track
		
		v.set_line_style(1);
        v.moveTo(segment.position(0.8*zstep/cosz));
        v.lineTo(p); 

		v.setColor("blue");
        v.moveTo(segment.position(0.));
        v.lineTo(segment.position(zstep/cosz));
    }
}

int KalTrack::numGaps() const
{
    int numGaps =0;
    if (numDataPoints() == 0) return numGaps;
    
    numGaps =1+ kplanelist.back().getIDPlane() -
        kplanelist.front().getIDPlane() - numDataPoints();
    
    return numGaps;
}

int KalTrack::compareFits(KalTrack& ktrack)
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

Point KalTrack::getHit(unsigned ipos)const
{
    if (ipos<kplanelist.size()){
        return kplanelist[ipos].getPoint(KalHit::MEAS);
    }
    //  else return Point(FLT_MAX,FLT_MAX,FLT_MAX);
    return Point(-999999.,-999999.,-999999.);
}

unsigned KalTrack::getHitIndex(unsigned ipos)const
{
    unsigned index= 0xFFFFFFFF;
    if (ipos<kplanelist.size()) 
        index=kplanelist[ipos].getIDHit();
    return index;
}

Point KalTrack::positionAtZ(double const z) const 
{
    if( kplanelist.empty() ) return Point(0.,0.,0.);
    const KalPlane& plane = *kplanelist.begin();
    double  x = plane.getHit(KalHit::SMOOTH).getPar().getXPosition();
    double  y = plane.getHit(KalHit::SMOOTH).getPar().getYPosition();
    double  z0  = plane.getZPlane()+0.01;
    double  tanX  = plane.getHit(KalHit::SMOOTH).getPar().getXSlope();
    double  tanY  = plane.getHit(KalHit::SMOOTH).getPar().getYSlope();
    return Point(x+(z-z0)*tanX, y+(z-z0)*tanY, z);
}

void KalTrack::setIniEnergy(double e) 
{
    m_energy0 = e;
    for (unsigned int iplane = 0; iplane < kplanelist.size(); iplane++)
        kplanelist[iplane].setEnergy(m_energy0);
} 

KalTrack::KalTrack()
: m_energy0(0) 
, m_chisq(1e6), m_chisqSmooth(1e6) 
, m_KalEnergy(0), m_KalThetaMS(0), m_rmsResid(0)
, m_numSegmentPoints(0), m_chisqSegment(0)
{}

double KalTrack::doFit()
{
    ini();
    
    int nplanes=kplanelist.size();
    if (nplanes<=4) return m_chisq;
      
    // Generate the initial hit to start the Kalman Filter
    //----------------------------------------------------
    KalHit hitf=generateFirstFitHit();
    if(hitf.getType() != KalHit::FIT) return m_chisq; // failure! 
    kplanelist[0].setHit(hitf); 

    m_chisq       = 0.;
    m_chisqSmooth = 0.;
    
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

void KalTrack::ini()
{
    //m_energy0      = 0.;
    m_x0           = Point(0., 0., 0.);
    m_dir          = Vector(0., 0., 0.);
    m_rmsResid     = 0.;
    m_KalEnergy    = 0.;
    m_chisq        = 1e6;
    m_chisqSmooth  = 1e6;
    m_KalThetaMS   = 0.;
    m_rmsResid     = 0.;
    m_numSegmentPoints = 0;
    m_chisqSegment = 1e6;
    
    unsigned int iplane =0;
    for (;iplane < kplanelist.size();iplane++) {
        kplanelist[iplane].clean();
    }
    
}

void KalTrack::finish()
{
    // Compute the fit variables  
    if (m_chisq>=0){
        int nplanes = kplanelist.size();
        double x =kplanelist[0].getHit(KalHit::SMOOTH).getPar().getXPosition();
        double y =kplanelist[0].getHit(KalHit::SMOOTH).getPar().getYPosition();
        double z =kplanelist[0].getZPlane() + .01; 
        m_x0 = Point(x,y,z);
        double x_slope =kplanelist[0].getHit(KalHit::SMOOTH).getPar().getXSlope();
        double y_slope =kplanelist[0].getHit(KalHit::SMOOTH).getPar().getYSlope();
        m_dir = Vector(-1.*x_slope,-1.*y_slope,-1.).unit();
        m_rmsResid=-1.;
        m_chisq=m_chisq/(1.*nplanes-4.); // 4 parameters in 3D fit
        m_chisqSmooth/=(1.*nplanes-4.);  
        m_rmsResid=0.;
        int iplane = 0; 
        double xm; 
        for (iplane=0;iplane<nplanes;iplane++){
            if(kplanelist[iplane].getProjection() == TkrCluster::X) {
                x =kplanelist[iplane].getHit(KalHit::SMOOTH).getPar().getXPosition();
                xm =kplanelist[iplane].getHit(KalHit::MEAS).getPar().getXPosition();
            }
            else {
                x =kplanelist[iplane].getHit(KalHit::SMOOTH).getPar().getYPosition();
                xm =kplanelist[iplane].getHit(KalHit::MEAS).getPar().getYPosition();
            }
            m_rmsResid+= (x-xm)*(x-xm);
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

void KalTrack::filterStep(int iplane) 
{
    KalHit hitp = kplanelist[iplane].predicted(kplanelist[iplane+1]);
    kplanelist[iplane+1].setHit(hitp);
    KalHit hitf1 = kplanelist[iplane+1].filter();
    kplanelist[iplane+1].setHit(hitf1);
    
    if (CONTROL_setDeltaEne) 
        kplanelist[iplane+1].setDeltaEne(kplanelist[iplane].getEnergy());
    
}

KalHit KalTrack::generateFirstFitHit()
{   
    int nplanes=kplanelist.size();
    if (nplanes<4) {
        std::cout << "ERROR - Kaltrack::generateFirstFitHit - too few planes" << '\n';
        return KalHit();
    }
    //Find first two x hits and first two y hits
    double x0,x1, y0,y1, zx0, zx1, zy0,zy1; 
    int nx = 0, ny=0;
    for(int i=0; i<8 && i< nplanes; i++) {
        if(kplanelist[i].getProjection() == TkrCluster::X && nx < 2) {
            nx++;
            if(nx == 1) {
                x0 =kplanelist[i].getHit(KalHit::MEAS).getPar().getXPosition();
                zx0 =kplanelist[i].getZPlane() + .01; 
            }
            else {
                x1 =kplanelist[i].getHit(KalHit::MEAS).getPar().getXPosition();
                zx1 =kplanelist[i].getZPlane() + .01;
            }
        }
        else if(kplanelist[i].getProjection() == TkrCluster::Y && ny < 2) {
            ny++;
            if(ny == 1) {
                y0 =kplanelist[i].getHit(KalHit::MEAS).getPar().getYPosition();
                zy0 =kplanelist[i].getZPlane() + .01; 
            }
            else {
                y1 =kplanelist[i].getHit(KalHit::MEAS).getPar().getYPosition();
                zy1 =kplanelist[i].getZPlane() + .01; 
            }
        }
        if(nx==2 && ny==2) break;
    }
    if(nx != 2 || ny!=2) {
        std::cout << "ERROR - Kaltrack::generateFirstFitHit: nx or ny != 2" << '\n';
        return KalHit();
    }
    double x_slope = (x1-x0)/(zx1-zx0);
    double y_slope = (y1-y0)/(zy1-zy0);
    m_dir = Vector(-1.*x_slope,-1.*y_slope,-1.).unit();
    double x_ini, y_ini, z_ini;
    if(zx0 > zy0) {// extrapolate the y co-ordinate back
        z_ini = zx0;
        x_ini = x0;
        y_ini = y0 + y_slope*(zx0-zy0);
    }
    else {         // ... extraoplate the x co-ord. back
        z_ini = zy0;
        y_ini = y0;
        x_ini = x0 + x_slope*(zy0-zx0);
    }
    m_x0 = Point(x_ini,y_ini,z_ini);
    
//    std::auto_ptr<IKalmanParticle> 
//        kalPart(TkrReconAlg::m_gismoSvc->kalmanParticle(m_x0, m_dir, fabs((zx0-zy0)/m_dir.z())));
    
    double energy = kplanelist[1].getEnergy();
    if (energy == 0.) energy = m_energy0;
//    double dist = kalPart->arcLength(); 
    KalMatrix m; // = kalPart->mScat_Covr(energy, dist); 
//    m(1,1) = GFcontrol::iniErrorPosition;
//    m(3,3) = GFcontrol::iniErrorPosition;
    m(2,2) = GFcontrol::iniErrorSlope * GFcontrol::iniErrorSlope;
    m(4,4) = GFcontrol::iniErrorSlope * GFcontrol::iniErrorSlope;
    //  The first error is arbitrary to a degree: choose 10*(ms + 1-plane hit precision)
   //m *= 10.; 
    
    KalPar parguess(Ray(m_x0, m_dir));
    KalHit hitf(KalHit::FIT,parguess, 
        (kplanelist[0].getHit(KalHit::MEAS)).getCov()+m);
    
    return hitf;
}

void KalTrack::eneDetermination()
{
    int nplanes = numDataPoints();
    int iplane = 2;
    double totalRad = 0.;
    double eneSum = 0.;
    double thetaSum = 0.;
    double x0 = kplanelist[1].getHit(KalHit::SMOOTH).getPar().getXPosition();
    double y0 = kplanelist[1].getHit(KalHit::SMOOTH).getPar().getYPosition();
    double z0 = kplanelist[1].getZPlane();
    Point x_ini(x0, y0, z0); 
    double slopeX = kplanelist[1].getHit(KalHit::SMOOTH).getPar().getXSlope();
    double slopeY = kplanelist[1].getHit(KalHit::SMOOTH).getPar().getYSlope();
    Vector dir_ini = Vector(-slopeX, -slopeY, -1.).unit();
    
    for (iplane = 2; iplane < nplanes; iplane++) {
        double chie = kplanelist[iplane].getDeltaChiEne(KalHit::PRED);
        double eta = ((chie-1.)*(chie-1.)*(chie-1.))/((chie+2.)*(chie+2.));
        eta = sqrt(fabs(eta));
        double z1 = kplanelist[iplane].getZPlane();
        std::auto_ptr<IKalmanParticle> 
            kalPart(TkrReconAlg::m_gismoSvc->kalmanParticle(m_x0, m_dir, (z1-z0)/dir_ini.z()));
        
        totalRad += kalPart->radLength();
        double factor = 1./(2.-exp(-1.*totalRad));
        // double factor = 1./(1.+totalRad);
        
        //       double sigma = sqrt(kplanelist[iplane].getHit(KalHit::MEAS).getCov().getsiga());
        //       double distance = abs(kplanelist[iplane].getZPlane()-kplanelist[iplane-1].getZPlane());
        
        // factor - energy loss factor
        // sigma/distance - dimensions
        // eta - undimesnionless parameters
        // cosX, etx - geometrical factors
        //       double theta = factor*(sigma/distance)*eta*(cosX*cosZ*sqrt(cosZ));
        //      theta /= (sqrt(radlen)*(1+0.038*log(radlen)));
        //      thetaSum += theta;
        
    }
    m_KalThetaMS = thetaSum/(nplanes-2.);
    m_KalEnergy = 0.0136/m_KalThetaMS;
    //    double radlen = KalPlane::radLen(kplanelist[0].getIDPlane());
    //    m_KalThetaMS *= sqrt(radlen)*(1+0.038*log(radlen));
}

double KalTrack::errorXPosition() const
{
    double error = 0.;
    if (m_chisqSmooth == 0) return error;
    error = sqrt(kplanelist[0].getHit(KalHit::SMOOTH).getCov().getcovX0X0());
    return error;
}

double KalTrack::errorXSlope() const
{
    double error = 0.;
    if (m_chisqSmooth == 0) return error;
    error = sqrt(kplanelist[0].getHit(KalHit::SMOOTH).getCov().getcovSxSx());
    return error;
}

double KalTrack::errorYPosition() const
{
    double error = 0.;
    if (m_chisqSmooth == 0) return error;
    error = sqrt(kplanelist[0].getHit(KalHit::SMOOTH).getCov().getcovY0Y0());
    return error;
}

double KalTrack::errorYSlope() const
{
    double error = 0.;
    if (m_chisqSmooth == 0) return error;
    error = sqrt(kplanelist[0].getHit(KalHit::SMOOTH).getCov().getcovSySy());
    return error;
}

int KalTrack::computeNumSegmentPoints(KalHit::TYPE typ)
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

double KalTrack::kink(int iplane) const
{
    double kink = 0.;
    if (iplane<0 || iplane > numDataPoints()-2) return kink;
    
    double slope0X = kplanelist[iplane].getHit(KalHit::SMOOTH).getPar().getXSlope();
    double slope1X= kplanelist[iplane+2].getHit(KalHit::SMOOTH).getPar().getXSlope();
    
    double slope0Y = kplanelist[iplane].getHit(KalHit::SMOOTH).getPar().getYSlope();
    double slope1Y= kplanelist[iplane+2].getHit(KalHit::SMOOTH).getPar().getYSlope();
    
    
    kink = (fabs(slope1X-slope0X) > fabs(slope1Y-slope0Y)) ? 
        (slope1X-slope0X) : (slope1Y-slope0Y);
    
    return kink;
}

double KalTrack::kinkNorma(int iplane) const
{
    double k = 0.;
    k = kink(iplane);
    
    if (iplane <0 || iplane >= numDataPoints()) return k;
    
    double errorXSQ = kplanelist[iplane+1].getHit(KalHit::SMOOTH).getCov().getcovSxSx();
    double errorYSQ = kplanelist[iplane+1].getHit(KalHit::SMOOTH).getCov().getcovSySy();
    errorXSQ += kplanelist[iplane+1].getQmaterial().getcovSxSx();
    errorYSQ += kplanelist[iplane+1].getQmaterial().getcovSySy();
    
    double error = sqrt(errorXSQ + errorYSQ);
    
    double kinkNorma =  0.;
    if (error != 0.) kinkNorma = fabs(k)/error;
    
    return kinkNorma;
    
}
