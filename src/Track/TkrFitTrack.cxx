//------------------------------------------------------------------------------
//
//     TkrFitTrack
//
//      Is a Kalman Track Follower class for GLAST
//
//      It uses a Ray direction as starting seed
//      and picks Si-Clusters inside a n-sigma region 
//      of the projected track into the Si Plane
//
//
//	Addapted from GFtrack by JA Hernando
//      W. B. Atwood, SCIPP/UCSC, Nov.,2001
//				
//------------------------------------------------------------------------------


#include "TkrRecon/Track/TkrFitTrack.h"
#include "TkrRecon/Track/GFtutor.h"
#include "TkrRecon/Track/GFcontrol.h"
//-----------------------------------------------------
//
//   TkrFitTrack
//
//-----------------------------------------------------

TkrFitTrack::TkrFitTrack(int ilyr, int itwr, double energy, const Ray& testRay)
                         : TkrBase(ilyr, itwr, energy, testRay.position(), testRay.direction()),
                          m_Q(-1e6), m_Xgaps(0), m_XistGaps(0), m_Ygaps(0), m_YistGaps(0)
{
    m_hits.clear();
}

TkrFitTrack::~TkrFitTrack()
{
    clear();
}

bool TkrFitTrack::empty() const
{
    bool empty = false;
    if (firstLayer()   < 0)                         empty = true;
    if (getNumHits()   < GFcontrol::minSegmentHits) empty = true;
    if (getChiSquare() < 0.)                        empty = true;

    return empty;
}

void TkrFitTrack::clear()
{   
    m_hits.clear();
    
    m_Xgaps        = m_Ygaps = 0;
    m_XistGaps     = m_YistGaps= 0;
    
    m_Q            = -1e6;
}

void TkrFitTrack::writeOut(MsgStream& log) const
{
    
    log << MSG::DEBUG << " --- TkrFitTrack::writeOut --- "            << endreq;
    log << MSG::DEBUG << " quality        = "   << getQuality()       << endreq;
    log << MSG::DEBUG << " num m_hits       = " << getNumHits()       << endreq;
    log << MSG::DEBUG << " num X Gaps     = "   << getNumXGaps()      << endreq;
    log << MSG::DEBUG << " num X 1st Gaps = "   << getNumXFirstGaps() << endreq;
    log << MSG::DEBUG << " num Y Gaps     = "   << getNumYGaps()      << endreq;
    log << MSG::DEBUG << " num Y 1st Gaps = "   << getNumYFirstGaps() << endreq;

    TkrBase::writeOut(log);    
}

void TkrFitTrack::draw(gui::DisplayRep& v) 
{
    v.markerAt(position());
    drawChiSq(v,TkrFitHit::SMOOTH);
    drawTrack(v,TkrFitHit::SMOOTH);
}


TkrFitPlane TkrFitTrack::getFirstPlane() const
{
    if (m_hits.size() == 0) {
        std::cout << "ERROR TkrFitTrack::thisKPlane " << endreq;
        return TkrFitPlane();
    }
    return m_hits.front();
}

TkrFitPlane TkrFitTrack::getLastPlane() const
{
    if (m_hits.size() == 0) {
        return TkrFitPlane();
    }
    return m_hits.back();
}

int TkrFitTrack::getNumGaps() const
{
    int numGaps =0;
    if (getNumHits() == 0) return numGaps;
    
    numGaps = 1 + m_hits.back().getIDPlane() - m_hits.front().getIDPlane() - getNumHits();
    
    return numGaps;
}

void TkrFitTrack::drawChiSq(gui::DisplayRep& v,  TkrFitHit::TYPE typ)
{
    TkrFitHit::TYPE fit = TkrFitHit::SMOOTH;
    int nplanes = m_hits.size();
    for (int iplane=0; iplane<nplanes; iplane++) {
        TkrFitPlane::AXIS prj = m_hits[iplane].getProjection();
        double x0, y0, z0, xl, xr, yl, yr;
        double delta= m_hits[iplane].getDeltaChiSq(typ);
        if(prj == TkrCluster::X){
            x0 = m_hits[iplane].getHit(typ).getPar().getXPosition();
            y0 = m_hits[iplane].getHit(fit).getPar().getYPosition(); 
            z0 = m_hits[iplane].getZPlane()+0.1;
            xl = x0-0.5*delta;
            xr = x0+0.5*delta;
            yl = y0;
            yr = y0;
        } 
        else {
            x0 = m_hits[iplane].getHit(fit).getPar().getXPosition();
            y0 = m_hits[iplane].getHit(typ).getPar().getYPosition(); 
            z0 = m_hits[iplane].getZPlane()+0.1;
            xl = x0;
            xr = x0;
            yl = y0-0.5*delta;
            yr = y0+0.5*delta;
        }		
        v.moveTo(Point(xl,yl,z0));
        v.lineTo(Point(xr,yr,z0));
    }
}

void TkrFitTrack::drawTrack(gui::DisplayRep& v, TkrFitHit::TYPE typ)
{
    TkrFitHit::TYPE fit = TkrFitHit::SMOOTH;
    int nplanes = m_hits.size();
    for (int iplane=0; iplane<nplanes-1; iplane++) {
        TkrFitPlane::AXIS prj = m_hits[iplane].getProjection();
        double x0, y0, z0;

		TkrFitHit::TYPE xtyp, ytyp;
		xtyp = (prj == TkrCluster::X ? typ : fit);
		ytyp = (prj == TkrCluster::X ? fit : typ);

		// this sets up the track segment to the next plane
    
        x0 = m_hits[iplane].getHit(xtyp).getPar().getXPosition();
        y0 = m_hits[iplane].getHit(ytyp).getPar().getYPosition(); 
        z0 = m_hits[iplane].getZPlane();
 
		double tanx = m_hits[iplane].getHit(typ).getPar().getXSlope();
        double tany = m_hits[iplane].getHit(typ).getPar().getYSlope();
        
        Point origin(x0,y0,z0);
        Vector dir = Vector(-1.*tanx,-1.*tany,-1.).unit();
        
        Ray segment(origin,dir);
        double zstep=m_hits[iplane+1].getZPlane()-z0;
        double cosz=dir.z();

        // this sets up the dotted line from the lower part of the extrapolated track
		//  to the next hit.


		prj = m_hits[iplane].getNextProj();

		xtyp = (prj == TkrCluster::X ? typ : fit);
		ytyp = (prj == TkrCluster::X ? fit : typ);

        x0 = m_hits[iplane+1].getHit(xtyp).getPar().getXPosition();
        y0 = m_hits[iplane+1].getHit(ytyp).getPar().getYPosition(); 
        z0 = m_hits[iplane+1].getZPlane();

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

double TkrFitTrack::getKink(int iplane) const
{
    double kink = 0.;
    if (iplane<0 || iplane > getNumHits()-2) return kink;
    
    double slope0X = m_hits[iplane  ].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
    double slope1X = m_hits[iplane+2].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
    
    double slope0Y = m_hits[iplane  ].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
    double slope1Y = m_hits[iplane+2].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
    
    
    kink = (fabs(slope1X-slope0X) > fabs(slope1Y-slope0Y)) ? (slope1X-slope0X) : (slope1Y-slope0Y);
    
    return kink;
}

double TkrFitTrack::getKinkNorma(int iplane) const
{
    double k = getKink(iplane);
    
    if (iplane <0 || iplane >= getNumHits()) return k;
    
    double errorXSQ = m_hits[iplane+1].getHit(TkrFitHit::SMOOTH).getCov().getcovSxSx();
    double errorYSQ = m_hits[iplane+1].getHit(TkrFitHit::SMOOTH).getCov().getcovSySy();

    errorXSQ += m_hits[iplane+1].getQmaterial().getcovSxSx();
    errorYSQ += m_hits[iplane+1].getQmaterial().getcovSySy();
    
    double error = sqrt(errorXSQ + errorYSQ);
    
    double kinkNorma =  0.;
    if (error != 0.) kinkNorma = fabs(k)/error;
    
    return kinkNorma;
}

double TkrFitTrack::getErrorXPosition() const
{
    double error = 0.;
    if (m_chisqSmooth == 0) return error;
    error = sqrt(m_hits[0].getHit(TkrFitHit::SMOOTH).getCov().getcovX0X0());
    return error;
}

double TkrFitTrack::getErrorXSlope() const
{
    double error = 0.;
    if (m_chisqSmooth == 0) return error;
    error = sqrt(m_hits[0].getHit(TkrFitHit::SMOOTH).getCov().getcovSxSx());
    return error;
}

double TkrFitTrack::getErrorYPosition() const
{
    double error = 0.;
    if (m_chisqSmooth == 0) return error;
    error = sqrt(m_hits[0].getHit(TkrFitHit::SMOOTH).getCov().getcovY0Y0());
    return error;
}

double TkrFitTrack::getErrorYSlope() const
{
    double error = 0.;
    if (m_chisqSmooth == 0) return error;
    error = sqrt(m_hits[0].getHit(TkrFitHit::SMOOTH).getCov().getcovSySy());
    return error;
}