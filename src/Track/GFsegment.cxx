

//--------------------------------------------------------
//
//   GFsegment 
//
//--------------------------------------------------------
#include "TkrRecon/Track/GFsegment.h"
#include "TkrRecon/Track/GFparticle.h"

//#########################################################
GFsegment::GFsegment(const GFtrack* _GFtrack)
: _mGFtrack(_GFtrack)
//#########################################################
{
    m_axis = _GFtrack->getAxis(); 
    clear();
}
//#########################################################
KalPlane GFsegment::getKPlane() const
//#########################################################
{   
    if (m_nextKplane.getIDPlane() < 0) {
        // std::cout << " error GFsegment::getKPlane " << "\n";
    }
    return m_nextKplane;
}

//########################################################
void GFsegment::clear()
//########################################################
{	
    m_indexhit = -1;
    m_statusHit = GFbase::EMPTY;
    m_nextKplane.clear();
    KalTrack::clear();
}
//########################################################
void GFsegment::next(int jplane)
//########################################################
{	
    clear();
    m_nextKplane.clear();
    
    KalPlane oriKplane = _mGFtrack->lastKPlane();
    
    doit(oriKplane,jplane,KalHit::FIT);
    
}
//########################################################
void GFsegment::previous(int jplane)
//########################################################
{	
    clear();
    m_nextKplane.clear();
    
    KalPlane oriKplane = _mGFtrack->firstKPlane();
    
    doit(oriKplane,jplane,KalHit::SMOOTH);
    
}

//#########################################################################
void GFsegment::best(int klayer)
//#########################################################################
{
    int indexhit1;
    std::vector<int> indexhit1list;
    indexhit1list.clear();
    
    double sigmac = _mGFtrack->sigmaCut();
    double best_chiSqhit = GFcontrol::minSegmentHits*sigmac*sigmac;
    GFsegment best_GFsegment(_mGFtrack);
    //unused:	int indexhit = -1;
    
    bool found = false;
    bool done = false;
    
    while (!done) {
        
        next(klayer);
        
        double chiSq = 1e6;
        if (!accept()) {
            done = true;
        } else {
            indexhit1 = m_indexhit;
            chiSq = chiGFSq();
            indexhit1list.push_back(indexhit1);
            
            GFtutor::_DATA->flagHit(m_axis,indexhit1);
        }
        
        if (accept() && chiSq < best_chiSqhit) {
            
            found = true;
            best_chiSqhit = chiSq;
            KalHit hitsmooth = getKPlane(klayer).getHit(KalHit::SMOOTH);
            KalHit hitfit = m_nextKplane.getHit(KalHit::FIT);
            KalHit hitNew(KalHit::FIT,hitsmooth.getPar(),hitfit.getCov());
            m_nextKplane.setHit(hitNew);
            best_GFsegment = *this;
        }
    }
    
    if (!found) clear();
    else *this = best_GFsegment;  
    
    for (unsigned int i=0; i<indexhit1list.size(); i++) 
        GFtutor::_DATA->unflagHit(m_axis, indexhit1list[i]);
    
}

//#########################################################
bool GFsegment::accept() const
//#########################################################
{
    bool accept = true;
    
    if (numDataPoints() < 3) {
        accept = false;
    }
    if (!accept) return accept;
    
    if (numGaps() > GFcontrol::maxGapsSegment) {
        accept = false;
    }
    if (!accept) return accept;
    
    if (chiGFSq() > GFcontrol::maxChiSqSegment) {
        accept = false;
    }
    if (!accept) return accept;
    
    return accept;
}

//########################################################
void GFsegment::flagUsedHits(int jplane)
//########################################################
{	
    unsigned i = 0;
    for ( ; i < kplanelist.size(); i++) {
        int kplane = kplanelist[i].getIDPlane();
        int indexhit = kplanelist[i].getIDHit();
        if (kplane >= jplane) {
            GFtutor::_DATA->flagHit(m_axis,indexhit);
        }
    }
    
}
//########################################################
void GFsegment::unFlagAllHits()
//########################################################
{	
    unsigned i = 0;
    for ( ; i < kplanelist.size(); i++) {
        //unused:		int kplane = kplanelist[i].getIDPlane();
        int indexhit = kplanelist[i].getIDHit();
        GFtutor::_DATA->unflagHit(m_axis,indexhit);
    }
    
}

//########################################################
double GFsegment::chiGFSq() const
//########################################################
{	
    double chiGF = 1e6;
    if (numDataPoints() < GFcontrol::minSegmentHits) return chiGF;
    
    double sigmaC = _mGFtrack->sigmaCut();
    chiGF =numGaps()*sigmaC*sigmaC/(1.*GFcontrol::minSegmentHits);
    chiGF += chiSquareSmooth();
    return chiGF;
}
//--------------------------------------------------------
// GFsegment - private
//--------------------------------------------------------

//#########################################################
KalPlane GFsegment::getKPlane(int kplane) const
//#########################################################
{
    
    KalPlane Kplane;
    
    //unused:	bool done = false;
    unsigned i = 0;
    for (; i< kplanelist.size(); i++) {
        if (kplanelist[i].getIDPlane() == kplane ) {
            Kplane =  kplanelist[i];
            break;
        }
    }
    
    return Kplane;
}
//#########################################################
KalPlane GFsegment::followingKPlane(int kplane) const
//#########################################################
{
    
    KalPlane Kplane;
    
    //unused:	bool done = false;
    unsigned i = 0;
    for (; i< kplanelist.size(); i++) {
        if (kplanelist[i].getIDPlane() > kplane ) {
            Kplane =  kplanelist[i];
            break;
        }
    }
    
    return Kplane;
    
}
//########################################################
void GFsegment::doit(KalPlane& oriKplane, int jplane, KalHit::TYPE type)
//########################################################
{	

    if (oriKplane.getIDPlane() != _mGFtrack->originalKPlane().getIDPlane()) 
        kplanelist.push_back(oriKplane);
    
    int kplane = jplane;
    int lstgaps = 0;
    int step_counter = 0; 
    GFbase::StatusHit statushit = GFbase::FOUND;
    while( -1 < kplane && kplane < GFtutor::numPlanes()) {

        step_counter++; 
        KalPlane prevKplane;
        KalPlane nextKplane;
        if (kplanelist.size() == 0) prevKplane = oriKplane; 
        else prevKplane = kplanelist.back();
        
        if (step_counter > 1) type = KalHit::FIT;
        GFbase::StatusHit statushit = nextKPlane(prevKplane, kplane, nextKplane, type); 
        if (kplane == jplane) m_statusHit = statushit;
        if (statushit != GFbase::FOUND) break; 

        lstgaps = nextKplane.getIDPlane() - prevKplane.getIDPlane(); 
        if (lstgaps == GFcontrol::maxConsecutiveGaps) break;
   
        kplanelist.push_back(nextKplane);
        if (kplanelist.size() >= 2) filterStep(kplanelist.size()-2);
        else {
             KalHit hitpred = nextKplane.getHit(KalHit::PRED);
             KalHit hitmeas = nextKplane.getHit(KalHit::MEAS);
             KalPar p(hitmeas.getPar().getXPosition(),hitpred.getPar().getXSlope(),
             hitmeas.getPar().getYPosition(),hitpred.getPar().getYSlope());
             KalHit hitfit(KalHit::FIT,p,hitpred.getCov());
             kplanelist[0].setHit(hitfit);
        }
     

        m_nextKplane = kplanelist.back();  // save the plane - the PRED changes afther doFit()
        m_indexhit   = m_nextKplane.getIDHit();
    }
    
    doFit();
}
//########################################################
GFbase::StatusHit GFsegment::nextKPlane(const KalPlane& previousKplane, 
                                        int kplane, KalPlane& nextKplane,
                                        KalHit::TYPE type)
//########################################################
{
    nextKplane = projectedKPlane(previousKplane, kplane, type);
    
    int indexhit = -1;
    double radius = 0.;
    GFbase::StatusHit statushit; 
    double sigma = sigmaFoundHit(previousKplane, nextKplane, indexhit, radius);
    if (sigma < _mGFtrack->m_sigmaCut) {
        statushit = GFbase::FOUND;
        incorporateFoundHit(nextKplane, indexhit);
    } else statushit = GFbase::EMPTY;
    
    return statushit;
}

//#########################################################################
KalPlane GFsegment::projectedKPlane(KalPlane prevKplane, int klayer, 
                                    KalHit::TYPE type) 
//########################################################################
{
    // The z end position of the next klayer plane
    int nlayers; 
    double zEnd, arc_min = 0.; 
    
    KalHit predhit = prevKplane.predicted(type, nlayers, klayer, zEnd, arc_min);
    KalPlane projectedKplane(prevKplane.getIDHit(),klayer+nlayers,prevKplane.getEnergy(), zEnd, 
        predhit, prevKplane.getNextProj());
    
    return projectedKplane;
}
//#########################################################################
void GFsegment::incorporateFoundHit(KalPlane& nextKplane, int indexhit)
//#########################################################################
{			
    Point nearHit = GFtutor::_DATA->position(m_axis, indexhit);
    double x0 = nearHit.x();
    double y0 = nearHit.y();
    double z0 = nearHit.z();
    
    KalPar measpar(x0,0.,y0,0.);
    
    double sigma = GFtutor::siResolution();
    double sigma_alt = 10.3; // 36/sqrt(12)
    double size  = GFtutor::_DATA->size(m_axis,indexhit);
    
    double factor = 1.;
    if (GFcontrol::sigmaCluster) factor = 1./size; // this must be the correct one But
    else factor = size*size;
    
    double cx, cy;
    if(m_axis == TkrCluster::X) {
        cx = sigma*sigma*factor;
        cy =10.3;
    }
    else {
        cx = 10.3;
        cy = sigma*sigma*factor;
    }
    KalMatrix meascov(1); 
    meascov(1,1) = cx;
    meascov(3,3) = cy;
    
    KalHit meashit(KalHit::MEAS, measpar, meascov);
    
    nextKplane.setIDHit(indexhit);
    nextKplane.setZPlane(z0);
    nextKplane.setHit(meashit);
}
//-------------------------------------------------------------
//   GFsegment -> Finding the Hit
//-------------------------------------------------------------
//#########################################################################
double GFsegment::sigmaFoundHit(const KalPlane& previousKplane, const KalPlane& nextKplane,
                                int& indexhit, double& radiushit)
//#########################################################################
{
    indexhit = -1;
    radiushit = 1e6;
    double sigmahit = 1e6;
    
    double DELTACHI2_CUT = _mGFtrack->m_sigmaCut;
    double MAX_RADIUS = GFtutor::trayWidth()/2.;
    double ERROR_ZPLANE= GFtutor::siThickness(); 
    
    KalHit hitp = nextKplane.getHit(KalHit::PRED);
    m_axis = nextKplane.getProjection();
    double tError = 0.;
    double zError = 0.; 
    if(m_axis == TkrCluster::X) {
        tError = sqrt(hitp.getCov().getcovX0X0());
        zError=ERROR_ZPLANE*hitp.getPar().getXSlope();
    }
    else {
        tError = sqrt(hitp.getCov().getcovY0Y0());
        zError=ERROR_ZPLANE*hitp.getPar().getYSlope();
    }
    double rError=sqrt(tError*tError+zError*zError);
    double outRadius=3.*DELTACHI2_CUT*rError; 
    if (outRadius > MAX_RADIUS ) outRadius = MAX_RADIUS;		
    
    double x0=hitp.getPar().getXPosition();
    double y0=hitp.getPar().getYPosition();
    double z0=nextKplane.getZPlane();
    Point center(x0,y0,z0);
    Point nearHit(0.,0.,z0);
    
    double inerRadius = -1.;
    int klayer = nextKplane.getIDPlane();
    
    // Must be inside Glast
    bool done = false;
    while (!done) {
        nearHit = GFtutor::_DATA->nearestHitOutside(m_axis, klayer, inerRadius, center, indexhit);
        done = foundHit(indexhit, inerRadius, outRadius, center, nearHit);
    }
    
    if (indexhit >= 0) {
        double dx_hit = nearHit.x() - center.x();
        double dy_hit = nearHit.y() - center.y();
        radiushit = sqrt(dx_hit*dx_hit + dy_hit*dy_hit);
        if (rError > 0.) sigmahit= radiushit/rError;
    }
    
    return sigmahit;
}

//#########################################################################
bool GFsegment::foundHit(int& indexhit, double& inerRadius, double outRadius,
                         const Point& centerX, const Point& nearHit)
//#########################################################################
{
    bool done = true;
    
    if (indexhit < 0) return done;
    double deltaX = (_mGFtrack->m_axis == TkrCluster::X? fabs(nearHit.x()-centerX.x()):
    fabs(nearHit.y()-centerX.y()));
    
    if (deltaX < outRadius) {
        if (GFtutor::_DATA->hitFlagged(_mGFtrack->m_axis, indexhit)) done = false;
    } else indexhit = -1; // outside region 
    
    // this condition is necessary for the large angle pattern recognition
    if (indexhit > 0 ) {
        if (GFtutor::okClusterSize(m_axis,indexhit,0.) == 0) done = false;
    }
    
    if (done == false) {
        indexhit = -1;
        inerRadius = deltaX + 0.1*GFtutor::siResolution();
    }
    
    return done;
}
//--------------------------------------------------------
// GFsegment - Utility functions
//--------------------------------------------------------
//########################################################
double GFsegment::getZklayer(enum TkrCluster::view axis, int klayer) const
//########################################################
{
    // Get the zEnd of the next plane	
    double zEnd =0;
    
    Point oneHit = GFtutor::_DATA->meanHit(axis, klayer);
    zEnd = oneHit.z();
    
    if (zEnd==0. && kplanelist.size()>0) zEnd = -50.; 
    
    return zEnd;
}	    
//########################################################
bool GFsegment::crack(const KalPlane& nextKplane) const
//########################################################
{
    bool crack = false;
    /*
    KalHit hitp = nextKplane.getHit(KalHit::PRED);
    double x0=0.;
    double y0=0.;
    if (m_axis==TkrCluster::X)  x0=hitp.getPar().getPosition();
    else  y0=hitp.getPar().getPosition();
    double z0=nextKplane.getZPlane();
    Point centerX(x0,y0,z0);
    
      float tower_width = GFtutor::side;
      + 2.*Tower::wall_thickness + Glast::wall_gap;
      float x_twr_bnd = fmod((centerX.x() + tower_width*2.5), tower_width);
      x_twr_bnd = .5*tower_width - fabs(x_twr_bnd -tower_width*.5);
      float y_twr_bnd = fmod((centerX.y() + tower_width*2.5), tower_width);
      y_twr_bnd = .5*tower_width - fabs(y_twr_bnd - tower_width*.5);
      float twr_bnd = (m_axis == TkrCluster::X) ? x_twr_bnd : y_twr_bnd;
      if(twr_bnd < .5) crack = true;
    */
    return crack;
}
