

//--------------------------------------------------------
//
//   GFsegment 
//
//--------------------------------------------------------
#include "TkrRecon/GFsegment.h"
#include "TkrRecon/GFparticle.h"

//Add 9/12/01 to put in tower info
#include "idents/TowerId.h"

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

    int trkSize        = _mGFtrack->kplanelist.size();

    if (trkSize > 0)
    {
        if (trkSize > 1)
        {
            KalPlane prevPlane = _mGFtrack->previousKPlane();
            kplanelist.push_back(prevPlane);
        }

        kplanelist.push_back(oriKplane);
    }
        
    doit(oriKplane,jplane,KalHit::FIT);
    
}
//########################################################
void GFsegment::previous(int jplane)
//########################################################
{	
    clear();
    m_nextKplane.clear();
    
    KalPlane oriKplane = _mGFtrack->firstKPlane();

    if (oriKplane.getIDPlane() != _mGFtrack->originalKPlane().getIDPlane()) 
	kplanelist.push_back(oriKplane);
    
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
    
    bool found      = false;
	bool done       = false;
    bool attemptFix  = true;
	
	while (!done) 
    {
	    next(klayer);
	    
	    double chiSq = 1e6;

	    if (!accept()) 
        {
            if (attemptFix && numDataPoints() == 2) 
            {
                int badHitIndex = kplanelist.back().getIDHit();
		        indexhit1list.push_back(badHitIndex);
                GFtutor::_DATA->flagHit(m_axis,badHitIndex);
                attemptFix = false;
            }
            else done = true;
	    } 
        else 
        {
		    indexhit1 = m_indexhit;
		    chiSq     = chiGFSq();

		    indexhit1list.push_back(indexhit1);

            GFtutor::_DATA->flagHit(m_axis,indexhit1);
	    }

	    if (accept() && chiSq < best_chiSqhit) 
        {
		    found            = true;
		    best_chiSqhit    = chiSq;
		    KalHit hitsmooth = getKPlane(klayer).getHit(KalHit::SMOOTH);
		    KalHit hitfit    = m_nextKplane.getHit(KalHit::FIT);
		    KalHit hitNew(KalHit::FIT,hitsmooth.getPar(),hitfit.getCov());
		    m_nextKplane.setHit(hitNew);
		    best_GFsegment   = *this;
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
void GFsegment::unFlagUsedHits(int jplane)
//########################################################
{	
    for (int i = 0 ; i < kplanelist.size(); i++) 
    {
	    int kplane = kplanelist[i].getIDPlane();
	    int indexhit = kplanelist[i].getIDHit();
        if (kplane >= jplane) 
        {
            GFtutor::_DATA->unflagHit(m_axis,indexhit);
        }
    }
    
}
//########################################################
void GFsegment::unFlagAllHits()
//########################################################
{	
    for (int i = 0; i < kplanelist.size(); i++) 
    {
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
    int iplane = oriKplane.getIDPlane();
    
    int up = ( iplane < jplane ? +1 : -1);

    //Loop over planes trying to form a "segment" with three hits
    int kplane = jplane;
    int lstgaps = 0;
    for (kplane = jplane ; -1 < kplane && kplane < GFtutor::numPlanes(); kplane+=up) 
    {
		KalPlane prevKplane;
		KalPlane nextKplane;

		if (kplanelist.size() == 0) prevKplane = oriKplane; 
		else                        prevKplane = kplanelist.back();
	
		if (kplane != jplane) type = KalHit::FIT;

		GFbase::StatusHit statushit = nextKPlane(prevKplane, kplane, nextKplane, type);
	
        //Found a candidate hit, include info here
        if (statushit == GFbase::FOUND) 
        { 
			lstgaps = 0;
            
			kplanelist.push_back(nextKplane);
			if (kplanelist.size() >= 2) 
            {
                filterStep(kplanelist.size()-2);
            }
			else 
            {
				KalHit hitpred = nextKplane.getHit(KalHit::PRED);
				KalHit hitmeas = nextKplane.getHit(KalHit::MEAS);
				KalPar p(hitmeas.getPar().getPosition(),hitpred.getPar().getSlope());
				KalHit hitfit(KalHit::FIT,p,hitpred.getCov());
				kplanelist[0].setHit(hitfit);
			}
		} 
        //No hit found, determine what to do
        else
        {
            int numNew = numDataPoints() - 1;
            int numOld = _mGFtrack->kplanelist.size();

            if (numNew + numOld < 3) break;

            lstgaps++;
        }
	
		if (kplane == jplane) 
        {
			m_statusHit = statushit;

			if (statushit != GFbase::FOUND) 
            {
                break;
            }

			m_nextKplane = kplanelist.back();  // save the plane - the PRED changes afther doFit()
			m_indexhit	 = m_nextKplane.getIDHit();
		}
	
		if (lstgaps == GFcontrol::maxConsecutiveGaps) 
        {
            break;
        }
		if (numDataPoints() == GFcontrol::minSegmentHits) 
        {
            break;
        }
    }
    
    doFit();
}
//########################################################
GFbase::StatusHit GFsegment::nextKPlane(const KalPlane& previousKplane, 
					int kplane, KalPlane& nextKplane,
					KalHit::TYPE type) const
//########################################################
{
    nextKplane = projectedKPlane(previousKplane, kplane, type);
    
    int indexhit = -1;
    double radius = 0.;
    GFbase::StatusHit statushit; 

    double sigma = sigmaFoundHit(previousKplane, nextKplane, indexhit, radius);

    if (sigma < _mGFtrack->m_sigmaCut) 
    {
	  statushit = GFbase::FOUND;
	  incorporateFoundHit(nextKplane, indexhit);
    }
    else statushit = GFbase::EMPTY;
    
    return statushit;
}

//#########################################################################
KalPlane GFsegment::projectedKPlane(KalPlane prevKplane, int klayer, 
				    KalHit::TYPE type) const
//########################################################################
{
    // The z end position of the next klayer plane
    // KalHit::TYPE type = KalHit::FIT;
    double zEnd = GFsegment::getZklayer(m_axis, klayer);
    
    KalHit predhit = prevKplane.predicted(type, zEnd, klayer);
    KalPlane projectedKplane(prevKplane.getIDHit(),klayer,prevKplane.getEnergy(), zEnd, 
	prevKplane.getOrthPar(), predhit);
    
    return projectedKplane;
}
//#########################################################################
void GFsegment::incorporateFoundHit(KalPlane& nextKplane, int indexhit) const
//#########################################################################
{			
    Point nearHit = GFtutor::_DATA->position(m_axis, indexhit);
    double x0 =0.;
    double z0 = nearHit.z();
    x0 = (_mGFtrack->m_axis == SiCluster::X ? nearHit.x() : nearHit.y());
    KalPar measpar(x0,0.);
    
    double  sigma   = GFtutor::siResolution();
    double  size    = GFtutor::_DATA->size(m_axis,indexhit);
    //int     towerId = GFtutor::_DATA->getHit(indexhit)->tower();

    idents::TowerId tower   = idents::TowerId(GFtutor::_DATA->getHit(indexhit)->tower());
    int     towerId = 10*(tower.iy()+1) + (tower.ix()+1);

	double factor = 1.;
    //if (GFcontrol::sigmaCluster) factor = 1./size; // this must be the correct one But
    if (GFcontrol::sigmaCluster) factor = size*size; 
	else factor = size*size;

	KalMatrix meascov(sigma*sigma*factor,0.,0.); 
    
    KalHit meashit(KalHit::MEAS, measpar, meascov);
    
    nextKplane.setIDHit(indexhit,towerId);
    nextKplane.setZPlane(z0);
    nextKplane.setHit(meashit);
}
//-------------------------------------------------------------
//   GFsegment -> Finding the Hit
//   9/18/01 Tracy Usher 
//   Have tried to patch this routine a bit. Try to add a check to make
//   sure hit is no further away than a neighboring tower... 
//-------------------------------------------------------------
//#########################################################################
double GFsegment::sigmaFoundHit(const KalPlane& previousKplane, const KalPlane& nextKplane,
				int& indexhit, double& radiushit) const
				//#########################################################################
{
    indexhit = -1;
    radiushit = 1e6;
    double sigmahit = 1e6;

    //Set up for current tower and tower tolerance
    double TowerWidth = 4.4 * GFtutor::trayWidth();
    int    oldIndex   = previousKplane.getIDHit();
    Point  oldCoords(0.,0.,0.);

    if (oldIndex > -1)
    {
        oldCoords  = GFtutor::_DATA->position(_mGFtrack->m_axis, oldIndex);
        TowerWidth = 1.1 * GFtutor::trayWidth();
    }

    
    double MAX_RADIUS = GFtutor::trayWidth()/2.;
    //double ERROR_ZPLANE= GFtutor::siThickness(); 
    
    KalHit hitp = nextKplane.getHit(KalHit::PRED);

    //Use predicted position to estimate the search region
    double xError = sqrt(hitp.getCov().getsiga());
    double outRadius = _mGFtrack->m_sigmaCut*xError;

    //The below commented out since I don't understand the zError part
    //I think the idea is to increase the search region for large angle 
    //tracks, but I am not sure why you would want to do that...
    //double zError=ERROR_ZPLANE*hitp.getPar().getSlope();
    //double rError=sqrt(xError*xError+zError*zError);
    //outRadius=3.*_mGFtrack->m_sigmaCut*rError;
    
    //Keep outside radius from getting too big
    if (outRadius > MAX_RADIUS ) outRadius = MAX_RADIUS;		
    
    double x0=0;
    double y0=0;
    if (_mGFtrack->m_axis==SiCluster::X)
    {
        x0 = hitp.getPar().getPosition();
        y0 = previousKplane.getOrthPar().getPosition();
        //y0 = oldCoords.y();
    }
    else
    {
        //x0 = oldCoords.x();
        x0 = previousKplane.getOrthPar().getPosition();
        y0 = hitp.getPar().getPosition();
    }

    double z0=nextKplane.getZPlane();

    Point centerX(x0,y0,z0);
    Point nearHit(0.,0.,z0);
    
    double inerRadius = -1.;
    int klayer = nextKplane.getIDPlane();
    double slope = hitp.getPar().getSlope();
    
    // Must be inside Glast
	bool done = false;
	while (!done) {
	    nearHit = GFtutor::_DATA->nearestHitOutside(m_axis, klayer, inerRadius, centerX, indexhit);
	    done = foundHit(indexhit, inerRadius, outRadius, TowerWidth, centerX, nearHit, slope);
	}
    
    if (indexhit >= 0) 
    {
	    double x_hit	= (_mGFtrack->m_axis == SiCluster::X? nearHit.x() : nearHit.y());
	    double x_center = (_mGFtrack->m_axis == SiCluster::X? centerX.x() : centerX.y());
	    radiushit = fabs(x_hit-x_center);
	    if (xError > 0.) sigmahit= radiushit/xError;
	    //if (rError > 0.) sigmahit= radiushit/rError;
    }
    
    return sigmahit;
}

//#########################################################################
bool GFsegment::foundHit(int& indexhit, double& inerRadius, double outRadius, double twrRadius,
			 const Point& centerX, const Point& nearHit, double slope) const
			 //#########################################################################
{
    bool done = true;
    
    if (indexhit < 0) return done;

    double hitDelta = fabs(nearHit.x() - centerX.x());
    double twrDelta = fabs(nearHit.y() - centerX.y());

    if (_mGFtrack->m_axis == SiCluster::Y)
    {
        hitDelta = fabs(nearHit.y() - centerX.y());
        twrDelta = fabs(nearHit.x() - centerX.x());
    }
    
    //if (hitDelta < outRadius && twrDelta < twrRadius)
    if (hitDelta < outRadius && twrDelta < twrRadius)
    {
	  if (GFtutor::_DATA->hitFlagged(_mGFtrack->m_axis, indexhit)) done = false;
    } 
    else indexhit = -1; // outside region 

    //double deltaX = (_mGFtrack->m_axis == SiCluster::X? fabs(nearHit.x()-centerX.x()):
    //fabs(nearHit.y()-centerX.y()));
    
    //if (deltaX < outRadius) {
	//if (GFtutor::_DATA->hitFlagged(_mGFtrack->m_axis, indexhit)) done = false;
    //} else indexhit = -1; // outside region 
    
    // this condition is necessary for the large angle pattern recognition
    if (indexhit > 0 ) 
    {
	  if (GFtutor::okClusterSize(m_axis,indexhit,slope) == 0) done = false;
    }
    
    if (done == false) 
    {
	  indexhit = -1;
	  //inerRadius = deltaX + 0.1*GFtutor::siResolution();
	  inerRadius = hitDelta + 0.1*GFtutor::siResolution();
    }
    
    return done;
}
//--------------------------------------------------------
// GFsegment - Utility functions
//--------------------------------------------------------
//########################################################
double GFsegment::getZklayer(enum SiCluster::view axis, int klayer) const
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
    if (m_axis==SiCluster::X)  x0=hitp.getPar().getPosition();
    else  y0=hitp.getPar().getPosition();
    double z0=nextKplane.getZPlane();
    Point centerX(x0,y0,z0);
    
    float tower_width = GFtutor::side;
	+ 2.*Tower::wall_thickness + Glast::wall_gap;
    float x_twr_bnd = fmod((centerX.x() + tower_width*2.5), tower_width);
    x_twr_bnd = .5*tower_width - fabs(x_twr_bnd -tower_width*.5);
    float y_twr_bnd = fmod((centerX.y() + tower_width*2.5), tower_width);
    y_twr_bnd = .5*tower_width - fabs(y_twr_bnd - tower_width*.5);
    float twr_bnd = (m_axis == SiCluster::X) ? x_twr_bnd : y_twr_bnd;
    if(twr_bnd < .5) crack = true;
    */
    return crack;
}
