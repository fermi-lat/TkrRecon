#include "TkrRecon/Display/TkrTracksRep.h"
#include "Event/TopLevel/EventModel.h"

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

using namespace Event;

//#############################################################################
TkrTracksRep::TkrTracksRep(IDataProviderSvc* dataProviderSvc)
//#############################################################################
{
	dps = dataProviderSvc;
}
//-------------------- private ----------------------
//##############################################
void TkrTracksRep::update()
//##############################################
{
    TkrFitTrackCol* pTracks = SmartDataPtr<TkrFitTrackCol>(dps,EventModel::TkrRecon::TkrFitTrackCol);

	//Now see if we can do the drawing
	if (pTracks)
	{
        gui::DisplayRep* pDisplay = this;

        pDisplay->setColor("blue");

        int numTracks = pTracks->size();

        if (numTracks > 0) 
        {
            int trkIdx = 0;
	        TkrFitTrackCol::const_iterator it = pTracks->begin();

            //for(it = pTracks->begin(); it != pTracks->end(); ++it)
            while(it < pTracks->end())
            {
                const TkrFitTrack& track = **it++;

                pDisplay->markerAt(track.getPosition());

                drawChiSq(track);
                drawTrack(track);
            }
        }
	}

    return;
}

void TkrTracksRep::drawChiSq(const TkrFitTrack& track)
{
    gui::DisplayRep* pDisplay = this;
    TkrFitHit::TYPE  fit      = TkrFitHit::SMOOTH;
    TkrFitHit::TYPE  typ      = TkrFitHit::SMOOTH;

    TkrFitPlaneConPtr hitIter = track.getHitIterBegin();

    while(hitIter < track.getHitIterEnd())
    {
        TkrFitPlane plane = *hitIter++;

        TkrFitPlane::AXIS prj = plane.getProjection();

        double x0, y0, z0, xl, xr, yl, yr;
        double delta= plane.getDeltaChiSq(typ);

        if(prj == TkrCluster::X){
            x0 = plane.getHit(typ).getPar().getXPosition();
            y0 = plane.getHit(fit).getPar().getYPosition(); 
            z0 = plane.getZPlane()+0.1;
            xl = x0-0.5*delta;
            xr = x0+0.5*delta;
            yl = y0;
            yr = y0;
        } 
        else {
            x0 = plane.getHit(fit).getPar().getXPosition();
            y0 = plane.getHit(typ).getPar().getYPosition(); 
            z0 = plane.getZPlane()+0.1;
            xl = x0;
            xr = x0;
            yl = y0-0.5*delta;
            yr = y0+0.5*delta;
        }		
        pDisplay->moveTo(Point(xl,yl,z0));
        pDisplay->lineTo(Point(xr,yr,z0));
    }
}

void TkrTracksRep::drawTrack(const TkrFitTrack& track)
{
    gui::DisplayRep* pDisplay = this;
    TkrFitHit::TYPE  fit      = TkrFitHit::SMOOTH;
    TkrFitHit::TYPE  typ      = TkrFitHit::SMOOTH;

    TkrFitPlaneConPtr hitIter = track.getHitIterBegin();

    while(hitIter < track.getHitIterEnd()-1)
    {
        TkrFitPlane plane     = *hitIter++;
        TkrFitPlane planeNext = *hitIter;

        TkrFitPlane::AXIS prj = plane.getProjection();

        double x0, y0, z0;

		TkrFitHit::TYPE xtyp, ytyp;
		xtyp = (prj == TkrCluster::X ? typ : fit);
		ytyp = (prj == TkrCluster::X ? fit : typ);

		// this sets up the track segment to the next plane
    
        x0 = plane.getHit(xtyp).getPar().getXPosition();
        y0 = plane.getHit(ytyp).getPar().getYPosition(); 
        z0 = plane.getZPlane();
 
		double tanx = plane.getHit(typ).getPar().getXSlope();
        double tany = plane.getHit(typ).getPar().getYSlope();
        
        Point origin(x0,y0,z0);
        Vector dir = Vector(-1.*tanx,-1.*tany,-1.).unit();
        
        Ray segment(origin,dir);
        double zstep=planeNext.getZPlane()-z0;
        double cosz=dir.z();

        // this sets up the dotted line from the lower part of the extrapolated track
		//  to the next hit.


		prj = plane.getNextProj();

		xtyp = (prj == TkrCluster::X ? typ : fit);
		ytyp = (prj == TkrCluster::X ? fit : typ);

        x0 = planeNext.getHit(xtyp).getPar().getXPosition();
        y0 = planeNext.getHit(ytyp).getPar().getYPosition(); 
        z0 = planeNext.getZPlane();

        Point p(x0, y0, z0);

		// do them in this order, so that the connection doesn't cover the track
		
		pDisplay->set_line_style(1);
        pDisplay->moveTo(segment.position(0.8*zstep/cosz));
        pDisplay->lineTo(p); 

		pDisplay->setColor("blue");
        pDisplay->moveTo(segment.position(0.));
        pDisplay->lineTo(segment.position(zstep/cosz));
    }
}
