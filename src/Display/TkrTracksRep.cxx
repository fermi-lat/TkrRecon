#include "src/Display/TkrTracksRep.h"
#include "Event/TopLevel/EventModel.h"
#include "geometry/Ray.h"

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
    TkrTrackCol* pTracks = SmartDataPtr<TkrTrackCol>(dps,EventModel::TkrRecon::TkrTrackCol);

    //Now see if we can do the drawing
    if (pTracks)
    {
        gui::DisplayRep* pDisplay = this;

        pDisplay->setColor("blue");

        int numTracks = pTracks->size();

        if (numTracks > 0) 
        {
      //int trkIdx = 0;
            TkrTrackCol::const_iterator it = pTracks->begin();

            while(it != pTracks->end())
            {
                const TkrTrack& track = **it++;

                pDisplay->markerAt(track.getInitialPosition());

                drawChiSq(track);
                drawTrack(track);
            }
        }
    }

    return;
}

void TkrTracksRep::drawChiSq(const TkrTrack& track)
{
    gui::DisplayRep* pDisplay = this;
    TkrTrackHit::ParamType  fit      = TkrTrackHit::SMOOTHED;
    TkrTrackHit::ParamType  typ      = TkrTrackHit::SMOOTHED;

    TkrTrackHitVecConItr hitIter = track.begin();

    while(hitIter < track.end())
    {
        SmartRef<TkrTrackHit> plane = *hitIter++;
        idents::TkrId tkrId = plane->getTkrId();

        int prj = tkrId.getView();

        double x0, y0, z0, xl, xr, yl, yr;
        double delta= plane->getChiSquareSmooth()*10.; //Scale factor! We're in mm now!
        Point pos = plane->getPoint(TkrTrackHit::SMOOTHED);
        x0 = pos.x();
        y0 = pos.y(); 
        z0 = plane->getZPlane()+0.1;

        if(prj == idents::TkrId::eMeasureX){
            xl = x0-0.5*delta;
            xr = x0+0.5*delta;
            yl = y0;
            yr = y0;
        } 
        else {
            xl = x0;
            xr = x0;
            yl = y0-0.5*delta;
            yr = y0+0.5*delta;
        }       
        pDisplay->moveTo(Point(xl,yl,z0));
        pDisplay->lineTo(Point(xr,yr,z0));
    }
}

void TkrTracksRep::drawTrack(const TkrTrack& track)
{
    gui::DisplayRep* pDisplay = this;
    TkrTrackHit::ParamType  fit      = TkrTrackHit::SMOOTHED;
    TkrTrackHit::ParamType  typ      = TkrTrackHit::SMOOTHED;

    TkrTrackHitVecConItr hitIter = track.begin();

    while(hitIter < track.end()-1)
    {
        SmartRef<TkrTrackHit> plane     = *hitIter++;
        SmartRef<TkrTrackHit> nextPlane = *hitIter;
        int prj = (plane->getTkrId()).getView();

        double x0, y0, z0;

        TkrTrackHit::ParamType xtyp, ytyp;
        xtyp = (prj == idents::TkrId::eMeasureX ? typ : fit);
        ytyp = (prj == idents::TkrId::eMeasureX ? fit : typ);

        // this sets up the track segment to the next plane
    
        Point pos = plane->getPoint(TkrTrackHit::SMOOTHED);
        x0 = pos.x();
        y0 = pos.y(); 
        z0 = pos.z();
 
        TkrTrackParams params = plane->getTrackParams(TkrTrackHit::SMOOTHED);
        double tanx = params.getxSlope();
        double tany = params.getySlope();
        
        Point origin(x0,y0,z0);
        Vector dir = Vector(-1.*tanx,-1.*tany,-1.).unit();
        
        Ray segment(origin,dir);
        double zstep=nextPlane->getZPlane()-z0;
        double cosz=dir.z();

        // this sets up the dotted line from the lower part of the extrapolated track
        //  to the next hit.

        prj = (nextPlane->getTkrId()).getView();

        xtyp = (prj == idents::TkrId::eMeasureX ? typ : fit);
        ytyp = (prj == idents::TkrId::eMeasureX ? fit : typ);

        pos = nextPlane->getPoint(TkrTrackHit::SMOOTHED);
        x0 = pos.x();
        y0 = pos.y(); 
        z0 = pos.z();

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
