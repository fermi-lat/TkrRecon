#include "Event/TopLevel/EventModel.h"
#include "TkrCandidate3DRep.h"

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

const char col3d_blue[]       = "blue";
const char col3d_violet[]     = "violet";
const char col3d_turquoise[]  = "turquoise";
const char col3d_orange[]     = "orange";
const char col3d_maroon[]     = "maroon";
const char col3d_aquamarine[] = "aquamarine";

const char* p3dColors[] = {col3d_blue,   col3d_violet, col3d_turquoise,
                           col3d_orange, col3d_maroon, col3d_aquamarine};

//#############################################################################
TkrCandidate3DRep::TkrCandidate3DRep(IDataProviderSvc* dataProviderSvc, ITkrGeometrySvc* tkrGeom)
//#############################################################################
{
    dps     = dataProviderSvc;
    m_tkrGeom = tkrGeom;
}
//-------------------- private ----------------------
//##############################################
void TkrCandidate3DRep::update()
//##############################################
{
    Event::TkrTrackCol* pTkrCandidates = SmartDataPtr<Event::TkrTrackCol>(dps,EventModel::TkrRecon::TkrTrackCol);

    //Now see if we can do the drawing
    if (pTkrCandidates)
    {
        /*
        int                     numCandTracks = pTkrCandidates->size();
        Event::TkrTrackColPtr cands         = pTkrCandidates->begin();
        int                     colorIdx      = 0;

    //        gui::DisplayRep* pDisplay = this;

        while(numCandTracks--)
        {
            Event::TkrTrack* pTkrCand = *cands++;

            //Put a marker at the start of the candidate
            Point  strtPoint = pTkrCand->getPosition();
            double x         = strtPoint.x();
            double y         = strtPoint.y();
            double zCur      = strtPoint.z() + 100.; //mm
            double zPrev     = strtPoint.z() + 1000.; //mm

            setColor(p3dColors[colorIdx]);
            markerAt(strtPoint);

            int                        numHits = pTkrCand->numPatCandHits();
            Event::CandHitVectorPtr hitPtr  = pTkrCand->getHitIterBegin();

            while(numHits--)
            {
                Event::TkrTrackHit* pHitCand = *hitPtr++;
                Point                 hitCoord = pHitCand->Position();

                if (pHitCand->View() == idents::TkrId::eMeasureX) x = hitCoord.x();
                else                                              y = hitCoord.y();
  
                zCur = hitCoord.z();

                if (fabs(zPrev-zCur) < 2.5) //mm
                {
                    Point newPoint(x,y,0.5*(zCur+zPrev));
                    moveTo(strtPoint);
                    lineTo(newPoint);

                    strtPoint = newPoint;
                }

                zPrev = zCur;
            }

            colorIdx = (colorIdx + 1) % 6;
        }

        setColor("blue");
    */
    }

    return;
}

