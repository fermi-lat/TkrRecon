#include "TkrPatCandRep.h"
#include "TkrRecon/MonteCarlo/TkrEventModel.h"
#include "Event/TopLevel/EventModel.h"

#include "TkrRecon/MonteCarlo/McEventStructure.h"
#include "TkrRecon/MonteCarlo/McLayerHit.h"

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

//#############################################################################
TkrPatCandRep::TkrPatCandRep(IDataProviderSvc* dataProviderSvc)
//#############################################################################
{
    dps = dataProviderSvc;
}
//-------------------- private ----------------------
//##############################################
void TkrPatCandRep::update()
//##############################################
{
    Event::TkrPatCandCol* pTkrCandidates = SmartDataPtr<Event::TkrPatCandCol>(dps,EventModel::TkrRecon::TkrPatCandCol);

    // If candidates exist, draw them!
    if (pTkrCandidates)
    {
        // Loop over the candidate tracks
        for(Event::TkrPatCandColPtr patCandIter = pTkrCandidates->begin(); patCandIter != pTkrCandidates->end(); patCandIter++)
        {
            // For now, draw all candidates in green
            drawTrack(*patCandIter, "green");
        }
    }

    return;
}

void TkrPatCandRep::drawTrack(Event::TkrPatCand* patCand, const std::string& color)
{
    gui::DisplayRep* pDisplay = this;

    // Set the starting point for drawing the candidate track
    Point x0        = patCand->getPosition();
    Point x1        = x0;
    int   lastLayer = -1;

    for(Event::CandHitVectorPtr hitIter = patCand->begin(); hitIter != patCand->end(); hitIter++)
    {
        // Pointer to candidate hit
        Event::TkrPatCandHit* candHit = *hitIter;
        // Point to move to
        Point xCur     = candHit->Position();
        int   curLayer = candHit->PlaneIndex();

        if (curLayer == lastLayer)
        {
            // Set the point to draw to
            if (candHit->View() == Event::TkrCluster::X) x1.setX(xCur.x());
            else                                         x1.setY(xCur.y());
            x1.setZ(0.5*(x1.z()+xCur.z()));

            // Draw a line between the two points
            pDisplay->setColor(color);
            pDisplay->moveTo(x0);
            pDisplay->lineTo(x1); 

            x0 = x1;
        }
        else
        {
            // Set the point to draw to
            if (candHit->View() == Event::TkrCluster::X) x1.setX(xCur.x());
            else                                         x1.setY(xCur.y());
            x1.setZ(xCur.z());

            lastLayer = curLayer;
        }
    }
}
