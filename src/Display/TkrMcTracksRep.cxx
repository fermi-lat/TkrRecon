#include "TkrRecon/Display/TkrMcTracksRep.h"
#include "TkrRecon/MonteCarlo/TkrEventModel.h"
#include "Event/TopLevel/EventModel.h"

#include "TkrRecon/MonteCarlo/McEventStructure.h"
#include "TkrRecon/MonteCarlo/McLayerHit.h"

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

//#############################################################################
TkrMcTracksRep::TkrMcTracksRep(IDataProviderSvc* dataProviderSvc)
//#############################################################################
{
    dps = dataProviderSvc;
}
//-------------------- private ----------------------
//##############################################
void TkrMcTracksRep::update()
//##############################################
{
    SmartDataPtr<Event::McEventStructure> mcEvent(dps,TkrEventModel::MC::McEventStructure);

    // no McEvent, no picture
    if (mcEvent)
    {
        // If the primary is charged we draw it first
        if (mcEvent->getClassificationBits() & Event::McEventStructure::CHARGED) drawTrack(mcEvent->getPrimaryParticle(),"blue");

        // Now draw the secondaries
        Event::McParticleRefVec::const_iterator partIter;

        for(partIter = mcEvent->beginSecondaries(); partIter != mcEvent->endSecondaries(); partIter++)
        {
            drawTrack(*partIter, "red");
        }

        for(partIter = mcEvent->beginAssociated(); partIter != mcEvent->endAssociated(); partIter++)
        {
            drawTrack(*partIter, "brown");
        }
    }

    return;
}

/*! A small class to use the sort algorithm */
class CompareTrackHits 
{
public:
    bool operator()(Event::McPartToHitRel *left, Event::McPartToHitRel *right)
    {
        bool                     leftTest   = false;

        const Event::McLayerHit* mcHitLeft  = left->getSecond();
        const Event::McLayerHit* mcHitRight = right->getSecond();

        // Find the McPositionHit associated with the McParticle
        const Event::McPositionHit* mcPosHitLeft  = findMcPosHit(mcHitLeft);
        const Event::McPositionHit* mcPosHitRight = findMcPosHit(mcHitRight);
        
        if (mcPosHitLeft && mcPosHitRight)
        {
            leftTest = mcPosHitLeft->timeOfFlight() < mcPosHitRight->timeOfFlight();
        }

        return leftTest;
    }
private:
    const Event::McPositionHit* findMcPosHit(const Event::McLayerHit* mcLayerHit)
    {
        const Event::McParticle*    mcPart = mcLayerHit->getMcParticle();
        const Event::McPositionHit* mcHit  = 0;

        const SmartRefVector<Event::McPositionHit>* mcPosHitVec = mcLayerHit->getMcPositionHitsVec();
        for(SmartRefVector<Event::McPositionHit>::const_iterator hitIter  = mcPosHitVec->begin();
                                                                 hitIter != mcPosHitVec->end(); hitIter++)
        {
            if ((*hitIter)->mcParticle() == mcPart)
            {
                mcHit = *hitIter;
                break;
            }
        }

        return mcHit;
    }
};

void TkrMcTracksRep::drawTrack(const Event::McParticle* mcPart, const std::string& color)
{
    gui::DisplayRep* pDisplay = this;

    // Retrieve the McParticle to hit relational table
    SmartDataPtr<Event::McPartToHitTabList> hitTable(dps,TkrEventModel::MC::McPartToHitTab);
    Event::McPartToHitTab mcPartToHitTab(hitTable);

    // Find the hits associated with this mcPart
    Event::McPartToHitVec hitVec = mcPartToHitTab.getRelByFirst(mcPart);
    std::sort(hitVec.begin(),hitVec.end(),CompareTrackHits());

    // Set line segment start and end points
    Point x0(mcPart->initialPosition().x(),mcPart->initialPosition().y(),mcPart->initialPosition().z());

    Event::McPartToHitVec::const_iterator hitIter;
    for(hitIter = hitVec.begin(); hitIter != hitVec.end(); hitIter++)
    {
        Event::McPartToHitRel*   mcHitRel = *hitIter;
        Event::McLayerHit*       lyrHit   =  mcHitRel->getSecond();

        const SmartRefVector<Event::McPositionHit>* mcPosHitVec = lyrHit->getMcPositionHitsVec();
        SmartRefVector<Event::McPositionHit>::const_iterator mcPosHitIter = mcPosHitVec->begin();

        const Event::McPositionHit* mcPosHit = *mcPosHitIter++;

        double sFac = 1.;
        double xPos = 0.5 * (mcPosHit->globalEntryPoint().x() + mcPosHit->globalExitPoint().x());
        double yPos = 0.5 * (mcPosHit->globalEntryPoint().y() + mcPosHit->globalExitPoint().y());
        double zPos = 0.5 * (mcPosHit->globalEntryPoint().z() + mcPosHit->globalExitPoint().z());

        for( ; mcPosHitIter != mcPosHitVec->end(); mcPosHitIter++)
        {
            xPos += 0.5 * (mcPosHit->globalEntryPoint().x() + mcPosHit->globalExitPoint().x());
            yPos += 0.5 * (mcPosHit->globalEntryPoint().y() + mcPosHit->globalExitPoint().y());
            zPos += 0.5 * (mcPosHit->globalEntryPoint().z() + mcPosHit->globalExitPoint().z());
            sFac += 1.;
        }

        Point  x1(xPos/sFac,yPos/sFac,zPos/sFac);

        // do them in this order, so that the connection doesn't cover the track
        
        pDisplay->setColor(color);
        pDisplay->moveTo(x0);
        pDisplay->lineTo(x1); 

        x0 = x1;
    }
}
