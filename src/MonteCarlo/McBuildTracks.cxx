 /**
 * @class McBuildTracks
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/MonteCarlo/McBuildTracks.cxx,v 1.5 2003/08/21 02:42:32 usher Exp $
 */
#include "McBuildTracks.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "Event/TopLevel/EventModel.h"
#include "TkrRecon/MonteCarlo/TkrEventModel.h"

#include "Event/MonteCarlo/McSiLayerHit.h"

Event::McBuildTracks::McBuildTracks(IDataProviderSvc* dataSvc)
{
    //Take care of insuring that data area has been created
    StatusCode  sc;

    //Define the tables needed for the MC tracks
    Event::McPartToHitTab  mcPartToHitTab;
    Event::ClusToLyrHitTab clusToLyrHitTab;
    Event::McLyrToHitTab   mcLyrToHitTab;

    //Initialize the tables
    mcPartToHitTab.init();
    clusToLyrHitTab.init();
    mcLyrToHitTab.init();

    //Store these objects in the TDS
    sc = dataSvc->registerObject(EventModel::MC::McPartToHitTab,mcPartToHitTab.getAllRelations());
    sc = dataSvc->registerObject(EventModel::MC::McClusToLyrHitTab,clusToLyrHitTab.getAllRelations());
    sc = dataSvc->registerObject(EventModel::MC::McLyrToHitTab,mcLyrToHitTab.getAllRelations());

    //Create the container for any McSiLayerHits we might create
    Event::McSiLayerHitCol* McSiLayerHits = new Event::McSiLayerHitCol();
    sc = dataSvc->registerObject(EventModel::MC::McSiLayerHitCol, McSiLayerHits);

    //Recover the McPositionHit to Cluster relational table
    SmartDataPtr<Event::ClusMcPosHitTabList> tkrTable(dataSvc,EventModel::Digi::TkrClusterHitTab);
    Event::ClusMcPosHitTab clusHitTab(tkrTable);

    //If no relations then no reason to continue
    if ((clusHitTab.getAllRelations())->size() == 0) return;
 
    // Get the MC info we want from the TDS
    SmartDataPtr<Event::McPositionHitVector> posHits(dataSvc, EventModel::MC::McPositionHitCol);

    Event::TkrCluster*       lastCluster = 0;
    Event::McSiLayerHit*     layerHit    = 0;
    const Event::McParticle* curParticle = 0;

    // Loop through McPositionHits to build McSiLayerHits 
    Event::McPositionHitVector::const_iterator hit;
    for (hit = posHits->begin(); hit != posHits->end(); hit++ ) 
    {
        const Event::McPositionHit*    mcPosHit = *hit;
        const Event::McParticle*       mcPart   = mcPosHit->mcParticle();
        const idents::VolumeIdentifier curVolId = mcPosHit->volumeID();

        // Only interested in volumes inside the tracker
        if (curVolId[0] != 0) continue;

        // Find the cluster associated with this McPositionHit
        std::vector<Event::ClusMcPosHitRel*> clusRelVector = clusHitTab.getRelBySecond(mcPosHit);

        // If no clusters associated to this McPositionHit, then skip
        int numRels = clusRelVector.size();
        //if (numRels == 0) continue;

        Event::TkrCluster* tkrCluster = numRels > 0 ? clusRelVector[0]->getFirst() : 0;

        // Can this happen?
        if (numRels > 1)
        {
            int mmm = -1;
        }

        // Just checking here
        int trayNum = curVolId[4];
        int botTop  = curVolId[6];
        int view    = curVolId[5];
        int layer   = 2*trayNum - 1 + botTop;

        // If an unseen cluster, then create a new McSiLayerHit, else use a previous one
        if (lastCluster != tkrCluster || curParticle != mcPart)
        {
            curParticle = mcPart;
            layerHit    = 0;

            // First step is to search for an McSiLayerHit associated with this mcPart 
            // which uses the current cluster
            Event::McPartTrack mcPartTrack = mcPartToHitTab.getRelByFirst(mcPart);
            if (mcPartTrack.size() > 0)
            {
                //Nightmare! Table is not ordered so must search for last layer...
                Event::McPartTrack::const_iterator tempIter;
                for(tempIter = mcPartTrack.begin(); tempIter != mcPartTrack.end(); tempIter++)
                {
                    Event::McSiLayerHit* lastHit = (*tempIter)->getSecond();
                    const Event::TkrCluster* vecClus = lastHit->getTkrCluster();

                    // Did we find one?
                    if (vecClus == tkrCluster)
                    {
                        layerHit = lastHit;
                        break;
                    }
                }
            }

            //No previous layerHit, so make a new one
            if (!layerHit) 
            {
                layerHit = new Event::McSiLayerHit(curParticle);

                layerHit->setTkrCluster(tkrCluster);

                McSiLayerHits->push_back(layerHit);

                // Relate the current McParticle to this McSiLayerHit
                Event::McPartToHitRel* partHitRel = new Event::McPartToHitRel(const_cast<Event::McParticle*>(mcPart), layerHit);
                mcPartToHitTab.addRelation(partHitRel);

                // Relate this McSiLayerHit to the TkrCluster and add to table
                Event::ClusToLyrHitRel* clusToLyrRel = new Event::ClusToLyrHitRel(tkrCluster, layerHit);
                clusToLyrHitTab.addRelation(clusToLyrRel);

                // Take the time to check if a cluster is shared by different McSiLayerHits
                Event::ClusToLyrHitVec clusHitVec = clusToLyrHitTab.getRelByFirst(tkrCluster);
                
                if (clusHitVec.size() > 1)
                {
                    //Go through and set the "shared" bits in the McSiLayerHits
                    for(Event::ClusToLyrHitVec::iterator clusHitVecIter = clusHitVec.begin();
                        clusHitVecIter != clusHitVec.end(); clusHitVecIter++)
                    {
                        (*clusHitVecIter)->getSecond()->setStatusBit(Event::McSiLayerHit::SHAREDCLUS);
                    }
                }
            }
        }

        // Add the McPositionHit 
        layerHit->addMcPositionHit(mcPosHit);
 
        // Relate this McSiLayerHit to McPositionHit and add to table
        Event::McLyrToHitRel* lyrToHitRel = new Event::McLyrToHitRel(layerHit, const_cast<Event::McPositionHit*>(mcPosHit));
        mcLyrToHitTab.addRelation(lyrToHitRel);

        lastCluster = tkrCluster;
    }

}

Event::McBuildTracks::~McBuildTracks()
{
    return;
}
