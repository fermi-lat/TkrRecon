/**
 * @class McBuildTracks
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/MonteCarlo/McBuildTracks.cxx,v 1.1 2003/08/04 20:11:24 usher Exp $
 */
#include "src/MonteCarlo/McBuildTracks.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "Event/TopLevel/EventModel.h"
#include "TkrRecon/MonteCarlo/TkrEventModel.h"

#include "TkrRecon/MonteCarlo/McLayerHit.h"

Event::McBuildTracks::McBuildTracks(IDataProviderSvc* dataSvc)
{
    //Take care of insuring that data area has been created
    //DataObject* pNode = 0;
    StatusCode  sc;
    //sc = dataSvc->retrieveObject(TkrEventModel::MC::Event, pNode);
    //if ( sc.isFailure() ) {
    //    sc = dataSvc->registerObject(TkrEventModel::MC::Event, new DataObject);
    //    if( sc.isFailure() ) {
    //     //   log << MSG::ERROR << "could not register /Event/tmp" << endreq;
    //        return;
    //    }
    //}

    //Define the tables needed for the MC tracks
    Event::McPartToHitTab  mcPartToHitTab;
    Event::ClusToLyrHitTab clusToLyrHitTab;
    Event::McLyrToHitTab   mcLyrToHitTab;

    //Initialize the tables
    mcPartToHitTab.init();
    clusToLyrHitTab.init();
    mcLyrToHitTab.init();

    //Store these objects in the TDS
    sc = dataSvc->registerObject(TkrEventModel::MC::McPartToHitTab,mcPartToHitTab.getAllRelations());
    sc = dataSvc->registerObject(TkrEventModel::MC::McClusToLyrHitTab,clusToLyrHitTab.getAllRelations());
    sc = dataSvc->registerObject(TkrEventModel::MC::McLyrToHitTab,mcLyrToHitTab.getAllRelations());

    //Create the container for any McLayerHits we might create
    Event::McLayerHitCol* mcLayerHits = new Event::McLayerHitCol();
    sc = dataSvc->registerObject(TkrEventModel::MC::McLayerHitCol, mcLayerHits);

    //Recover the McPositionHit to Cluster relational table
    SmartDataPtr<Event::ClusMcPosHitTabList> tkrTable(dataSvc,EventModel::Digi::TkrClusterHitTab);
    Event::ClusMcPosHitTab clusHitTab(tkrTable);

    //If no relations then no reason to continue
    if ((clusHitTab.getAllRelations())->size() == 0) return;
 
    // Get the MC info we want from the TDS
    SmartDataPtr<Event::McPositionHitVector> posHits(dataSvc, EventModel::MC::McPositionHitCol);

    Event::TkrCluster*       lastCluster = 0;
    Event::McLayerHit*       layerHit    = 0;
    const Event::McParticle* curParticle = 0;

    // Loop through McPositionHits to build McLayerHits 
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

        // If an unseen cluster, then create a new McLayerHit, else use a previous one
        if (lastCluster != tkrCluster || curParticle != mcPart)
        {
            curParticle = mcPart;
            layerHit    = 0;

            // First step is to search for an McLayerHit associated with this mcPart 
            // which uses the current cluster
            Event::McPartTrack mcPartTrack = mcPartToHitTab.getRelByFirst(mcPart);
            if (mcPartTrack.size() > 0)
            {
                //Nightmare! Table is not ordered so must search for last layer...
                Event::McPartTrack::const_iterator tempIter;
                for(tempIter = mcPartTrack.begin(); tempIter != mcPartTrack.end(); tempIter++)
                {
                    Event::McLayerHit* lastHit = (*tempIter)->getSecond();
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
                layerHit = new Event::McLayerHit(curParticle);

                layerHit->setTkrCluster(tkrCluster);

                mcLayerHits->push_back(layerHit);

                Event::McPartToHitRel* partHitRel = new Event::McPartToHitRel(const_cast<Event::McParticle*>(mcPart), layerHit);
                mcPartToHitTab.addRelation(partHitRel);

                // Relate this McLayerHit to the TkrCluster and add to table
                Event::ClusToLyrHitRel* clusToLyrRel = new Event::ClusToLyrHitRel(tkrCluster, layerHit);
                clusToLyrHitTab.addRelation(clusToLyrRel);
            }
        }

        // Add the McPositionHit 
        layerHit->addMcPositionHit(mcPosHit);
 
        // Relate this McLayerHit to McPositionHit and add to table
        Event::McLyrToHitRel* lyrToHitRel = new Event::McLyrToHitRel(layerHit, const_cast<Event::McPositionHit*>(mcPosHit));
        mcLyrToHitTab.addRelation(lyrToHitRel);

        lastCluster = tkrCluster;
    }

}

Event::McBuildTracks::~McBuildTracks()
{
    return;
}
