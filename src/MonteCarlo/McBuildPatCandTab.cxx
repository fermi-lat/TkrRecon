/**
 * @class McBuildPatCandTab
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/MonteCarlo/McBuildPatCandTab.cxx,v 1.1 2003/08/04 20:11:24 usher Exp $
 */
#include "McBuildPatCandTab.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"

#include "TkrRecon/MonteCarlo/McLayerHit.h"
#include "TkrRecon/MonteCarlo/McPatCand.h"

Event::McBuildPatCandTab::McBuildPatCandTab(DataSvc* dataSvc)
{
    StatusCode sc = StatusCode::SUCCESS;
 
    // Look up the Pattern Track collection from the TDS
    Event::TkrPatCandCol* pTkrCands = SmartDataPtr<Event::TkrPatCandCol>(dataSvc,EventModel::TkrRecon::TkrPatCandCol);
    Event::TkrClusterCol* pTkrClus = SmartDataPtr<Event::TkrClusterCol>(dataSvc,EventModel::TkrRecon::TkrClusterCol); 

    // Look up the Monte Carlo track tables
    SmartDataPtr<Event::ClusToLyrHitTabList> tkrTable(dataSvc,"/Event/tmp/ClusToLyrHit");
    Event::ClusToLyrHitTab clusToLyrHitTab(tkrTable);
    SmartDataPtr<Event::McPartToHitTabList> partTable(dataSvc,"/Event/tmp/McPartToHit");
    Event::McPartToHitTab mcPartToHitTab(partTable);

    //Create the pattern track relational tables
    Event::PatHitToLyrHitTab  patHitToLyrHitTab;
    Event::PatCandToMcCandTab patCandToMcCandTab;

    patHitToLyrHitTab.init();
    patCandToMcCandTab.init();

    //Store them in the TDS
    sc = dataSvc->registerObject("/Event/tmp/PatHitToLyrHit",patHitToLyrHitTab.getAllRelations());
    sc = dataSvc->registerObject("/Event/tmp/PatCandToMcCand",patCandToMcCandTab.getAllRelations());

    //McPatCandCol defined here, stored in TDS to make sure it gets cleaned up
    Event::McPatCandCol* mcPatCands = new Event::McPatCandCol();
    sc = dataSvc->registerObject("/Event/tmp/McPatCandCol", mcPatCands);

    // Loop over them
    for(Event::TkrPatCandColPtr cands = pTkrCands->begin(); cands != pTkrCands->end(); cands++)
    {
        Event::TkrPatCand* patCand   = *cands;
        Event::McPatCand*  mcPatCand = new Event::McPatCand();

        mcPatCands->push_back(mcPatCand);

        Event::PatCandToMcCandRel* patMcCandRel = new Event::PatCandToMcCandRel(patCand,mcPatCand);
        patCandToMcCandTab.addRelation(patMcCandRel);

        int                     numHits = patCand->numPatCandHits();
        Event::CandHitVectorPtr candPtr = patCand->getHitIterBegin();
        while(numHits--)
        {
            ////Event::TkrPatCandHit candHit = *candPtr++;
            Event::TkrPatCandHit*  candHit = *candPtr++;
            int                    clusIdx = candHit->HitIndex();
            Event::TkrCluster*     cluster = pTkrClus->getHit(clusIdx);
            Event::ClusToLyrHitVec lyrHits = clusToLyrHitTab.getRelByFirst(cluster);

            int                    numLyrHits = lyrHits.size();
            if (numLyrHits > 1)
            {
                numLyrHits = 1;
            }

            for(Event::ClusToLyrHitVec::iterator lyrHitIter = lyrHits.begin(); lyrHitIter != lyrHits.end(); lyrHitIter++)
            {
                Event::ClusToLyrHitRel*  lyrHitRel = *lyrHitIter;
                Event::McLayerHit*       lyrHit    = lyrHitRel->getSecond();
                const Event::McParticle* mcPart    = lyrHit->getMcParticle();

                Event::PatHitToLyrHitRel* patHitRel = new Event::PatHitToLyrHitRel(candHit, lyrHit);
                patHitToLyrHitTab.addRelation(patHitRel);

                Event::PatHitToMcPartRel* patHitMcPartRel = new Event::PatHitToMcPartRel(candHit,const_cast<Event::McParticle*>(mcPart));
                mcPatCand->addRelation(patHitMcPartRel);
            }
        }

        // Finish buiding the McPatCand class
        mcPatCand->fillMcParticleVec();

        // Check number of hits per track
        ///int numMcParts = mcPatCand->getNumMcParticles();
        ///for(std::vector<Event::McParticle*>::iterator mcPartIter  = mcPatCand->begin();
        ///                                              mcPartIter != mcPatCand->end();    mcPartIter++)
        ///{
        ///    Event::McParticle* mcPart = *mcPartIter;
        ///    int numHits = mcPatCand->getNumHits(mcPart);
        ///    int chkum   = 0;
        ///}
    }

    //Quick cross check
    for(Event::TkrPatCandColPtr cands = pTkrCands->begin(); cands != pTkrCands->end(); cands++)
    {
        Event::TkrPatCand*        patCand         = *cands;
        Event::PatCandToMcCandVec patCandToMcVec  = patCandToMcCandTab.getRelByFirst(patCand);
        Event::McPatCand*         mcPatCand       = (patCandToMcVec.front())->getSecond();

        for(std::vector<Event::McParticle*>::iterator mcPartIter  = mcPatCand->begin();
                                                      mcPartIter != mcPatCand->end();    mcPartIter++)
        {
            Event::McParticle*         mcPart  = *mcPartIter;
            Event::McPartTrack         mcTrack = mcPartToHitTab.getRelByFirst(mcPart);

            int numMcHits   = mcTrack.size();
            int numCandHits = mcPatCand->getNumHits(mcPart);
            int nonthng = 0;
        }

        int                     numHits = patCand->numPatCandHits();
        Event::CandHitVectorPtr candPtr = patCand->getHitIterBegin();
        while(numHits--)
        {
            Event::TkrPatCandHit*    candHit      = *candPtr++;
            Event::PatHitToLyrHitVec patHitRelVec = patHitToLyrHitTab.getRelByFirst(candHit);

            int                      numLyrHits   = patHitRelVec.size();
            int mnnn=0;

        }
    }
    return;    
}

Event::McBuildPatCandTab::~McBuildPatCandTab()
{
    return;
}
