/** 
 * @class TkrBuildMcRelationsAlg
 *
 * @brief TkrRecon Gaudi Algorithm for building the Relational tables between the recon information
 *        and the Monte Carlo "truth" information. 
 *        NOTE: This algorithm will only operate in "pruneCal" or "full" McParticle pruning modes
 *
 * Created 6-Aug-2003
 * 
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/GaudiAlg/TkrBuildMcRelationsAlg.h,v 1.7 2002/08/28 22:55:46 usher Exp $
 */

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"

#include "TkrRecon/MonteCarlo/TkrEventModel.h"
#include "TkrRecon/MonteCarlo/McEventStructure.h"
#include "TkrRecon/MonteCarlo/McLayerHit.h"
#include "TkrRecon/MonteCarlo/McPatCand.h"
#include "src/MonteCarlo/McBuildTracks.h"

class TkrBuildMcRelationsAlg : public Algorithm
{
public:
    // Standard Gaudi Algorithm constructor format
    TkrBuildMcRelationsAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrBuildMcRelationsAlg() {}

    // The thee phases in the life of a Gaudi Algorithm
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    void       buildPatCandRelations();
};

static const AlgFactory<TkrBuildMcRelationsAlg>  Factory;
const IAlgFactory& TkrBuildMcRelationsAlgFactory = Factory;

TkrBuildMcRelationsAlg::TkrBuildMcRelationsAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)  
{ 
}

StatusCode TkrBuildMcRelationsAlg::initialize()
{
    // Purpose and Method: Initialization method for the Tracker MC relations algorithm
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: none
    // Restrictions and Caveats:  None

    MsgStream log(msgSvc(), name());

    setProperties();
   
    return StatusCode::SUCCESS;
}

StatusCode TkrBuildMcRelationsAlg::execute()
{
    // Purpose and Method: Called each event to build first the tables relating McParticles, McPositionHits
    //                     and Tkr Clusters to form MC tracks, then to relate this information to the 
    //                     Tracker Recon output.
    // Inputs:  None
    // Outputs:  Relations tables in the Tracker TDS, StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // Retrieve the pointer to the McEventStructure
    SmartDataPtr<Event::McEventStructure> mcEvent(eventSvc(),TkrEventModel::MC::McEventStructure);

    // If it doesn't exist then we need to build the MC structure
    if (mcEvent == 0)
    {
        //stuff
        IParticlePropertySvc* ppsvc;
        if( (sc = service("ParticlePropertySvc", ppsvc)).isFailure() ) 
        {
            throw GaudiException("Service [ParticlePropertySvc] not found", name(), sc);
        }

        DataObject* pNode;
        sc = eventSvc()->retrieveObject(TkrEventModel::MC::Event, pNode);
        if ( sc.isFailure() ) 
        {
            sc = eventSvc()->registerObject(TkrEventModel::MC::Event, new DataObject);
            if( sc.isFailure() ) 
            {
                throw GaudiException("Unable to establish the Tracker MC TDS section", name(), sc);
            }
        }

        //This builds the Monte Carlo event structure - basically a description of the event
        mcEvent = new Event::McEventStructure(eventSvc(), ppsvc);

        //Store in the Tkr TDS 
        sc = eventSvc()->registerObject(TkrEventModel::MC::McEventStructure, mcEvent);
        if (sc.isFailure())
        {
            throw GaudiException("Cannot store the McEventStructure in the TDS", name(), sc);
        }

        //This builds the Monte Carlo tracks
        Event::McBuildTracks eventTracks(eventSvc());
    }

    //Have the Monte Carlo base info, now relate to the recon info
    buildPatCandRelations();
        
    return sc;
}


StatusCode TkrBuildMcRelationsAlg::finalize()
{   
    return StatusCode::SUCCESS;
}

void TkrBuildMcRelationsAlg::buildPatCandRelations()
{
    StatusCode sc = StatusCode::SUCCESS;

    // Look up the Pattern Track collection from the TDS
    Event::TkrPatCandCol* pTkrCands = SmartDataPtr<Event::TkrPatCandCol>(eventSvc(),EventModel::TkrRecon::TkrPatCandCol);
    Event::TkrClusterCol* pTkrClus  = SmartDataPtr<Event::TkrClusterCol>(eventSvc(),EventModel::TkrRecon::TkrClusterCol); 

    // Look up the Monte Carlo track tables
    SmartDataPtr<Event::ClusToLyrHitTabList> tkrTable(eventSvc(),TkrEventModel::MC::McClusToLyrHitTab);
    Event::ClusToLyrHitTab clusToLyrHitTab(tkrTable);
    SmartDataPtr<Event::McPartToHitTabList> partTable(eventSvc(),TkrEventModel::MC::McPartToHitTab);
    Event::McPartToHitTab mcPartToHitTab(partTable);

    //Create the pattern track relational tables
    Event::PatHitToLyrHitTab  patHitToLyrHitTab;
    Event::PatCandToMcCandTab patCandToMcCandTab;

    patHitToLyrHitTab.init();
    patCandToMcCandTab.init();

    //Store them in the TDS
    sc = eventSvc()->registerObject(TkrEventModel::MC::PatHitToLyrHit,patHitToLyrHitTab.getAllRelations());
    if( sc.isFailure() ) 
    {
        throw GaudiException("Cannot store PatHitToLyrHit relations in the TDS!", name(), sc);
    }

    sc = eventSvc()->registerObject(TkrEventModel::MC::PatCandToMcCand,patCandToMcCandTab.getAllRelations());
    if( sc.isFailure() ) 
    {
        throw GaudiException("Cannot store PatCandToMcCand relations in the TDS!", name(), sc);
    }

    //McPatCandCol defined here, stored in TDS to make sure it gets cleaned up
    Event::McPatCandCol* mcPatCands = new Event::McPatCandCol();
    sc = eventSvc()->registerObject(TkrEventModel::MC::McPatCandCol, mcPatCands);
    if( sc.isFailure() ) 
    {
        throw GaudiException("Cannot store the McPatCandCol in the TDS!", name(), sc);
    }

    // Loop over them
    for(Event::TkrPatCandColPtr cands = pTkrCands->begin(); cands != pTkrCands->end(); cands++)
    {
        Event::TkrPatCand* patCand   = *cands;
        Event::McPatCand*  mcPatCand = new Event::McPatCand();

        mcPatCands->push_back(mcPatCand);

        Event::PatCandToMcCandRel* patMcCandRel = new Event::PatCandToMcCandRel(patCand,mcPatCand);
        patCandToMcCandTab.addRelation(patMcCandRel);

        for(Event::CandHitVectorPtr candPtr = patCand->begin(); candPtr != patCand->end(); candPtr++)
        {
            Event::TkrPatCandHit*  candHit = *candPtr;
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
}

