/**
 * @class MonteCarloFindTrackTool
 *
 * @brief Implements a Gaudi Tool for find track candidates. This particular tool uses the 
 *        "MonteCarlo" method, described in detail in TkrMonteCarloPatRec, which is very similar to 
 *        the method employed by Glast up to the time of the PDR.
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/MonteCarlo/MonteCarloFindTrackTool.h,v 1.7 2003/07/29 15:08:01 cohen Exp $
 */

#ifndef MONTECARLOFINDTRACKTOOL_H
#define MONTECARLOFINDTRACKTOOL_H

#include "GaudiKernel/IParticlePropertySvc.h"
#include "src/PatRec/PatRecBaseTool.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"

class MonteCarloFindTrackTool : public PatRecBaseTool 
{
public:
    /// Standard Gaudi Tool interface constructor
    MonteCarloFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~MonteCarloFindTrackTool() {}

    /// @brief Method to find candidate tracks. Will retrieve the necessary information from
    ///        the TDS, including calorimeter energy, and then use TkrMonteCarloPatRec to find all
    ///        possible track candidates. The resulting track candidate collection is then 
    ///        stored in the TDS for the next stage.
	
    /// put actual init stuff here
    StatusCode initialize();
    /// does the work
    StatusCode findTracks();

private:
    Event::TkrPatCand* buildTrack(const Event::McParticle* mcPart);

    IParticlePropertySvc* m_ppsvc;
};

#endif
