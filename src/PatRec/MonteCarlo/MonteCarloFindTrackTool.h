/**
 * @class MonteCarloFindTrackTool
 *
 * @brief Implements a Gaudi Tool for find track candidates. This particular tool uses the 
 *        "MonteCarlo" method. Here, McPositionHits are related to their originating McParticles to
 *        form Monte Carlo tracks. The McPositionHits are then related to their corresponding Tracker
 *        cluster hits and these are then used to form Pattern Candidate tracks. 
 *        The aim of this method is to allow downstream testing of fitting, vertexing and analysis 
 *        assuming "perfect" knowledge of the pattern recognition.
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/MonteCarlo/MonteCarloFindTrackTool.h,v 1.1 2003/08/04 20:17:35 usher Exp $
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
	
    /// @brief Intialization of the tool
    StatusCode initialize();
    /// @brief Method to association the Monte Carlo hits into Pattern Candidate tracks
    StatusCode findTracks();

private:
    /// private method to build an individual Monte Carlo track
    Event::TkrPatCand* buildTrack(const Event::McParticle* mcPart);

    IParticlePropertySvc* m_ppsvc;
};

#endif
