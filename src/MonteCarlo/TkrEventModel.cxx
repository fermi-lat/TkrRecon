// File and Version Information:
// $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/MonteCarlo/TkrEventModel.cxx,v 1.2 2003/08/06 21:53:21 usher Exp $

#define _TkrEventModel_CPP_


#include "TkrRecon/MonteCarlo/TkrEventModel.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/ClassID.h"

    
/** @class EvModel
 *  @brief Event Model: Definition of logical paths and class identifiers
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/MonteCarlo/TkrEventModel.cxx,v 1.2 2003/08/06 21:53:21 usher Exp $
 */
class TkrEvModel {
        
public:
    
    TkrEvModel() {
        // Access to GLAST event
        TkrEventModel::EventHeader           = "/Event";
            
        // Monte Carlo 
        TkrEventModel::MC::Event             = TkrEventModel::EventHeader + "/tmp";
//        TkrEventModel::MC::McEventStructure  = TkrEventModel::MC::Event  + "/McEventStructure";

//        TkrEventModel::MC::McPartToHitTab    = TkrEventModel::MC::Event  + "/McPartToHitTab";
//        TkrEventModel::MC::McClusToLyrHitTab = TkrEventModel::MC::Event  + "/McClusToLyrHitTab";
//        TkrEventModel::MC::McLyrToHitTab     = TkrEventModel::MC::Event  + "/McLyrToHitTab";
//        TkrEventModel::MC::McSiLayerHitCol   = TkrEventModel::MC::Event  + "/McSiLayerHitCol";

        TkrEventModel::MC::PatHitToLyrHit    = TkrEventModel::MC::Event  + "/PatHitToLyrHit";
        TkrEventModel::MC::PatCandToMcCand   = TkrEventModel::MC::Event  + "/PatCandToMcCand";
        TkrEventModel::MC::McPatCandCol      = TkrEventModel::MC::Event  + "/McPatCandCol";
    }
};
    
    
static TkrEvModel TkrMod;    // where  used? has file scope     
    
