// File and Version Information:
// $Header: /nfs/slac/g/glast/ground/cvs/Event/src/TopLevel/EventModel.cpp,v 1.55 2003/03/12 23:21:44 usher Exp $

#define _TkrEventModel_CPP_


#include "TkrRecon/MonteCarlo/TkrEventModel.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/ClassID.h"

    
/** @class EvModel
 *  @brief Event Model: Definition of logical paths and class identifiers
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Event/src/TopLevel/EventModel.cpp,v 1.55 2003/03/12 23:21:44 usher Exp $
 */
class TkrEvModel {
        
public:
    
    TkrEvModel() {
        // Access to GLAST event
        TkrEventModel::EventHeader           = "/Event";
            
        // Monte Carlo 
        TkrEventModel::MC::Event             = TkrEventModel::EventHeader + "/tmp";
        TkrEventModel::MC::McEventStructure  = TkrEventModel::MC::Event  + "/McEventStructure";

        TkrEventModel::MC::McPartToHitTab    = TkrEventModel::MC::Event  + "/McPartToHitTab";
        TkrEventModel::MC::McClusToLyrHitTab = TkrEventModel::MC::Event  + "/McClusToLyrHitTab";
        TkrEventModel::MC::McLyrToHitTab     = TkrEventModel::MC::Event  + "/McLyrToHitTab";
        TkrEventModel::MC::McLayerHitCol     = TkrEventModel::MC::Event  + "/McLayerHitCol";
    }
};
    
    
static TkrEvModel TkrMod;    // where  used? has file scope     
    
