//====================================================================
//  GlastSvc_load.cpp
//--------------------------------------------------------------------
//
//  Package    : Gaudi/System
//
//  Description: Implementation of <Package>_load routine.
//               This routine is needed for forcing the linker
//               to load all the components of the library. 
//
//====================================================================

#include "GaudiKernel/DeclareFactoryEntries.h"


DECLARE_FACTORY_ENTRIES(TkrRecon) 
{
    DECLARE_ALGORITHM( TkrReconAlg             );
    DECLARE_ALGORITHM( TkrClusterAlg           );
    DECLARE_ALGORITHM( TkrFindAlg              );
    DECLARE_ALGORITHM( TkrTrackFitAlg          );
    DECLARE_ALGORITHM( TkrVertexAlg            );
    DECLARE_ALGORITHM( TkrDisplayAlg           );

    //DECLARE_SERVICE(   TkrGeometrySvc          );
    DECLARE_SERVICE(   TkrInitSvc              );
    //DECLARE_SERVICE(   TkrBadStripsSvc         );

    DECLARE_TOOL(      VtxSingleTrkTool        );
    DECLARE_TOOL(      VtxKalFitTool           );
    DECLARE_TOOL(      ComboVtxTool            );
    DECLARE_TOOL(      TkrComboFitTool         );
    DECLARE_TOOL(      TkrLinkAndTreeFitTool   );
    DECLARE_TOOL(      TkrNeuralNetFitTool     );
    DECLARE_TOOL(      ComboFindTrackTool      );
    DECLARE_TOOL(      LinkAndTreeFindTrackTool);
    DECLARE_TOOL(      NeuralNetFindTrackTool  );
    DECLARE_TOOL(      TkrMeritTool            );
} 

