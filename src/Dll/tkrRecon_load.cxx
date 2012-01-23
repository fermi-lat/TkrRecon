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
    DECLARE_ALGORITHM( TkrFilterAlg            );
    DECLARE_ALGORITHM( TkrFindAlg              );
    DECLARE_ALGORITHM( TkrTrackFitAlg          );
    DECLARE_ALGORITHM( TkrVertexAlg            );
    DECLARE_ALGORITHM( TkrDisplayAlg           );

    DECLARE_SERVICE(   TkrInitSvc              );

    DECLARE_TOOL(      VtxSingleTrkTool        );
    DECLARE_TOOL(      VtxKalFitTool           );
    DECLARE_TOOL(      ComboVtxTool            );
    DECLARE_TOOL(      KalmanTrackFitTool      );
    DECLARE_TOOL(      ComboFindTrackTool      );
    DECLARE_TOOL(      LinkAndTreeFindTrackTool);
    DECLARE_TOOL(      NeuralNetFindTrackTool  );
    DECLARE_TOOL(      MonteCarloFindTrackTool );
    DECLARE_TOOL(      VectorLinksTool         );
    DECLARE_TOOL(      TkrTrackEnergyTool      );
    DECLARE_TOOL(      TkrEnergySplitTool      );
    DECLARE_TOOL(      TkrAlignHitsTool        );
    DECLARE_TOOL(      FindTrackHitsTool       );
    DECLARE_TOOL(      TkrCalFilterTool        );
    DECLARE_TOOL(      TkrFilterTool           );
    DECLARE_TOOL(      TkrVecPointsFilterTool  );
    DECLARE_TOOL(      TkrHoughFilterTool      );
    DECLARE_TOOL(      TreeBasedTool           );   // added 27-Oct-2010 LSR

//    DECLARE_TOOL(      TkrComboFitTool         );
//    DECLARE_TOOL(      TkrLinkAndTreeFitTool   );
//    DECLARE_TOOL(      TkrNeuralNetFitTool     );

//LSR moved to TkrUtil 07-Nov-08, the version in TkrRecon is obsolete!
//    DECLARE_TOOL(      TkrGhostTool            ); 

} 

