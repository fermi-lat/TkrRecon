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

#include "GaudiKernel/ICnvFactory.h"
#include "GaudiKernel/ISvcFactory.h"
#include "GaudiKernel/IAlgFactory.h"
#include "GaudiKernel/IToolFactory.h"


#define DLL_DECL_SERVICE(x)    extern const ISvcFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_CONVERTER(x)  extern const ICnvFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_ALGORITHM(x)  extern const IAlgFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_OBJECT(x)     extern const IFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_TOOL(x)       extern const IToolFactory& x##Factory; x##Factory.addRef();

//! Load all  services: 
void TkrRecon_load() {
    DLL_DECL_SERVICE(   TkrGeometrySvc     );
    DLL_DECL_SERVICE(   TkrInitSvc         );
    DLL_DECL_SERVICE(   TkrBadStripsSvc    );
    DLL_DECL_ALGORITHM( TkrReconAlg        );
    DLL_DECL_ALGORITHM( TkrClusterAlg      );
    DLL_DECL_ALGORITHM( TkrFindAlg         );
    DLL_DECL_ALGORITHM( TkrTrackFitAlg     );
    DLL_DECL_ALGORITHM( TkrVertexAlg       );
    DLL_DECL_ALGORITHM( TkrDisplayAlg      );
//    DLL_DECL_ALGORITHM( McReconAlg       );
//    DLL_DECL_ALGORITHM( TkrNtupleAlg     );
//    DLL_DECL_ALGORITHM( RecNtupleAlg     );
    DLL_DECL_TOOL( VtxSingleTrkTool        );
    DLL_DECL_TOOL( VtxKalFitTool           );
    DLL_DECL_TOOL( ComboVtxTool            );
    DLL_DECL_TOOL( TkrComboFitTool         );
    DLL_DECL_TOOL( TkrLinkAndTreeFitTool   );
    DLL_DECL_TOOL( TkrNeuralNetFitTool     );
    DLL_DECL_TOOL( ComboFindTrackTool      );
    DLL_DECL_TOOL( LinkAndTreeFindTrackTool);
    DLL_DECL_TOOL( NeuralNetFindTrackTool  );

} 

extern "C" void TkrRecon_loadRef()    {
  TkrRecon_load();
}

