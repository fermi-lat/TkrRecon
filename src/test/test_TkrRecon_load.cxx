// $Header: /nfs/slac/g/glast/ground/cvs/TkrDigi/src/test/test_TkrDigi_load.cxx,v 1.1 2002/08/25 17:13:12 lsrea Exp $
//====================================================================
//
//  Description: Implementation of <Package>_load routine.
//               This routine is needed for forcing the linker
//               to load all the components of the library. 
//
//====================================================================

#include "GaudiKernel/ICnvFactory.h"
#include "GaudiKernel/ISvcFactory.h"
#include "GaudiKernel/IAlgFactory.h"


#define DLL_DECL_SERVICE(x)    extern const ISvcFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_CONVERTER(x)  extern const ICnvFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_ALGORITHM(x)  extern const IAlgFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_OBJECT(x)     extern const IFactory& x##Factory; x##Factory.addRef();

//! Load all  services: 
void test_TkrRecon_load() {
    DLL_DECL_ALGORITHM( test_TkrRecon );
} 

extern "C" void test_TkrRecon_loadRef()    {
    test_TkrRecon_load();
}

