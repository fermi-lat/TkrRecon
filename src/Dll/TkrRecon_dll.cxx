/** @file TkrRecon_dll.cxx
@brief dll installation routine

$Header$
*/

//====================================================================
//  TkrRecon_dll.cxx
//--------------------------------------------------------------------
//
//  Package    : GlastSvc
//
//  Description: Implementation of DllMain routine.
//               The DLL initialisation must be done seperately for 
//               each DLL. 
//
//  Author     : H. Gillespie
//  Created    : 1 Aug 00
//  Changes    : 
//               Version copied from Glast, and addopted for GlastSvc
//
//====================================================================


#include "GaudiKernel/LoadFactoryEntries.h"
LOAD_FACTORY_ENTRIES(TkrRecon)
