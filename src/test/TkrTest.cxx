
// Include files
#include "GaudiKernel/SmartIF.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IAppMgrUI.h"
#include "GaudiKernel/IProperty.h"
#include "GaudiKernel/Property.h"
#include "GaudiKernel/System.h"


//------------------------------------------------------------------------------
//
//  Package    : CalRecon
//
//  Description: Test Main Program
//
//------------------------------------------------------------------------------
#include <iostream>
// this needed for instrument.lib??
#ifndef WIN32
void FATAL(const char* text){std::cerr << text << std::endl;}
#endif // WIN32



int main( int argn, char** argc) {
    


    // get the path to this package from its root environment variable: if not there,
    // assume that we are in the root
    const char * local_path = ::getenv("TKRRECONROOT");
    std::string joboptions_file = std::string(local_path? local_path: "");

    joboptions_file +=  std::string("/src/test/jobOptions.txt");
    
    // Create an instance of an application manager
    IInterface* iface = Gaudi::createApplicationMgr();
    
    SmartIF<IProperty>     propMgr ( IID_IProperty, iface );
    SmartIF<IAppMgrUI>     appMgr  ( IID_IAppMgrUI, iface );
    
    // Set properties of algorithms and services
    StatusCode status = StatusCode::SUCCESS;
    if ( iface && propMgr == iface )    {
        status = propMgr->setProperty( StringProperty("JobOptionsPath", joboptions_file) );
    }
    else  {
        exit(-1);
    }

    
    // Run the application manager and process events
    if ( appMgr )   {
        status = appMgr->run();
    }
    else  {
        return 0;
    }
    
    // All done - exit
    return 0;
    
}
