
// Include files
#include "TkrRecon/TkrRecoAlg.h"

#include "Gaudi/MessageSvc/MsgStream.h"
#include "Gaudi/Kernel/AlgFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "Gaudi/DataSvc/SmartDataPtr.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "xml/IFile.h"

#include "GlastEvent/data/TdSiData.h"
#include "GlastEvent/data/TdGlastData.h"
;
static const AlgFactory<TkrRecoAlg>  Factory;
const IAlgFactory& TkrRecoAlgFactory = Factory;

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.
TkrRecoAlg::TkrRecoAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator), m_detSvc(0) {
}


//------------------------------------------------------------------------------
/*! The "functional" part of the class: For the EmptyAlgorithm example they do
nothing apart from print out info messages.
NB in the initialize method: you must explicitly initialize the base class
before using any services (message service, event data service etc.) otherwise 
the behaviour will be unpredictable at best.
*/
StatusCode TkrRecoAlg::initialize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    // now try to find the GlastDevSvc service
    StatusCode sc = service("GlastDetSvc", m_detSvc);

    if (!sc.isSuccess ()){
        log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
    }

    // test: get a constant from the ini file
    m_ini = const_cast<xml::IFile*>(m_detSvc->iniFile()); //OOPS!
    int nx = m_ini->getInt("glast", "xNum");


    // define the tuple
    return sc;
}


//------------------------------------------------------------------------------
StatusCode TkrRecoAlg::execute() {
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << "execute" << endreq;
        
    // get the CsiData object from the TDS by a converter
    SmartDataPtr<TdSiData> si(eventSvc(), "/Event/Raw/TdSiDatas");
    SmartDataPtr<TdCsIData> csi(eventSvc(), "/Event/Raw/TdCsIDatas");

    //Now build the TdGlastData object;

    // see what is there

    return sc;
}


//------------------------------------------------------------------------------
StatusCode TkrRecoAlg::finalize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize" << endreq;
    
    return StatusCode::SUCCESS;
}






