#include "Gaudi/MessageSvc/MsgStream.h"
#include "Gaudi/Kernel/AlgFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "Gaudi/DataSvc/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "TkrRecon/SiLayersIRFAlg.h"
#include "TkrRecon/SiLayers.h"
#include "TkrRecon/trackerGeo.h"
#include "GlastEvent/data/TdGlastData.h"

static const AlgFactory<SiLayersIRFAlg>  Factory;
const IAlgFactory& SiLayersIRFAlgFactory = Factory;

//#############################################################################
SiLayersIRFAlg::SiLayersIRFAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
//#############################################################################
{	
}

//##############################################
StatusCode SiLayersIRFAlg::initialize()
//##############################################
{
	m_SiLayers = 0;
	m_SiData = 0;
	return StatusCode::SUCCESS;
}
//##############################################
StatusCode SiLayersIRFAlg::execute()
//##############################################
{
	StatusCode sc = retrieve();
	if (sc != StatusCode::SUCCESS) return sc;
	int ToT = -1;

	for (int iview = 0; iview < trackerGeo::numViews(); iview++) {
		SiData::Axis a = SiData::X;
		if (iview == 1) a = SiData::Y;
		for (int iplane = 0; iplane < trackerGeo::numPlanes(); iplane++) {
			int nstrips = m_SiData->nHits(a,iplane);
			if (nstrips >0) {
				SiLayer* slayer = new SiLayer(iplane,iview,ToT);
				for (int istrip = 0; istrip < nstrips; istrip++) {
					int idstrip = m_SiData->hitId(a,iplane,istrip);
					// most likely there is a conversion from idstrip to real stripid
					// crosscheck this code
					slayer->addStrip(idstrip);
				}
				m_SiLayers->add(slayer);
			}
		}
	}
	
        m_SiLayers->writeOut();
	return StatusCode::SUCCESS;

}
//##############################################
StatusCode SiLayersIRFAlg::finalize()
//##############################################
{
	return StatusCode::SUCCESS;
}

//##############################################
StatusCode SiLayersIRFAlg::retrieve()
//##############################################
{
    MsgStream log(msgSvc(), name());

    StatusCode sc;

    m_SiLayers = new SiLayers();
    m_SiLayers->clear();
    
    /*! Check to see if we can get the subdirectory. If not create it
    */
    DataObject* pnode=0;

    sc = eventSvc()->retrieveObject( "/Event/TkrRecon", pnode );
    
    if( sc.isFailure() ) {
        sc = eventSvc()->registerObject("/Event/TkrRecon",new DataObject);
        if( sc.isFailure() ) {
            
            log << MSG::ERROR << "Could not create Raw directory" << endreq;
            return sc;
        }
    }

    sc = eventSvc()->registerObject( "/Event/TkrRecon/SiLayers", m_SiLayers );


    /*! Instead of loading the object directly we load a TdGlastData object
        and get the part important to the Tkr.
    */

    SmartDataPtr<TdGlastData> glastData(eventSvc(),"/Event/TdGlastData");
    // get the rdSiData object from the TDS by a converter
    m_SiData   = glastData->getSiData();
    
    if (m_SiLayers == 0 || m_SiData == 0) sc = StatusCode::FAILURE;
    return sc;
}