#include "Gaudi/MessageSvc/MsgStream.h"
#include "Gaudi/Kernel/AlgFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "Gaudi/DataSvc/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "TkrRecon/SiLayersIRFAlg.h"
#include "TkrRecon/SiLayers.h"
#include "TkrRecon/trackerGeo.h"
#include "GlastEvent/Raw/TdSiData.h"

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
	
	return StatusCode::SUCCESS;
}
//##############################################
StatusCode SiLayersIRFAlg::finalize()
//##############################################
{
//	m_SiLayers->writeOut();
	return StatusCode::SUCCESS;
}

//##############################################
StatusCode SiLayersIRFAlg::retrieve()
//##############################################
{
	m_SiLayers = new SiLayers();
	m_SiLayers->clear();
	StatusCode sc = eventSvc()->registerObject("/Event/Recon/TkrRecon","SiLayers",m_SiLayers);

     // get the rdSiData object from the TDS by a converter
     m_SiData   = SmartDataPtr<TdSiData>(eventSvc(),"/Event/Raw/TdSiDatas");

	if (m_SiLayers == 0 || m_SiData == 0) sc = StatusCode::FAILURE;
	return sc;
}