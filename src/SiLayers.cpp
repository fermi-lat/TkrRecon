
#include "TkrRecon/SiLayers.h"

#include "Gaudi/MessageSvc/MsgStream.h"
//#include "Event/messageManager.h"
//#include "instrument/trackerGeo.h"

//----------------- SiLayer ------------------

//################################################
SiLayer::SiLayer(int ilayer, int iview, int ToT):
m_layer(ilayer),m_view(iview),m_ToT(ToT)
//################################################
{
	clear();
}
//################################################
void SiLayer::addStrip(int istrip)
//################################################
{
//	if (validStrip(istrip)) 
		m_stripList.push_back(istrip);
}
//################################################
void SiLayer::writeOut() const
//################################################
{
      std::cout << "Layer " << m_layer << " xy " << m_view << " tot " << m_ToT <<" nstrips " << m_stripList.size() << " ";
      for (int i = 0; i < m_stripList.size();i++) {
          std::cout << m_stripList[i] << " ";
    }
      std::cout <<"\n";
}
//----------------- private -----------------------
/*
//################################################
bool SiLayer::validStrip(int istrip)
//################################################
{
	bool valid = true;
	
	//  check if the strip is physical
	int nladders = trackerGeo::nLadders(m_layer,detGeo::makeAxis(m_view));
	int iladder = istrip/trackerGeo::ladderNStrips();
	if (iladder >= nladders) valid = false;
	if (!valid) {
		messageManager::instance()->message(" SiLayer: strip not valid",istrip,"DEBUG");
	}
	return valid; 
}
*/
//----------------- SiLayers -----------------
//################################################
void SiLayers::clear()
//################################################
{
	int nlayers = num();
	for (int ilayer = 0; ilayer < nlayers; ilayer++) {
		delete m_SiLayersList[ilayer];
	}
	m_SiLayersList.clear();
}

//################################################
void SiLayers::ini()
//################################################
{
	m_SiLayersList.clear();
}
//################################################
void SiLayers::writeOut() const
//################################################
{
	if (m_SiLayersList.size() <=0) return;

        std::cout << " SiLayers " << m_SiLayersList.size() <<"\n";
	for (int i = 0; i < m_SiLayersList.size();i++) {
		m_SiLayersList[i]->writeOut();
	}
}
//################################################
int SiLayers::numStrips()
//################################################
{
	int nStrips = 0;
	for (int i = 0; i < m_SiLayersList.size();i++) {
		nStrips += m_SiLayersList[i]->nstrips();
	}
	return nStrips;
}
//################################################
SiLayer* SiLayers::getSiLayer(int klayer, int iview)
//################################################
{
	SiLayer* sd = 0;
	int nlayers = m_SiLayersList.size();
	for (int ilayer = 0; ilayer < nlayers; ilayer++){
		if (m_SiLayersList[ilayer]->layer() == klayer &&
			m_SiLayersList[ilayer]->view()  == iview) {
			sd = m_SiLayersList[ilayer];
			break;
		}
	}
	return sd;
}
