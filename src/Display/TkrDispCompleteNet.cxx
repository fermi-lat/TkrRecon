//------------------------------------------------------------------------------
// TkrDispCompleteNet for TkrNeurlPR Implementation
//
// Used for displaying the entire neural network. 
// Mainly for debugging.  Copied from Tracy's display routines.
//
// b. allgood and w. atwood, 3/02  
//------------------------------------------------------------------------------

#include "TkrDispCompleteNet.h"
#include "Event/TopLevel/EventModel.h"
#include "src/PatRec/NeuralNet/TkrNeuralNet.h"
#include "src/PatRec/NeuralNet/TkrNeuron.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

TkrDispCompleteNet::TkrDispCompleteNet(IDataProviderSvc* dataProviderSvc, 
                                       ITkrGeometrySvc* tkrGeom)
{
    dps     = dataProviderSvc;
    m_tkrGeom = tkrGeom;
}

//-------------------- private ----------------------

void TkrDispCompleteNet::update()
{
    Event::TkrTrackCol* pTkrTrackCol = 
      SmartDataPtr<Event::TkrTrackCol>(dps,EventModel::TkrRecon::TkrTrackCol);

    //Now see if we can do the drawing
    if (pTkrTrackCol)
    {  
      DataObject* dataObj;
      std::string s("/Event/NeuralNet");
      StatusCode sc = dps->retrieveObject(s,dataObj);
      if(sc.isFailure())
	{
	  std::cout<<"WARNING: NeuralNet could not be retrieved from TDS"<<std::endl;
	  return;
	}

      TkrNeuralNet* pTkrNeuralNet = dynamic_cast<TkrNeuralNet*>(dataObj);

      const char violet[]     = "violet";
      setColor(violet);
      
      TkrNeuronList tmpList = pTkrNeuralNet->neurons();
      TkrNeuronList::const_iterator hypo;
      
      for(hypo  = tmpList.begin(); hypo != tmpList.end();   hypo++){
	
	Point point0 = (*hypo).getPnt(top);
	Point point1 = (*hypo).getPnt(bottom);
	
	moveTo(point0);
	lineTo(point1);
	
      }
      
      setColor("blue");
    }
    
    return;
}
