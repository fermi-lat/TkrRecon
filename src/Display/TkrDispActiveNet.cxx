//------------------------------------------------------------------------------
// TkrDispActiveNet for TkrNeurlPR Implementation
//
// Used for displaying the active neural network. 
// Mainly for debugging.  Copied from Tracy's display routines.
//
// b. allgood and w. atwood, 3/02  
//------------------------------------------------------------------------------

#include "TkrRecon/Display/TkrDispActiveNet.h"

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

const char col_blue[]       = "blue";
const char col_red[]        = "red";
const char col_yellow[]     = "yellow";
const char col_violet[]     = "violet";
const char col_turquoise[]  = "turquoise";
const char col_orange[]     = "orange";
const char col_maroon[]     = "maroon";
const char col_aquamarine[] = "aquamarine";

const char* pNNColors[] = {col_red, col_orange, col_yellow, col_aquamarine,
                           col_blue, col_violet};


TkrDispActiveNet::TkrDispActiveNet(IDataProviderSvc* dataProviderSvc, 
                                   ITkrGeometrySvc* pTkrGeometry)
{
    dps     = dataProviderSvc;
    pTkrGeo = pTkrGeometry;
}

//-------------------- private ----------------------

void TkrDispActiveNet::update()
{
    TkrPatCandCol* pTkrPatCandCol = SmartDataPtr<TkrPatCandCol>(dps,"/Event/TkrRecon/TkrPatCandCol");

	//Now see if we can do the drawing
	if (pTkrPatCandCol)
	{

		TkrNeuralNet* pTkrNeuralNet = dynamic_cast<TkrNeuralNet*>(pTkrPatCandCol);

        int numDispNeurons = pTkrNeuralNet->numNeurons();
        int colorIdx      = 1;
		setColor(pNNColors[colorIdx]);

        gui::DisplayRep* pDisplay = this;

		TkrNeuralNet::TkrNeuronList tmpList = pTkrNeuralNet->neurons();
		TkrNeuralNet::TkrNeuronList::const_iterator hypo;

		for(hypo  = tmpList.begin(); 
		hypo != tmpList.end();   hypo++){
			
			if((*hypo).getActivity() >= 0.9) colorIdx = 0;
			else continue;

            setColor(pNNColors[colorIdx]);    
			Point point0 = (*hypo).getPnt(top);
			Point point1 = (*hypo).getPnt(bottom);

			moveTo(point0);
			lineTo(point1);

		}

	    setColor("blue");
    }

    return;
}