#include "src/Vertex/Combo/TkrComboVtxRep.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

//#############################################################################
TkrComboVtxRep::TkrComboVtxRep(IDataProviderSvc* dataProviderSvc, ITkrGeometrySvc* pTkrGeometry)
//#############################################################################
{
    dps     = dataProviderSvc;
    pTkrGeo = pTkrGeometry;
}
//-------------------- private ----------------------
//##############################################
void TkrComboVtxRep::update()
//##############################################
{
    Event::TkrVertexCol* pVertices = SmartDataPtr<Event::TkrVertexCol>(dps,"/Event/TkrRecon/TkrVertexCol");

    //Now see if we can do the drawing
    if (pVertices)
    {
        gui::DisplayRep* pDisplay = this;
      
        Event::TkrVertexCol::const_iterator iter;

        for(iter = pVertices->begin(); iter != pVertices->end(); ++iter)
        {
	        const Event::TkrVertex& pVertex = **iter;
	  
	        Point startPoint = Point(pVertex.getPosition());
	  
	        // draw reconstructed gamma
	        setColor("yellow");
	        markerAt(startPoint);
	        moveTo(startPoint);
	        lineTo(startPoint - 300.*pVertex.getDirection());
	        setColor("black");
	    }
    }
  
    return;
}

