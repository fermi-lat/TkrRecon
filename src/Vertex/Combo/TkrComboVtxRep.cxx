#include "src/Vertex/Combo/TkrComboVtxRep.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

//#############################################################################
TkrComboVtxRep::TkrComboVtxRep(IDataProviderSvc* dataProviderSvc, ITkrGeometrySvc* tkrGeom)
//#############################################################################
{
    dps     = dataProviderSvc;
    m_tkrGeom = tkrGeom;
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
      //gui::DisplayRep* pDisplay = this;
      
        Event::TkrVertexCol::const_iterator iter;

        for(iter = pVertices->begin() + 1; iter < pVertices->end(); ++iter)
        {
            const Event::TkrVertex& pVertex = **iter;
      
            Point startPoint = Point(pVertex.getPosition());
      
            // draw reconstructed gamma
            setColor("brown");
            markerAt(startPoint);
            moveTo(startPoint);
            lineTo(startPoint - 1000.*pVertex.getDirection());
            setColor("black");
        }
    }
  
    return;
}

