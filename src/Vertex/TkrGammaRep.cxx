#include "src/Vertex/TkrGammaRep.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

//#############################################################################
TkrGammaRep::TkrGammaRep(IDataProviderSvc* dataProviderSvc, ITkrGeometrySvc* pTkrGeometry)
//#############################################################################
{
    dps     = dataProviderSvc;
    pTkrGeo = pTkrGeometry;
}
//-------------------- private ----------------------
//##############################################
void TkrGammaRep::update()
//##############################################
{
    Event::TkrVertexCol* pVertices = SmartDataPtr<Event::TkrVertexCol>(dps,"/Event/TkrRecon/TkrVertexCol");

    //Now see if we can do the drawing
    if (pVertices)
    {
        if (pVertices->size() > 0)
        {
            Event::TkrVertexCol::const_iterator iter = pVertices->begin();

            const Event::TkrVertex& pVertex = **iter;
      
            Point startPoint = Point(pVertex.getPosition());
      
            // draw reconstructed gamma
            setColor("yellow");
            markerAt(startPoint);
            moveTo(startPoint);
            lineTo(startPoint - 1000.*pVertex.getDirection());
            setColor("black");
        }
    }
  
    return;
}

