#include "src/Vertex/Combo/TkrComboVtxRep.h"

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
        int numVertices = pVertices->getNumVertices();

        gui::DisplayRep* pDisplay = this;

        while(numVertices--)
        {
            Event::TkrVertex* pVertex = pVertices->getVertex(numVertices);

            Point startPoint = Point(pVertex->getPosition());
            Point endPoint   = startPoint;

            endPoint -= 300.*pVertex->getDirection();

            // draw reconstructed gamma
            setColor("yellow");
            markerAt(startPoint);
            moveTo(startPoint);
            lineTo(endPoint);
            setColor("black");
        }
    }

    return;
}

