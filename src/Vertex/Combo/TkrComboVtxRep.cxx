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
    TkrRecon::TkrVertexCol* pVertices = SmartDataPtr<TkrRecon::TkrVertexCol>(dps,"/Event/TkrRecon/TkrVertexCol");

	//Now see if we can do the drawing
	if (pVertices)
	{
        int numVertices = pVertices->getNumVertices();

        gui::DisplayRep* pDisplay = this;

        while(numVertices--)
        {
            TkrRecon::TkrVertex* pVertex = pVertices->getVertex(numVertices);

            Point startPoint = Point(pVertex->position());
            Point endPoint   = startPoint;

            endPoint -= 300.*pVertex->direction();

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

