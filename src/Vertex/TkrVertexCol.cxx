/*
	Code to implement the TkrVertexCol class
	Tracy Usher March 1, 2002
*/

#include "TkrRecon/Vertex/TkrVertexCol.h"

TkrVertexCol::TkrVertexCol(ITkrGeometrySvc* pTkrGeo, TkrClusters* pTkrClus)
{
    m_Vertices.clear();

	return;
}

TkrVertexCol::~TkrVertexCol()
{
    int             numVertices = getNumVertices();
    TkrVertexVecPtr vtxs        = getVertexPtr();

    while(numVertices--) delete *vtxs++;
}


void TkrVertexCol::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG << " --- TkrVertexCol::writeOut --- " << endreq;
    log << MSG::DEBUG << "     Number of vertices: " << getNumVertices() << endreq;

    int vertexNo = 0;

    while(vertexNo < getNumVertices())
    {
        log << MSG::DEBUG << "     Vertex Number: " << vertexNo << endreq;
        m_Vertices[vertexNo++]->writeOut(log);
    }

    return;
}

void TkrVertexCol::addVertex(TkrVertex* vertex)
{
    m_Vertices.push_back(vertex);

    return;
}