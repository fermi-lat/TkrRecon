
#ifndef __TKRVERTEXALG_H
#define __TKRVERTEXALG_H 1

#include "TkrRecon/Vertex/TkrFindVertex.h"

#include "GaudiKernel/Algorithm.h"


/** 
* @class TkrVertexAlg
*
* @brief Algorithm to construct TkrVertexCol and its contents
*
* 01-03-2002 
*
* Adapted and augmented from vertex finding code in Atwood/Hernando code
*
* @author Tracy Usher, Leon Rochester
*
* $Header$
*/



class TkrVertexAlg : public Algorithm
{
public:
    TkrVertexAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrVertexAlg() {}
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    
    TkrFindVertex* pFindVertex;
};

#endif // __TKRVERTEXALG_H
