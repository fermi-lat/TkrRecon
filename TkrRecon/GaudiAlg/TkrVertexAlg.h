
#ifndef __TKRVERTEXALG_H
#define __TKRVERTEXALG_H 1

#include "GaudiKernel/Algorithm.h"
#include "src/Vertex/IVtxBaseTool.h"


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
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/GaudiAlg/TkrVertexAlg.h,v 1.5 2002/07/22 13:45:53 cohen Exp $
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
    std::string m_VertexerType;
    IVtxBaseTool* m_VtxTool;
};

#endif // __TKRVERTEXALG_H
