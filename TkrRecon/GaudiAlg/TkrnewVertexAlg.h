
#ifndef __TKRNEWVERTEXALG_H
#define __TKRNEWVERTEXALG_H 1

#include "GaudiKernel/Algorithm.h"
#include "src/Vertex/DocaVtx/DocaVtxAlg.h"


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
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/GaudiAlg/TkrVertexAlg.h,v 1.4 2002/05/12 05:52:58 usher Exp $
*/



class TkrnewVertexAlg : public Algorithm
{
public:
    TkrnewVertexAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrnewVertexAlg() {}
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    
    Algorithm* m_docaVtxAlg;
};

#endif // __TKRNEWVERTEXALG_H
