
#ifndef TkrVertexAlg_H
#define TkrVertexAlg_H 1

#include "TkrRecon/Vertex/TkrFindVertex.h"

#include "GaudiKernel/Algorithm.h"

//----------------------------------------------
//
//   TkrVertexAlg
//
//   Algorithm Data constructor of TkrVertexCol
//----------------------------------------------
//   Tracy Usher 03/01/02
//----------------------------------------------
//##########################################################
class TkrVertexAlg : public Algorithm
//##########################################################
{
public:
    //! Constructor of this form must be provided
    TkrVertexAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrVertexAlg() {}
    //! mandatory
    StatusCode initialize();
    //! mandatory
    StatusCode execute();
    //! mandatory
    StatusCode finalize();
    
private:
    
    TkrFindVertex* pFindVertex;
};

#endif
