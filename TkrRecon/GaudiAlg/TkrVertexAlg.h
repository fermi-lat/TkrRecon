#ifndef __TKRVERTEXALG_H
#define __TKRVERTEXALG_H 1
/*
#include "GaudiKernel/Algorithm.h"
#include "src/Vertex/IVtxBaseTool.h"
*/
/** 
 * @class TkrVertexAlg
 *
 * @brief TkrRecon Gaudi Algorithm for controlling the search for and fitting
 *        of all possible vertices from the collection of fit tracks. 
 *        Gaudi Tools are used to implement a particular type of vertex algorithm, 
 *        and this algorithm controls their creation and use. 
 *        The algorithm depends upon input from the track finding and fitting
 *        stages of TkrRecon. Results are output to the TDS class TkrVertex
 *
 * 01-03-2002 
 *
 * Adapted and augmented from vertex finding code in Atwood/Hernando code
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/TkrRecon/GaudiAlg/TkrVertexAlg.h,v 1.2 2004/09/08 15:32:41 usher Exp $
 */
/*
class TkrVertexAlg : public Algorithm
{
public:
    // Standard Gaudi Algorithm constructor format
    TkrVertexAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrVertexAlg() {}

    // The thee phases in the life of a Gaudi Algorithm
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    // Type of vertexing algorithm to run
    std::string   m_VertexerType;

    // Yet another fine tool from Sears
    IVtxBaseTool* m_VtxTool;
};
*/
#endif // __TKRVERTEXALG_H
