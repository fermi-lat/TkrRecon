#ifndef DocaVtxAlg_H
#define TkrKalVtxAlg_H

#include "GaudiKernel/Algorithm.h"
#include "Event/Recon/TkrRecon/TkrVertexCol.h"


using namespace Event;


class DocaVtxAlg : public Algorithm
{
 public:
  
  DocaVtxAlg(const std::string& name, ISvcLocator* pSvcLoc);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

  //need to retrieve the VtxCol to fulfill TkrVertexAlg design
  TkrVertexCol* getTkrVertexCol();

 private:
  TkrVertexCol* m_VtxCol;
  double m_docaLimit; // mm
};

#endif
