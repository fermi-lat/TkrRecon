
#ifndef _VtxKalFitTool_H
#define _VtxKalFitTool_H 1

/**
 * @class VtxKalFitTool
 *
 * @brief Vertexing algorithm based on Kalman filter equations.
 * see:
 * Luchsinger and Grab, Comp. Phys. Comm. 76 (1993) 263-280
 * for more information.
 * 
 *
 * @author Johann Cohen-Tanugi
 */

#include "VtxBaseTool.h"

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

class VtxKalFitTool : public VtxBaseTool
{
 public:
  // Constructor
  VtxKalFitTool( const std::string& type,
		 const std::string& name,
		 const IInterface* parent);

  // Standard Destructor
  virtual ~VtxKalFitTool() { }

  StatusCode initialize();
  StatusCode initVertex(Event::TkrFitTrackCol&);

  StatusCode doVtxFit(Event::TkrVertexCol&);

  StatusCode doVtx();

  double const getChi2() {return m_chi2;}

  HepVector computeQatZref(const Event::TkrFitTrack& theTrack);

 private:

  void initVertex();

  //properties
  double m_chi2max;

  //Reference z-plane.
  double m_Zref;

  //computation helpers:
  HepVector computeVectorH(const HepVector, const HepVector);
  HepMatrix computeMatrixA(const HepVector, const HepVector);
  HepMatrix computeMatrixB(const HepVector, const HepVector);

  double m_chi2;

  std::vector<HepVector>    m_VtxEstimates;
  std::vector<HepSymMatrix> m_VtxCovEstimates;

  std::vector<double>       m_chi2f;
  std::vector<HepVector>    m_c0;
  std::vector<HepMatrix>    m_A;
  std::vector<HepMatrix>    m_B;
  std::vector<HepSymMatrix> m_C;
  std::vector<HepSymMatrix> m_D;
  std::vector<HepMatrix>    m_E;
  std::vector<HepMatrix>    m_G;
  std::vector<HepMatrix>    m_W;

  //Final matrices:
  HepSymMatrix m_CovXX;
  std::vector<HepVector>    m_m;
  std::vector<HepVector>    m_Q;
  std::vector<HepMatrix>    m_CovQQ;
  std::vector<HepMatrix>    m_CovXQ;
  std::vector<HepMatrix>    m_Skew;
  
};
#endif
