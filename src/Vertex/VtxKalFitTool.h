
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
  VtxKalFitTool( const std::string& type,
		 const std::string& name,
		 const IInterface* parent);

  virtual ~VtxKalFitTool() { }

  ///base tool overwritten method
  StatusCode initialize();

  ///@brief Initialization: first estimate of Vertex is first hit of best track
  ///best track is assumed to be first in the TkrFitTrack list
  StatusCode initVertex(Event::TkrFitTrackCol&);

  ///main method: implements the filter
  StatusCode doVtxFit(Event::TkrVertexCol&);

  double const getChi2() {return m_chi2;}

  //bring momentum (Sx,Sy) to Zref. currently does nothing
  HepVector computeQatZref(const Event::TkrFitTrack& theTrack);

  //@brief propagate Cov matrix from first hit location to Zref plane
  //method not yet implemented
  Event::TkrFitMatrix propagCovToVtx(const Event::TkrFitMatrix, 
				     const HepVector);

 private:

  ///max value of single track chi2, above which track is rejected
  double m_chi2max;

  ///Reference z-plane where linearization occurs.
  double m_Zref;

  ///compute measurement matrix
  HepVector computeVectorH(const HepVector, const HepVector);
  ///compute derivate of measurement matrix wrt vertex
  HepMatrix computeMatrixA(const HepVector, const HepVector);
  ///compute derivate of measurement matrix wrt momentum
  HepMatrix computeMatrixB(const HepVector, const HepVector);
  ///returns TkrFitMatrix as an HepSymMatrix
  HepSymMatrix getHepSymCov(const Event::TkrFitMatrix& );


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
  std::vector<HepSymMatrix>    m_W;

  //Final matrices:
  HepSymMatrix m_CovXX;
  std::vector<HepVector>    m_m;
  std::vector<HepVector>    m_Q;
  std::vector<HepMatrix>    m_CovQQ;
  std::vector<HepMatrix>    m_CovXQ;
  std::vector<HepMatrix>    m_Skew;
  
};
#endif
