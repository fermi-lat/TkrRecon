
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

  ///@brief bring geometrical momentum (Sx,Sy,E) close to current vertex estimate.
  ///In theory this method should return the geometrical
  ///momentum at POCA to vtx (conventional but reasonnable choice, advocated by the author of the paper). 
  ///In our case only the energy could possibly be changed by that.
  ///This is not yet implemented.
  HepVector computeQatVtx(const Event::TkrFitTrack& /*theTrack*/,const HepVector /*theVertex*/);

  ///@brief Get the weight matrix G = Cov(m)^-1 with m the track parameters 
  ///Cov(m) is first propagated back to the vertex estimate, before being inverted.
  HepSymMatrix computeWeightMatrix(const Event::TkrFitTrack& theTrack,const HepVector Vtx);

  ///@brief propagate Cov matrix from first hit location to Zref plane
  ///method not yet implemented
  Event::TkrFitMatrix propagCovToVtx(const Event::TkrFitMatrix, 
				     const HepVector);

  //@brief (X,Sx,Y,Sy) are returned as an HepVector.
  //(X,Sx,Y,Sy) are the track fitted parameters. They play for
  //the Kalman vertexer the role of measurement vector. 
  HepVector getTkrParVec(const Event::TkrFitTrack& /*theTrack*/);

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

  ///Compute Transformation matrix (Sx,Sy)->(ux,uy,uz)
  HepMatrix SlopeToDir(HepVector /*Q*/);

  ///Vector of successive estimates.
  std::vector<HepVector>    m_VtxEstimates;
  std::vector<HepSymMatrix> m_VtxCovEstimates;

};
#endif
