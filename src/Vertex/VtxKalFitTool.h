#ifndef _VtxKalFitTool_H
#define _VtxKalFitTool_H 1

#include "VtxBaseTool.h"

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

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
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Vertex/VtxKalFitTool.h,v 1.12 2004/12/16 05:04:24 usher Exp $
 */
class VtxKalFitTool : public VtxBaseTool
{
 public:
  VtxKalFitTool( const std::string& type,
                 const std::string& name,
                 const IInterface* parent);

  virtual ~VtxKalFitTool() { }

  ///base tool overwritten method: need to set local properties
  StatusCode initialize();

  ///@brief Finds start estimate of vertex and its Cov matrix
  /// First estimate is first hit of first track on the list, which is 
  /// the best fitted track.
  StatusCode initVertex(Event::TkrTrackCol&);

  ///main method: implements the filter
  StatusCode doVtxFit(Event::TkrVertexCol&);

  ///@brief bring geometrical momentum (Sx,Sy,E) close to current vtx estimate.
  ///In theory this method should return the geometrical momentum at POCA to
  ///current vertex estimate
  ///(arbitrary yet reasonnable choice, advocated by Luchsinger et al.). 
  ///In our case only the energy might be changed by this.
  CLHEP::HepVector computeQatVtx(const Event::TkrTrack& /*theTrack*/,
                                 const CLHEP::HepVector /*theVertex*/);

  ///@brief Get the weight matrix G = Cov(m)^-1 with m the track parameters 
  ///Cov(m) is first propagated back to the vertex estimate, before being 
  ///inverted.
  ///@param theTrack: current fitted track;
  ///@param Vtx:      current vertex estimates;
  ///@return The weight matrix is returned as an HepSymMatrix object.
  CLHEP::HepSymMatrix computeWeightMatrix(const Event::TkrTrack& theTrack,
                                          const CLHEP::HepVector Vtx);

  ///@brief propagate Cov from first hit location to current Vtx estimate.
  ///@param CovMatrix TkrFitMatrix object to be propagated;
  ///@param Vtx Current estimate of the vertex; 
  Event::TkrTrackParams propagCovToVtx(const Event::TkrTrackParams, 
                                       const CLHEP::HepVector);

  ///@brief returns (X,Sx,Y,Sy,E) as an HepVector.
  ///(X,Sx,Y,Sy) are the track fitted parameters, E is its estimated energy. 
  ///These play for the Kalman vertexer the role of measurement vector, 
  ///the measurement error being their Cov. matrix.
  CLHEP::HepVector getTkrParVec(const Event::TkrTrack& /*theTrack*/);

 private:

  ///max value of single track chi2, above which track is rejected
  double m_chi2max;

  ///Reference z-plane where linearization occurs.
  double m_Zref;

  ///compute measurement matrix
  CLHEP::HepVector computeVectorH(const CLHEP::HepVector, const CLHEP::HepVector);

  ///compute derivative of measurement matrix wrt vertex
  CLHEP::HepMatrix computeMatrixA(const CLHEP::HepVector, const CLHEP::HepVector);

  ///compute derivative of measurement matrix wrt momentum
  CLHEP::HepMatrix computeMatrixB(const CLHEP::HepVector, const CLHEP::HepVector);

  ///translate a  TkrFitMatrix object into an HepSymMatrix object
  CLHEP::HepSymMatrix getHepSymCov(const Event::TkrTrackParams& );

  ///Compute derivative matrix of (Sx,Sy,E)->(Eux,Euy,Euz) transformation
  ///@param Q HepVector (Sx,Sy,E)
  ///@param sign_uz Sign(uz) needed to remove Sx and Sy sign ambiguity
  CLHEP::HepMatrix SlopeToDir(CLHEP::HepVector /*Q*/, int /*sign_uz*/);

  ///Vector of successive estimates.
  std::vector<CLHEP::HepVector>    m_VtxEstimates;
  std::vector<CLHEP::HepSymMatrix> m_VtxCovEstimates;

};
#endif

