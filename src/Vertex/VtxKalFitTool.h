#ifndef _VtxKalFitTool_H
#define _VtxKalFitTool_H 1

#include "VtxBaseTool.h"

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
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
 * $Header: /home/cvs/SLAC/TkrRecon/src/Vertex/VtxKalFitTool.h,v 1.8 2002/09/02 19:46:15 cohen Exp $
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
  StatusCode initVertex(Event::TkrFitTrackCol&);

  ///main method: implements the filter
  StatusCode doVtxFit(Event::TkrVertexCol&);

  ///@brief bring geometrical momentum (Sx,Sy,E) close to current vtx estimate.
  ///In theory this method should return the geometrical momentum at POCA to
  ///current vertex estimate
  ///(arbitrary yet reasonnable choice, advocated by Luchsinger et al.). 
  ///In our case only the energy might be changed by this.
  HepVector computeQatVtx(const Event::TkrFitTrack& /*theTrack*/,
                          const HepVector /*theVertex*/);

  ///@brief Get the weight matrix G = Cov(m)^-1 with m the track parameters 
  ///Cov(m) is first propagated back to the vertex estimate, before being 
  ///inverted.
  ///@param theTrack: current fitted track;
  ///@param Vtx:      current vertex estimates;
  ///@return The weight matrix is returned as an HepSymMatrix object.
  HepSymMatrix computeWeightMatrix(const Event::TkrFitTrack& theTrack,
                                   const HepVector Vtx);

  ///@brief propagate Cov from first hit location to current Vtx estimate.
  ///@param CovMatrix TkrFitMatrix object to be propagated;
  ///@param Vtx Current estimate of the vertex; 
  Event::TkrFitMatrix propagCovToVtx(const Event::TkrFitMatrix, 
                                     const HepVector);

  ///@brief returns (X,Sx,Y,Sy,E) as an HepVector.
  ///(X,Sx,Y,Sy) are the track fitted parameters, E is its estimated energy. 
  ///These play for the Kalman vertexer the role of measurement vector, 
  ///the measurement error being their Cov. matrix.
  HepVector getTkrParVec(const Event::TkrFitTrack& /*theTrack*/);

 private:

  ///max value of single track chi2, above which track is rejected
  double m_chi2max;

  ///Reference z-plane where linearization occurs.
  double m_Zref;

  ///compute measurement matrix
  HepVector computeVectorH(const HepVector, const HepVector);

  ///compute derivative of measurement matrix wrt vertex
  HepMatrix computeMatrixA(const HepVector, const HepVector);

  ///compute derivative of measurement matrix wrt momentum
  HepMatrix computeMatrixB(const HepVector, const HepVector);

  ///translate a  TkrFitMatrix object into an HepSymMatrix object
  HepSymMatrix getHepSymCov(const Event::TkrFitMatrix& );

  ///Compute derivative matrix of (Sx,Sy,E)->(Eux,Euy,Euz) transformation
  ///@param Q HepVector (Sx,Sy,E)
  ///@param sign_uz Sign(uz) needed to remove Sx and Sy sign ambiguity
  HepMatrix SlopeToDir(HepVector /*Q*/, int /*sign_uz*/);

  ///Vector of successive estimates.
  std::vector<HepVector>    m_VtxEstimates;
  std::vector<HepSymMatrix> m_VtxCovEstimates;

};
#endif

