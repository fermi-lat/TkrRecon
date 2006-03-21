// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Vertex/VtxKalFitTool.cxx,v 1.26 2005/02/11 07:14:55 lsrea Exp $
// Description:                                                  
//      Implementation of the Kalman vertexer
//
//
// Author
//      Johann Cohen-Tanugi


#include "VtxKalFitTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"

//In the base class once and for all?
#include "GaudiKernel/SmartDataPtr.h"
#include "Event/TopLevel/EventModel.h"

static ToolFactory<VtxKalFitTool> s_factory;
const IToolFactory& VtxKalFitToolFactory = s_factory;


VtxKalFitTool::VtxKalFitTool(const std::string& type, 
                             const std::string& name, 
                             const IInterface* parent)                        
  : VtxBaseTool(type, name, parent)
{
  //max tolerated chi2 increase by track
  declareProperty("maxChi2Contrib", m_chi2max = 10); 
}


// FROM GaudiExamples, it seems that setProperties is not needed....
StatusCode VtxKalFitTool::initialize()
{
  VtxBaseTool::initialize();
  setProperties();
  return StatusCode::SUCCESS;
}



StatusCode VtxKalFitTool::doVtxFit(Event::TkrVertexCol& /*VtxCol*/)
{/*
  m_VtxCovEstimates.erase(m_VtxCovEstimates.begin(),m_VtxCovEstimates.end());
  m_VtxEstimates.erase(m_VtxEstimates.begin(),m_VtxEstimates.end());

  std::vector<HepSymMatrix>  CovQQ;
  std::vector<CLHEP::HepMatrix>     CovXQ;
  std::vector<CLHEP::HepVector>     fQ_list;
  std::vector<CLHEP::HepVector>     m_list;
  std::vector<CLHEP::HepVector>     sQ_list;
  std::vector<CLHEP::HepMatrix>     Tmp_list;

  int    ifail      = 0; //flag for matrix inversion check
  double totalChi2  = 0; 

  MsgStream log(msgSvc(), name());
 
  Event::TkrFitTrackCol* theTracks = SmartDataPtr<Event::TkrFitTrackCol>
    (m_evtSvc,EventModel::TkrRecon::TkrFitTrackCol);


  //set initial values for vertex and its Cov.
  StatusCode sc = initVertex(*theTracks);
  if(sc.isFailure())
    {
      log << MSG::ERROR <<"Vertex not Properly initialized"<<endreq;
      return sc;
    }

  //Vector that indexes the tracks kept to build the vertex
  // (ie the daughters)
  std::vector<Event::TkrFitTrack*> usedTracks;

  log<<MSG::DEBUG;
  if (log.isActive() ) {
      log <<"NUMBER OF TRACKS TO VERTEX: "<<theTracks->size();
  }
  log <<endreq;

  //FILTER STEP:


  sc = StatusCode::FAILURE;

  Event::TkrFitColPtr tkrIter = theTracks->begin();
  while(tkrIter != theTracks->end())
    { 
      Event::TkrFitTrack* theTrack  = dynamic_cast<Event::TkrFitTrack*>(*tkrIter++);

      if(theTrack->getQuality()<0) continue;

      theTrack->writeOut(log);


      //Step 1: define new point for linearization:
      //-------------------------------------------

      //current estimations of vertex and its Cov. Matrix:
      CLHEP::HepVector    Vtx = m_VtxEstimates.back();
      HepSymMatrix invC = m_VtxCovEstimates.back();
      invC.invert(ifail);
      if(ifail) 
        {
          log << MSG::ERROR <<"ERROR INVERTING C!"<<endreq;
          return sc;
        }

      // Track parameters (X,Sx,Y,Sy,E) as Measurement Vector m
      CLHEP::HepVector m = getTkrParVec(*theTrack);

      // get "momentum" vector (Sx,Sy,E) of track at
      // POCA to vertex (conventional)
      CLHEP::HepVector    q = computeQatVtx(*theTrack,Vtx); 

      CLHEP::HepMatrix A  = computeMatrixA(Vtx,q);
      CLHEP::HepMatrix B  = computeMatrixB(Vtx,q);
      CLHEP::HepVector h0 = computeVectorH(Vtx,q);
      CLHEP::HepVector c0 = h0 - A*Vtx - B*q;

     //Step 2: Compute Covariance matrices
      //-----------------------------------
      HepSymMatrix G = computeWeightMatrix(*theTrack,Vtx);
      HepSymMatrix W = G.similarityT(B); //W = B^T*invG*B
      W.invert(ifail);
      if(ifail) 
        {
          log << MSG::ERROR <<"ERROR INVERTING W!"<<endreq;
          return StatusCode::FAILURE;
        }

      HepSymMatrix Gb = G - W.similarity(G*B);

      //Update vertex Covariance matrix:
      // (Ck+1)^-1 = Ck^-1 + A^T*Gb*A
      HepSymMatrix newC = invC + Gb.similarityT(A); 
      newC.invert(ifail);
      if(ifail)
        {
          log << MSG::ERROR <<"UPDATED MATRIX C NOT INVERTIBLE!"<<endreq;
          return StatusCode::FAILURE;
        }
      
      // D = Cov(Q,Q)     E = cov(x,Q)
      HepSymMatrix D = W + newC.similarity(W*B.T()*G*A);      
      CLHEP::HepMatrix    E = -newC*A.T()*G*B*W;

      
      //Step 3: Update vertex estimation and check chi2
      //-----------------------------------------------
      //new vertex and momentum at vertex for current track:
      CLHEP::HepVector newVtx = newC*(invC*Vtx + A.T()*Gb*(m-c0));
      CLHEP::HepVector newQ   = W*B.T()*G*(m-c0-A*newVtx);      


      CLHEP::HepVector r = m - c0 - A*newVtx - B*newQ;
      CLHEP::HepVector dx = newVtx - Vtx;
      // chi2f = r.T()*G*r + dx.T()*C*dx SCALAR
      double chi2f = G.similarity(r) + invC.similarity(dx);
      
      log<<MSG::DEBUG;
      if (log.isActive() ) {
          log <<"chi2f= "<<chi2f<<", and chi2max= "<<m_chi2max;
      }
      log << endreq;
      if(chi2f < m_chi2max) 
        {
          totalChi2 += chi2f;
          usedTracks.push_back(theTrack);
          //save for smoother step
          fQ_list.push_back(newQ);
          m_list.push_back(m);
          //keep new updates
          m_VtxEstimates.push_back(newVtx);
          m_VtxCovEstimates.push_back(newC);
          sc = StatusCode::SUCCESS;
        }
    }//end of track while loop
  
  if(sc.isFailure()) 
    {
      log<<MSG::ERROR<<"KALMAN Filter FAILED!"<<endreq;
      return sc;
    }

      
  log<<MSG::DEBUG;
  if (log.isActive() ) {
      log <<"Number of tracks used: "<<usedTracks.size();
  }
  log << endreq;

  //SMOOTHER STEP:


  //Final estimates:
  CLHEP::HepVector Vtx      = m_VtxEstimates.back();
  HepSymMatrix CovXX = m_VtxCovEstimates.back();

  std::vector<Event::TkrFitTrack*>::iterator usedIter = usedTracks.begin();
  int i=0;
  while(usedIter != usedTracks.end())
    { 
      Event::TkrFitTrack* theUsedTrack = *usedIter++;

      CLHEP::HepMatrix A  = computeMatrixA(Vtx, fQ_list[i]);
      CLHEP::HepMatrix B  = computeMatrixB(Vtx, fQ_list[i]);
      CLHEP::HepVector h0 = computeVectorH(Vtx, fQ_list[i]);
      CLHEP::HepVector c0 = h0 - A*Vtx - B*fQ_list[i];
      
      Event::TkrFitMatrix measCov = 
        propagCovToVtx(theUsedTrack->getTrackCov(), Vtx);
      measCov.invert(ifail);
      HepSymMatrix G = computeWeightMatrix(*theUsedTrack,Vtx);
      HepSymMatrix W = G.similarityT(B); //W^{-1} = B^T*G*B
      W.invert(ifail);

      //Final Q=(Sx,Sy,E)
      CLHEP::HepVector    Q = W*B.T()*G*(m_list[i] - c0 - A*Vtx);
      //Final Cov(Q)
      HepSymMatrix D = W + CovXX.similarity(W*B.T()*G*A);
      //Final Cov(Vtx,Q) 
      CLHEP::HepMatrix    E = -CovXX*A.T()*G*B*W;
      //help for cov(Q_i,Q_j) = Tmp*CovXX*Tmp.T()
      CLHEP::HepMatrix Tmp = W*B.T()*G*A;

      sQ_list.push_back(Q);
      CovQQ.push_back(D);
      CovXQ.push_back(E);
      Tmp_list.push_back(Tmp);

      i++;
    }


  // Physical Momentum and Cov

  Vector totP;
  HepSymMatrix totCovP(3,0);
  CLHEP::HepMatrix totCovXP(3,3,0);

  usedIter = usedTracks.begin();
  i=0;
  while(usedIter != usedTracks.end())
    { 
      Event::TkrFitTrack* theUsedTrack = *usedIter++;
      double uz = theUsedTrack->getDirection().z();
      int sgn = uz>0?+1:-1;

      CLHEP::HepVector Qi = sQ_list[i]; 

      //Physical Momentum:
      //-----------------
      Vector momentum = Vector(sgn*Qi[0],sgn*Qi[1],sgn*1).unit();
      momentum.setMag(Qi[2]);
      totP += momentum;

      //Cov Matrix of Physical Momentum:
      //-------------------------------
      CLHEP::HepMatrix Ti = SlopeToDir(Qi,sgn);
      totCovP += CovQQ[i].similarity(Ti);

            std::cout<<CovQQ[i]<<Ti<<totCovP<<std::endl;

      HepSymMatrix Cov_ij(3,0);
      int j;
      int numTracks = usedTracks.size();
      for(j=i+1;j<numTracks;j++)
        {
          int sgn2 = theUsedTrack->getDirection().z()>0?+1:-1;
          CLHEP::HepVector Qj = sQ_list[j];
          CLHEP::HepMatrix Tj = SlopeToDir(Qj,sgn2);
          
          CLHEP::HepMatrix Qij = Ti*Tmp_list[i]*CovXX*Tmp_list[j].T()*Tj.T();
          HepSymMatrix tmp;
          tmp.assign(Qij + Qij.T()); 
          Cov_ij += tmp;

          totCovP += Cov_ij;
        }

      totCovXP += CovXQ[i]*Ti;
      i++;
    }

    log << MSG::DEBUG;
    if (log.isActive()) {
        log.stream() <<  "totE " << totP.magnitude();
        log.stream() <<  ", totP " << totP; 
        log << endreq;
        log.stream() <<  "totCovP " << totCovP;
        log.stream() <<  ", totCovXP " << totCovXP;
    }
    log << endreq;
  

  //Creating the TkrVertex instance:

  Point Pt(Vtx[0],Vtx[1],Vtx[2]);
  double totE = totP.magnitude();
  Ray theRay(Pt,totP);
  Event::TkrVertex* 
    vertex = new Event::TkrVertex(usedTracks.front()->getLayer(),
                                  usedTracks.front()->getTower(),
                                  totE, 0.,
                                  theRay
                                  );

  usedIter = usedTracks.begin();
  // Comment out for now (changing from TkrFitTrackBase to TkrTrack - 8/26/04 TU)
  //while(usedIter != usedTracks.end()) vertex->addTrack(*usedIter++);

  vertex->writeOut(log);

  VtxCol.push_back(vertex); 
*/
  return StatusCode::SUCCESS;
}


StatusCode  VtxKalFitTool::initVertex(Event::TkrTrackCol& /*theTracks*/)
{
  // Purpose and Method: defines starting vertex estimate as first hit on 
  //                     first (aka best) track
  //                     and starting vertex Cov. estimate as Cov. matrix at 
  //                     the latter.
  //                     Also defines reference surface for computation of 
  //                     measurement equation:
  //                     here we take the horizontal plane located at 
  //                     z = bestTrack->firstHit->z 
  // Inputs: TkrFitTrackCol object
  // Output: StatusCode upon completion. Starting vertex and Cov., and 
  //         reference plane are stored in member data
  // Dependencies: None
  //
  // Restrictions and Caveats: First track in TkrFitTrackCol needs to be the 
  //                           best track.
  //                           Error on Z-coordinate needs to be iterated and 
  //                           should come from TkrGeoSvc.
/*
  Point p0 = theTracks.front()->getPosition();
  CLHEP::HepVector vtx0(3);
  vtx0[0] = p0.x();
  vtx0[1] = p0.y();
  vtx0[2] = p0.z();

  m_VtxEstimates.push_back(vtx0);

  m_Zref = theTracks.front()->getTrackParZ();

  //Now the Cov matrix:
  Event::TkrFitMatrix  trkCov = theTracks.front()->getTrackCov();
  HepSymMatrix Cov0(3,0); 
  Cov0(1,1) = trkCov.getcovX0X0();
  Cov0(2,2) = trkCov.getcovY0Y0();
  Cov0(1,2) = trkCov.getcovX0Y0();
  Cov0(2,1) = trkCov.getcovX0Y0();
  Cov0(3,3) = 0.4/sqrt(12.0); //waffer thickness = 400um

  m_VtxCovEstimates.push_back(Cov0);
*/
  return StatusCode::SUCCESS;
}

CLHEP::HepVector VtxKalFitTool::computeVectorH(const CLHEP::HepVector x, const CLHEP::HepVector q)
{
  // Purpose and Method: 
  //Compute measurement vector H(x,q)
  //Kalman vertexer technique relies on the definition of a
  //measurement equation: m = H(x,q) + v, where 
  // - m is the measurement vector of the track: m=(X,Sx,Y,Sy);
  // - v is the measurement noise (assumed to be unbasied, 
  //   of finite variance and independent on the other tracks)
  // - x is the vertex position;
  // - q is the geometrical momentum of the track at vertex.
  // In our case q = (Sx,Sy,E) is constant along the track (apart maybe for E).
  // To obtain H, one must choose a reference surface where m is assumed
  // to have been measured. This is chosen to be the horizontal 
  // plane containing the first hit of the best track
  // (it should be close enough to the vertex), as implemented 
  // in method VtxKalFitTool::initVertex. As a result, the measurement vector 
  // H(v,q) is simply (with v=(x,y,z) and q=(Sx,Sy,E)):
  //  H_0 = x + Sx*(Zref - z)
  //  H_1 = Sx
  //  H_2 = y + Sy*(Zred - z)
  //  H_3 = Sy
  //  H_4 = E
  //                     
  // Inputs: In practice x is the point of linearization, 
  // and q is the geometrical momentum of current track close to x.
  // Output: Measurement vector H (dimension 5 as (X,Sx,Y,Sy,E).)
  // Dependencies: None
  //
  // Restrictions and Caveats: We might need to generalize 
  //                           w.r.t sign(uz), and if E can be correlated...
  
  CLHEP::HepVector H(5,0);
  H[0] = x[0] + q[0]*(m_Zref-x[2]);        // X
  H[1] = q[0];                             // SX
  H[2] = x[1] + q[1]*(m_Zref-x[2]);        // Y
  H[3] = q[1];                             // SY
  H[4] = q[2];                             // E
  return H;
}


CLHEP::HepMatrix VtxKalFitTool::computeMatrixA(const CLHEP::HepVector /*x*/, const CLHEP::HepVector q)
{
  // Purpose and Method:  matrix A is simply {partial H}/{partial x}
  // Inputs: position and momentum at linearization point
  // Output: derivative matrix A
  // Dependencies: None
  //
  // Restrictions and Caveats: None
  
  CLHEP::HepMatrix A(5,3,0);
  A(1,1) = 1;
  A(1,3) = -q[0];
  A(3,2) = 1;
  A(3,3) = -q[1];
  return A;
}

CLHEP::HepMatrix VtxKalFitTool::computeMatrixB(const CLHEP::HepVector x, const CLHEP::HepVector /*q*/)
{
  // Purpose and Method:  matrix B is simply {partial H}/{partial q}
  // Inputs: position x and momentum q at linearization point
  // Output: derivative matrix B
  // Dependencies: None
  //
  // Restrictions and Caveats: None

  CLHEP::HepMatrix B(5,3,0);
  B(1,1) = m_Zref - x[2];
  B(2,1) = 1;
  B(3,2) = m_Zref - x[2];
  B(4,2) = 1;
  B(5,3) = 1;
  return B;
}


CLHEP::HepVector VtxKalFitTool::getTkrParVec(const Event::TkrTrack& /*theTrack*/)
{
  // Purpose and Method: Simple translation of TkrFitPar (+ energy) 
  //                     into an CLHEP::HepVector
  // Inputs: fitted Track
  // Output: CLHEP::HepVector (X,Sy,Y,Sy,E)
  // Dependencies: None
  //
  // Restrictions and Caveats: None

///  Event::TkrFitPar par0 = theTrack.getTrackPar();
  CLHEP::HepVector m(5,0);
///  m(1) = par0.getXPosition();
///  m(2) = par0.getXSlope();
///  m(3) = par0.getYPosition();
///  m(4) = par0.getYSlope();
///  m(5) = theTrack.getEnergy();
  return m;
}


CLHEP::HepVector VtxKalFitTool::computeQatVtx(const Event::TkrTrack& /*theTrack*/,
                                       const CLHEP::HepVector /*theVertex*/)
{
  // Purpose and Method: Simple building of the geometrical momentum (Sx,Sy,E).
  //                     In general it should be computed at the POCA to the 
  //                     vertex estimate.
  //                     In practice Sx,Sy don't change! E might need 
  //                     reevaluation...
  // Inputs: TkrFitTrack
  // Output: geometrical momentum (Sx,Sy,E)
  // Dependencies: None
  //
  // Restrictions and Caveats: No attempt to play with energy, maybe not even 
  // necessary
  CLHEP::HepVector q(3);
///  q[0] = theTrack.getTrackPar().getXSlope();
///  q[1] = theTrack.getTrackPar().getYSlope();
///  q[2] = theTrack.getEnergy();
  return q;
}


//This is just a reformatting of TkrFitMatrix into HepSymMatrix
CLHEP::HepSymMatrix VtxKalFitTool::getHepSymCov(const Event::TkrTrackParams& measCov)
{
  CLHEP::HepSymMatrix G(4,0);
  G(1,1) = measCov(1,1);
  G(2,2) = measCov(2,2);
  G(1,2) = measCov(1,2);
  G(2,1) = measCov(2,1);
  
  G(3,3) = measCov(3,3);
  G(4,4) = measCov(4,4);
  G(3,4) = measCov(3,4);
  G(4,3) = measCov(4,3);
  
  G(1,3) = measCov(1,3);
  G(1,4) = measCov(1,4);
  G(4,1) = measCov(4,1);
  
  G(3,1) = measCov(3,1);
  G(3,2) = measCov(3,2);
  G(2,3) = measCov(2,3);
  
  G(2,4) = measCov(2,4);
  G(4,2) = measCov(4,2);
  return G;
}


CLHEP::HepMatrix VtxKalFitTool::SlopeToDir(CLHEP::HepVector Q, int sign_uz)
{
  // Purpose and Method: Transformation Matrix T (Sx,Sy,E)->(Eux,Euy,Euz)
  //                     newCov = T*oldCov*T.T() and T_{ij} = dy_i/dx_j   
  //                     (i for line and j for column of matrix)
  //                     where x and y are respectively the vector of old and 
  //                     new parameters.
  // Inputs: geometrical momentum Q=(Sx,Sy,E), sign_uz = uz/|uz|
  // Output: Matrix T
  // Dependencies: None
  //
  // Restrictions and Caveats: This is an approximation working only if the 
  //                           terms left out in the Taylor expansion
  //                           are small compared to the cov. matrix elements.
      double sx = Q[0];
      double sy = Q[1];
      double  E = Q[2];

      CLHEP::HepMatrix T(3,3,0);
      // partial{Euz} / partial{Sy} up to norm. factor:
      T(1,1) = E * (1 + sy*sy);
      T(1,2) = - E * sx*sy;
      T(1,3) = sx * (1 + sx*sx + sy*sy);
      // partial{Euy} / partial{Sy} up to norm. factor:
      T(2,1) = - E * sx*sy;
      T(2,2) = E * (1 + sx*sx);
      T(2,3) = sy * (1 + sx*sx + sy*sy);
      // partial{Euz} / partial{Sy} up to norm. factor:
      T(3,1) = - E * sx;
      T(3,2) = - E * sy;
      T(3,3) = (1 + sx*sx + sy*sy);
      // norm. factor:
      T *= (sign_uz) / pow(1 + sx*sx + sy*sy,1.5);
      
      return T;
}

Event::TkrTrackParams 
VtxKalFitTool::propagCovToVtx(const Event::TkrTrackParams Cov, 
                              const CLHEP::HepVector /*Vtx*/)
{
  // Purpose and Method: Propagate Cov matrix to vicinity of current 
  //                     vertex estimate
  // Inputs: Cov Matrix at first hit location
  // Output: HepSymMatrix translation of Cov matrix, propagated to POCA 
  //         of track w.r.t. vertex
  // Dependencies: None
  //
  // Restrictions and Caveats:  NOT YET IMPLEMENTED!
  return Cov;
}


CLHEP::HepSymMatrix 
VtxKalFitTool::computeWeightMatrix(const Event::TkrTrack& theTrack,
                                   const CLHEP::HepVector /*Vtx*/)
{
  // Purpose and Method: Computation of the weight matrix in 3 steps:
  //                     1) Propagate Cov. matrix to the Poca w.r.t Vtx;
  //                     2) Invert it;
  //                     3) estimate energy error and returns final Cov;
  // Inputs: TkrFitTrack and current estimate of the vertex
  // Output: Final Cov matrix in HepSymMatrix format
  // Dependencies: None
  //
  // Restrictions and Caveats: Energy error might need reviewing
  
  MsgStream log(msgSvc(), name());

  //first bring Cov(X,Sx,Y,Sy) close to current vertex
///  Event::TkrFitMatrix measCov = propagCovToVtx(theTrack.getTrackCov(), Vtx);
  //Then get the weight matrix of this part:
///  measCov.invert(ifail);
///  if(ifail)
///    log << MSG::ERROR <<"ERROR INVERTING measCov!"<<endreq;
 
  //We clearly assumes here that Energy is independant:
  CLHEP::HepSymMatrix G(5,0);
///  G.sub(1,getHepSymCov(measCov));
   
  //Now add the energy error:
  double E = theTrack.getInitialEnergy();
  if(E>1000000) 
    {
      log << MSG::WARNING << "Track Energy Not determined: will put huge errors" << endreq;
      G(5,5) = 0.001; //arbitrary but should be enough to kill it
    }
  else
    {
      int nHits = theTrack.getNumHits();
      double delta_E = E/sqrt((double)nHits/2);
      G(5,5) = 1/delta_E;
    }
  return G;
}
