
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


StatusCode VtxKalFitTool::doVtxFit(Event::TkrVertexCol& VtxCol)
{
  m_VtxCovEstimates.erase(m_VtxCovEstimates.begin(),m_VtxCovEstimates.end());
  m_VtxEstimates.erase(m_VtxEstimates.begin(),m_VtxEstimates.end());

  std::vector<HepSymMatrix>  CovQQ;
  std::vector<HepMatrix>     CovXQ;
  std::vector<HepVector>     fQ_list;
  std::vector<HepVector>     m_list;
  std::vector<HepVector>     sQ_list;
  std::vector<HepMatrix>     Tmp_list;

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

  log<<MSG::DEBUG<<"NUMBER OF TRACKS TO VERTEX: "<<theTracks->size()<<endreq;
  //***********************
  //FILTER STEP:
  //***********************

  sc = StatusCode::FAILURE;

  Event::TkrFitConPtr tkrIter = theTracks->begin();
  while(tkrIter != theTracks->end())
    { 
      Event::TkrFitTrack* theTrack  = *tkrIter++;

      if(theTrack->getQuality()<0) continue;

      std::cout<<"INFO OF NEW TRACK: "<<endl;
      std::cout<<"Direction "<<theTrack->getDirection()<<endl;
      std::cout<<"Position " <<theTrack->getPosition()<<endl;
      std::cout<<"Energy "   <<theTrack->getEnergy()<<endl;
      std::cout<<"# HITS "   <<theTrack->getNumHits()<<endl;
      std::cout<<"Quality "  <<theTrack->getQuality()<<endl;


      //Step 1: define new point for linearization:
      //-------------------------------------------

      //current estimations of vertex and its Cov. Matrix:
      HepVector    Vtx = m_VtxEstimates.back();
      HepSymMatrix invC = m_VtxCovEstimates.back();
      invC.invert(ifail);
      if(ifail) 
        {
          log << MSG::ERROR <<"ERROR INVERTING C!"<<endreq;
          return sc;
        }

      // Track parameters (X,Sx,Y,Sy,E) as Measurement Vector m
      HepVector m = getTkrParVec(*theTrack);

      // get "momentum" vector (Sx,Sy,E) of track at
      // POCA to vertex (conventional)
      HepVector    q = computeQatVtx(*theTrack,Vtx); 

      HepMatrix A  = computeMatrixA(Vtx,q);
      HepMatrix B  = computeMatrixB(Vtx,q);
      HepVector h0 = computeVectorH(Vtx,q);
      HepVector c0 = h0 - A*Vtx - B*q;

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
      HepSymMatrix newC = invC + Gb.similarityT(A); // (Ck+1)^-1 = Ck^-1 + A^T*Gb*A
      newC.invert(ifail);
      if(ifail)
        {
          log << MSG::ERROR <<"UPDATED MATRIX C NOT INVERTIBLE!"<<endreq;
          return StatusCode::FAILURE;
        }
      
      // D = Cov(Q,Q)     E = cov(x,Q)
      HepSymMatrix D = W + newC.similarity(W*B.T()*G*A);      
      HepMatrix    E = -newC*A.T()*G*B*W;

      
      //Step 3: Update vertex estimation and check chi2
      //-----------------------------------------------
      //new vertex and momentum at vertex for current track:
      HepVector newVtx = newC*(invC*Vtx + A.T()*Gb*(m-c0));
      HepVector newQ   = W*B.T()*G*(m-c0-A*newVtx);      


      HepVector r = m - c0 - A*newVtx - B*newQ;
      HepVector dx = newVtx - Vtx;
      double chi2f = G.similarity(r) + invC.similarity(dx);// returns r.T()*G*r + dx.T()*C*dx SCALAR
      
      std::cout<<"chi2f= "<<chi2f<<", and chi2max= "<<m_chi2max<<endl;
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

      
  log<<MSG::DEBUG<<"Number of tracks used: "<<usedTracks.size()<<endreq;
  //**************************
  //SMOOTHER STEP:
  //**************************

  //Final estimates:
  HepVector Vtx      = m_VtxEstimates.back();
  HepSymMatrix CovXX = m_VtxCovEstimates.back();

  Event::TkrFitConPtr usedIter = usedTracks.begin();
  int i=0;
  while(usedIter != usedTracks.end())
    { 
      Event::TkrFitTrack* theUsedTrack = *usedIter++;

      HepMatrix A  = computeMatrixA(Vtx, fQ_list[i]);
      HepMatrix B  = computeMatrixB(Vtx, fQ_list[i]);
      HepVector h0 = computeVectorH(Vtx, fQ_list[i]);
      HepVector c0 = h0 - A*Vtx - B*fQ_list[i];
      
      Event::TkrFitMatrix measCov = propagCovToVtx(theUsedTrack->getTrackCov(), Vtx);
      measCov.invert(ifail);
      HepSymMatrix G = computeWeightMatrix(*theUsedTrack,Vtx);
      HepSymMatrix W = G.similarityT(B); //W^{-1} = B^T*G*B
      W.invert(ifail);

      //Final Q=(Sx,Sy)
      HepVector    Q = W*B.T()*G*(m_list[i] - c0 - A*Vtx);
      //Final Cov(Q)
      HepSymMatrix D = W + CovXX.similarity(W*B.T()*G*A);
      //Final Cov(Vtx,Q) 
      HepMatrix    E = -CovXX*A.T()*G*B*W;
      //help for cov(Q_i,Q_j) = Tmp*CovXX*Tmp.T()
      HepMatrix Tmp = W*B.T()*G*A;

      sQ_list.push_back(Q);
      CovQQ.push_back(D);
      CovXQ.push_back(E);
      Tmp_list.push_back(Tmp);

      i++;
    }

  //**************************
  // Physical Momentum and Cov
  //**************************
  Vector totP;
  HepSymMatrix totCovP(3,0);
  HepMatrix totCovXP(3,3,0);

  usedIter = usedTracks.begin();
  i=0;
  while(usedIter != usedTracks.end())
    { 
      Event::TkrFitTrack* theUsedTrack = *usedIter++;
      double sgn = theUsedTrack->getDirection().z();
      int e = sgn>0?+1:-1;

      HepVector Qi = sQ_list[i]; 

      //Physical Momentum:
      //-----------------
      Vector momentum = Vector(e*Qi[0],e*Qi[1],e*1).unit();
      momentum.setMag(theUsedTrack->getEnergy());
      totP += momentum;

      //Cov Matrix of Physical Momentum:
      //-------------------------------
      HepMatrix Ti = SlopeToDir(Qi);
      totCovP += CovQQ[i].similarityT(Ti);

      //std::cout<<"CovQQ "<<CovQQ[i]<<endl;
      //std::cout<<"Matrix T "<<Ti<<endl;
      //std::cout<<"totCovP "<<totCovP<<endl;
      
      
      HepSymMatrix Cov_ij(3,0);
      int j;
      for(j=i+1;j<usedTracks.size();j++)
	{
	  HepVector Qj = sQ_list[j];
	  HepMatrix Tj = SlopeToDir(Qj);
	  
	  HepMatrix Qij = Ti.T()*Tmp_list[i]*CovXX*Tmp_list[j].T()*Tj;
	  //	  std::cout<<Qj<<Tj<<Qij<<endl;    
	  HepSymMatrix tmp;
	  tmp.assign(Qij + Qij.T()); 
	  Cov_ij += tmp;

	  totCovP += Cov_ij;
	}
      //      std::cout<<totCovXP<<CovXQ[i]<<endl;

      totCovXP += CovXQ[i]*Ti;
      i++;
    }
  //  std::cout<<"totP "<<totP<<endl;
  //  std::cout<<"totCovP "<<totCovP<<endl;
  

  //********************************
  //Creating the TkrVertex instance:
  //********************************
  Point Pt(Vtx[0],Vtx[1],Vtx[2]);
  double totE = totP.magnitude();
  Ray theRay(Pt,totP);
  Event::TkrVertex* vertex = new Event::TkrVertex(usedTracks.front()->getLayer(),
						  usedTracks.front()->getTower(),
						  totE, 0.,
						  theRay
						  );

  tkrIter = usedTracks.begin();
  while(tkrIter != usedTracks.end()) vertex->addTrack(*tkrIter++);

  vertex->writeOut(log);

  VtxCol.push_back(vertex); 

  return StatusCode::SUCCESS;
}


StatusCode  VtxKalFitTool::initVertex(Event::TkrFitTrackCol& theTracks)
{
  //use the first hit of the first (aka best) track as initial Vtx estimate
  Point p0 = theTracks.front()->getPosition();
  HepVector vtx0(3);
  vtx0[0] = p0.x();
  vtx0[1] = p0.y();
  vtx0[2] = p0.z();

  m_VtxEstimates.push_back(vtx0);

  m_Zref = theTracks.front()->getTrackParZ();

  //Now the Cov matrix:
  Event::TkrFitMatrix  trkCov = theTracks.front()->getTrackCov();
  HepSymMatrix Cov0(3,0); //start with 0 matrix
  Cov0(1,1) = trkCov.getcovX0X0();
  Cov0(2,2) = trkCov.getcovY0Y0();
  Cov0(1,2) = trkCov.getcovX0Y0();
  Cov0(2,1) = trkCov.getcovX0Y0();
  Cov0(3,3) = 0.4/sqrt(12); //waffer thickness = 400um

  m_VtxCovEstimates.push_back(Cov0);

  return StatusCode::SUCCESS;
}

HepVector VtxKalFitTool::computeVectorH(const HepVector x, const HepVector q)
{
  HepVector H(5,0);
  H[0] = x[0] + q[0]*(m_Zref-x[2]);        // X
  H[1] = q[0];                             // SX
  H[2] = x[1] + q[1]*(m_Zref-x[2]);        // Y
  H[3] = q[1];                             // SY
  H[4] = q[2];                             // E
  return H;
}

HepMatrix VtxKalFitTool::computeMatrixA(const HepVector x, const HepVector q)
{
  HepMatrix A(5,3,0.);
  A(1,1) = 1;
  A(1,3) = -q[0];
  A(3,2) = 1;
  A(3,3) = -q[1];
  return A;
}


HepMatrix VtxKalFitTool::computeMatrixB(const HepVector x, const HepVector q)
{
  HepMatrix B(5,3,0.);
  B(1,1) = m_Zref - x[2];
  B(2,1) = 1;
  B(3,2) = m_Zref - x[2];
  B(4,2) = 1;
  B(5,3) = 1;
  return B;
}


HepVector VtxKalFitTool::getTkrParVec(const Event::TkrFitTrack& theTrack)
{
  // par0 = (X,SX,Y,SY,E) 
  Event::TkrFitPar par0 = theTrack.getTrackPar();
  HepVector m(5,0);
  m(1) = par0.getXPosition();
  m(2) = par0.getXSlope();
  m(3) = par0.getYPosition();
  m(4) = par0.getYSlope();
  m(5) = theTrack.getEnergy();
  return m;
}


HepVector VtxKalFitTool::computeQatVtx(const Event::TkrFitTrack& theTrack,const HepVector theVertex)
{
  HepVector q(3);
  q[0] = theTrack.getTrackPar().getXSlope();
  q[1] = theTrack.getTrackPar().getYSlope();
  q[2] = theTrack.getEnergy();
  return q;
}


HepSymMatrix VtxKalFitTool::getHepSymCov(const Event::TkrFitMatrix& measCov)
{
  HepSymMatrix G(4,0);
  G(1,1) = measCov.getcovX0X0();
  G(2,2) = measCov.getcovSxSx();
  G(1,2) = measCov.getcovX0Sx();
  G(2,1) = measCov.getcovSxX0();
  
  G(3,3) = measCov.getcovY0Y0();
  G(4,4) = measCov.getcovSySy();
  G(3,4) = measCov.getcovY0Sy();
  G(4,3) = measCov.getcovSyY0();
  
  G(1,3) = measCov.getcovX0Y0();
  G(1,4) = measCov.getcovX0Sy();
  G(4,1) = measCov.getcovSyX0();
  
  G(3,1) = measCov.getcovY0X0();
  G(3,2) = measCov.getcovY0Sx();
  G(2,3) = measCov.getcovSxY0();
  
  G(2,4) = measCov.getcovSxSy();
  G(4,2) = measCov.getcovSySx();
  return G;
}


HepMatrix VtxKalFitTool::SlopeToDir(HepVector Q)
{
      double sx = Q[0];
      double sy = Q[1];
      double  E = Q[2];

      //Transformation: (Sx,Sy,E) -> (Eux,Euy,Euz)
      //sgn_uz discarded for now: should not matter,
      //even when transforming skew matrices....
      HepMatrix T(3,3,0);
      T(1,1) = 1+sy*sy;
      T(1,2) = -sx*sy;
      T(1,3) = -sx;
      T(2,1) = -sx*sy;
      T(2,2) = 1+sx*sx;
      T(2,3) = -sy;
      T *= E;
      T(3,1) = sx*(1+sx*sx+sy*sy);
      T(3,2) = sy*(1+sx*sx+sy*sy);
      T(3,3) =    (1+sx*sx+sy*sy);
      T /= pow(1+sx*sx+sy*sy,1.5);
      
      return T;
}


Event::TkrFitMatrix VtxKalFitTool::propagCovToVtx(const Event::TkrFitMatrix Cov, 
						  const HepVector Vtx)
{
//not implemented yet...
  return Cov;
}

HepSymMatrix VtxKalFitTool::computeWeightMatrix(const Event::TkrFitTrack& theTrack,const HepVector Vtx)
{
  int ifail;
  MsgStream log(msgSvc(), name());

  //first bring Cov(X,Sx,Y,Sy) close to current vertex
  Event::TkrFitMatrix measCov = propagCovToVtx(theTrack.getTrackCov(), Vtx);
  //Then get the weight matrix of this part:
  measCov.invert(ifail);
  if(ifail)
    log << MSG::ERROR <<"ERROR INVERTING measCov!"<<endreq;
 
  //We clearly assumes here that Energy is independant:
  HepSymMatrix G(5,0);
  G.sub(1,getHepSymCov(measCov));
   
  //Now add the energy error:
  double E = theTrack.getEnergy();
  if(E>1000000) 
    {
      log<<MSG::WARNING<<"Track Energy Not determined: will put huge errors"<<endreq;
      G(5,5) = 0.001; //arbitrary but should be enough to kill it
    }
  else
    {
      int nHits = theTrack.getNumHits();
      double delta_E = E/sqrt(nHits/2);
      G(5,5) = 1/delta_E;
    }
  return G;
}
