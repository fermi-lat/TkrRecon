
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
  : VtxBaseTool(type, name, parent), m_chi2(0.0)
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
  m_chi2f.erase(m_chi2f.begin(),m_chi2f.end());
  m_c0.erase(m_c0.begin(),m_c0.end());
  m_m.erase(m_m.begin(),m_m.end());
  m_Q.erase(m_Q.begin(),m_Q.end());
  m_CovQQ.erase(m_CovQQ.begin(),m_CovQQ.end());
  m_CovXQ.erase(m_CovXQ.begin(),m_CovXQ.end());
  m_Skew.erase(m_Skew.begin(),m_Skew.end());
  m_A.erase(m_A.begin(),m_A.end());
  m_B.erase(m_B.begin(),m_B.end());
  m_C.erase(m_C.begin(),m_C.end());
  m_D.erase(m_D.begin(),m_D.end());
  m_E.erase(m_E.begin(),m_E.end());
  m_G.erase(m_G.begin(),m_G.end());
  m_W.erase(m_W.begin(),m_W.end());

  int    ifail      = 0; //flag for matrix inversion check

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
      
      //Step 1: define new point for linearization:
      //paragraph before eq(4): last estimated vertex and POCA of track to the former are good
      //new starting points for linearisation of the measurement equation:
      
      //current estimations of vertex and its Cov. Matrix:
      HepVector    x = m_VtxEstimates.back();
      HepSymMatrix C = m_VtxCovEstimates.back();

      //(X,SX,Y,SY) track parameters:
      Event::TkrFitPar par0 = theTrack->getTrackPar();

      HepVector par(4,0);
      par(1) = par0.getXPosition();
      par(2) = par0.getXSlope();
      par(3) = par0.getYPosition();
      par(4) = par0.getYSlope();

      //std::cout<<"INITIAL (UX,UY,UZ) OF TRACK "<<theTrack->getDirection()<<endl;
      // get "momentum" vector of track at ZRef
      HepVector    q = computeQatZref(*theTrack); 

      //Only C^(-1) is needed from now on:
      C.invert(ifail);
      if(ifail) 
        {
          log << MSG::ERROR <<"ERROR INVERTING C!"<<endreq;
          return sc;
        }

      // Now computes A, B and c0 at linearization point:
      HepMatrix A  = computeMatrixA(x,q);
      HepMatrix B  = computeMatrixB(x,q);
      HepVector h0 = computeVectorH(x,q);
      HepVector c0 = h0 - A*x - B*q;


      //Step 2: Compute Covariance matrices
      Event::TkrFitMatrix measCov = propagCovToVtx(theTrack->getTrackCov(),x);
      measCov.invert(ifail);
      if(ifail)
        {
          log << MSG::ERROR <<"ERROR INVERTING measCov!"<<endreq;
          return sc;
        }

      HepSymMatrix invG = getHepSymCov(measCov);
      HepSymMatrix W = invG.similarityT(B); //W = B^T*invG*B
      W.invert(ifail);
      if(ifail) 
        {
          log << MSG::ERROR <<"ERROR INVERTING W!"<<endreq;
          return StatusCode::FAILURE;
        }

      HepSymMatrix Gb = invG - W.similarity(invG*B);

      //Update vertex Covariance matrix:
      HepSymMatrix newC = C + Gb.similarityT(A); // (Ck+1)^-1 = Ck^-1 + A^T*Gb*A
      newC.invert(ifail);
      if(ifail)
        {
          log << MSG::ERROR <<"UPDATED MATRIX C NOT INVERTIBLE!"<<endreq;
          return StatusCode::FAILURE;
        }

      // D = Cov(Q,Q)     E = cov(x,Q)
      HepSymMatrix D = W + newC.similarity(W*B.T()*invG*A);      
      HepMatrix    E = -newC*A.T()*invG*B*W;


      //Step 3: Update vertex estimation and check chi2

      //new vertex and momentum at vertex for current track:
      HepVector newX = newC*(C*x + A.T()*Gb*(par-c0));
      HepVector newQ   = W*B.T()*invG*(par-c0-A*newX);      


      HepVector r = par - c0 - A*newX - B*newQ;
      HepVector dx = newX - x;
      double chi2f = invG.similarity(r) + C.similarity(dx);// returns r.T()*G*r + dx.T()*C*dx SCALAR

      //      std::cout<<"NEW VERTEX POSITION"<< newX<<endl;
      //      std::cout<<"NEW PARAM of current track"<< newQ<<endl;
      //      std::cout<<"CHI2f "<<chi2f<<endl;
      if(chi2f < m_chi2max) 
        {
	  usedTracks.push_back(theTrack);
          m_chi2 += chi2f;
	  m_A.push_back(A);
	  m_m.push_back(par);
	  m_c0.push_back(c0);
	  m_B.push_back(B);
	  m_G.push_back(invG);
	  m_W.push_back(W);

	  //keep new updates
          m_VtxEstimates.push_back(newX);
          m_VtxCovEstimates.push_back(newC);
	  sc = StatusCode::SUCCESS;
	}
    }//end of track while loop
  
  if(sc.isFailure()) 
    {
      log<<MSG::ERROR<<"KALMAN Filter FAILED!"<<endreq;
      return sc;
    }



  //**************************
  //SMOOTHER STEP:
  //**************************

  //Final estimates:
  HepVector    Vtx   = m_VtxEstimates.back();
  HepSymMatrix CovXX = m_VtxCovEstimates.back();

  Vector totQ;
  HepSymMatrix totCovQ(3,0);

  Event::TkrFitConPtr usedIter = usedTracks.begin();
  int i=0;
  while(usedIter != usedTracks.end())
    { 
      Event::TkrFitTrack* theUsedTrack = *usedIter++;
      
      double sgn = theUsedTrack->getDirection().z();
      
      //Final momentum (sx,sy) of track i at vertex:
      HepVector       Q = m_W[i]*m_B[i].T()*m_G[i]*(m_m[i]-m_c0[i]-m_A[i]*Vtx);
      //Final Cov Matrix of momentum of current track:
      HepSymMatrix newD = m_W[i] + CovXX.similarity(m_W[i]*m_B[i].T()*m_G[i]*m_A[i]);
      //Final Cov Matrix (Vtx,Q) between final vertex and momentum of current track
      HepMatrix    newE = -CovXX*m_A[i].T()*m_G[i]*m_B[i]*m_W[i];
      
      int j;
      HepSymMatrix covSkew(2,0);
      for(j=i+1;j<usedTracks.size();j++)
	{
	  HepMatrix Skew = m_W[i]*m_B[i].T()*m_G[i]*m_A[i]*CovXX*m_A[j].T()*m_G[j]*m_B[j]*m_W[j];
	  HepSymMatrix tmp;
	  tmp.assign(Skew+Skew.T()); //first need to change representation to (ux,uy,uz)
	  covSkew += tmp;
	}
      m_Skew.push_back(covSkew);
      
      //Physical Momentum:
      //-----------------
      int e = sgn>0?+1:-1;
      Vector momentum = Vector(e*Q[0],e*Q[1],e*1).unit();
      momentum.setMag(theUsedTrack->getEnergy());
      
      totQ += momentum;
      
      //Cov Matrix of Physical Momentum:
      //-------------------------------
      double sx = Q[0];
      double sy = Q[1];
      
      //Transformation: (Sx,Sy,E) -> (Eux,Euy,Euz)
      HepMatrix T(2,3,0);
      T(1,1) = 1+sy*sy;
      T(1,2) = sx*sy;
      T(1,3) = sx;
      T(2,1) = sx*sy;
      T(2,2) = 1+sx*sx;
      T(2,3) = sy;
      
      HepSymMatrix tmp = newD.similarityT(T);
      tmp /= (1+sx*sx+sy*sy);
            
      totCovQ += tmp;
      
      //all this and m_Skew should be fed to TkrVertex and the tracks
      m_Q.push_back(Q);
      m_CovQQ.push_back(newD);
      m_CovXQ.push_back(newE);
      
      //	  std::cout<<"FINAL (ux,uy,uz) OF TRACK "<<i<<":"<<momentum<<endl;
      //	  std::cout<<"FINAL ENERGY OF TRACK "<<i<<": "<<theUsedTrack->getEnergy()<<endl;
      i++;
    }
  
  
  //  std::cout<<"GAMMA DIR ()UX,UY,UZ COV MATRIX"<<totCovQ<<endl;
  

  //********************************
  //Creating the TkrVertex instance:
  //********************************
  Point Pt(Vtx[0],Vtx[1],Vtx[2]);
  double totE = totQ.magnitude();
  Ray theRay(Pt,totQ);
  Event::TkrVertex* vertex = new Event::TkrVertex(usedTracks.front()->getLayer(),
						  usedTracks.front()->getTower(),
						  totE, 0.,
						  theRay
						  );

  tkrIter = usedTracks.begin();
  while(tkrIter != theTracks->end())
    { 
      vertex->addTrack(*tkrIter++);
    }
  
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
  HepVector H(4,0);
  H[0] = x[0] + q[0]*(m_Zref-x[2]);        // X
  H[1] = q[0];                             // SX
  H[2] = x[1] + q[1]*(m_Zref-x[2]);        // Y
  H[3] = q[1];                             // SY

  return H;
}

HepMatrix VtxKalFitTool::computeMatrixA(const HepVector x, const HepVector q)
{
  HepMatrix A(4,3,0.);
  A(1,1) = 1;
  A(1,3) = -q[0];
  A(3,2) = 1;
  A(3,3) = -q[1];
  return A;
}


HepMatrix VtxKalFitTool::computeMatrixB(const HepVector x, const HepVector q)
{
  HepMatrix B(4,2,0.);
  B(1,1) = m_Zref - x[2];
  B(2,1) = 1;
  B(3,2) = m_Zref - x[2];
  B(4,2) = 1;
  return B;
}

HepVector VtxKalFitTool::computeQatZref(const Event::TkrFitTrack& theTrack)
{
  //This is the most involving part: need to bring parameters and 
  //errors at the reference plane m_Zref, but it might be absurd
  //if a less good track has a hit above the best track...

  //dummy for now:
  HepVector q(2);
  q[0] = theTrack.getTrackPar().getXSlope();
  q[1] = theTrack.getTrackPar().getYSlope();
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


Event::TkrFitMatrix VtxKalFitTool::propagCovToVtx(const Event::TkrFitMatrix Cov, 
						  const HepVector Vtx)
{
//not implemented yet...
  return Cov;
}
