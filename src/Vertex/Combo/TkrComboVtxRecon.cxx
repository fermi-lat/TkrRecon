/*
    Code to implement the Combo vertex finding class

    Tracy Usher March 1, 2002
*/

#include "src/Vertex/Combo/TkrComboVtxRecon.h"
#include "src/Vertex/Combo/RayDoca.h"

TkrComboVtxRecon::TkrComboVtxRecon(ITkrGeometrySvc* /*pTkrGeo*/, TkrVertexCol* vertexCol, TkrFitTrackCol* pTracks, TkrPatCandCol* /*pCandTracks*/)
{
       //Define a vector to contain a list of "isolated" tracks
    int    numTracks = pTracks->size();
    bool*  unused    = new bool[numTracks];
    
    while(numTracks--) unused[numTracks] = true;
    
    //Track counter
    int   trk1Idx = 0;
    
    //Loop over the number of Fit tracks
    TkrFitConPtr pTrack1 = pTracks->begin();
    
    while(pTrack1 != pTracks->end())
    {
        TkrFitTrack* track1  = *pTrack1++;
        TkrFitConPtr pTrack2 = pTrack1;
        int          trk2Idx = trk1Idx +1;
        
        while(pTrack2 < pTracks->end())
        {
            TkrFitTrack* track2 = *pTrack2++;
            

            RayDoca doca    = RayDoca(Ray(track1->getPosition(),track1->getDirection()),
                                      Ray(track2->getPosition(),track2->getDirection()));
            double  dist    = doca.docaRay1Ray2();
            
            double cost1t2 = track1->getDirection()*track2->getDirection();
            double t1t2 = acos(cost1t2);  
            
            //Check that the DOCA is not too big
            if ((dist < 5. && doca.arcLenRay1() <= 15. && doca.arcLenRay2() <= 15.) ||
                (dist < 1.5 && t1t2 < .07)) {

                Point  gamPos;
                Vector gamDir;

                double gamEne = track1->getEnergy() + track2->getEnergy();
                
                Vector trk1Dir = track1->getDirection();
                Vector trk2Dir = track2->getDirection();
                
                trk1Dir.setMag(track1->getEnergy());
                trk2Dir.setMag(track2->getEnergy());
                
                double qual_1 = track1->getQuality();
                double qual_2 = track2->getQuality();
                qual_1 *= qual_1;
                qual_2 *= qual_2; 
                
                double w1 = qual_1/(qual_1+qual_2);
                double w2 = qual_2/(qual_1+qual_2);
                
                gamDir  = w1*trk1Dir + w2*trk2Dir;
                
                gamDir.setMag(1.);
/*               
       //Use tracking covariances & energies to weight track directions
                TkrFitMatrix cov_t1 = track1->getTrackCov();
                TkrFitMatrix cov_t2 = track2->getTrackCov();
                TkrFitPar    par_t1 = track1->getTrackPar();
                double sx_1= par_t1.getXSlope();
                double sy_1= par_t1.getYSlope();
                TkrFitPar    par_t2 = track2->getTrackPar();
                double sx_2= par_t2.getXSlope();
                double sy_2= par_t2.getYSlope();
                double e_t1 =  track1->getEnergy();
                double e_t2 =  track2->getEnergy();
                double gamEne = e_t1 + e_t2;

                double norm_1 = sqrt(1.+ sx_1*sx_1 + sy_1*sy_1); 
      
                double norm_2 = sqrt(1.+ sx_2*sx_2 + sy_2*sy_2);  
                
                double tx_1 = -sx_1/norm_1;
                double ty_1 = -sy_1/norm_1;
                double tx_2 = -sx_2/norm_2;
                double ty_2 = -sy_2/norm_2;

                double dtx_1_sq = ((1.+sy_1*sy_1)*(1.+sy_1*sy_1)*cov_t1.getcovSxSx() +
                                       sx_1*sx_1*sy_1*sy_1*cov_t1.getcovSySy())/pow(norm_1,6);
                double dty_1_sq = ((1.+sx_1*sx_1)*(1.+sx_1*sx_1)*cov_t1.getcovSySy() +
                                       sx_1*sx_1*sy_1*sy_1*cov_t1.getcovSxSx())/pow(norm_1,6);
                double dtx_2_sq = ((1.+sy_2*sy_2)*(1.+sy_2*sy_2)*cov_t2.getcovSxSx() +
                                       sx_2*sx_2*sy_2*sy_2*cov_t2.getcovSySy())/pow(norm_2,6);
                double dty_2_sq = ((1.+sx_2*sx_2)*(1.+sx_2*sx_2)*cov_t2.getcovSySy() +
                                       sx_2*sx_2*sy_2*sy_2*cov_t2.getcovSxSx())/pow(norm_2,6);

                double wt_1x  = (e_t1/dtx_1_sq) / (e_t1/dtx_1_sq + e_t2/dtx_2_sq);
                double wt_2x  = (e_t2/dtx_2_sq) / (e_t1/dtx_1_sq + e_t2/dtx_2_sq);
                double wt_1y  = (e_t1/dty_1_sq) / (e_t1/dty_1_sq + e_t2/dty_2_sq);
                double wt_2y  = (e_t2/dty_2_sq) / (e_t1/dty_1_sq + e_t2/dty_2_sq);;
                   
                double gam_tx = tx_1*wt_1x + tx_2*wt_2x;
                double gam_ty = ty_1*wt_1y + ty_2*wt_2y;

                double gam_tz = sqrt(1. - gam_tx*gam_tx - gam_ty*gam_ty);

                gamDir = Vector(gam_tx, gam_ty, -gam_tz).unit();


 

              // Second try using full cov. matrices... 
                TkrFitMatrix cov_t1 = track1->getTrackCov();
                TkrFitMatrix cov_t2 = track2->getTrackCov();
                TkrFitPar    par_t1 = track1->getTrackPar();
                TkrFitPar    par_t2 = track2->getTrackPar();
                double e_t1 =  track1->getEnergy();
                double e_t2 =  track2->getEnergy();

              // Recast matrices and parameters to elliminate postions
                int ierror; 
                TkrFitMatrix cs_t1; 
                cs_t1(1,1) = cov_t1(2,2);
                cs_t1(1,2) = cs_t1(2,1) = cov_t1(2,4);
                cs_t1(2,2) = cov_t1(4,4); 
                cs_t1(3,3) = cs_t1(4,4) = 1.;
                double nrm_1   = sqrt(1.+sx_1*sx_1 + sy_1*sy_1); 
                double nrm_1_3 = nrm_1*nrm_1*nrm_1;
                TkrFitPar pt_t1(e_t1*track1->getDirection().x(), 
                                e_t1*track1->getDirection().y(),
                                e_t1*track1->getDirection().z(),1.);

                TkrFitMatrix T1;
            
                T1(1,1) = e_t1*(1.+sy_1*sy_1)/nrm_1_3;
                T1(1,2) = T1(2,1) = -e_t1*sx_1*sy_1/nrm_1_3;
                T1(1,3) = e_t1*sx_1*(sx_1*sx_1 + sy_1*sy_1)/nrm_1_3;
                T1(2,2) = e_t1*(1.+sx_1*sx_1)/nrm_1_3;
                T1(2,3) = e_t1*sy_1*(sx_1*sx_1 + sy_1*sy_1)/nrm_1_3;
                T1(3,3) = e_t1/nrm_1; 
                T1(4,4) = 1; 

                TkrFitMatrix wt_t1 = T1 * cs_t1 * T1.T();
                wt_t1.inverse(ierror); 

                TkrFitMatrix cs_t2; 
                cs_t2(1,1) = cov_t2(2,2);
                cs_t2(1,2) = cs_t2(2,1) = cov_t2(2,4);
                cs_t2(2,2) = cov_t2(4,4); 
                cs_t2(3,3) = cs_t2(4,4) = 1.;
                double sx_2= par_t2.getXSlope();
                double sy_2= par_t2.getYSlope();
                double nrm_2   = sqrt(1.+sx_2*sx_2 + sy_1*sy_2); 
                double nrm_2_3 = nrm_2*nrm_2*nrm_1;
                TkrFitPar pt_t2(e_t2*track2->getDirection().x(), 
                                e_t2*track2->getDirection().y(),
                                e_t2*track2->getDirection().z(),1.);

                TkrFitMatrix T2;
            
                T2(1,1) = e_t2*(1.+sy_2*sy_2)/nrm_2_3;
                T2(1,2) = T2(2,1) = -e_t2*sx_2*sy_2/nrm_2_3;
                T2(1,3) = e_t2*sx_2*(sx_2*sx_2 + sy_2*sy_2)/nrm_2_3;
                T2(2,2) = e_t2*(1.+sx_2*sx_2)/nrm_2_3;
                T2(2,3) = e_t2*sy_2*(sx_2*sx_2 + sy_2*sy_2)/nrm_2_3;
                T2(3,3) = e_t2/nrm_2; 
                T2(4,4) = 1.; 

                TkrFitMatrix wt_t2 = T2 * cs_t2 * T2.T();
                wt_t2.inverse(ierror); 

                TkrFitMatrix cov_gam = wt_t1+wt_t2;
                cov_gam.inverse(ierror);

                TkrFitPar p_gam = cov_gam * (wt_t1*pt_t1 + wt_t2*pt_t2); 

               gamDir = Vector(p_gam.getXPosition(), p_gam.getXSlope(), p_gam.getYPosition()).unit();
*/
//                double gamEne = e_t1 + e_t2;
              
                gamPos  = doca.docaPointRay1();
                gamPos += doca.docaPointRay2();
                gamPos *= 0.5;
                
                Ray        gamma  = Ray(gamPos,gamDir);
                TkrVertex* vertex = new TkrVertex(track1->getLayer(),track1->getTower(),gamEne,dist,gamma);
//                vertex->setDist1(doca.arcLenRay1());
//                vertex->setDist2(doca.arcLenRay2())
//                vertex->setAngle(t1t2); 
                vertex->addTrack(track1);
                vertex->addTrack(track2);
                
                vertexCol->push_back(vertex); // addVertex(vertex);
                
                unused[trk1Idx] = false;
                unused[trk2Idx] = false;
            }
        }
        trk1Idx++;
    }
    
    
    //Go through unused list looking for isolated tracks
    TkrFitConPtr pTrack = pTracks->begin();
    int          trkIdx = 0;
    
    while(pTrack != pTracks->end())
    {
        TkrFitTrack* track1 = *pTrack++;

        if (unused[trkIdx++])
        {
            TkrVertex* vertex = new TkrVertex(track1->getLayer(),track1->getTower(),track1->getEnergy(),0.,Ray(track1->getPosition(),track1->getDirection()));

            vertex->addTrack(track1);

            vertexCol->push_back(vertex); //addVertex(vertex);
        }
    }

    //Don't leave anything dangling
    delete unused;

    return;
}


TkrComboVtxRecon::~TkrComboVtxRecon()
{

    return;
}

