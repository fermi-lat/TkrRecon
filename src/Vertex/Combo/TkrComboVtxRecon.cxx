/*
	Code to implement the Combo vertex finding class

	Tracy Usher March 1, 2002
*/

#include "src/Vertex/Combo/TkrComboVtxRecon.h"
#include "src/Vertex/Combo/RayDoca.h"

TkrComboVtxRecon::TkrComboVtxRecon(ITkrGeometrySvc* pTkrGeo, TkrFitTrackCol* pTracks, TkrPatCandCol* pCandTracks)
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
            if ((dist < 5. && doca.arcLenRay1() <= 10. && doca.arcLenRay2() <= 10.) ||
                (dist < 1. && t1t2 < .005)) {

                Point  gamPos;
                Vector gamDir;

/*              double gamEne = track1->getEnergy() + track2->getEnergy();
                
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
 */               
       //Use tracking covariances & energies to weight track directions
                TkrFitMatrix cov_t1 = track1->getTrackCov();
                TkrFitMatrix cov_t2 = track2->getTrackCov();
                TkrFitPar    par_t1 = track1->getTrackPar();
                TkrFitPar    par_t2 = track2->getTrackPar();
                double e_t1 =  track1->getEnergy();
                double e_t2 =  track2->getEnergy();
                double gamEne = e_t1 + e_t2;

                double norm_1 = sqrt(1.+par_t1.getXSlope()*par_t1.getXSlope() +
                                        par_t1.getYSlope()*par_t1.getYSlope());      
                double norm_2 = sqrt(1.+par_t2.getXSlope()*par_t2.getXSlope() +
                                        par_t2.getYSlope()*par_t2.getYSlope());  
                
                double tx_1 = -par_t1.getXSlope()/norm_1;
                double ty_1 = -par_t1.getYSlope()/norm_1;
                double tx_2 = -par_t2.getXSlope()/norm_2;
                double ty_2 = -par_t2.getYSlope()/norm_2;

                double dtx_1_sq = ((1.+par_t1.getYSlope()*par_t1.getYSlope())*
                                                         cov_t1.getcovSxSx() +
                                    par_t1.getXSlope()*par_t1.getYSlope()*
                                                         cov_t1.getcovSySy())/
                                                         pow(norm_1,3);
                double dty_1_sq =  ((1.+par_t1.getXSlope()*par_t1.getXSlope())*
                                                         cov_t1.getcovSySy() +
                                    par_t1.getXSlope()*par_t1.getYSlope()*
                                                         cov_t1.getcovSxSx())/
                                                         pow(norm_1,3);
                double dtx_2_sq = ((1.+par_t2.getYSlope()*par_t2.getYSlope())*
                                                         cov_t2.getcovSxSx() +
                                    par_t2.getXSlope()*par_t2.getYSlope()*
                                                         cov_t2.getcovSySy())/
                                                         pow(norm_2,3);
                double dty_2_sq =  ((1.+par_t2.getXSlope()*par_t2.getXSlope())*
                                                         cov_t2.getcovSySy() +
                                    par_t2.getXSlope()*par_t2.getYSlope()*
                                                         cov_t2.getcovSxSx())/
                                                         pow(norm_2,3);

                double wt_1x  = (e_t1/dtx_1_sq) / (e_t1/dtx_1_sq + e_t2/dtx_2_sq);
                double wt_2x  = (e_t2/dtx_2_sq) / (e_t1/dtx_1_sq + e_t2/dtx_2_sq);
                double wt_1y  = (e_t1/dty_1_sq) / (e_t1/dty_1_sq + e_t2/dty_2_sq);
                double wt_2y  = (e_t2/dty_2_sq) / (e_t1/dty_1_sq + e_t2/dty_2_sq);;
                   
                double gam_tx = tx_1*wt_1x + tx_2*wt_2x;
                double gam_ty = ty_1*wt_1y + ty_2*wt_2y;

                double gam_tz = sqrt(1. - gam_tx*gam_tx - gam_ty*gam_ty);

                gamDir = Vector(gam_tx, gam_ty, -gam_tz).unit();
                /*
                int i_error;
                cov_t1.invert(i_error);
                TkrFitMatrix wgh_1 = cov_t1;
  
                cov_t2.invert(i_error);
                TkrFitMatrix wgh_2 = cov_t2;

                TkrFitMatrix wgh_12 = wgh_1 + wgh_2;
                wgh_12.invert(i_error);
                TkrFitMatrix cov_12 = wgh_12; 
                TkrFitPar par_1(0.,
                                par_t1.getXSlope()*e_t1,
                                0.,
                                par_t1.getYSlope()*e_t1);
                TkrFitPar par_2(0.,
                                par_t2.getXSlope()*e_t2,
                                0.,
                                par_t2.getYSlope()*e_t2);

                TkrFitPar gam_par = cov_12*(wgh_1*par_1 + wgh_2*par_2);
                TkrFitMatrix gam_cov = cov_12; 
                
                gamDir = Vector(-gam_par.getXSlope()/gamEne, 
                                -gam_par.getYSlope()/gamEne, -1.).unit();
              */  
                gamPos  = doca.docaPointRay1();
                gamPos += doca.docaPointRay2();
                gamPos *= 0.5;
                
                Ray        gamma  = Ray(gamPos,gamDir);
                TkrVertex* vertex = new TkrVertex(track1->getLayer(),track1->getTower(),gamEne,dist,gamma);
 //               vertex->setDist1(doca.arcLenRay1());
 //               vertex->setDist2(doca.arcLenRay2())
//                vertex->setAngle(t1t2); 
                vertex->addTrack(track1);
                vertex->addTrack(track2);
                
                push_back(vertex); // addVertex(vertex);
                
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

            push_back(vertex); //addVertex(vertex);
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

