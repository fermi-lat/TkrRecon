// Implementation file for combining the first
// two tracks to estimate the gamma direction
// Bill Atwood Spring, 2002

#include "src/Vertex/Combo/TkrComboVtxRecon.h"
#include "src/Vertex/Combo/RayDoca.h"

double erf(double x) {
    double t = 1./(1.+.47047*x);
    double results = 1. -(.34802 - (.09587 -.74785*t)*t)*t*exp(-x*x);
    return results;
}
double thrshold(double x) {
    if(x < 0) return (.5*(1. + erf(-x)));
    else      return (.5*(1. - erf( x)));
}

TkrComboVtxRecon::TkrComboVtxRecon(ITkrGeometrySvc* /*pTkrGeo*/, TkrVertexCol* vertexCol, TkrFitTrackCol* pTracks, TkrPatCandCol* /*pCandTracks*/)
{
    //Define a vector to contain a list of "isolated" tracks
    int    numTracks = pTracks->size();
    std::vector<bool> unused(numTracks);
    
    while(numTracks--) unused[numTracks] = true;
    
    //Track counter
    int   trk1Idx = 0;
    
    double gamEne = 0.;
    
    
    TkrFitConPtr pTrack1 = pTracks->begin();
    TkrFitConPtr pTrack2 = pTrack1;
    pTrack2++;
    
    int          trk2Idx = trk1Idx+1;
    int      bst_trk2Idx = trk2Idx; 
    
    TkrFitTrack* track1  = *pTrack1;
    TkrFitTrack* best_track2 = 0;
    
    double e_t1 =  track1->getEnergy();
    Point  gamPos; 
    double best_doca; 
    
    double max_wgt = 0.; 
    const double root2 = sqrt(2.);
    
    //Loop over the number of Fit tracks to find the best 2nd to pair with #1
    //while(pTrack2 < pTracks->end())
    for (; pTrack2!=pTracks->end(); pTrack2++, trk2Idx++) 
    {
        TkrFitTrack* track2 = *pTrack2;
        
        if(gamEne == 0.) gamEne = e_t1+track2->getEnergy();
        
        RayDoca doca    = RayDoca(Ray(track1->getPosition(),track1->getDirection()),
                                  Ray(track2->getPosition(),track2->getDirection()));
        double dist = doca.docaRay1Ray2();
        double s1   = doca.arcLenRay1();
        double s2   = doca.arcLenRay2(); 
        
        double cost1t2 = track1->getDirection()*track2->getDirection();
        double t1t2    = acos(cost1t2);
        
        double q2       = track2->getQuality();   
        double wgt_doca = thrshold((dist-100./gamEne)/root2);  
        double wgt_s1   = thrshold((s1+.5)/(root2*1.5 * gamEne/100.));     
        double wgt_t1t2 = thrshold((t1t2-200./gamEne)/(root2*200./gamEne));
        double wgt_q2   = thrshold((40. - q2)/(root2*20.)); 
        
        double total_wgt = wgt_doca*wgt_s1*wgt_q2*wgt_t1t2; 
         if(total_wgt > max_wgt) {
            max_wgt = total_wgt;
            best_track2 = track2;
            gamPos  = doca.docaPointRay1();
            gamPos += doca.docaPointRay2();
            gamPos *= 0.5;
            best_doca = dist; 
            bst_trk2Idx = trk2Idx; 
        }
    }  
    
    //Check the overall weight
    if (best_track2 != 0 && max_wgt > .0005) {
        

        Vector gamDir;
 
        // Begin Method 1
        //Use tracking covariances & energies to weight track directions
        TkrFitMatrix cov_t1 = track1->getTrackCov();
        TkrFitMatrix cov_t2 = best_track2->getTrackCov();
        TkrFitPar    par_t1 = track1->getTrackPar();
        double sx_1= par_t1.getXSlope();
        double sy_1= par_t1.getYSlope();
        TkrFitPar    par_t2 = best_track2->getTrackPar();
        double sx_2= par_t2.getXSlope();
        double sy_2= par_t2.getYSlope();
        double e_t2 =  best_track2->getEnergy();
        
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
        

        double transistion = 1.; //thrshold((gamEne - 400.)/150.);
        
        double wt_x1 = 1./(dtx_1_sq*transistion + (1. - transistion));
        double wt_x2 = 1./(dtx_2_sq*transistion + (1. - transistion));       
        double wt_1x  = (e_t1*wt_x1) / (e_t1*wt_x1 + e_t2*wt_x2);
        double wt_2x  = (e_t2*wt_x2) / (e_t1*wt_x1 + e_t2*wt_x2);
        double wt_y1 = 1./(dty_1_sq*transistion + (1. - transistion));
        double wt_y2 = 1./(dty_2_sq*transistion + (1. - transistion));  
        double wt_1y  = (e_t1*wt_y1) / (e_t1*wt_y1 + e_t2*wt_y2);
        double wt_2y  = (e_t2*wt_y2) / (e_t1*wt_y1 + e_t2*wt_y2);
        
        double gam_tx = tx_1*wt_1x + tx_2*wt_2x;
        double gam_ty = ty_1*wt_1y + ty_2*wt_2y;     
        double gam_tz = sqrt(1. - gam_tx*gam_tx - gam_ty*gam_ty);
        
        gamDir = Vector(gam_tx, gam_ty, -gam_tz).unit();
        
        Ray        gamma  = Ray(gamPos,gamDir);
        TkrVertex* vertex = new TkrVertex(track1->getLayer(),track1->getTower(),gamEne,best_doca,gamma);
        //                vertex->setDist1(doca.arcLenRay1());
        //                vertex->setDist2(doca.arcLenRay2())
        //                vertex->setAngle(t1t2); 
        vertex->addTrack(track1);
        vertex->addTrack(best_track2);
        
        vertexCol->push_back(vertex); // addVertex(vertex);
        
        unused[trk1Idx] = false;
        unused[bst_trk2Idx] = false;
    }
            
            
            //Go through unused list looking for isolated tracks
    TkrFitConPtr pTrack = pTracks->begin();
    int          trkIdx = 0;
            
    //while(pTrack != pTracks->end())
    for (; pTrack != pTracks->end(); pTrack++, trkIdx++)
    {
         TkrFitTrack* track1 = *pTrack;
                
         if (unused[trkIdx])
         {
            TkrVertex* vertex = new TkrVertex(track1->getLayer(),track1->getTower(),track1->getEnergy(),0.,Ray(track1->getPosition(),track1->getDirection()));
                    
            vertex->addTrack(track1);
                    
            vertexCol->push_back(vertex); //addVertex(vertex);
         }
    }            
    return;
}


TkrComboVtxRecon::~TkrComboVtxRecon()
{   
    return;
}

