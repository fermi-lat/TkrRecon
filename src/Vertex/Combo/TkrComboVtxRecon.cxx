// Implementation file for combining the first
// two tracks to estimate the gamma direction
// Bill Atwood Spring, 2002

#include "src/Vertex/Combo/TkrComboVtxRecon.h"
#include "src/Vertex/Combo/RayDoca.h"
#include "Event/Recon/TkrRecon/TkrVertexTab.h"

#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"

double fast_erf(double x) {
    double t = 1./(1.+.47047*x);
    double results = 1. -(.34802 - (.09587 -.74785*t)*t)*t*exp(-x*x);
    return results;
}
double thrshold(double x) {
    if(x < 0) return (.5*(1. + fast_erf(-x)));
    else      return (.5*(1. - fast_erf( x)));
}

using namespace Event;

TkrComboVtxRecon::TkrComboVtxRecon(ITkrGeometrySvc* /*tkrGeom*/, 
                                   Event::TkrVertexCol* vertexCol, 
                                   Event::TkrTrackCol* pTracks, 
                                   Event::TkrPatCandCol* /*pCandTracks*/,
                                   Event::TkrVertexTrackTab* vertexRelTab)
{
    //Define a vector to contain a list of "isolated" tracks
    int    numTracks = pTracks->size();
    std::vector<bool> unused(numTracks);
    
    while(numTracks--) unused[numTracks] = true;
    
    //Track counter
    int   trk1Idx = 0;
    
    double gamEne = 0.;
    
    
    TkrTrackColPtr pTrack1 = pTracks->begin();
    TkrTrackColPtr pTrack2 = pTrack1;
    pTrack2++;
    
    int          trk2Idx = trk1Idx+1;
    int      bst_trk2Idx = trk2Idx; 
    
    TkrTrack* track1      = *pTrack1;
    TkrTrack* best_track2 = 0;

    /// --* What do we put here? Use initial energy for now
    double e_t1 =  track1->getInitialEnergy();
    Point  gamPos; 
    double best_doca; 
    
    double max_wgt = 0.; 
    const double root2 = sqrt(2.);
    
    //Loop over the number of Fit tracks to find the best 2nd to pair with #1
    //while(pTrack2 < pTracks->end())
    for (; pTrack2!=pTracks->end(); pTrack2++, trk2Idx++) 
    {
        TkrTrack* track2 = *pTrack2;
        
        /// --* Same here, see above comment
        if(gamEne == 0.) gamEne = e_t1+track2->getInitialEnergy();

        Point   trk1Pos = track1->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
        Vector  trk1Dir = track1->front()->getDirection(Event::TkrTrackHit::SMOOTHED);
        Point   trk2Pos = track2->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
        Vector  trk2Dir = track2->front()->getDirection(Event::TkrTrackHit::SMOOTHED);

        RayDoca doca    = RayDoca(Ray(trk1Pos, trk1Dir), Ray(trk2Pos, trk2Dir));
        double dist = doca.docaRay1Ray2();
//        double s1   = doca.arcLenRay1();
//        double s2   = doca.arcLenRay2(); 
        
        double cost1t2 = trk1Dir * trk2Dir;
        double t1t2    = acos(cost1t2);
        
        double head_sep = (trk1Pos - trk2Pos).magnitude();
        if(head_sep > 40.) continue; 

        double q2       = track2->getQuality();   
 //       double wgt_doca = thrshold((dist-100./gamEne)/root2);  
 //       double wgt_s1   = thrshold((s1+.5)/(root2*1.5 * gamEne/100.));
  //      VTX_DOCA_Wgt = thrshold((dist-50./gamEne)/root2);  
  //      VTX_S1_Wgt   = thrshold((s1-2.5)/(root2*1.5 * gamEne/100.));     
  //      VTX_T12_Wgt  = thrshold((.01*t1t2*gamEne -.2)/(root2*.2));
  //      VTX_T2Q_Wgt  = thrshold((40. - q2)/(root2*20.)); 
  //    double wgt_t1t2 = thrshold((t1t2-200./gamEne)/(root2*200./gamEne));

        double wgt_t1t2 = thrshold((.01*t1t2*gamEne -.2)/(root2*.2));
        double wgt_q2   = thrshold((40. - q2)/(root2*20.));
        double wgt_hs   = thrshold((head_sep - 2.)/(root2*2.));
        
        double total_wgt = wgt_hs*wgt_q2*wgt_t1t2; 
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
    if (best_track2 != 0 && max_wgt > 0.) {
        

        Vector gamDir;
 
        // Begin Method 1
        //Use tracking covariances & energies to weight track directions
        const TkrTrackParams& par_t1 = track1->front()->getTrackParams(Event::TkrTrackHit::SMOOTHED);
        const TkrTrackParams& par_t2 = best_track2->front()->getTrackParams(Event::TkrTrackHit::SMOOTHED);

        double sx_1 = par_t1.getxSlope();
        double sy_1 = par_t1.getySlope();

        double sx_2 = par_t2.getxSlope();
        double sy_2 = par_t2.getySlope();

        //double e_t2 = best_track2->getI

        /// --* see above question on best energy to use here
        double e_t2 =  best_track2->getInitialEnergy();
        
        double norm_1 = sqrt(1.+ sx_1*sx_1 + sy_1*sy_1); 
        
        double norm_2 = sqrt(1.+ sx_2*sx_2 + sy_2*sy_2);  
        
        double tx_1 = -sx_1/norm_1;
        double ty_1 = -sy_1/norm_1;
        double tx_2 = -sx_2/norm_2;
        double ty_2 = -sy_2/norm_2;
        double tz_1 = -1./norm_1;
        double tz_2 = -1./norm_2;
        
        double dtx_1_sq = ((1.+sy_1*sy_1)*(1.+sy_1*sy_1)*par_t1.getxSlpxSlp() +
                           sx_1*sx_1*sy_1*sy_1*par_t1.getySlpySlp())/pow(norm_1,6);
        double dty_1_sq = ((1.+sx_1*sx_1)*(1.+sx_1*sx_1)*par_t1.getySlpySlp() +
                           sx_1*sx_1*sy_1*sy_1*par_t1.getxSlpxSlp())/pow(norm_1,6);
        double dtz_1_sq = (tx_1*tx_1/(tz_1*tz_1))*dtx_1_sq + 
                          (ty_1*ty_1/(tz_1*tz_1))*dty_1_sq;
        double dtx_2_sq = ((1.+sy_2*sy_2)*(1.+sy_2*sy_2)*par_t2.getxSlpxSlp() +
                           sx_2*sx_2*sy_2*sy_2*par_t2.getySlpySlp())/pow(norm_2,6);
        double dty_2_sq = ((1.+sx_2*sx_2)*(1.+sx_2*sx_2)*par_t2.getySlpySlp() +
                           sx_2*sx_2*sy_2*sy_2*par_t2.getxSlpxSlp())/pow(norm_2,6);
        double dtz_2_sq = (tx_2*tx_2/(tz_2*tz_2))*dtx_2_sq + 
                          (ty_2*ty_2/(tz_2*tz_2))*dty_2_sq;
        

        double transistion = 1.; //thrshold((gamEne - 400.)/150.);
        
        double wt_x1 = 1./(dtx_1_sq*transistion + (1. - transistion));
        double wt_x2 = 1./(dtx_2_sq*transistion + (1. - transistion));       
        double wt_1x  = (e_t1*wt_x1) / (e_t1*wt_x1 + e_t2*wt_x2);
        double wt_2x  = (e_t2*wt_x2) / (e_t1*wt_x1 + e_t2*wt_x2);
        double wt_y1 = 1./(dty_1_sq*transistion + (1. - transistion));
        double wt_y2 = 1./(dty_2_sq*transistion + (1. - transistion));  
        double wt_1y  = (e_t1*wt_y1) / (e_t1*wt_y1 + e_t2*wt_y2);
        double wt_2y  = (e_t2*wt_y2) / (e_t1*wt_y1 + e_t2*wt_y2);
        double wt_z1 = 1./(dtz_1_sq*transistion + (1. - transistion));
        double wt_z2 = 1./(dtz_2_sq*transistion + (1. - transistion));  
        double wt_1z  = (e_t1*wt_z1) / (e_t1*wt_z1 + e_t2*wt_z2);
        double wt_2z  = (e_t2*wt_z2) / (e_t1*wt_z1 + e_t2*wt_z2);
        
        double gam_tx = tx_1*wt_1x + tx_2*wt_2x;
        double gam_ty = ty_1*wt_1y + ty_2*wt_2y;     
        double gam_tz = tz_1*wt_1z + tz_2*wt_2z;
        
        gamDir = Vector(gam_tx, gam_ty, gam_tz).unit();
        
        Ray        gamma  = Ray(gamPos,gamDir);
        int        layer  = 2 * track1->front()->getTkrId().getTray() - 1 + track1->front()->getTkrId().getBotTop();
        int        tower  = 4 * track1->front()->getTkrId().getTowerX() + track1->front()->getTkrId().getTowerY(); 
        TkrVertex* vertex = new TkrVertex(layer,tower,gamEne,best_doca,gamma);
        //                vertex->setDist1(doca.arcLenRay1());
        //                vertex->setDist2(doca.arcLenRay2())
        //                vertex->setAngle(t1t2); 
        vertex->addTrack(track1);
        vertex->addTrack(best_track2);
        
        vertexCol->push_back(vertex); // addVertex(vertex);

        // Add track/vertex to relational table
        Event::TkrVertexTrackRel* rel1 = new Event::TkrVertexTrackRel(vertex, track1);
        Event::TkrVertexTrackRel* rel2 = new Event::TkrVertexTrackRel(vertex, best_track2);
        vertexRelTab->addRelation(rel1);
        vertexRelTab->addRelation(rel2);

        unused[trk1Idx] = false;
        unused[bst_trk2Idx] = false;
    }
            
            
     //Go through unused list looking for isolated tracks
    TkrTrackColPtr pTrack = pTracks->begin();
    int            trkIdx = 0;
            
    //while(pTrack != pTracks->end())
    for (; pTrack != pTracks->end(); pTrack++, trkIdx++)
    {
         TkrTrack* track1 = *pTrack;
                
         if (unused[trkIdx])
         {
            Point   trk1Pos = track1->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
            Vector  trk1Dir = track1->front()->getDirection(Event::TkrTrackHit::SMOOTHED);
            int     layer   = 35 - 2 * track1->front()->getTkrId().getTray() - 1 + track1->front()->getTkrId().getBotTop();
            int     tower   = 4 * track1->front()->getTkrId().getTowerX() + track1->front()->getTkrId().getTowerY(); 

            TkrVertex* vertex = new TkrVertex(layer,tower,track1->getInitialEnergy(),0.,Ray(trk1Pos,trk1Dir));
                    
            vertex->addTrack(track1);
                    
            vertexCol->push_back(vertex); //addVertex(vertex);
 
            // Add track/vertex to relational table
            Event::TkrVertexTrackRel* rel = new Event::TkrVertexTrackRel(vertex, track1);
            vertexRelTab->addRelation(rel);
        }
    }            
    return;
}


TkrComboVtxRecon::~TkrComboVtxRecon()
{   
    return;
}

