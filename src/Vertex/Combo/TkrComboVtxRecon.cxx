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
	double docaLimit = 10.0;  //mm

    while(numTracks--) unused[numTracks] = true;

    //Track counter
    int    trk1Idx = 0;

    //Loop over the number of Fit tracks
    TkrFitConPtr pTrack1 = pTracks->begin();

    while(pTrack1 != pTracks->end() && unused[trk1Idx])
    {
        TkrFitTrack* track1  = *pTrack1++;
        TkrFitConPtr pTrack2 = pTrack1;
        int          trk2Idx = trk1Idx + 1;

        while(pTrack2 != pTracks->end() && unused[trk2Idx])
        {
            TkrFitTrack* track2 = *pTrack2++;

            RayDoca doca    = RayDoca(Ray(track1->getPosition(),track1->getDirection()),
                                      Ray(track2->getPosition(),track2->getDirection()));
            double  dist    = doca.docaRay1Ray2();

            //Arc length to doca from track start
            double  arcLen1 = doca.arcLenRay1();
            double  arcLen2 = doca.arcLenRay2();

            //Length of the track
            double  trkLen1 = (track1->getPosition(TkrRecInfo::Start) - track1->getPosition(TkrRecInfo::End)).mag();
            double  trkLen2 = (track2->getPosition(TkrRecInfo::Start) - track2->getPosition(TkrRecInfo::End)).mag();

            //Keep the candidate vertex if 1) doca is not too big
            bool    keepIt  = dist < docaLimit;

            // 2) the doca occurs along track1 or not too far past it
            keepIt = keepIt && (-20. < arcLen1 && arcLen1 < trkLen1);

            // 3) ditto for track 2
            keepIt = keepIt && (-20. < arcLen2 && arcLen2 < trkLen2);

            //Keep it?
            if (keepIt)
            {
                Point  gamPos;
                Vector gamDir;
                double gamEne = track1->getEnergy() + track2->getEnergy();

                //Ok, check if first track vertex's with second before first hit
                if      (arcLen1 > 0. && arcLen2 < 1.)
                {
                    //Ok, make sure end of track 2 close to track 1
                    if (doca.docaRay1Point2() > docaLimit) continue; 

                    gamPos = track1->getPosition();
                    gamDir = track1->getDirection();
                }
                //Same check for track 2
                else if (arcLen2 > 0. && arcLen1 < 1.)
                {
                    //Make sure end of track 1 is close to track 2
                    if (doca.docaPoint1Ray2() > docaLimit) continue; 

                    gamPos = track2->getPosition();
                    gamDir = track2->getDirection();
                }
                //Else both tracks make good vertex
                else
                {
                    Vector trk1Dir = track1->getDirection();
                    Vector trk2Dir = track2->getDirection();

                    trk1Dir.setMag(track1->getEnergy());
                    trk2Dir.setMag(track2->getEnergy());

                    gamDir  = trk1Dir + trk2Dir;

                    gamDir.setMag(1.);

                    gamPos  = doca.docaPointRay1();
                    gamPos += doca.docaPointRay2();
                    gamPos *= 0.5;
                }

                Ray        gamma  = Ray(gamPos,gamDir);
                TkrVertex* vertex = new TkrVertex(track1->getLayer(),track1->getTower(),gamEne,dist,gamma);

                vertex->addTrack(track1);
                vertex->addTrack(track2);

                push_back(vertex); // addVertex(vertex);

                unused[trk1Idx] = false;
                unused[trk2Idx] = false;

                trk2Idx++;
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

