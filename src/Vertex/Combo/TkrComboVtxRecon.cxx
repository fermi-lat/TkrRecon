/*
	Code to implement the Combo vertex finding class

	Tracy Usher March 1, 2002
*/

#include "src/Vertex/Combo/TkrComboVtxRecon.h"
#include "src/Vertex/Combo/RayDoca.h"

TkrComboVtxRecon::TkrComboVtxRecon(ITkrGeometrySvc* pTkrGeo, TkrFitTrackCol* pTracks, TkrPatCandCol* pCandTracks)
{
    //Define a vector to contain a list of "isolated" tracks
    int   numTracks = pTracks->getNumTracks();
    bool* unused    = new bool[numTracks];
	double docaLimit = 10.0;  //mm

    while(numTracks--) unused[numTracks] = true;

    //Track counter
    int   trk1Idx = 0;

    //Loop over the number of Fit tracks
    TkrFitColPtr pTrack1 = pTracks->getTrackPtr();

    while(pTrack1 < pTracks->getTrackEnd() && unused[trk1Idx])
    {
        TkrFitTrack* track1  = *pTrack1++;
        TkrFitColPtr pTrack2 = pTrack1;
        int          trk2Idx = trk1Idx + 1;

        while(pTrack2 < pTracks->getTrackEnd() && unused[trk2Idx])
        {
            TkrFitTrack* track2 = *pTrack2++;

            RayDoca doca    = RayDoca(Ray(track1->position(),track1->direction()),
                                      Ray(track2->position(),track2->direction()));
            double  dist    = doca.docaRay1Ray2();

            //Arc length to doca from track start
            double  arcLen1 = doca.arcLenRay1();
            double  arcLen2 = doca.arcLenRay2();

            //Length of the track
            double  trkLen1 = (track1->position(TkrRecInfo::Start) - track1->position(TkrRecInfo::End)).mag();
            double  trkLen2 = (track2->position(TkrRecInfo::Start) - track2->position(TkrRecInfo::End)).mag();

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
                double gamEne = track1->energy() + track2->energy();

                //Ok, check if first track vertex's with second before first hit
                if      (arcLen1 > 0. && arcLen2 < 1.)
                {
                    //Ok, make sure end of track 2 close to track 1
                    if (doca.docaRay1Point2() > docaLimit) continue; 

                    gamPos = track1->position();
                    gamDir = track1->direction();
                }
                //Same check for track 2
                else if (arcLen2 > 0. && arcLen1 < 1.)
                {
                    //Make sure end of track 1 is close to track 2
                    if (doca.docaPoint1Ray2() > docaLimit) continue; 

                    gamPos = track2->position();
                    gamDir = track2->direction();
                }
                //Else both tracks make good vertex
                else
                {
                    Vector trk1Dir = track1->direction();
                    Vector trk2Dir = track2->direction();

                    trk1Dir.setMag(track1->energy());
                    trk2Dir.setMag(track2->energy());

                    gamDir  = trk1Dir + trk2Dir;

                    gamDir.setMag(1.);

                    gamPos  = doca.docaPointRay1();
                    gamPos += doca.docaPointRay2();
                    gamPos *= 0.5;
                }

                Ray        gamma  = Ray(gamPos,gamDir);
                TkrVertex* vertex = new TkrVertex(track1->layer(),track1->tower(),gamEne,dist,gamma);

                vertex->addTrack(track1);
                vertex->addTrack(track2);

                addVertex(vertex);

                unused[trk1Idx] = false;
                unused[trk2Idx] = false;

                trk2Idx++;
            }
        }

        trk1Idx++;
    }


    //Go through unused list looking for isolated tracks
    numTracks = pTracks->getNumTracks();

    while(numTracks--)
    {
        if (unused[numTracks])
        {
            TkrFitTrack* track1 = pTracks->getTrack(numTracks);

            TkrVertex* vertex = new TkrVertex(track1->layer(),track1->tower(),track1->energy(),0.,Ray(track1->position(),track1->direction()));

            vertex->addTrack(track1);

            addVertex(vertex);
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

