/*
	Code to implement the Combo vertex finding class

	Tracy Usher March 1, 2002
*/

#include "src/Vertex/Combo/TkrComboVtxRecon.h"
#include "src/Vertex/Combo/RayDoca.h"

TkrComboVtxRecon::TkrComboVtxRecon(ITkrGeometrySvc* pTkrGeo, TkrTracks* pTracks, TkrCandidates* pCandTracks)
{
    //Define a vector to contain a list of "isolated" tracks
    int   numTracks = pTracks->getNumTracks();
    bool* unused    = new bool[numTracks];
	double docaLimit = 10.0;  //mm

    while(numTracks--) unused[numTracks] = true;

    //Track counter
    int   trk1Idx = 0;

    //Loop over the number of Fit tracks
    TkrVectorPtr pTrack1 = pTracks->getTrackPtr();

    while(pTrack1 < pTracks->getTrackEnd())
    {
        TkrFitTrack* track1  = *pTrack1++;
        TkrVectorPtr pTrack2 = pTrack1;
        int          trk2Idx = trk1Idx;

        while(pTrack2 < pTracks->getTrackEnd())
        {
            TkrFitTrack* track2 = *pTrack2++;

            RayDoca doca = RayDoca(Ray(track1->position(),track1->direction()),
                                   Ray(track2->position(),track2->direction()));

            trk2Idx++;

            //Check that doca not too big and that vertex starts before or at first hit
            if (doca.docaRay1Ray2() < docaLimit && (doca.arcLenRay1() <= 0. || doca.arcLenRay2() <= 0.))
            {
                Point  gamPos;
                Vector gamDir;
                double gamEne = track1->energy() + track2->energy();

                //Ok, check if first track vertex's with second before first hit
                if      (doca.arcLenRay1() > 0.)
                {
                    if (doca.docaRay1Point2() > docaLimit) continue; 

                    gamPos = track1->position();
                    gamDir = track1->direction();
                }
                //Same check for track 2
                else if (doca.arcLenRay2() > 0.)
                {
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
                TkrVertex* vertex = new TkrVertex(track1->firstLayer(),track1->tower(),gamEne,gamma);

                vertex->addTrack(track1);
                vertex->addTrack(track2);

                addVertex(vertex);

                unused[trk1Idx] = false;
                unused[trk2Idx] = false;
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

            TkrVertex* vertex = new TkrVertex(track1->firstLayer(),track1->tower(),track1->energy(),Ray(track1->position(),track1->direction()));

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

