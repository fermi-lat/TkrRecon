#include "TkrRecon/Display/TkrClustersRep.h"
#include "TkrRecon/ITkrGeometrySvc.h"


//------------------------------------------------------------------------------
//  Code to implement the TkrClustersRep class. 
//  This class does the display of clusters in the Silicon
//  Tracker. 
//
//#############################################################################
//  Constructor for the class
TkrClustersRep::TkrClustersRep(IDataProviderSvc* dataProviderSvc)
//#############################################################################
{
    dps = dataProviderSvc;
}

//##############################################
//  Display routine for the class
void TkrClustersRep::update()
//##############################################
{
    TkrRecon::TkrClusterCol* pClusters = SmartDataPtr<TkrRecon::TkrClusterCol>(dps,"/Event/TkrRecon/TkrClusterCol");

    if (pClusters)
    {
		int    nHits      = pClusters->nHits();
		double stripPitch = pTkrGeo->siStripPitch();
		double towerPitch = pTkrGeo->towerPitch();
		
		setColor("green");

		//Loop over all cluster hits in the SiClusters vector
		while(nHits--)
		{
			TkrRecon::TkrCluster* pCluster = pClusters->getHit(nHits);
			Point                 clusPos  = pCluster->position();
	 
			double      x        = clusPos.x();
			double      y        = clusPos.y();
			double      z        = clusPos.z();
	
			double      delta    =  20.*stripPitch;
			double      Offset   = -0.5 * towerPitch;

			//Draw a cross at the position of the cluster center
			if (pCluster->v() == TkrRecon::TkrCluster::X)
			{
                Offset += y;
                moveTo(Point(x - delta, Offset, z - delta));
				lineTo(Point(x + delta, Offset, z + delta));
				moveTo(Point(x - delta, Offset, z + delta));
				lineTo(Point(x + delta, Offset, z - delta));
			}
			else
			{
                Offset += x;
				moveTo(Point(Offset, y - delta, z - delta));
				lineTo(Point(Offset, y + delta, z + delta));
				moveTo(Point(Offset, y - delta, z + delta));
				lineTo(Point(Offset, y + delta, z - delta));
			}

			//Now draw the width of the cluster (if more than one)
			if (pCluster->size() > 1)
			{
				double stripWid = 0.5 * pCluster->size() * stripPitch;

				if (pCluster->v() == TkrRecon::TkrCluster::X)
				{
					moveTo(Point(x - stripWid, Offset, z));
					lineTo(Point(x + stripWid, Offset, z));
				}
				else
				{
					moveTo(Point(Offset, y - stripWid, z));
					lineTo(Point(Offset, y + stripWid, z));
				}
			}
		}
	}

    return;
}

ITkrGeometrySvc * pTkrGeo;

