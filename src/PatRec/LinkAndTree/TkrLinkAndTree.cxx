/*
	Code to implement the PatRecTracks class
	Tracy Usher Nov 27, 2000
*/

#include "src/PatRec/LinkAndTree/TkrLinkAndTree.h"

TkrLinkAndTree::TkrLinkAndTree(ITkrGeometrySvc* pTkrGeo, TkrClusters* pTkrClus)
{
	ini();

	//How many clusters are we dealing with?
	setNumClusters(pTkrClus->nHits());

	//Build the links in the X view
	setLinkList(new TkrClusterLinkList(pTkrClus, X), X);

	//Build the trees from this
	setForest(new TkrLinkForest(getLinkList(X)), X);
	
	//Build the links in the Y view
	setLinkList(new TkrClusterLinkList(pTkrClus, Y), Y);

	//Build the trees from this
	setForest(new TkrLinkForest(getLinkList(Y)), Y);

    //Now build the 3-D track candidates from this
    buildCand3D();


	return;
}


TkrLinkAndTree::~TkrLinkAndTree()
{
	clear();

	return;
}

void TkrLinkAndTree::setLinkList(TkrClusterLinkList* pLinks, TkrPlaneType plane)
{
	if (plane == X)
	{
		pLinkListX = pLinks;
		numLinksX  = pLinks->getNumLinksTotal();
	}
	else
	{
		pLinkListY = pLinks;
		numLinksY  = pLinks->getNumLinksTotal();
	}

	return;
}

TkrClusterLinkList* TkrLinkAndTree::getLinkList(TkrPlaneType plane)
{
	return plane == X ? pLinkListX : pLinkListY;
}

int TkrLinkAndTree::getNumLinks(TkrPlaneType plane)
{
	return plane == X ? numLinksX : numLinksY;
}

void TkrLinkAndTree::setForest(TkrLinkForest* pForest, TkrPlaneType plane)
{
	if (plane == X)
	{
		pForestX   = pForest;
		numTracksX = pForest->getNumTrees();
	}
	else
	{
		pForestY   = pForest;
		numTracksY = pForest->getNumTrees();
	}

	return;
}

TkrLinkForest* TkrLinkAndTree::getForest(TkrPlaneType plane)
{
	return plane == X ? pForestX : pForestY;
}

void TkrLinkAndTree::buildCand3D()
{
    int numXcands = pForestX->getNumTrees();
    int numYcands = pForestY->getNumTrees();

    //This routine blindly assumes that, since the trees have been presorted
    //by longest and straightest, that the two projections already match up.
    //We will proceed (for now!) by just tossing the info from the two into 
    //the pattern recognition output object
    if (numXcands > 0 && numYcands > 0)
    {
        int           num3Dtrks = numXcands < numYcands ? numXcands : numYcands;
        treeListPtr   xIter     = pForestX->getListStart();
        TkrLinkTree*  pTreeX    = &(*xIter);
        BestNodeList* nodeListX = pTreeX->getBestNodeList(0);
        treeListPtr   yIter     = pForestY->getListStart();
        TkrLinkTree*  pTreeY    = &(*yIter);
        BestNodeList* nodeListY = pTreeY->getBestNodeList(0);
        int           nTwrTries = 0;
        int           maxLength = 0;

        //Ok, make a pass through to find the maximum length...
        while(xIter != pForestX->getListEnd() && yIter != pForestY->getListEnd())
        {
            pTreeX    = &(*xIter++);
            pTreeY    = &(*yIter++);
            nodeListX = pTreeX->getBestNodeList(0);
            nodeListY = pTreeY->getBestNodeList(0);

            int             nNodesX      = nodeListX->getNumNodes();
            int             nNodesY      = nodeListY->getNumNodes();

            if (maxLength < nNodesX) maxLength = nNodesX;
            if (maxLength < nNodesY) maxLength = nNodesY;

            if (maxLength >= 17-pTreeX->getFirstLayer()) break;
            if (maxLength >= 17-pTreeY->getFirstLayer()) break;
        }

        int minLength = maxLength > 11 ? maxLength / 3 : 3;

        xIter = pForestX->getListStart();
        yIter = pForestY->getListStart();

        //Loop over possible number of 3D tracks
        //boy is this ugly...
        while(xIter != pForestX->getListEnd() && yIter != pForestY->getListEnd())
        {
            pTreeX    = &(*xIter);
            pTreeY    = &(*yIter);
            nodeListX = pTreeX->getBestNodeList(0);
            nodeListY = pTreeY->getBestNodeList(0);

            int             nNodesX      = nodeListX->getNumNodes();
            int             nNodesY      = nodeListY->getNumNodes();

            linkNodeListPtr nodePtrX     = nodeListX->getNodeList();
            linkNodeListPtr nodePtrY     = nodeListY->getNodeList();

            TkrLinkNode*    TkrNodeX     = dynamic_cast<TkrLinkNode*>(*nodePtrX);
            TkrLinkNode*    TkrNodeY     = dynamic_cast<TkrLinkNode*>(*nodePtrY);

	        TkrClusterLink* pClusLinkX   = dynamic_cast<TkrClusterLink*>(TkrNodeX->getThisLayerLink());
            TkrCluster*     pTopClusterX = pClusLinkX->pTopClus();
            TkrCluster*     pBotClusterX = pClusLinkX->pBotClus();

	        TkrClusterLink* pClusLinkY   = dynamic_cast<TkrClusterLink*>(TkrNodeY->getThisLayerLink());
            TkrCluster*     pTopClusterY = pClusLinkY->pTopClus();
            TkrCluster*     pBotClusterY = pClusLinkY->pBotClus();

            //Starting point must be the same tower
            if (pTopClusterX->tower() != pTopClusterY->tower())
            {
                //Ok, try going to next x combo, if that doesn't work try the y...
                if (nTwrTries++ < 2)
                {
                    if (nTwrTries == 1) {xIter++;}
                    else                {xIter--;yIter++;}
                }
                else
                {
                    xIter++;
                    nTwrTries = 0;
                }

                continue;
            }

            //ok, this is to check how long the candidate track is
            int longest = nNodesX > nNodesY ? nNodesX : nNodesY;

            //Must be long enough
            if (longest < minLength)
            {
                xIter++;
                yIter++;
                continue;
            }

            // Both should start within a layer of each other
            int strtDiff = pTopClusterX->plane()-pTopClusterY->plane();
            int nodeDiff = nNodesX - nNodesY;

            if (abs(strtDiff) > 1)
            {
                //Ok, bump the track which starts highest up...
                if (strtDiff > 0) yIter++;
                else              xIter++;
                continue;
            }

            //Check the lengths... same start layer
            //if (strtDiff == 0 && abs(nodeDiff) > 0 && (nNodesX < 4 || nNodesY < 4))
            if (strtDiff == 0 && abs(nodeDiff) > 0 && (nNodesX < 3 || nNodesY < 3))
            {
                if (nodeDiff > 0) xIter++;
                else              yIter++;
                continue;
            }

            //Check the lengths... different start layer
            //if (strtDiff != 0 && (nNodesX < 4 || nNodesY < 4 || abs(nodeDiff) > 2))
            if (strtDiff != 0 && (nNodesX < 3 || nNodesY < 3 || abs(nodeDiff) > 3))
            {
                if (strtDiff > 0) yIter++;
                else              xIter++;
                continue;
            }

            //Ok, its a keeper
	        double          x            = pTopClusterX->position().x();
	        double          y            = pTopClusterY->position().y();
            double          z            = pTopClusterX->position().z() > pTopClusterY->position().z()
                                         ? pTopClusterX->position().z()
                                         : pTopClusterY->position().z();

            Point  start(x,y,z);
            Vector vDir(pBotClusterX->position().x()-x,
                        pBotClusterY->position().y()-y,
                        pBotClusterX->position().z()-z);
            Ray    initDir(start,vDir.unit());
            double quality = TkrNodeX->getQuality()->getAngleRMS()
                           + TkrNodeY->getQuality()->getAngleRMS();

            TkrPatCand* newTrack = new TkrPatCand(pTopClusterX->plane(),pTopClusterX->tower(),0.,quality,initDir);

            newTrack->addCandHit(pTopClusterX);
            newTrack->addCandHit(pTopClusterY);
           
            while(nNodesX--)
            {
                TkrNodeX     = dynamic_cast<TkrLinkNode*>(*nodePtrX++);
	            pClusLinkX   = dynamic_cast<TkrClusterLink*>(TkrNodeX->getThisLayerLink());
                
                newTrack->addCandHit(pClusLinkX->pBotClus());
            }
           
            while(nNodesY--)
            {
                TkrNodeY     = dynamic_cast<TkrLinkNode*>(*nodePtrY++);
	            pClusLinkY   = dynamic_cast<TkrClusterLink*>(TkrNodeY->getThisLayerLink());
                
                newTrack->addCandHit(pClusLinkY->pBotClus());
            }

            addTrack(newTrack);

            xIter++;
            yIter++;
        }
    }

    return;
}

int TkrLinkAndTree::getNumTrees(TkrPlaneType plane)
{
	return plane == X ? numTracksX : numTracksY;
}

void TkrLinkAndTree::ini()
{
	numClusters = 0;
	pLinkListX  = 0;
	numLinksX   = 0;
	pLinkListY  = 0;
	numLinksY   = 0;
	pForestX    = 0;
	numTracksX  = 0;
	pForestY    = 0;
	numTracksY  = 0;

	return;
}

void TkrLinkAndTree::clear()
{
	//Make sure we clean things up before departing
	if (pForestX)   delete pForestX;
	if (pForestY)   delete pForestY;
	if (pLinkListX) delete pLinkListX;
	if (pLinkListY) delete pLinkListY;

	numLinksX  = 0;
	numLinksY  = 0;
	numTracksX = 0;
	numTracksY = 0;

	return;
}

/*
//########################################################
void TkrLinkAndTree::writeOut() const
//########################################################
{
	if (numTracksX <=0 && numTracksY <=0) return;
	
	std::ostream& out = messageManager::instance()->out();
	if (!messageManager::instance()->acceptLevel("DEBUG")) return;

	out << " --- TkrLinkAndTree ---- " << "\n";
	out << " num Tracks X = " << numTracksX << "\n";
	out << " num Tracks Y = " << numTracksY << "\n";
//	if (numGammas()>0) {
//		for (int ig=0; ig < numGammas(); ig++) {
//			m_GFgammaList[ig]->GFdata::writeOut(out);
//		}
//	}
//	out << " num Particles = " << numParticles() << "\n";
//	if (numParticles()>0) {
//		for (int ip=0; ip < numParticles(); ip++) { 
//			m_GFparticleList[ip]->GFdata::writeOut(out);
//		}
//	}
	return;
}
*/
/*
//########################################################
void TkrLinkAndTree::update(GraphicsRep& v)
//########################################################
{
	//Draw the links for the "best" track
	if (numTracksX > 0)
	{
		v.setColor("green");
		pForestX->draw(v);
	}

	//Draw the links for the "best" track
	if (numTracksY > 0)
	{
		v.setColor("aquamarine");
		pForestY->draw(v);
	}

	v.setColor("blue");

	return;
}
*/