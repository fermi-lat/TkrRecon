#include "TkrRecon/Display/TkrBestCandRep.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

//#############################################################################
TkrBestCandRep::TkrBestCandRep(IDataProviderSvc* dataProviderSvc, ITkrGeometrySvc* pTkrGeometry)
//#############################################################################
{
    dps     = dataProviderSvc;
    pTkrGeo = pTkrGeometry;
}
//-------------------- private ----------------------
//##############################################
void TkrBestCandRep::update()
//##############################################
{
    TkrPatCandCol* pTkrCandidates = SmartDataPtr<TkrPatCandCol>(dps,"/Event/TkrRecon/TkrPatCandCol");

    if (pTkrCandidates)
    {
        TkrLinkAndTree* pTkrCands = dynamic_cast<TkrLinkAndTree*>(pTkrCandidates);

	    //Now see if we can do the drawing
	    if (pTkrCands)
	    {
		    gui::DisplayRep* pDisplay = this;
	        if (pTkrCands->getNumTrees(X) > 0)
            {
		        setColor("green");
		        TkrDrawBestCand(pTkrCands, X);
            }

	        //Draw the links for the "best" track
	        if (pTkrCands->getNumTrees(Y) > 0)
            {
		        setColor("aquamarine");
		        TkrDrawBestCand(pTkrCands, Y);
            }

	        setColor("blue");
        }
    }

    return;
}


const char bstcol_blue[]       = "blue";
const char bstcol_violet[]     = "violet";
const char bstcol_turquoise[]  = "turquoise";
const char bstcol_orange[]     = "orange";
const char bstcol_maroon[]     = "maroon";
const char bstcol_aquamarine[] = "aquamarine";

const char* pBstColors[] = {bstcol_blue,   bstcol_violet, bstcol_turquoise,
                            bstcol_orange, bstcol_maroon, bstcol_aquamarine};

void TkrBestCandRep::TkrDrawBestCand(TkrPatCandCol* pTkrCandidates, TkrPlaneType plane)
{
    TkrLinkAndTree* pTkrCands      = dynamic_cast<TkrLinkAndTree*>(pTkrCandidates);

    //Draw the candidate tracks
    TkrLinkForest* pForest  = pTkrCands->getForest(plane);
    int            nTrees   = pForest->getNumTrees();
	treeListPtr    treePtr  = pForest->getListStart();
	bool           fullTree = true;
    int            colorIdx = 0;

	//messageManager::instance()->message(" ****************************************");
	//messageManager::instance()->message(" Number of trees found:", numTrees);

	//Only plot six longest, straightest trees
	if (nTrees > 6) nTrees = 6;

	while(nTrees--)
    {
        TkrLinkTree* pTree = &(*treePtr++);

		//messageManager::instance()->message(" Tree Height:", pTree->getTreeHeight());

	    //Now come back and do the "best" set of links in the tree
        int nBest = pTree->getNumBestVectors();
        bestNodeVecPtr nodeVecPtr = pTree->getNodeVectorPtr();

        while(nBest--)
        {

            BestNodeList*   pNodeList = *nodeVecPtr++;
     	    int             nNodes    = pNodeList->getNumNodes();
	        linkNodeListPtr nodePtr   = pNodeList->getNodeList();

	        setColor(pBstColors[colorIdx]);

	        while(nNodes--)
            {
                LayerLinkNode* pNode   = *nodePtr++;
		        TkrLinkNode*   pTkrNode = dynamic_cast<TkrLinkNode*>(pNode);

		        drawLinkNode(pTkrNode);
            }
        }

	    //Put a marker at the tree head
	    TkrLinkNode*    pTkrNode  = dynamic_cast<TkrLinkNode*>(pTree->getHeadNode());
	    LayerLink*      pLink     = pTkrNode->getThisLayerLink();
	    TkrClusterLink* pTkrLink  = dynamic_cast<TkrClusterLink*>(pLink);
	    TkrCluster*     pClus     = pTkrLink->pTopClus();
    	double          x         = pClus->position().x();
	    double          y         = pClus->position().y();
	    double          z         = pClus->position().z();
	    //double          offset    = -0.5*trackerGeo::trayWidth();

	    if (pClus->v() == TkrCluster::X) y += 0.5 * pTkrGeo->towerPitch();
	    else                             x += 0.5 * pTkrGeo->towerPitch();

	    setColor(pBstColors[colorIdx]);
	    markerAt(Point(x,y,z));

		//fullTree = false;

		colorIdx = ++colorIdx % 6;
    }

    return;
}

void TkrBestCandRep::drawLinkNode(TkrLinkNode* pTkrNode)
{
	TkrClusterLink* pClusLink = dynamic_cast<TkrClusterLink*>(pTkrNode->getThisLayerLink());

    TkrCluster* pTopCluster = pClusLink->pTopClus();

	double x      = pTopCluster->position().x();
	double y      = pTopCluster->position().y();
	double z      = pTopCluster->position().z();
	//double offset = -0.5*trackerGeo::trayWidth();

	if (pTopCluster->v() == TkrCluster::X) y += 0.5 * pTkrGeo->towerPitch();
	else                                         x += 0.5 * pTkrGeo->towerPitch();

	moveTo(Point(x,y,z));

    TkrCluster* pBotCluster = pClusLink->pBotClus();

	x = pBotCluster->position().x();
	y = pBotCluster->position().y();
	z = pBotCluster->position().z();

	if (pBotCluster->v() == TkrCluster::X) y += 0.5 * pTkrGeo->towerPitch();
	else                                         x += 0.5 * pTkrGeo->towerPitch();

	lineTo(Point(x,y,z));


	return;
}
