#include "Event/TopLevel/EventModel.h"
#include "TkrCandidatesRep.h"

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

//#############################################################################
TkrCandidatesRep::TkrCandidatesRep(IDataProviderSvc* dataProviderSvc, ITkrGeometrySvc* pTkrGeometry)
//#############################################################################
{
    dps     = dataProviderSvc;
    pTkrGeo = pTkrGeometry;
}
//-------------------- private ----------------------
//##############################################
void TkrCandidatesRep::update()
//##############################################
{
    TkrPatCandCol* pTkrCandidates = SmartDataPtr<TkrPatCandCol>(dps,EventModel::TkrRecon::TkrPatCandCol);

    //Now see if we can do the drawing
    if (pTkrCandidates)
    {
        TkrLinkAndTree* pTkrCands = dynamic_cast<TkrLinkAndTree*>(pTkrCandidates);
    //  gui::DisplayRep*  pDisplay  = this;

        if (pTkrCands->getNumTrees(X) > 0)
        {
            setColor("green");
            TkrDrawCandidates(pTkrCands, X);
        }

        //Draw the links for the "best" track
        if (pTkrCands->getNumTrees(Y) > 0)
        {
            setColor("aquamarine");
            TkrDrawCandidates(pTkrCands, Y);
        }

        setColor("blue");
    }

    return;
}


const char color_blue[]       = "blue";
const char color_violet[]     = "violet";
const char color_turquoise[]  = "turquoise";
const char color_orange[]     = "orange";
const char color_maroon[]     = "maroon";
const char color_aquamarine[] = "aquamarine";

const char* pColors[] = {color_blue,   color_violet, color_turquoise,
                         color_orange, color_maroon, color_aquamarine};

void TkrCandidatesRep::TkrDrawCandidates(TkrPatCandCol* pTkrCandidates, TkrPlaneType plane)
{
    TkrLinkAndTree* pTkrCands      = dynamic_cast<TkrLinkAndTree*>(pTkrCandidates);

    //Draw the candidate tracks
    TkrLinkForest* pForest  = pTkrCands->getForest(plane);
    int            nTrees   = pForest->getNumTrees();
    treeListPtr    treePtr  = pForest->getListStart();
    int            colorIdx = 0;

    //messageManager::instance()->message(" ****************************************");
    //messageManager::instance()->message(" Number of trees found:", numTrees);

    //Only plot two longest, straightest trees
    if (nTrees > 6) nTrees = 6;

    while(nTrees--)
    {
        TkrLinkTree* pTree = &(*treePtr++);

        //Draw the full tree
        setColor(pColors[colorIdx]);
        drawFullTree(pTree->getHeadNode());

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

        setColor(pColors[colorIdx]);
        markerAt(Point(x,y,z));

        colorIdx = ++colorIdx % 6;
    }

    return;
}

void TkrCandidatesRep::drawFullTree(LayerLinkNode* pNode)
{
    //First, draw the children of this node
    int numChildren = pNode->getNumKids();

    while(numChildren--) drawFullTree(pNode->getThisKid(numChildren));
    
    //Ok, now draw this link
    TkrLinkNode* pTkrNode = dynamic_cast<TkrLinkNode*>(pNode);

    drawLinkNode(pTkrNode);

    return;
}

void TkrCandidatesRep::drawLinkNode(TkrLinkNode* pTkrNode)
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
