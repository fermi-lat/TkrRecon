/*
    Code to implement the TkrLinkNode class
    Tracy Usher Dec 4, 2000
*/

#include <math.h>
#include "src/PatRec/LinkAndTree/TkrLinkNode.h"
//#include "instrument/calorimeterGeo.h"

//Null creator
TkrLinkNode::TkrLinkNode()
{
    pLayerLink    = 0;
    pParent       = 0;
    nKidsTotal    = 0;
    bestPathsMask = 0;
    pChildren.clear();

    pQuality      = new TkrNodeQuality();

    keepNode      = true;

    return;
}

//Creates a link node with pointer to a link
TkrLinkNode::TkrLinkNode(LayerLink* pLink)
{
    pLayerLink    = pLink;
    pParent       = 0;
    nKidsTotal    = 0;
    bestPathsMask = 0;
    pChildren.clear();

    pQuality      = new TkrNodeQuality();

    keepNode      = true;

    if (pLink->getLinkNode() == 0) pLink->setLinkNode(this);

    return;
}

//Creates a link node with pointer to a link
//with angles initialized by node passed in
TkrLinkNode::TkrLinkNode(LayerLink* pLink, TkrNodeQuality* pQual, double angle) 
{
    pLayerLink    = pLink;
    pParent       = 0;
    nKidsTotal    = 0;
    bestPathsMask = 0;
    pChildren.clear();

    pQuality      = new TkrNodeQuality(pQual, angle);

    if (pLink->getLinkNode() == 0) pLink->setLinkNode(this);

    return;
}

TkrLinkNode::~TkrLinkNode()
{
    if (pLayerLink) 
    {
        if (pLayerLink->getLinkNode() == this) pLayerLink->setLinkNode(0);
    }

    if (pQuality)   delete pQuality;

    return;
}


//This destructor makes sure all child nodes are deleted
LayerLinkNode::~LayerLinkNode()
{
    return;
}

void TkrLinkNode::remThisKid(LayerLinkNode* pNode)
{
    int               numKidsInVector = pChildren.size();
    linkNodeVectorPtr pVector         = pChildren.begin();

    //This removes the node we want to get rid of
    while(numKidsInVector--)
    {
        if (pNode == *pVector)
        {
            pChildren.erase(pVector);
            break;
        }
        pVector++;
    }

    //Reset the number of kids counter (which really exists for debugging)
    nKidsTotal = getNumKids();

    //Now set the next best node as the primary path node
    //(helps prevent it from getting stolen in arbitration)
    if (nKidsTotal)
    {
        TkrLinkNode*    pBest     = dynamic_cast<TkrLinkNode*>(getBestKid());
        TkrNodeQuality* pBestQual = pBest->getQuality();

        pBestQual->setNumSince(pBestQual->getNumNodes());
    }

    return;
}


//Return pointer to the best child in our list
LayerLinkNode* TkrLinkNode::getBestKid()
{
    int               nKids     = pChildren.size() - 1;
    LayerLinkNode*    pNode     = pChildren[nKids];
    TkrLinkNode*      pBest     = dynamic_cast<TkrLinkNode*>(pNode);
    TkrNodeQuality*   pBestQual = pBest->getQuality();

    while(nKids--)
    {
        pNode = pChildren[nKids];
        TkrLinkNode* pNext = dynamic_cast<TkrLinkNode*>(pNode);

        if (!pBestQual->isNodeWorseRMS(pNext)) pBest = pNext;
    }

    return pBest;
}


//Return pointer to the worst child in our list
LayerLinkNode* TkrLinkNode::getWorstKid()
{
    int               nKids      = pChildren.size() - 1;
    LayerLinkNode*    pNode      = pChildren[nKids];
    TkrLinkNode*      pWorse     = dynamic_cast<TkrLinkNode*>(pNode);
    TkrNodeQuality*   pWorseQual = pWorse->getQuality();

    while(nKids--)
    {
        pNode = pChildren[nKids];
        TkrLinkNode* pNext = dynamic_cast<TkrLinkNode*>(pNode);

        if (pWorseQual->isNodeWorseRMS(pNext)) pWorse = pNext;
    }

    return pWorse;
}

//Delete the worst kid in our list
void TkrLinkNode::delWorstKid()
{
    int               nKids      = pChildren.size() - 1;
    linkNodeVectorPtr pNodeVec   = pChildren.begin();
    linkNodeVectorPtr pWorseVec  = pNodeVec;
    LayerLinkNode*    pNode      = *pNodeVec++;
    TkrLinkNode*      pWorse     = dynamic_cast<TkrLinkNode*>(pNode);
    TkrNodeQuality*   pWorseQual = pWorse->getQuality();

    while(nKids--)
    {
        pNode = *pNodeVec;
        TkrLinkNode* pNext = dynamic_cast<TkrLinkNode*>(pNode);

        if (pWorseQual->isNodeWorseRMS(pNext)) 
        {
            pWorse    = pNext;
            pWorseVec = pNodeVec;
        }

        pNodeVec++;
    }

    delete pWorse;
    pChildren.erase(pWorseVec);

    //Reset the number of kids counter (which really exists for debugging)
    nKidsTotal = getNumKids();

    //Now set the next best node as the primary path node
    //(helps prevent it from getting stolen in arbitration)
    if (nKidsTotal)
    {
        TkrLinkNode*    pBest     = dynamic_cast<TkrLinkNode*>(getBestKid());
        TkrNodeQuality* pBestQual = pBest->getQuality();

        pBestQual->setNumSince(pBestQual->getNumNodes());
    }

    return;
}

bool TkrLinkNode::isOldPathStraighter(LayerLinkNode* pOld)
{
    bool         isBetter     = false;
    TkrLinkNode* pOldLinkNode = dynamic_cast<TkrLinkNode*>(pOld);
    int          depthDiff    = pQuality->nodeDepthDiff(pOldLinkNode);

    if (depthDiff >= 0)
    {
        if (depthDiff > 0 || !pQuality->isNodeWorseRMS(pOldLinkNode)) isBetter = true;
    }

    return isBetter;
}

bool TkrLinkNode::isOldPathBetter(LayerLinkNode* pOld)
{
    bool         isBetter     = false;
    TkrLinkNode* pOldLinkNode = dynamic_cast<TkrLinkNode*>(pOld);
    //int          depthDiff    = pQuality->branchDepthDiff(pOldLinkNode);
    int          depthDiff    = pQuality->nodeDepthDiff(pOldLinkNode);

    if (depthDiff >= 0)
    {
        if (depthDiff > 0 || !pQuality->isNodeWorseRMS(pOldLinkNode)) isBetter = true;
    }

    return isBetter;
}

bool TkrLinkNode::isOldNodeBetter(LayerLinkNode* pOld)
{
    bool         isBetter     = false;
    TkrLinkNode* pOldLinkNode = dynamic_cast<TkrLinkNode*>(pOld);

    //Get the difference in length to this node. If depthDiff is negative then
    //the new node (ie, this one) is longer
    int          depthDiff = pQuality->nodeDepthDiff(pOldLinkNode);

    //If old node is longer then its better
    if (depthDiff >= 0)
    {
        if (depthDiff > 0 || !pQuality->isNodeWorseAVE(pOldLinkNode)) isBetter = true;
    }

    return isBetter;
}


bool TkrLinkNode::isOldAngleSmaller(LayerLinkNode* pOld)
{
    bool            isBetter     = false;
    TkrLinkNode*    pOldLinkNode = dynamic_cast<TkrLinkNode*>(pOld);
    TkrNodeQuality* pOldQuality  = pOldLinkNode->getQuality();

    if (fabs(pQuality->getLastAngle()) > fabs(pOldQuality->getLastAngle())) isBetter = true;

    return isBetter;
}

//Checks to see if a link matches the current link
//If reasonable match then create a candidate new node with this link
LayerLinkNode* TkrLinkNode::doLinksMatch(LayerLink* pLinkToMatch)
{
    TkrLinkNode*    pNewLinkNode  = 0;

    TkrClusterLink* pTkrLinkToMatch = dynamic_cast<TkrClusterLink*>(pLinkToMatch);

    //If link in use by a previous tree then skip
    if (!pTkrLinkToMatch->linkInUse())
    {
        TkrClusterLink* pCurLink = dynamic_cast<TkrClusterLink*>(pLayerLink);

        //Require exact match between clusters
        if (pCurLink->pBotClus() == pTkrLinkToMatch->pTopClus())
        {
            double newAngle = pLayerLink->angleWith(pLinkToMatch);

            //Link must lie within a reasonable angle of current link
            if (fabs(newAngle) <= pQuality->getAngleMAX()) pNewLinkNode = new TkrLinkNode(pLinkToMatch, pQuality, newAngle);
        }
    }

    return pNewLinkNode;
}

//Provides a more restrictive test of a candidate new node
bool TkrLinkNode::keepNewNode(LayerLinkNode* pNewNode)
{
    bool            keepNode      = false;
    TkrLinkNode*    pNewLinkNode  = dynamic_cast<TkrLinkNode*>(pNewNode);
    LayerLink*      pNewLinkToAdd = pNewLinkNode->getThisLayerLink();
    double          newAngle      = pLayerLink->angleWith(pNewLinkToAdd);

    //More restrictive test on angle to see if we want to keep this link
    if (pQuality->keepNewLink(newAngle))
    {
        keepNode = true;
        //If link is already in use then recover pointer to that node
        if (pNewNode->isOldLinkBetter()) keepNode = false;
    }

    return keepNode;
}

//A candidate link may already be in use. Recover that node and 
//check if adding to the new tree will be "better"
bool TkrLinkNode::isOldLinkBetter()
{
    bool            oldLinkBetter = false;

    //If link is already in use then recover pointer to that node
    LayerLinkNode*  pCurLinkNode = pLayerLink->getLinkNode();

    //If this link previously used then we need to arbitrate
    if (pCurLinkNode != 0 && pCurLinkNode != this)
    {
        //If the previous node is better then mark current node
        //to be deleted if it doesn't have any other offspring
        //if (pNewLinkNode->isOldNodeBetter(pCurLinkNode))
        if (isOldPathStraighter(pCurLinkNode))
        {
            if (getNumKids() == 0) clearKeepNode();

            oldLinkBetter = true;
        }
        //Otherwise life is a bit more complicated as we have to
        //move the link to this node and clear it from the other's 
        //parent
        else
        {
            //Get the pointer to the parent 
            LayerLinkNode* pOldParent = pCurLinkNode->getParentNode();

            //Get the pointer to the old parent node and clear this link
            TkrLinkNode*   pOldOne    = dynamic_cast<TkrLinkNode*>(pOldParent);

            pOldOne->remThisKid(pCurLinkNode);

            //If old parent has no children left then mark it for deletion
            if (pOldOne->getNumKids() == 0) pOldOne->clearKeepNode();

            //Delete all the nodes from this link and beyond
            delNodes(pCurLinkNode);
        }
    }

    return oldLinkBetter;
}


//Do not allow to links to end on the same cluster
//This will arbitrate between two links that wish to do so
bool TkrLinkNode::isOldClusBetter(layerLinkVector* pLayerLinkVector)
{
    bool               oldClusBetter = false;

    //Set up to loop through all vectors in this layer
    layerLinkVectorPtr pLinkVector   = pLayerLinkVector->begin();
    int                nLinks        =  pLayerLinkVector->size();

    //Loop over all the links in this layer (yawn...)
    while(nLinks--)
    {
        LayerLink* pOldLink = *pLinkVector++;

        //Look for a link which shares the bottom cluster (but is not the same one!)
        if (pLayerLink != pOldLink && pLayerLink->sameBotCluster(pOldLink))
        {
            LayerLinkNode* pOldLinkNode = pOldLink->getLinkNode();

            //If this link previously used then we need to arbitrate
            if (pOldLinkNode != 0)
            {
                //If the previous node is better then mark current node
                //to be deleted if it doesn't have any other offspring
                if (isOldPathBetter(pOldLinkNode))
                {
                    if (getNumKids() == 0) clearKeepNode();

                    oldClusBetter = true;

                    break;
                }
                //Otherwise life is a bit more complicated as we have to
                //move the link to this node and clear it from the other's 
                //parent
                else
                {
                    //Get the pointer to the parent 
                    LayerLinkNode* pOldParent = pOldLinkNode->getParentNode();

                    //Get the pointer to the old parent node and clear this link
                    TkrLinkNode*   pOldOne    = dynamic_cast<TkrLinkNode*>(pOldParent);

                    pOldOne->remThisKid(pOldLinkNode);

                    //If old parent has no children left then mark it for deletion
                    if (pOldOne->getNumKids() == 0) pOldOne->clearKeepNode();

                    //Delete all the nodes from this link and beyond
                    delNodes(pOldLinkNode);

                    break;
                }
            }
        }
    }

    return oldClusBetter;
}



//Adds a the new node to the list of children of this parent
void TkrLinkNode::addNewNode(LayerLinkNode* pNewNode, int branchIdx)
{
    TkrLinkNode*    pNewLinkNode  = dynamic_cast<TkrLinkNode*>(pNewNode);
    LayerLink*      pNewLinkToAdd = pNewLinkNode->getThisLayerLink();

    pNewLinkNode->setParentNode(this);
    pNewLinkNode->setKeepNode();
    setNewChildNode(pNewLinkNode);
    setKeepNode();

    if (branchIdx > 0) (pNewLinkNode->getQuality())->setNumSince(0);

    pNewLinkToAdd->setLinkNode(pNewNode);

    nKidsTotal = getNumKids();
 
    return;
}

void TkrLinkNode::delNodes(LayerLinkNode* pNode)
{
    if (pNode)
    {
        //Make sure link is made available again
        LayerLink* pLink = pNode->getThisLayerLink();
        pLink->setLinkNode(0);

        //Go through and get rid of the kids
        int nKids = pNode->getNumKids();

        while(nKids--) delNodes(pNode->getThisKid(nKids));

        delete pNode;
    }
    return;
}

bool TkrLinkNode::keepThisNode()
{
    bool returnVal = false;

    if (keepNode) returnVal = true;

    return returnVal;
}

void TkrLinkNode::setKeepNode()
{
    keepNode = true;
    return;
}

void TkrLinkNode::clearKeepNode()
{
    keepNode = false;
    return;
}

void TkrLinkNode::setLinkUsedList(int id)
{
    if (id >= 0 && id <= 31) bestPathsMask |= 1 << (id - 1);
    return;
}

int  TkrLinkNode::getNodeDepth()
{
    return pQuality->getNumNodes();
}
