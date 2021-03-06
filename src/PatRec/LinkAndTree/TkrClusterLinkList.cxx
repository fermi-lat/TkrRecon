/*
Code to implement the SiLinkList class
Tracy Usher Nov 30, 2000
*/

#include "src/PatRec/LinkAndTree/TkrClusterLinkList.h"

//Ugliness that will need to be dealt with
#define NLAYERS     18  //fix
#define NLINKLAYERS 17  //fix


TkrClusterLinkList::TkrClusterLinkList()
{
    layerView   = UNDEFINED;

    return;
}

TkrClusterLinkList::TkrClusterLinkList(ITkrQueryClustersTool* clusTool, TkrPlaneType plane)
{
    layerView       = plane;
    numLinksTotal   =  0;
    mostLinksLayer  = -1;
    numLinksInLayer =  0;

    //Do nothing if not a legal plane
    if (plane == X || plane == Y)
    {
        //      int layerNum = plane == X ? NLAYERS : NLAYERS - 1;
        int layerNum = NLAYERS - 1;

        //Loop is over "link" layers - incorporate two silicon layers
        int ilayer = 0;
        for (;ilayer<NLAYERS; ++ ilayer) {
            //while(layerNum--)
            //std::vector<TkrCluster*> pClus = pClusters->getHits((TkrCluster::view)plane, layerNum);
            Event::TkrClusterVec pClus = clusTool->getClusters(plane, ilayer);

            if (pClus.size() > 0)
            {
                TkrClusterLinkVector* pClusLinks = new TkrClusterLinkVector(clusTool, ilayer, plane);
                LayerLinkVector*      pVector    = pClusLinks;
                int                   numVectors = pVector->size();

                numLinksTotal += numVectors;
                if (numVectors > numLinksInLayer)
                {
                    numLinksInLayer = numVectors;
                    mostLinksLayer  = layerNum;
                }

                push_front(pVector);
            }
        }

        //If we have some links in our list then create dummy terminating link
        if (size())
        {
            TkrClusterLinkVector*  pClusLinks = new TkrClusterLinkVector();
            LayerLinkVector*      pVector    = pClusLinks;

            push_back(pVector);
        }
    }

    return;
}

//TkrClusterLinkList the destroyer
TkrClusterLinkList::~TkrClusterLinkList()
{
    //Loop over elements in the list and call their destructors
    int nListElems = size();

    if (nListElems)
    {
        layerLinkListPtr pList = begin();

        while(nListElems--) delete *pList++;

        clear();
    }

    return;
}

//TkrClusterLinkList::~TkrClusterLinkList()
LayerLinkList::~LayerLinkList()
{
    return;
}


