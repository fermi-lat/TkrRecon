/*
	Code to implement the PatRecTracks class
	Tracy Usher Nov 27, 2000
*/

#include "TkrRecon/PatRec/TkrCandidates.h"

TkrCandidates::TkrCandidates(ITkrGeometrySvc* pTkrGeo, TkrClusters* pTkrClus)
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


	return;
}


TkrCandidates::~TkrCandidates()
{
	clear();

	return;
}

void TkrCandidates::setLinkList(TkrClusterLinkList* pLinks, TkrPlaneType plane)
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

TkrClusterLinkList* TkrCandidates::getLinkList(TkrPlaneType plane)
{
	return plane == X ? pLinkListX : pLinkListY;
}

int TkrCandidates::getNumLinks(TkrPlaneType plane)
{
	return plane == X ? numLinksX : numLinksY;
}

void TkrCandidates::setForest(TkrLinkForest* pForest, TkrPlaneType plane)
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

TkrLinkForest* TkrCandidates::getForest(TkrPlaneType plane)
{
	return plane == X ? pForestX : pForestY;
}

int TkrCandidates::getNumTrees(TkrPlaneType plane)
{
	return plane == X ? numTracksX : numTracksY;
}

void TkrCandidates::ini()
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

void TkrCandidates::clear()
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
void TkrCandidates::writeOut() const
//########################################################
{
	if (numTracksX <=0 && numTracksY <=0) return;
	
	std::ostream& out = messageManager::instance()->out();
	if (!messageManager::instance()->acceptLevel("DEBUG")) return;

	out << " --- TkrCandidates ---- " << "\n";
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
void TkrCandidates::update(GraphicsRep& v)
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