//------------------------------------------------------------------------------
// TkrNeuralNet Implementation
//
// Search for track candidates using a Neural Network.  Based on the 
// work of Stimpfl-Abele & Garrido Comp. Phys. Commun. 64 (1991) 46-50
// 
// b. allgood and w. atwood, 3/02 
//------------------------------------------------------------------------------

#include "src/PatRec/NeuralNet/TkrNeuralNet.h"
#include "src/Track/TkrControl.h"
#include "src/Utilities/TkrPoints.h"
#include "src/Utilities/TkrPoint.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include <algorithm>


// I would like these to be passed from the joboptions file.
#define MAX_LAYER_SEPERATION 3
#define MAX_PITCH            0.2
#define NEURON_LAMBDA    5     // tunable parameter for neuron angle
#define NEURON_MU        2     // tunable parameter for layer seperation 
#define NEURON_ALPHA_UP  5     // tunable parameter for up connections
#define NEURON_ALPHA_DO  5     // tunable parameter for down connections
#define NEURON_BIAS      0.2   // tunable parameter for stimulation
#define NEURON_GAMMA     10    // tunable parameter for weight
#define TEMP             1     // tempurature for relaxation process

// Right now it doesn't use calEne, but in the future it might.
TkrNeuralNet::TkrNeuralNet(ITkrGeometrySvc* pTkrGeo, TkrClusterCol* pClusters, 
                           double calEne, Point calHit): m_Pcal(calHit), 
                           m_energy(calEne), m_tkrGeo(pTkrGeo), 
                           m_clusters(pClusters)
{ 
    // make the neurons
    m_numNeurons = generateNeurons(); 

    // make the neural net
    buildNet();

    if(m_numNeurons > 0){

        // relax the neurons
        relax();

        // build the candidate tracks
        buildCand();
    }

    return;
}

/*-------------------Private-------------------*/

// generateNeurons() 
// This function first gathers the cluster pairs and stores them in 
// m_pointList.  It then constructs neurons and adds them to m_neuronList
// if they meet the criteria of layer seperation and max pitch angle.
// Each neuron is assigned a random bias using the standard rand() C++ call.
// We may want to have this use randSvc() from Gaudi?  This returns
// the total number of neurons created.
unsigned int TkrNeuralNet::generateNeurons()
{
    // This gets all of the TkrPoints in to a vector
    for (int ilayer = 0 ; ilayer < 18; ilayer++)
    {
        TkrPoints tempTkrPoints(ilayer, m_clusters);
        if(!tempTkrPoints.finished()){
            TkrPointList tmpList = tempTkrPoints.getAllLayerPoints();
            m_pointList.insert(m_pointList.end(),tmpList.begin(),tmpList.end());
        }
    }

    m_numPoints = m_pointList.size();

    // This creates the neuron list.  Because the point list is in order of 
    // layer (0-17) the search only looks for pairs which are further down 
    // the list.  Because a neuron with two points in the same layer makes no
    // sense we exclude those and because neurons which consist of points 
    // seperated by more than MAX_LAYER_SEPERATION are very unlikely we exclude
    // them.  Finally, neurons with an angle greater than MAX_PITCH are very 
    // unlikely we exclude them as well.

    srand((unsigned int)(time(NULL)));

    for (unsigned int iTop = 0; iTop < m_numPoints; iTop++)
    {
        for(unsigned int iBottom = iTop+1; iBottom < m_numPoints; iBottom++)
        {
            TkrNeuron neuron(&m_pointList[iTop],&m_pointList[iBottom],
                             (float) rand()/(RAND_MAX+1.0));

            if (neuron.getLayerDiff() == 0) continue;
            if (neuron.getLayerDiff() > MAX_LAYER_SEPERATION) break;
            if((neuron.getDirection()).z() < MAX_PITCH) continue;

            neuron.setBias((float) NEURON_BIAS);
            m_neuronList.push_back(neuron);
        }
    }

    return (unsigned int) m_neuronList.size();
}

// buildNet() 
// This function constructs the synapse lists for all of the neurons.  The
// synapse lists contains pointers to other neurons sharing the same end
// points and the wieghts associated with that connection.
void TkrNeuralNet::buildNet()
{
    //float weight=0.0;

    for(unsigned int i = 0; i < m_numNeurons; i++)
    {
        for(unsigned int j = i+1; j < m_numNeurons; j++)
        {
            // do they share the same top point?
            if(m_neuronList[i].getPnt(top) == m_neuronList[j].getPnt(top))
            {
                m_neuronList[i].addConnection(&TkrNeuralNet::calcWieght, 
                    this, top, &m_neuronList[j]);
                m_neuronList[j].addConnection(&TkrNeuralNet::calcWieght, 
                    this, top, &m_neuronList[i]);

            }
            // do they share the same bottom point? 
            else if(m_neuronList[i].getPnt(bottom) == m_neuronList[j].getPnt(bottom))
            {
                m_neuronList[i].addConnection(&TkrNeuralNet::calcWieght, 
                    this, bottom, &m_neuronList[j]);
                m_neuronList[j].addConnection(&TkrNeuralNet::calcWieght, 
                    this, bottom, &m_neuronList[i]);
            }
            // are they consecutive neurons?
            else if(m_neuronList[i].getPnt(bottom) == m_neuronList[j].getPnt(top))
            {
                m_neuronList[i].addConnection(&TkrNeuralNet::calcWieght, 
                    this, bottom, &m_neuronList[j]);
                m_neuronList[j].addConnection(&TkrNeuralNet::calcWieght, 
                    this, top, &m_neuronList[i]);
            }
        }
    }

    return;
}

// relax()
// This function evolves the neural net to the ground state.  This is done
// by shuffling the neurons and then updating the neurons using the update 
// function (explained below).  This is repeated until the normalized 
// difference in activity of the last two update rounds is less then some
// convergence value.
void TkrNeuralNet::relax()
{
    double prevActivity = 0.0;
    double cumActivityDiff = 0.0;
    double converge = 0.0005;
    std::vector<TkrNeuron*> tmpList;

    unsigned int i;
    for(i=0;i < m_numNeurons; i++) tmpList.push_back(&m_neuronList[i]);

    do{

        // shuffle the list
        std::random_shuffle(tmpList.begin(),tmpList.end());  //just for now

        cumActivityDiff = 0.0;

        // update neurons from the shuffled list
        for(i=0;i < m_numNeurons; i++){
            prevActivity = tmpList[i]->getActivity();
            updateNeuron(tmpList[i],(float) TEMP);
            cumActivityDiff += fabs(prevActivity - tmpList[i]->getActivity());
        }

        float tmp = 0.0;

        // readjust the biases in each neuron to help speed up convergence.
        for(i=0;i < m_numNeurons; i++){
            tmp = 0.0;
	    unsigned int j;
            for(j = 0; j < tmpList[i]->numSynapse(top); j++) 
                tmp += (tmpList[i]->getNextNeuron(top,j))->getActivity();
	    unsigned int k;
            for(k = 0; k < tmpList[i]->numSynapse(bottom); k++) 
                tmp += (tmpList[i]->getNextNeuron(bottom,k))->getActivity();
            tmp *= 4/((float)(j+k));
            if(tmp <= 0.4) tmpList[i]->setBias(tmp);
            else tmpList[i]->setBias((float) NEURON_BIAS);
        }

    }while((cumActivityDiff/((double)m_numNeurons)) > converge);

}

// buildCand()
// This is function is based on one taken from the TkrCombo PatRec routine.
// This funciton takes all neurons with an activity above 0.9 and places
// them in m_candidates.  It then proceed to do some preliminary fitting
// of the track.  This will be changed soon.
void TkrNeuralNet::buildCand()
{

    // take neurons with activity above 0.9.
    for(unsigned int i = 0;i < m_numNeurons; i++){
        if(m_neuronList[i].getActivity() >= 0.9){
            Candidate trial = Candidate( m_neuronList[i].getLayer(top), m_neuronList[i].getTower(top), .03, 
                m_neuronList[i].getPnt(top),m_neuronList[i].getDirection());
            m_candidates.push_back(trial);
        }
    }

        TkrControl * control = TkrControl::getPtr(); 
    if (m_candidates.size() > 0) {
        
        TkrNeuralNet::const_iterator hypo;
            
            for(hypo  = m_candidates.begin(); 
                hypo != m_candidates.end();   hypo++){
                
                int   iniLayer = (*hypo).firstLayer();
                int   iniTower = (*hypo).tower();
                Ray   testRay  = (*hypo).ray();
                float energy   = (*hypo).energy();
                

                KalFitTrack* _track = new KalFitTrack(m_clusters, m_tkrGeo, iniLayer, iniTower, 
                                           control->getSigmaCut(), energy, testRay); 

                _track->findHits();
                _track->doFit();

                if (!_track->empty(control->getMinSegmentHits())) 
                {
                    //Keep pointer to the track temporarily
                    m_tracks.push_back(_track);

                    //Keep this track (but as a candidate)
                    TkrPatCand* newTrack = new TkrPatCand(_track->getLayer(),
                        _track->getTower(),energy,1.,_track->getQuality(),_track->getRay());

                    push_back(newTrack);

                    //_track->flagAllHits();

                } 
                else delete _track;
            }
    }

    // Ok, go through all the attempted track fits and unflag the hits for the 
    // real fit
    if (m_tracks.size())
    {
        TkrFitTrackCol::iterator iter = m_tracks.begin();

        while(iter != m_tracks.end())
        {
            KalFitTrack* pTrack = (KalFitTrack*)(*iter++);

            pTrack->unFlagAllHits();
        }
    }

    return;
}

// calcWieght()
// This function returns the weight value for the two neurons' connection.
float TkrNeuralNet::calcWieght(TkrNeuron* neuron1, TkrNeuron* neuron2)
{
    // If the neurons are joined at the bottom or at the top the weight is zero.
    if((neuron1->getPnt(top) == neuron2->getPnt(top)) || 
        (neuron1->getPnt(bottom) == neuron2->getPnt(bottom))) return 0.0;

    float tmp = pow(neuron1->getDirection().dot(neuron2->getDirection()),
        (double) NEURON_LAMBDA);

    // to avoid negative weights, retrun 0 if agnle is > 90 degrees.
    if(tmp < 0.8) return 0.0;  

    tmp /= pow((double) neuron1->getLayerDiff(),(double) NEURON_MU) + 
        pow((double) neuron2->getLayerDiff(),(double) NEURON_MU);

    return ((float) tmp);
}

// updateNeurons()
// This function updates a given neurons activity level based on the update
// rule.  Neurons which are incoming or outgoing reinforce the neuron.
// Neurons which are competing detract from the neurons activity.
void TkrNeuralNet::updateNeuron(TkrNeuron* neuron, float temp)
{

    float tmp1 = 0.0;  // for reinforcing from the neurons activity.
    float tmp2 = 0.0;  // for detracting from the neurons activity.
    float tmp3 = 0.0;  // for detracting from the neurons activity.

    unsigned int i;
    unsigned int j;
    for(i=0, j=0;i < neuron->numSynapse(bottom); i++){
        if(neuron->getNextNeuron(bottom, i)->getPnt(top) == neuron->getPnt(bottom)) 
            tmp1 +=(neuron->getWieght(bottom, i)) * 
                   (neuron->getNextNeuron(bottom, i)->getActivity());
        else{
            tmp2 += neuron->getNextNeuron(bottom, i)->getActivity();
            j++;
        }
    }

    for(i = 0;i < neuron->numSynapse(top); i++){
        if(neuron->getNextNeuron(top, i)->getPnt(bottom) == neuron->getPnt(top))
            tmp1 +=(neuron->getWieght(top, i)) * 
                   (neuron->getNextNeuron(top, i)->getActivity());
        else{
            tmp3 += neuron->getNextNeuron(top, i)->getActivity();
        }
    }

    tmp1 *= (float) NEURON_GAMMA;
    tmp2 *= (float) NEURON_ALPHA_UP;
    tmp3 *= (float) NEURON_ALPHA_DO;

    // update the neuron
    neuron->setActivity((0.5*(1+tanh((1/temp)*(tmp1-tmp2-tmp3+(neuron->getBias()))))));

    return;

}
