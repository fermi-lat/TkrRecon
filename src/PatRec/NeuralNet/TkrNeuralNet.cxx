//------------------------------------------------------------------------------
// TkrNeuralNet Implementation
//
// Search for track candidates using a Neural Network.  Based on the 
// work of Stimpfl-Abele & Garrido Comp. Phys. Commun. 64 (1991) 46-50
// 
// b. allgood and w. atwood, 3/02 
//------------------------------------------------------------------------------

#include "src/PatRec/NeuralNet/TkrNeuralNet.h"

#include "src/Utilities/TkrPoints.h"
#include "src/Utilities/TkrPoint.h"

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <algorithm>


// Right now it doesn't use calEne, but in the future it might.
TkrNeuralNet::TkrNeuralNet(Event::TkrClusterCol* pClusters, std::map<const char*,double,ltstr>& params,
                           double calEne, Point calHit) : 
  m_numNeurons(0), m_Pcal(calHit), m_energy(calEne), 
  m_clusters(pClusters), m_params(params)
{ 
  // make the neurons
  m_numNeurons = generateNeurons(); 
  
  // make the neural net
  buildNet();
  
  if(m_numNeurons > 0) relax();
  
  return;
}

/*-------------------Private-------------------*/

// generateNeurons() 
// This function first gathers the cluster pairs and stores them in 
// m_pointList.  It then constructs neurons and adds them to m_neuronList
// if they meet the criteria of layer separation and max pitch angle.
// Each neuron is assigned a random bias using the standard rand() C++ call.
// We may want to have this use randSvc() from Gaudi?  This returns
// the total number of neurons created.
unsigned int TkrNeuralNet::generateNeurons()
{
  // This gets all of the TkrPoints into a vector
  for (int ilayer = 0 ; ilayer < 18; ilayer++)
    {
      TkrPoints tempTkrPoints(ilayer, m_clusters);
      if(!tempTkrPoints.finished())
	{
	  TkrPointList tmpList = tempTkrPoints.getAllLayerPoints();
	  m_pointList.insert(m_pointList.end(),tmpList.begin(),tmpList.end());
	}
    }
  
  unsigned int numPoints = m_pointList.size();
  
  // This creates the neuron list.  Because the point list is in order of 
  // layer (0-17) the search only looks for pairs which are further down 
  // the list.  Because a neuron with two points in the same layer makes no
  // sense we exclude those and because neurons which consist of points 
  // separated by more than MAX_LAYER_SEPARATION are very unlikely we exclude
  // them.  Finally, neurons with an angle greater than MAX_PITCH are very 
  // unlikely we exclude them as well.
  
  srand((unsigned int)(time(NULL)));
  
  double maxLayerDiff = m_params["maxLayerDiff"];
  double maxPitch     = m_params["maxPitch"];
  double bias         = m_params["Bias"];
  
  for (unsigned int iTop = 0; iTop < numPoints; iTop++)
    {
      for(unsigned int iBottom = iTop+1; iBottom < numPoints; iBottom++)
        {
	  TkrNeuron neuron(&m_pointList[iTop],&m_pointList[iBottom],
			   (float) rand()/(RAND_MAX+1.0));
	  
	  if (neuron.getLayerDiff() == 0) continue;
	  if (neuron.getLayerDiff() > maxLayerDiff) break;
	  if((neuron.getDirection()).z() < maxPitch) continue;
	  
	  neuron.setBias((float) bias);
	  m_neuronList.push_back(neuron);
        }
    }
  
  return (unsigned int) m_neuronList.size();
}


// buildNet() 
// This function constructs the synapse lists for all of the neurons.  The
// synapse lists contains pointers to other neurons sharing the same end
// points and the weights associated with that connection.
void TkrNeuralNet::buildNet()
{
  
  double lambda = m_params["Lambda"];
  double mu     = m_params["Mu"];

  for(unsigned int i = 0; i < m_numNeurons; i++)
    {
      for(unsigned int j = i+1; j < m_numNeurons; j++)
        {
	  // do they share the same top point?
	  if(m_neuronList[i].getPnt(top) == m_neuronList[j].getPnt(top))
            {
	      m_neuronList[i].addConnection(top, &m_neuronList[j], lambda, mu);
	      m_neuronList[j].addConnection(top, &m_neuronList[i], lambda, mu);
            }
	  // do they share the same bottom point? 
	  else if(m_neuronList[i].getPnt(bottom) == m_neuronList[j].getPnt(bottom))
            {
	      m_neuronList[i].addConnection(bottom, &m_neuronList[j], lambda, mu);
	      m_neuronList[j].addConnection(bottom, &m_neuronList[i], lambda, mu);
            }
	  // are they consecutive neurons?
	  else if(m_neuronList[i].getPnt(bottom) == m_neuronList[j].getPnt(top))
	    {
	      m_neuronList[i].addConnection(bottom, &m_neuronList[j], lambda, mu);
	      m_neuronList[j].addConnection(top, &m_neuronList[i],    lambda, mu);
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
  
  double bias       = m_params["Bias"];
  double temp       = m_params["temperature"];
  double gamma      = m_params["Gamma"];
  double alpha_up   = m_params["AlphaUP"];
  double alpha_down = m_params["AlphaDOWN"];


  unsigned int i;
  for(i=0;i < m_numNeurons; i++) tmpList.push_back(&m_neuronList[i]);
  
  do{
    // shuffle the list
    std::random_shuffle(tmpList.begin(),tmpList.end());  //just for now
    
    cumActivityDiff = 0.0;
    
    // update neurons from the shuffled list
    for(i=0; i < m_numNeurons; i++)
      {
	prevActivity = tmpList[i]->getActivity();
	
	tmpList[i]->update((float) temp, gamma, alpha_up, alpha_down);
	
	cumActivityDiff += fabs(prevActivity - tmpList[i]->getActivity());
      }
    
    float tmp = 0.0;
    
    // readjust the biases in each neuron to help speed up convergence.
    for(i=0;i < m_numNeurons; i++)
      {
	tmp = 0.0;
	unsigned int j;
	for(j = 0; j < tmpList[i]->numSynapse(top); j++) 
	  tmp += (tmpList[i]->getNextNeuron(top,j))->getActivity();
	
	unsigned int k;
	for(k = 0; k < tmpList[i]->numSynapse(bottom); k++) 
	  tmp += (tmpList[i]->getNextNeuron(bottom,k))->getActivity();
	
	tmp *= 4/((float)(j+k));
	
	if(tmp <= 0.4) 
	  tmpList[i]->setBias(tmp);
	else 
	  tmpList[i]->setBias((float) bias);
      }
    
  }while((cumActivityDiff/((double)m_numNeurons)) > converge);
}

