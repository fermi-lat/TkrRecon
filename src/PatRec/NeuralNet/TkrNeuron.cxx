//------------------------------------------------------------------------------
// TkrNeuron for TkrNeuralNet Implementation
//
// Container class for Neurons to be used in the Neural Net PR routine.
// This class contains information about possible subtracks (neurons)
// which are used to determine track candidates by the Nerual Net PR
// routine for the Kalman filter.
//  
// b. allgood and w. atwood 3/02
//------------------------------------------------------------------------------

#include "src/PatRec/NeuralNet/TkrNeuron.h"

// constructor
TkrNeuron::TkrNeuron(TkrPoint* pnt0, TkrPoint* pnt1, float act=0.0):
  m_pnt0(pnt0),m_pnt1(pnt1)
{
    m_l = pnt1->getLayer() - pnt0->getLayer();
    m_direction = (pnt0->getPoint() - pnt1->getPoint()).unit();
    m_bias = 0.0;
    m_activity = act;

}

// assignment oper.
TkrNeuron& TkrNeuron::operator=(const TkrNeuron& neuron)
{
    m_pnt0         = neuron.m_pnt0;
    m_pnt1         = neuron.m_pnt1;
    m_bias         = neuron.m_bias;
    m_activity     = neuron.m_activity;
    m_l            = neuron.m_l;
    m_direction    = neuron.m_direction;
    m_synapseList0.insert(m_synapseList0.end(),neuron.m_synapseList0.begin(),
                          neuron.m_synapseList0.end());
    m_synapseList1.insert(m_synapseList1.end(),neuron.m_synapseList1.begin(),
                          neuron.m_synapseList1.end());
    return *this;
}

// getNextNeuron 
// param: the neuron pointed to by the synapse is != NULL.
TkrNeuron* TkrNeuron::getNextNeuron(position pos, unsigned int syn) const
{
    TkrNeuron* tmp;
    tmp=NULL;

    if(pos == top){
        if(m_synapseList0.size() > syn) tmp=m_synapseList0[syn]->getNeuron();
    }else{
        if(m_synapseList1.size() > syn) tmp=m_synapseList1[syn]->getNeuron();
    }

    return tmp;
}

// pre: syn is not larger than the size of list.
float TkrNeuron::getWeight(position pos, unsigned int syn) const
{
    if(pos == top){
        return m_synapseList0[syn]->getWeight();
    }else{
        return m_synapseList1[syn]->getWeight();
    }
}

// calcWeight()
// This function returns the weight value for the two neurons' connection.
float TkrNeuron::calcWeight(TkrNeuron* neuron2, double lambda, double mu)
{
    // If the neurons are joined at the bottom or at the top the weight is zero.
    if(((this)->getPnt(top) == neuron2->getPnt(top)) || 
        ((this)->getPnt(bottom) == neuron2->getPnt(bottom))) return 0.0;

    float tmp = pow((this)->getDirection().dot(neuron2->getDirection()),lambda);

    // to avoid negative weights, retrun 0 if agnle is > 90 degrees.
    if(tmp < 0.8) return 0.0;  

    tmp /= pow((double) (this)->getLayerDiff(), mu) + 
        pow((double) neuron2->getLayerDiff(), mu);

    return ((float) tmp);
}


void TkrNeuron::addConnection(position pos, TkrNeuron* neuron,double lambda,double mu)
{

    float weight = calcWeight(neuron, lambda, mu);

    Synapse* connection = new Synapse(neuron,weight);

    if(pos == top){
        m_synapseList0.push_back(connection);
    }else{
        m_synapseList1.push_back(connection);
    }

    return;
}


inline bool TkrNeuron::operator==(const TkrNeuron& neuron) const
{
    return (getPnt(top) == neuron.getPnt(top) && 
	    getPnt(bottom) == neuron.getPnt(bottom));
}


// update()
// This function updates a given neurons activity level based on the update
// rule.  Neurons which are incoming or outgoing reinforce the neuron.
// Neurons which are competing detract from the neurons activity.
void TkrNeuron::update(float temp, double gamma, double alpha_up, double alpha_do)
{
  float tmp1 = 0.0;  // for reinforcing from the neurons activity.
  float tmp2 = 0.0;  // for detracting from the neurons activity.
  float tmp3 = 0.0;  // for detracting from the neurons activity.
  
  unsigned int i;
  unsigned int j;
  for(i=0, j=0;i < numSynapse(bottom); i++){
    if(getNextNeuron(bottom, i)->getPnt(top) == getPnt(bottom)) 
      tmp1 +=(getWeight(bottom, i)) * 
	(getNextNeuron(bottom, i)->getActivity());
    else{
      tmp2 += getNextNeuron(bottom, i)->getActivity();
      j++;
    }
  }
  
  for(i = 0;i < numSynapse(top); i++){
    if(getNextNeuron(top, i)->getPnt(bottom) == getPnt(top))
      tmp1 +=(getWeight(top, i)) * 
	(getNextNeuron(top, i)->getActivity());
    else{
      tmp3 += getNextNeuron(top, i)->getActivity();
    }
  }
  
  tmp1 *= (float) gamma;
  tmp2 *= (float) alpha_up;
  tmp3 *= (float) alpha_do;
  
  // update the neuron
  setActivity((0.5*(1+tanh((1/temp)*(tmp1-tmp2-tmp3+(getBias()))))));
  
  return;
}
