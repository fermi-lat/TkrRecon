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
#include "src/PatRec/NeuralNet/TkrNeuralNet.h"

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
float TkrNeuron::getWieght(position pos, unsigned int syn) const
{

    if(pos == top){
        return m_synapseList0[syn]->getWieght();
    }else{
        return m_synapseList1[syn]->getWieght();
    }

}

void TkrNeuron::addConnection(float (TkrNeuralNet::*func)(TkrNeuron*,TkrNeuron*),
                              TkrNeuralNet* const net, position pos, 
                              TkrNeuron* neuron)
{

    float wieght = (net->*func)(this,neuron);

    Synapse* connection = new Synapse(neuron,wieght);

    if(pos == top){
        m_synapseList0.push_back(connection);
    }else{
        m_synapseList1.push_back(connection);
    }

    return;
}

inline bool TkrNeuron::operator==(const TkrNeuron& neuron) const
{
    return(getPnt(top) == neuron.getPnt(top) && 
        getPnt(bottom) == neuron.getPnt(bottom));
}
