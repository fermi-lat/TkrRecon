/**
 * @class NeuralNetFindTrackTool
 *
 * @brief Implements a Gaudi Tool for find track candidates. This particular tool uses the 
 *        Neural Net method, described in detail in TkrNeuralNet.h. 
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/NeuralNet/NeuralNetFindTrackTool.h,v 1.5 2003/07/18 21:26:59 lsrea Exp $
 */

#ifndef NEURALNETFINDTRACKTOOL_H
#define NEURALNETFINDTRACKTOOL_H

#include "GaudiKernel/DataSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "src/PatRec/NeuralNet/TkrNeuron.h"
#include "src/PatRec/PatRecBaseTool.h"

class NeuralNetFindTrackTool : public PatRecBaseTool
{
 public:

  /// Standard Gaudi Tool interface constructor
  NeuralNetFindTrackTool(const std::string& type, 
			 const std::string& name, 
			 const IInterface* parent);

  virtual ~NeuralNetFindTrackTool() {}
  
  /// initialize method allows to load the properties
  StatusCode initialize();

  /// @brief Method to find candidate tracks. Will retrieve the necessary information from
  ///        the TDS, including calorimeter energy, and then use TkrNeuralNet to find all
  ///        possible track candidates. The resulting track candidate collection is then   
  ///        stored in the TDS for the next stage.
  StatusCode findTracks();


  /// Instantiation and fake fit of the TkrPatCand candidate tracks.
  void buildCand( Event::TkrPatCandCol&, 
		  const TkrNeuronList&, 
		  Event::TkrClusterCol* );
  

 private:
  /// The properties to be passed to TkrNeuralNet
  double m_MaxLayerDiff;
  double m_MaxPitch;
  double m_Lambda;
  double m_Mu;
  double m_AlphaUP;
  double m_AlphaDOWN;
  double m_Bias;
  double m_Gamma;
  double m_temperature;
};

#endif
