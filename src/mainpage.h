// Mainpage for doxygen

/** @mainpage TkrRecon
 *
 * @section introduction Introduction
 *
 * TkrRecon contains the code that reconstructs tracks from
 * hit strips in the TKR.
 * 
 * It's organized a a series of algorithms that act
 * successively, taking input from, and sending output to
 * the Transient Data Store.
 * 
 * The major driver algorithms are written to support multiple 
 * implementations. 
 * For example, TkrFindAlg currently has three implementations
 * in various stages of development:
 * ComboPR, LinkAndTreePR and NeuralNetPR.
 *
 * A diagram of the program flow can be seen
 * <A href= "../images/Tkr_diagram.gif>here</A>.
 * 
 * @section algorithms Major algorithms
 * 
 * The algorithms are:
 *
 * TkrClusterAlg: Groups hits into clusters. Takes as input 
 * the digitized hits in TkrDigi, and fills TkrClusters
 *
 * TkrFindAlg: Finds possible tracks. Takes the list of clusters as input 
 * and fills TkrCandidates. 
 *
 * TkrReconAlg: (bad name! It's called TkrFitAlg in the diagram.) 
 * Fits the candidate tracks and fills TkrTracks.
 *
 * TkrVertexAlg: Tries to find vertices among the TkrTracks, 
 * to produce TkrVertices. A simple implementation based on 
 * the original Atwood/Hernando code as been installed.
 *
 * @section services Services
 * 
 * TkrGeometrySvc: an interface to GlastDetSvc, it supplies the geometry constants. 
 *
 * TkrBadStripsSvc: manages the list of bad strips.
 * This service currently interacts with TkrClusterAlg, but 
 * ultimately will also be used by other algorithms.
 *
 * TkrInitSvc: allows for loading of code required by the
 * various specific implementations of the algorithms.
 * This service may go away if we can teach
 * the algorithms about sub-algorithms.
 *
 * @authors Tracy Usher, Bill Atwood, Leon Rochester
 * <hr>
 * @section requirements CMT Requirements
 * @verbinclude requirements
 * <hr>
 * @section notes Release notes
 * release.notes
 * <hr>
 * @todo
 * 
 */
