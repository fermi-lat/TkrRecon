// Mainpage for doxygen

/** @mainpage TkrRecon
 *
 * @section introduction Introduction
 *
 * TkrRecon contains the code that reconstructs tracks from
 * hit strips in the TKR. 
 * 
 * It's organized as a series of algorithms that act
 * successively, taking input from, and sending output to
 * the Transient Data Store. Starting from TkrDigi objects, it generates
 * TkrCluster's, then TkrPatCand's, TkrFitTrack's and finally TkrVertex objects.
 * 
 * The major driver algorithms are written to support multiple 
 * implementations. 
 * For example, TkrFindAlg currently has three implementations
 * in various stages of development:
 * TkrComboPatRec, TkrLinkAndTree and TkrNeuralNet.
 *
 * A diagram of the program flow can be seen
 * <A href= "../images/Tkr_diagram.gif>here</A>.
 * 
 * @section algorithms Major algorithms
 * 
 * The algorithms themselves don't contain the "algorithms." They merely handle the Gaudi manipulations
 * and then pass the work on to other classes, which in turn, operate on the TDS classes. The main algorithms
 * in the chain are:
 *
 * TkrClusterAlg: The work is done by TkrMakeClusters, which groups hits into clusters. Takes as input 
 * the digitized hits in TkrDigi, and fills TkrClusters with TkrCluster objects.
 *
 * TkrFindAlg: The work is done by TkrComboPatRec, TkrLinkAndTree, and TkrNeuralNet. 
 * All these find possible tracks. Takes the list of clusters as input 
 * and fills TkrCandidates with TkrPatCand objects. 
 *
 * TkrReconAlg: (bad name! It's called TkrFitAlg in the diagram.)
 * The candidate tracks are fit, and TkrTracks is filled with TkrFitTrack objects.
 * Here the work is done by one of TkrLinkAndTreeFit, TkrComboFit or TkrNeuralNetFit.
 * All three of these currently pass patrec candidates to the same Kalman Filter, KalFitTrack, 
 * but because each patrec produces slightly different information, we still need
 * three separate classes. We need to fix this.
 *
 * TkrVertexAlg: Finds vertices among the TkrTracks, 
 * to produce TkrVertexCol, the collection of TkrVertex objects.. 
 * The work is done in TkrComboVtxRecon, a simple implementation meant to be a placeholder and
 * minimal baseline.
 *
 * @section access Access to the "interesting" quantities
 *
 * The abstract class TkrRecInfo provides the user of TkrRecon TDS objects with a common
 * interface to the basic information about TkrPatCand, TkrFitTrack and TkrVertex objects.
 * The class will live in GlastEvent, so that other subsystems may access this information
 * without having to depend on TkrRecon.
 *
 * The information includes the track (or vertex) parameters and error matrices 
 * at both ends of the track (or one "end" for the vertex), the energy and quality.
 * 
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
 * @authors Tracy Usher, Leon Rochester (SLAC) <br>
 * Bill Atwood, Brandon Allgood(UC Santa Cruz) <br>
 * Michael Kuss, Johann Cohen-Tanugi (INFN Pisa) <br>
 * Jose Angel Hernando Morata (Santiago de Compostela)
 * <hr>
 * @section requirements CMT Requirements
 * @verbinclude requirements
 * <hr>
 * @section notes Release notes
 * release.notes
 * <hr>
 * @todo
 * Develop mechanism to handle constants required by code. Current plan is to put the constants in an xml file,
 * and to make them properties as well, so that they can be over-written from the jobOptions.
 *
 * Consolidate references to geometry - remove from GFtutor and other places
 *
 * Root out last hardwired numbers
 *
 * Respond to DSTF review
 *
 * Doxygenate remainder of code
 */
