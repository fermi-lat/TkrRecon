# -*- python -*-
# $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/SConscript,v 1.35 2011/01/12 00:23:59 lsrea Exp $ 
# Authors: Leon Rochester <lsrea@slac.stanford.edu>, Tracy Usher <usher@slac.stanford.edu>
# Version: TkrRecon-10-23-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('addLinkDeps', package='TkrRecon', toBuild='component')
TkrRecon=libEnv.SharedLibrary('TkrRecon',
                              listFiles(['src/Dll/*.cxx','src/GaudiAlg/*.cxx',
                                         'src/Display/*.cxx','src/Services/*.cxx',
                                         'src/Cluster/*.cxx','src/Filter/*.cxx',
                                         'src/PatRec/*.cxx',
                                         'src/PatRec/LinkAndTree/*.cxx',
                                         'src/PatRec/Combo/*.cxx',
                                         'src/PatRec/NeuralNet/*.cxx',
                                         'src/PatRec/VectorLinks/*.cxx',
                                         'src/PatRec/VectorLinks/StdRelTable/*.cxx',
                                         'src/PatRec/TreeBased/*.cxx',
                                         'src/PatRec/MonteCarlo/*.cxx',
                                         'src/PatRec/KalFitTrack/*.cxx'
                                         'src/PatRec/Utilities/*.cxx',
                                         'src/Track/*.cxx','src/TrackFit/*.cxx',
                                         'src/TrackFit/KalmanFilter/*.cxx',
                                         'src/TrackFit/KalmanFilterUtils/*.cxx',
                                         'src/TrackFit/KalmanFilterFit/*.cxx',
                                         'src/TrackFit/KalmanFilterFit/FitMatrices/*.cxx',
                                         'src/TrackFit/KalmanFilterFit/HitErrors/*.cxx',
                                         'src/TrackFit/KalmanFilterFit/TrackEnergy/*.cxx',
                                         'src/TrackFit/KalFitTrack/*.cxx',
                                         'src/TrackFit/LineFit2D/*.cxx',
                                         'src/TrackFit/LineFit3D/*.cxx',
                                         'src/Utilities/*.cxx','src/Vertex/*.cxx',
                                         'src/Vertex/Combo/*.cxx',
                                         'src/Vertex/DocaVtx/*.cxx']))
progEnv.Tool('TkrReconLib')

test_TkrRecon = progEnv.GaudiProgram('test_TkrRecon',
                                     listFiles(['src/test/*.cxx']),
                                     test = 1, package='TkrRecon')

progEnv.Tool('registerTargets', package='TkrRecon',
             libraryCxts = [[TkrRecon,libEnv]],
             testAppCxts = [[test_TkrRecon, progEnv]],
             includes = listFiles(['TkrRecon/*'], recursive = 1),
             jo = ['src/test/jobOptions.txt'])




