import sys
import rootpy.ROOT as ROOT
import numpy as np
ROOT.gSystem.Load("libDelphes")
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
from include.utils import pTof, pTetaphiIDof, get_process_ID_VV3body

num_particle_IDs = 6
num_jet_IDs = 8

inputFile = sys.argv[1]
outputFile = sys.argv[2]
chain = ROOT.TChain("Delphes")

try:
    chain.Add(inputFile)
except:
    print("Unable to load file", inputFile)
    
treeReader = ROOT.ExRootTreeReader(chain)
numevents = treeReader.GetEntries()

branchParticle = treeReader.UseBranch("Particle")
branchFatJet = treeReader.UseBranch("FatJet10")
branchJet = treeReader.UseBranch("Jet")
branchEFlowTrack = treeReader.UseBranch("EFlowTrack")
branchEFlowECal = treeReader.UseBranch("EFlowPhoton")
branchEFlowHCal = treeReader.UseBranch("EFlowNeutralHadron")
branchMuon = treeReader.UseBranch("Muon")
branchElectron = treeReader.UseBranch("Electron")
branchPhoton = treeReader.UseBranch("Photon")


# Figure out how many good events we have
numpass = 0
for event_i in range(numevents):

    treeReader.ReadEntry(event_i)
    
    numjets = branchFatJet.GetEntriesFast()
    if numjets == 0:
        continue
    jet = branchFatJet.At(0)
    if jet.PT >= 500:
        numpass += 1

print("Number of jets: ", numpass, " from ", numevents, " events.")

# Fill out the numpy array

jets_np = np.zeros((numpass,200,3 + num_particle_IDs))
jet_IDs = np.zeros((numpass,num_jet_IDs))
index = 0
j = 0

try:
    for event_i in range(numevents):
        if event_i%1000 is 0:
            print("Working on event", event_i)
        treeReader.ReadEntry(event_i)
        numjets = branchFatJet.GetEntriesFast()
        if numjets == 0:
            continue
        jet = branchFatJet.At(0)
        
        if jet.PT < 500:
            continue
        
        constituents = jet.Constituents
        particles_np = np.zeros((200,3 + num_particle_IDs))
        
        for j, particle in enumerate(constituents):
            particles_np[j] = pTetaphiIDof(particle)
            order = np.flip(np.argsort(particles_np[:,0]))
            particles_np = particles_np[order]
            
        jets_np[index] = particles_np
        jet_IDs[index] = get_process_ID_VV3body(branchParticle)

        index += 1
except:
    print("\n\nFailed to identify particle.")
    print("Event", event_i)
    print("Particle", j)

try:
    np.savez(outputFile,
             constituents=jets_np,
             jet_ID=jet_IDs)
except:
    print("Unable to save", outputFile)
