import sys
import rootpy.ROOT as ROOT
import numpy as np
ROOT.gSystem.Load("libDelphes")
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')

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

def pTof(obj):
    if isinstance(obj,ROOT.Track):
        return obj.PT
    elif isinstance(obj,ROOT.Muon):
        return obj.PT
    elif isinstance(obj,ROOT.Electron):
        return obj.PT
    elif isinstance(obj,ROOT.Photon):
        return obj.PT
    elif isinstance(obj,ROOT.Tower):
        return obj.ET
    else:
        raise

num_IDs = 6
def pTetaphiIDof(obj):
    part_ID = 0
    if isinstance(obj,ROOT.Tower):
        if obj.Ehad > obj.Eem:
            part_ID = 0
        else:
            part_ID = 1
    elif isinstance(obj,ROOT.Track):
        part_ID = 2
    elif isinstance(obj,ROOT.Muon):
        part_ID = 3
    elif isinstance(obj,ROOT.Electron):
        part_ID = 4
    elif isinstance(obj,ROOT.Photon):
        part_ID = 5
        
    part_ID_onehot = np.zeros(6)
    part_ID_onehot[part_ID] = 1
    
    
    try:
        return np.append(np.array([pTof(obj),obj.Eta,obj.Phi]),
                         part_ID_onehot)
    except:
        raise


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

jets_np = np.zeros((numpass,200,3+num_IDs))
index = 0

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
        particles_np = np.zeros((200,3+num_IDs))
        
        for j, particle in enumerate(constituents):
            particles_np[j] = pTetaphiIDof(particle)
            order = np.flip(np.argsort(particles_np[:,0]))
            particles_np = particles_np[order]
            
        jets_np[index] = particles_np
        index += 1
except:
    print("\n\nFailed to identify particle.")
    print("Event", event_i)
    print("Particle", j)

try:
    np.save(outputFile, jets_np)
except:
    print("Unable to save", outputFile)
