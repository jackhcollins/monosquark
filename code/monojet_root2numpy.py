import sys
import rootpy.ROOT as ROOT
import numpy as np
ROOT.gSystem.Load("libDelphes")
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
from include.utils import pTof, pTetaphiIDof, get_process_ID_VV3body
from include.utils import get_process_ID_iterate, get_first_ino_ID, get_first_squark_ID, get_squark_ino_ID, get_squark_q_ID     # For ino-jet

num_particle_IDs = 6
num_jet_IDs = 8
debug = 0

inputFile = sys.argv[1]
outputFile = sys.argv[2]
R = "12"
event_type = ""
if len(sys.argv) > 3:
    R = sys.argv[3]
if len(sys.argv) > 4:
    event_type = sys.argv[4]
if len(sys.argv) > 5:
    debug = int(sys.argv[5])
    print("Debug mode is",debug)

chain = ROOT.TChain("Delphes")

try:
    chain.Add(inputFile)
except:
    print("Unable to load file", inputFile, flush=True)
    
treeReader = ROOT.ExRootTreeReader(chain)
numevents = treeReader.GetEntries()

branchParticle = treeReader.UseBranch("Particle")
branchFatJet = treeReader.UseBranch("FatJet" + R)
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

print("Number of jets: ", numpass, " from ", numevents, " events.", flush=True)

# Fill out the numpy array


if event_type == "VV3body":
    print("Event type is VV3body")
    jet_IDs = np.zeros((numpass,num_jet_IDs))
elif event_type == "Zq":
    print("Event type is Zq")
    jet_IDs = np.ones(numpass)
elif event_type == "Zg":
    print("Event type is Zg")
    jet_IDs = np.zeros(numpass)
elif event_type == "ino":
    print("Event type is ino")
    jet_IDs = np.array(np.zeros(numpass),dtype=[('nW','int32'),
                             ('nZ','int32'),
                             ('nh','int32'),
                             ('nq','int32'),
                             ('nb','int32'),
                             ('nnu','int32'),
                             ('nl','int32'),
                             ('ntau','int32')])
elif event_type == "squark-wino":
    print("Event type is squark-wino")
    jet_IDs = np.array(np.zeros(numpass),dtype=[('nW','int32'),
                             ('nZ','int32'),
                             ('nh','int32'),
                             ('nq','int32'),
                             ('nb','int32'),
                             ('nnu','int32'),
                             ('nl','int32'),
                             ('ntau','int32'),
                             ('ass-nW','int32'),
                             ('ass-nZ','int32'),
                             ('ass-nh','int32'),
                             ('ass-nq','int32'),
                             ('ass-nb','int32'),
                             ('ass-nnu','int32'),
                             ('ass-nl','int32'),
                             ('ass-ntau','int32'),
                             ('ass-wino-pT','f8'),
                             ('ass-wino-eta','f8'),
                             ('ass-wino-phi','f8'),
                             ('boosted-wino-pT','f8'),
                             ('boosted-wino-eta','f8'),
                             ('boosted-wino-phi','f8'),
                             ('q-pT','f8'),
                             ('q-eta','f8'),
                             ('q-phi','f8')])  
else:
    print("Event type is unknown")
    jet_IDs = np.zeros(numpass)
    
jets_np = np.zeros((numpass,200,3 + num_particle_IDs))
    
index = 0
j = 0

try:
    for event_i in range(numevents):
        if event_i%1000 is 0:
            print("Working on event", event_i, flush=True)
        treeReader.ReadEntry(event_i)
        numjets = branchFatJet.GetEntriesFast()
        if numjets == 0:
            continue
        jet = branchFatJet.At(0)
        
        if jet.PT < 500:
            continue
        
        constituents = jet.Constituents
        particles_np = np.zeros((constituents.GetEntriesFast(),3 + num_particle_IDs))
        
        for j, particle in enumerate(constituents):
            if j > 199:
                break
            particles_np[j] = pTetaphiIDof(particle)
            order = np.flip(np.argsort(particles_np[:,0]))
            particles_np = particles_np[order]

        particles_np_padded = np.zeros((200,3 + num_particle_IDs))
        if len(particles_np) >= 200:
            particles_np_padded = particles_np[:200]
        else:
            particles_np_padded[:len(particles_np)] = particles_np
            
        jets_np[index] = particles_np_padded

        if event_type == "VV3body":
            jet_IDs[index] = get_process_ID_VV3body(branchParticle)
            
        if event_type == "ino":
            jet_IDs[index] = get_process_ID_iterate(get_first_ino_ID(branchParticle),jet_IDs[index],branchParticle)
            
        if event_type == "squark-wino":
            if(debug > 1):
                print("Getting ass_wino_idx")
            ass_wino_idx = get_first_ino_ID(branchParticle)
            if(debug > 1):
                print("Getting squark_idx")
            squark_idx = get_first_squark_ID(branchParticle, debug=debug)
            if(debug > 1):
                print("squark_idx =",squark_idx)
                print("squark PID =",branchParticle.At(squark_idx).PID)
            if(debug > 1):
                print("Getting boosted_wino_idx")
            boosted_wino_idx = get_squark_ino_ID(branchParticle, squark_ID = squark_idx,debug=debug)
            if(debug > 1):
                print("Getting squark_q_idx")
            squark_q_idx = get_squark_q_ID(branchParticle, squark_ID = squark_idx,debug=debug)
            if(debug > 1):
                print("squark_q_idx =",squark_idx)
                print("q PID =",branchParticle.At(squark_q_idx).PID)
            
            jet_IDs['ass-wino-pT'][index] = branchParticle.At(ass_wino_idx).PT
            jet_IDs['ass-wino-eta'][index] = branchParticle.At(ass_wino_idx).Eta
            jet_IDs['ass-wino-phi'][index] = branchParticle.At(ass_wino_idx).Phi
            
            jet_IDs['boosted-wino-pT'][index] = branchParticle.At(boosted_wino_idx).PT
            jet_IDs['boosted-wino-eta'][index] = branchParticle.At(boosted_wino_idx).Eta
            jet_IDs['boosted-wino-phi'][index] = branchParticle.At(boosted_wino_idx).Phi
            
            jet_IDs['q-pT'][index] = branchParticle.At(squark_q_idx).PT
            jet_IDs['q-eta'][index] = branchParticle.At(squark_q_idx).Eta
            jet_IDs['q-phi'][index] = branchParticle.At(squark_q_idx).Phi
            
            empty_event_ID = np.array(np.zeros(1),dtype=[('nW','int32'),
                             ('nZ','int32'),
                             ('nh','int32'),
                             ('nq','int32'),
                             ('nb','int32'),
                             ('nnu','int32'),
                             ('nl','int32'),
                             ('ntau','int32')])
            
            if (debug > 0):
                print("Getting boosted wino process ID") 
                print(boosted_wino_idx)
            boosted_wino_ID = get_process_ID_iterate(boosted_wino_idx,empty_event_ID,branchParticle,debug=debug)
            
            if (debug > 0):
                print("Getting ass-wino process ID")
            ass_wino_ID = get_process_ID_iterate(ass_wino_idx,empty_event_ID,branchParticle,debug=debug)

            jet_IDs['nW'][index] = boosted_wino_ID['nW'][0]
            jet_IDs['nZ'][index] = boosted_wino_ID['nZ'][0]
            jet_IDs['nh'][index] = boosted_wino_ID['nh'][0]
            jet_IDs['nq'][index] = boosted_wino_ID['nq'][0]
            jet_IDs['nb'][index] = boosted_wino_ID['nb'][0]
            jet_IDs['nnu'][index] = boosted_wino_ID['nnu'][0]
            jet_IDs['nl'][index] = boosted_wino_ID['nl'][0]
            jet_IDs['ntau'][index] = boosted_wino_ID['ntau'][0]
            
            jet_IDs['ass-nW'][index] = ass_wino_ID['nW'][0]
            jet_IDs['ass-nZ'][index] = ass_wino_ID['nZ'][0]
            jet_IDs['ass-nh'][index] = ass_wino_ID['nh'][0]
            jet_IDs['ass-nq'][index] = ass_wino_ID['nq'][0]
            jet_IDs['ass-nb'][index] = ass_wino_ID['nb'][0]
            jet_IDs['ass-nnu'][index] = ass_wino_ID['nnu'][0]
            jet_IDs['ass-nl'][index] = ass_wino_ID['nl'][0]
            jet_IDs['ass-ntau'][index] = ass_wino_ID['ntau'][0]

        index += 1
except:
    print("\n\nError: exception raised")
    print("Event", event_i)
    print("Particle", j)
    print(particle.ClassName(), flush=True)
    
try:
    np.savez(outputFile,
             constituents=jets_np,
             jet_ID=jet_IDs)
except:
    print("Unable to save", outputFile, flush=True)
