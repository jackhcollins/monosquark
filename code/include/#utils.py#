import numpy as np
import rootpy.ROOT as ROOT
ROOT.gSystem.Load("libDelphes")
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')

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


def find_last(startID,branchParticle):
    particle = branchParticle.At(startID)
    next_entry = startID
    
    for i in range(10000):
        if particle.D1 == particle.D2:
            next_entry = particle.D1
            particle = branchParticle.At(next_entry)
        else:
            D1 = branchParticle.At(particle.D1)
            D2 = branchParticle.At(particle.D2)
            if D1.PID == particle.PID:
                next_entry = particle.D1
                particle = D1
            elif D2.PID == particle.PID:
                next_entry = particle.D2
                particle = D2
            else:
                return next_entry
    raise

def get_process_ID_VV3body(branchParticle):

    quarks = np.array([-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,21])
    bottom = np.array([5,-5])
    leptons = np.array([11,13,-11,-13])
    neutrinos = np.array([12,14,16,-12,-14,-16])
    tau = np.array([15,-15])
    
    inos = np.append(np.arange(1000023,1000100),np.arange(-1000100,-1000022))
    LSPs = np.array([1000022,-1000022])
    squarks = np.concatenate((np.arange(1000001,1000006),
                              np.arange(2000001,2000006),
                              np.arange(-1000005,-1000000),
                              np.arange(-2000005,-2000000)))

    ino = None
    LSP = None
    current_particle = None
    nW = 0
    nZ = 0
    nh = 0
    nl = 0
    nq = 0
    nb = 0
    ntau = 0
    nnu = 0
    i_lastwino=0
    
    # Find first ino entry
    for i, particle in enumerate(branchParticle):
        if np.any(particle.PID == inos):
            ino = particle
            ino_ID = i
            break
        
    #Find last ino entry
    i_lastwino = find_last(ino_ID,branchParticle)
    ino = branchParticle.At(i_lastwino)
    
    #Find three daughters of ino
    daughters = np.array([None, None])
    daughter_IDs = [None,None]
    current_daughter = 0
    for i in range(i_lastwino+1,branchParticle.GetEntriesFast()):
        part = branchParticle.At(i)
        if part.M1 == i_lastwino:
            if(np.any(part.PID == LSPs)):
                LSP = part
            else:
                daughters[current_daughter] = part
                daughter_IDs[current_daughter] = i
                current_daughter += 1
            if current_daughter > 1 and LSP is not None:
                break
            
    if abs(daughters[0].PID) == 24:
        nW += 1
    if abs(daughters[1].PID) == 24:
        nW += 1
    if daughters[0].PID == 25:
        nh += 1
    if daughters[1].PID == 25:
        nh += 1
    if daughters[0].PID == 23:
        nZ += 1
    if daughters[1].PID == 23:
        nZ += 1
    
    for i_daughter, daughter in enumerate(daughters):
        for i in range(10000):
            daughter_ID = find_last(daughter_IDs[i_daughter],branchParticle)
            daughter = branchParticle.At(daughter_ID)
            
            p12 = [branchParticle.At(daughter.D1),branchParticle.At(daughter.D2)]
            for particle in p12:
                if np.any(particle.PID == quarks):
                    nq += 1
                if np.any(particle.PID == bottom):
                    nb += 1
                if np.any(particle.PID == leptons):
                    nl += 1
                if np.any(particle.PID == tau):
                    ntau += 1
                if np.any(particle.PID == neutrinos):
                    nnu += 1
            break
        
    return np.array([nW,nZ,nh,nq,nl,nnu,ntau,nb])
