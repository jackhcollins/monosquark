import numpy as np
import rootpy.ROOT as ROOT
ROOT.gSystem.Load("libDelphes")
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


quarks = np.array([-6,-5,-4,-3,-2,-1,1,2,3,4,5,6])
bottom = np.array([5,-5])
leptons = np.array([11,13,-11,-13])
neutrinos = np.array([12,14,16,-12,-14,-16])
tau = np.array([15,-15])
photon = np.array([22])

PIDs = {}
PIDs['q'] = np.array([-6,-5,-4,-3,-2,-1,1,2,3,4,5,6])
PIDs['g'] = np.array([21])
PIDs['j'] = np.concatenate((PIDs['q'],PIDs['g']))
PIDs['b'] = np.array([5,-5])
PIDs['l'] = leptons
PIDs['nu'] = neutrinos
PIDs['tau'] = tau
PIDs['W'] = np.array([24,-24])
PIDs['h'] = np.array([25])
PIDs['Z'] = np.array([23])
PIDs['V'] = np.concatenate((PIDs['W'],
                           PIDs['h'],
                           PIDs['Z']))
PIDs['photon'] = np.array([22])


inos = np.append(np.arange(1000023,1000100),np.arange(-1000100,-1000022))
LSPs = np.array([1000022,-1000022])
squarks = np.concatenate((np.arange(1000001,1000006),
                          np.arange(2000001,2000006),
                          np.arange(-1000005,-1000000),
                          np.arange(-2000005,-2000000)))

PIDs['squark'] = squarks
PIDs['LSP'] = LSPs
PIDs['ino'] = inos

PIDs['final'] = np.concatenate((PIDs['j'],
                                PIDs['b'],
                                PIDs['l'],
                                PIDs['nu'],
                                PIDs['tau'],
                                PIDs['LSP'],
                                PIDs['photon']))


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

    
def get_all_daughters(startID,branchParticle,stopnum = -1,debug=0):
    daughters = []
    daughter_IDs = []

    for i in range(startID+1,branchParticle.GetEntriesFast()):
        part = branchParticle.At(i)
        if part.M1 == startID:
            daughters.append(part)
            daughter_IDs.append(i)
        if len(daughters) == stopnum:
            break
         
    return daughter_IDs, daughters


def update_event_ID(event_ID, PID):
    if np.any(PID == PIDs['W']):
        event_ID['nW'] += 1
    if np.any(PID == PIDs['Z']):
        event_ID['nZ'] += 1
    if np.any(PID == PIDs['h']):
        event_ID['nh'] += 1
    if np.any(PID == PIDs['q']):
        event_ID['nq'] += 1
    if np.any(PID == PIDs['b']):
        event_ID['nb'] += 1
    if np.any(PID == PIDs['l']):
        event_ID['nl'] += 1
    if np.any(PID == PIDs['tau']):
        event_ID['ntau'] += 1
    if np.any(PID == PIDs['nu']):
        event_ID['nnu'] += 1
                
    return event_ID


def get_process_ID_iterate(winoID,event_ID,branchParticle,debug = False, twoprong = False):
   
    if np.any(branchParticle.At(winoID).PID == PIDs['LSP']):
        return event_ID

    last_wino_ID = find_last(winoID,branchParticle)
    last_wino = branchParticle.At(last_wino_ID)
    

    
    if last_wino.D1 == last_wino.D2 or last_wino.D2 == -1 or last_wino.D1 == -1:
        print("Invalid daughters")
        raise
       
    if twoprong:
        daughter_IDs = [last_wino.D1,last_wino.D2]
        daughters = [branchParticle.At(daughter_IDs[0]), branchParticle.At(daughter_IDs[1])]
    else:
        if debug > 2:
            print("Finding daughters")
        daughter_IDs, daughters = get_all_daughters(last_wino_ID,branchParticle,debug=debug)
        
        if len(daughter_IDs) is 0:
            print("No daughters found")
            raise
    
    for i, daughter in enumerate(daughters):
        D_PID = daughter.PID
        if debug > 1:
            print(D_PID)
        event_ID = update_event_ID(event_ID,D_PID)
        if not np.any(D_PID == PIDs['final']):
            event_ID = get_process_ID_iterate(daughter_IDs[i],event_ID,branchParticle,debug=debug,twoprong=twoprong)

    return event_ID


def get_first_ino_ID(branchParticle):
    
    for i in range(branchParticle.GetEntriesFast()):
        if np.any(branchParticle.At(i).PID == PIDs["ino"]):
            return i
        
    print("Error: Failed to find a -ino")
    raise
    
    
def get_first_squark_ID(branchParticle, debug = False):
    
    for i in range(branchParticle.GetEntriesFast()):
        if np.any(branchParticle.At(i).PID == PIDs["squark"]):
            return i
        
    print("Error: Failed to find a squark")
    raise
    

def get_squark_ino_ID(branchParticle, squark_ID = None,debug = 0):
    
    if squark_ID is None:
        squark_ID = find_last(get_first_squark_ID,branchParticle)
    else:
        squark_ID = find_last(squark_ID,branchParticle)
    
    for i in range(squark_ID+1,branchParticle.GetEntriesFast()):
        part = branchParticle.At(i)
        if debug > 2:
            print(i, part.PID, part.M1, part.M2)
        if part.M1 == squark_ID and (np.any(part.PID == PIDs['ino']) or np.any(part.PID == PIDs['LSP'])):
            if debug > 0:
                print("squark_ino_ID found:",i)
            return i

    print("Failed to find a squark-ino")
    raise
    
def get_squark_q_ID(branchParticle, squark_ID = None,debug=0):
    
    if squark_ID is None:
        squark_ID = find_last(get_first_squark_ID,branchParticle)
    else:
        squark_ID = find_last(squark_ID,branchParticle)
    
    for i in range(squark_ID+1,branchParticle.GetEntriesFast()):
        part = branchParticle.At(i)
        if debug > 2:
            print(i, part.PID, part.M1, part.M2, np.any(part.PID == PIDs['q']))
        if part.M1 == squark_ID and np.any(part.PID == PIDs['q']):
            return i

    print("Failed to find a squark-q")
    raise
   
    

def get_process_ID_VV3body(branchParticle):

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


def PtEtaPhiMvec_to_pxpypzEvec(PtEtaPhiMvec):
    thetavec = 2*np.arctan(np.exp(-PtEtaPhiMvec[:,1]))
    phivec = PtEtaPhiMvec[:,2]
    pzvec = PtEtaPhiMvec[:,0]/np.tan(thetavec)
    pxvec = PtEtaPhiMvec[:,0]*np.cos(phivec)
    pyvec = PtEtaPhiMvec[:,0]*np.sin(phivec)
    Evec = np.sqrt(PtEtaPhiMvec[:,3]*PtEtaPhiMvec[:,3]+
                  pxvec*pxvec+
                  pyvec*pyvec+
                  pzvec*pzvec)
    return np.stack((pxvec,pyvec,pzvec,Evec),axis=-1)


def PtEtaPhiMeventlist_to_pxpypzEeventlist(PtEtaPhiMvec):
    thetavec = 2*np.arctan(np.exp(-PtEtaPhiMvec[:,:,1]))
    phivec = PtEtaPhiMvec[:,:,2]
    pzvec = PtEtaPhiMvec[:,:,0]/np.tan(thetavec)
    pxvec = PtEtaPhiMvec[:,:,0]*np.cos(phivec)
    pyvec = PtEtaPhiMvec[:,:,0]*np.sin(phivec)
    Evec = np.sqrt(PtEtaPhiMvec[:,:,3]*PtEtaPhiMvec[:,:,3]+
                  pxvec*pxvec+
                  pyvec*pyvec+
                  pzvec*pzvec)
    return np.stack((pxvec,pyvec,pzvec,Evec),axis=-1)


def pxpypzE_to_M(pxpypzE):
    return np.sqrt(pxpypzE[3]*pxpypzE[3] - np.power(np.linalg.norm(pxpypzE[:3]),2))


def PtEtaPhiMvec_to_M(PtEtaPhiMvec):
    jet4vec = np.sum(PtEtaPhiMvec_to_pxpypzEvec(PtEtaPhiMvec),axis=0)
    return pxpypzE_to_M(jet4vec)
    

def event_list_to_PtEtaPhiM_list(event_list):
    numevents, numparticles = event_list.shape[:2]
    return np.append(event_list[:,:,:3],np.zeros((numevents,numparticles,1)),axis=-1)