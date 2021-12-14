#!/usr/bin/env python
import sys
import ROOT
from collections import OrderedDict
from ROOT import TLorentzVector
from array import array
import numpy as np
import argparse
from math import sqrt,pow

class TreeProducer:
    def __init__(self, debug):

         # flat tree branches
         self.debug = debug

         self.t = ROOT.TTree( "mytree","TestTree" )
         self.maxn = 9999

         # declare arrays
         self.evt_size = array( 'i', [ 0 ] )

         self.vtx_size         = array( 'i', [ 0 ] )
         self.vtx_pt2          = array( 'f', self.maxn*[0. ] )
         self.vtx_x            = array( 'f', self.maxn*[0. ])  
         self.vtx_y            = array( 'f', self.maxn*[0. ])  
         self.vtx_z            = array( 'f', self.maxn*[0. ])  
         
         ## put dummy value for now
         self.true_int          = array( 'i', [ -1 ] )

         ## main MC weight
         self.genweight         = array( 'f', [ 0 ] )
         
         ## variation event weights (from "LHEEventProduct")
         self.lheweight_size   = array( 'i', [ 0 ] )
         self.lheweight_val    = array( 'f', self.maxn*[ 0 ] )

         self.elec_size        = array( 'i', [ 0 ] )
         self.elec_pt          = array( 'f', self.maxn*[ 0. ] )
         self.elec_eta         = array( 'f', self.maxn*[ 0. ] )
         self.elec_phi         = array( 'f', self.maxn*[ 0. ] )
         self.elec_mass        = array( 'f', self.maxn*[ 0. ] )
         self.elec_charge      = array( 'i', self.maxn*[ 0  ] )
         self.elec_idvar       = array( 'f', self.maxn*[ 0. ] )
         self.elec_reliso      = array( 'f', self.maxn*[ 0. ] )
         self.elec_idpass      = array( 'i', self.maxn*[ 0  ] )
         self.elec_isopass     = array( 'i', self.maxn*[ 0  ] )

         self.muon_size        = array( 'i', [ 0 ] )
         self.muon_pt          = array( 'f', self.maxn*[ 0. ] )
         self.muon_eta         = array( 'f', self.maxn*[ 0. ] )
         self.muon_phi         = array( 'f', self.maxn*[ 0. ] )
         self.muon_mass        = array( 'f', self.maxn*[ 0. ] )
         self.muon_charge      = array( 'i', self.maxn*[ 0  ] )
         self.muon_idvar       = array( 'f', self.maxn*[ 0. ] )
         self.muon_reliso      = array( 'f', self.maxn*[ 0. ] )
         self.muon_idpass      = array( 'i', self.maxn*[ 0  ] )
         self.muon_isopass     = array( 'i', self.maxn*[ 0  ] )

         self.jetpuppi_size         = array( 'i', [ 0 ] )
         self.jetpuppi_pt           = array( 'f', self.maxn*[ 0. ] )
         self.jetpuppi_eta          = array( 'f', self.maxn*[ 0. ] )
         self.jetpuppi_phi          = array( 'f', self.maxn*[ 0. ] )
         self.jetpuppi_mass         = array( 'f', self.maxn*[ 0. ] )
         self.jetpuppi_idpass       = array( 'i', self.maxn*[ 0  ] )
         self.jetpuppi_DeepJET      = array( 'f', self.maxn*[ 0. ] )
         self.jetpuppi_btag         = array( 'i', self.maxn*[ 0  ] )
         self.jetpuppi_Nmedium         = array( 'i', [ 0 ] )

         self.metpuppi_size         = array( 'i', [ 0 ] )
         self.metpuppi_pt           = array( 'f', self.maxn*[ 0. ] )
         self.metpuppi_phi          = array( 'f', self.maxn*[ 0. ] )

         self.metpf_size         = array( 'i', [ 0 ] )
         self.metpf_pt           = array( 'f', self.maxn*[ 0. ] )
         self.metpf_phi          = array( 'f', self.maxn*[ 0. ] )

         # declare tree branches
         self.t.Branch( "evt_size",self.evt_size, "evt_size/I")

         self.t.Branch( "vtx_size",self.vtx_size, "vtx_size/I")
         self.t.Branch( "trueInteractions",self.true_int, "trueInteractions/I")
         self.t.Branch( "npuVertices",self.vtx_size, "npuVertices/I")
         self.t.Branch( "vtx_pt2",self.vtx_pt2, "vtx_pt2[vtx_size]/F")
         self.t.Branch( "vtx_x", self.vtx_x, "vtx_x[vtx_size]/F") 
         self.t.Branch( "vtx_y", self.vtx_y, "vtx_y[vtx_size]/F") 
         self.t.Branch( "vtx_z", self.vtx_z, "vtx_z[vtx_size]/F") 

         self.t.Branch( "elec_size",self.elec_size, "elec_size/I")
         self.t.Branch( "elec_pt",self.elec_pt, "elec_pt[elec_size]/F")
         self.t.Branch( "elec_eta",self.elec_eta, "elec_eta[elec_size]/F")
         self.t.Branch( "elec_phi",self.elec_phi, "elec_phi[elec_size]/F")
         self.t.Branch( "elec_mass",self.elec_mass, "elec_mass[elec_size]/F")
         self.t.Branch( "elec_charge",self.elec_charge, "elec_charge[elec_size]/I")
         self.t.Branch( "elec_idvar",self. elec_idvar, "elec_idvar[elec_size]/F")
         self.t.Branch( "elec_reliso",self.elec_reliso, "elec_reliso[elec_size]/F")
         self.t.Branch( "elec_idpass",self.elec_idpass, "elec_idpass[elec_size]/i")
         self.t.Branch( "elec_isopass",self. elec_isopass, "elec_isopass[elec_size]/i")

         self.t.Branch( "muon_size",self.muon_size, "muon_size/I")
         self.t.Branch( "muon_pt",self.muon_pt, "muon_pt[muon_size]/F")
         self.t.Branch( "muon_eta",self.muon_eta, "muon_eta[muon_size]/F")
         self.t.Branch( "muon_phi",self.muon_phi, "muon_phi[muon_size]/F")
         self.t.Branch( "muon_mass",self.muon_mass, "muon_mass[muon_size]/F")
         self.t.Branch( "muon_charge",self.muon_charge, "muon_charge[muon_size]/I")
         self.t.Branch( "muon_idvar",self. muon_idvar, "muon_idvar[muon_size]/F")
         self.t.Branch( "muon_reliso",self.muon_reliso, "muon_reliso[muon_size]/F")
         self.t.Branch( "muon_idpass",self. muon_idpass, "muon_idpass[muon_size]/i")
         self.t.Branch( "muon_isopass",self. muon_isopass, "muon_isopass[muon_size]/i")

         self.t.Branch( "jetpuppi_size",self.jetpuppi_size, "jetpuppi_size/I")
         self.t.Branch( "jetpuppi_pt",self.jetpuppi_pt, "jetpuppi_pt[jetpuppi_size]/F")
         self.t.Branch( "jetpuppi_eta",self.jetpuppi_eta, "jetpuppi_eta[jetpuppi_size]/F")
         self.t.Branch( "jetpuppi_phi",self.jetpuppi_phi, "jetpuppi_phi[jetpuppi_size]/F")
         self.t.Branch( "jetpuppi_mass",self.jetpuppi_mass, "jetpuppi_mass[jetpuppi_size]/F")
         self.t.Branch( "jetpuppi_idpass",self. jetpuppi_idpass, "jetpuppi_idpass[jetpuppi_size]/i")
         #self.t.Branch( "jetpuppi_DeepJET",self.jetpuppi_DeepJET,"jetpuppi_DeepJET[jetpuppi_size]/F")
         self.t.Branch( "jetpuppi_btag",self.jetpuppi_btag,"jetpuppi_btag[jetpuppi_size]/I")
         self.t.Branch( "jetpuppi_Nmedium",self.jetpuppi_Nmedium,"jetpuppi_Nmedium/I")


         self.t.Branch( "metpuppi_size",self.metpuppi_size, "metpuppi_size/I")
         self.t.Branch( "metpuppi_pt",self. metpuppi_pt, "metpuppi_pt[metpuppi_size]/F")
         self.t.Branch( "metpuppi_phi",self.metpuppi_phi, "metpuppi_phi[metpuppi_size]/F")


    #___________________________________________
    def processEvent(self, entry):
        self.evt_size[0] = entry

    #___________________________________________
    def processVertices(self, vertices):
        i = 0
        for item in vertices:
            self.vtx_pt2[i] = item.SumPT2
            self.vtx_x[i] = item.X   #Gamze
            self.vtx_y[i] = item.Y   #Gamze
            self.vtx_z[i] = item.Z   #Gamze
            i += 1
        self.vtx_size[0] = i

    #___________________________________________
    def processWeights(self, event, weights):
        i = 0
        
        self.genweight[0] = event[0].Weight
        for item in weights:
            self.lheweight_val[i] = item.Weight
            i += 1
        self.lheweight_size[0] = i


        '''
        for item in vertices:
            self.vtx_pt2[i] = item.SumPT2
            i += 1
        self.vtx_size[0] = i
        '''
    #___________________________________________
    def processElectrons(self, electrons, electrons_loose, electrons_medium, electrons_tight):
        i = 0

        for item in electrons:

            if(item.PT<15 or abs(item.Eta)>2.4): continue

            self.elec_pt      [i] = item.PT
            self.elec_eta     [i] = item.Eta
            self.elec_phi     [i] = item.Phi
            self.elec_mass    [i] = item.P4().M()
            self.elec_charge  [i] = item.Charge
            self.elec_idvar   [i] = 0.            # DUMMY
            self.elec_reliso  [i] = item.IsolationVar
            self.elec_idpass  [i] = 0             # DUMMY 
            self.elec_isopass [i] = 0             # DUMMY

            if self.elec_reliso[i] < 0.06:
               self.elec_isopass[i] |= 1 << 2

            if self.elec_reliso[i] < 0.2:
               self.elec_isopass[i] |= 1 << 1

            if self.elec_reliso[i] < 0.3:
               self.elec_isopass[i] |= 1 << 0

            # loop over ID collections
            for loose in electrons_loose:
                if dr_match(item,loose,0.005):
                    self.elec_idpass[i] |= 1 << 0

            for medium in electrons_medium:
                if dr_match(item,medium,0.005):
                    self.elec_idpass[i] |= 1 << 1

            for tight in electrons_tight:
                if dr_match(item,tight,0.005):
                    self.elec_idpass[i] |= 1 << 2

            if(self.elec_idpass[i]<4 or self.elec_isopass[i]<4): continue

            i += 1

        '''        
        ## now add fake rates
        for item in electrons_loose:
            if(item.PT<15 or abs(item.Eta)>2.4): continue
            if item.IsolationVar<0.:
                self.elec_pt      [i] = item.PT
                self.elec_eta     [i] = item.Eta
                self.elec_phi     [i] = item.Phi
                self.elec_mass    [i] = 0.
                self.elec_idvar   [i] = 0.            # DUMMY
                self.elec_reliso  [i] = item.IsolationVar
                self.elec_isopass [i] = 7             # DUMMY
                self.elec_idpass  [i] = 0
                self.elec_idpass[i] |= 1 << 0
                if(self.elec_idpass[i]<4 or self.elec_isopass[i]<4): continue
                i+=1

        for item in electrons_medium:
            if(item.PT<15 or abs(item.Eta)>2.4): continue
            if item.IsolationVar<0.:
                self.elec_pt      [i] = item.PT
                self.elec_eta     [i] = item.Eta
                self.elec_phi     [i] = item.Phi
                self.elec_mass    [i] = 0.
                self.elec_idvar   [i] = 0.            # DUMMY
                self.elec_reliso  [i] = item.IsolationVar
                self.elec_isopass [i] = 7             # DUMMY
                self.elec_idpass  [i] = 0
                self.elec_idpass[i] |= 1 << 1
                if(self.elec_idpass[i]<4 or self.elec_isopass[i]<4): continue
                i+=1
                
        for item in electrons_tight:
            if(item.PT<15 or abs(item.Eta)>2.4): continue
            if item.IsolationVar<0.:
                self.elec_pt      [i] = item.PT
                self.elec_eta     [i] = item.Eta
                self.elec_phi     [i] = item.Phi
                self.elec_mass    [i] = 0.
                self.elec_idvar   [i] = 0.            # DUMMY
                self.elec_reliso  [i] = item.IsolationVar
                self.elec_isopass [i] = 7             # DUMMY
                self.elec_idpass  [i] = 0
                self.elec_idpass[i] |= 1 << 2
                if(self.elec_idpass[i]<4 or self.elec_isopass[i]<4): continue
                i+=1
        '''

        self.elec_size[0] = i


    #___________________________________________
    def processMuons(self, muons, muons_loose, muons_medium, muons_tight):
        i = 0
        for item in muons:

            if(item.PT<15 or abs(item.Eta)>2.4): continue

            self.muon_pt      [i] = item.PT
            self.muon_eta     [i] = item.Eta
            self.muon_phi     [i] = item.Phi
            self.muon_mass    [i] = item.P4().M()
            self.muon_charge  [i] = item.Charge
            self.muon_idvar   [i] = 0.            # DUMMY
            self.muon_reliso  [i] = item.IsolationVar
            self.muon_idpass  [i] = 0             # DUMMY 
            self.muon_isopass [i] = 0             # DUMMY

            if self.muon_reliso[i] < 0.15:
               self.muon_isopass[i] |= 1 << 2

            if self.muon_reliso[i] < 0.20:
               self.muon_isopass[i] |= 1 << 1

            if self.muon_reliso[i] < 0.25:
              self.muon_isopass[i] |= 1 << 0

            #if self.muon_reliso[i] < 0.4:
            #   self.muon_isopass[i] |= 1 << 0

            # loop over ID collections
            for loose in muons_loose:
                if dr_match(item,loose,0.005):
                    self.muon_idpass[i] |= 1 << 0
                  
            for medium in muons_medium:
                if dr_match(item,medium,0.005):
                    self.muon_idpass[i] |= 1 << 1
        
            for tight in muons_tight:
                if dr_match(item,tight,0.005):
                    self.muon_idpass[i] |= 1 << 2
    
            if(self.muon_isopass[i]<1 or self.muon_idpass[i]<4): continue

            i += 1

        '''
        ## now add fake rates
        for item in muons_loose:
            if(item.PT<15 or abs(item.Eta)>2.4): continue
            if item.IsolationVar<0.:
                self.muon_pt      [i] = item.PT
                self.muon_eta     [i] = item.Eta
                self.muon_phi     [i] = item.Phi
                self.muon_mass    [i] = 0.
                self.muon_idvar   [i] = 0.            # DUMMY
                self.muon_reliso  [i] = item.IsolationVar
                self.muon_isopass [i] = 7             # DUMMY
                self.muon_idpass  [i] = 0
                self.muon_idpass[i] |= 1 << 0
                if(self.muon_isopass[i]<1 or self.muon_idpass[i]<4): continue
                i+=1

        for item in muons_medium:
            if(item.PT<15 or abs(item.Eta)>2.4): continue
            if item.IsolationVar<0.:
                self.muon_pt      [i] = item.PT
                self.muon_eta     [i] = item.Eta
                self.muon_phi     [i] = item.Phi
                self.muon_mass    [i] = 0.
                self.muon_idvar   [i] = 0.            # DUMMY
                self.muon_reliso  [i] = item.IsolationVar
                self.muon_isopass [i] = 7             # DUMMY
                self.muon_idpass  [i] = 0
                self.muon_idpass[i] |= 1 << 1
                if(self.muon_isopass[i]<1 or self.muon_idpass[i]<4): continue            
                i+=1
                
        for item in muons_tight:
            if(item.PT<15 or abs(item.Eta)>2.4): continue            
            if item.IsolationVar<0.:
                self.muon_pt      [i] = item.PT
                self.muon_eta     [i] = item.Eta
                self.muon_phi     [i] = item.Phi
                self.muon_mass    [i] = 0.
                self.muon_idvar   [i] = 0.            # DUMMY
                self.muon_reliso  [i] = item.IsolationVar
                self.muon_isopass [i] = 7             # DUMMY
                self.muon_idpass  [i] = 0
                self.muon_idpass[i] |= 1 << 2
                if(self.muon_isopass[i]<1 or self.muon_idpass[i]<4): continue
                i+=1
        '''

        self.muon_size[0] = i

    #___________________________________________
    def processPuppiJets(self, jets, jets_loose, jets_tight,el_s,el_pt,el_eta,el_phi,el_M,mu_s,mu_pt,mu_eta,mu_phi,mu_M):

        i = 0
        self.jetpuppi_Nmedium[0] = 0 

        for item in jets:
            jetp4 = item.P4()
            
            if(jetp4.Pt()<30 or abs(jetp4.Eta())>2.4): continue

            self.jetpuppi_pt      [i] = jetp4.Pt()
            self.jetpuppi_eta     [i] = jetp4.Eta()
            self.jetpuppi_phi     [i] = jetp4.Phi()
            self.jetpuppi_mass    [i] = jetp4.M()
            self.jetpuppi_idpass  [i] = 0             # DUMMY 

            # loop over ID collections
            for loose in jets_loose:
                if dr_match(item,loose,0.1):
                    self.jetpuppi_idpass[i] |= 1 << 0
                    self.jetpuppi_idpass[i] |= 1 << 1

            for tight in jets_tight:
                if dr_match(item,tight,0.1):
                    self.jetpuppi_idpass[i] |= 1 << 2
                
            if(self.jetpuppi_idpass[i]<4): continue
    
            ### jet-lepton cleaning  
            flag = False  
            j1,j = TLorentzVector(), TLorentzVector()
            j.SetPtEtaPhiM(jetp4.Pt(),jetp4.Eta(),jetp4.Phi(),jetp4.M())
            for indx in range(0,el_s[0]):
                j1.SetPtEtaPhiM(el_pt[indx],el_eta[indx],el_phi[indx],el_M[indx])
                dr=j.DeltaR(j1)
                if (dr<0.4): 
                    flag = True
            for indx in range(0,mu_s[0]):
                j1.SetPtEtaPhiM(mu_pt[indx],mu_eta[indx],mu_phi[indx],mu_M[indx])
                dr=j.DeltaR(j1)
                if (dr<0.4): 
                    flag = True            
            if(flag): continue
            
            #### BTagging

            self.jetpuppi_DeepJET [i] = 0.  ## some dummy value
            self.jetpuppi_btag[i] = item.BTag

            #if(self.jetpuppi_btag[i]==3 or self.jetpuppi_btag[i]==4 or self.jetpuppi_btag[i]==5 or self.jetpuppi_btag[i]==6 or self.jetpuppi_btag[i]==7):
            if(self.jetpuppi_btag[i]>1):
                self.jetpuppi_Nmedium[0] +=1  
            #print  jetp4.Pt(), jetp4.Eta(), jetp4.Phi(), self.jetpuppi_btag[i], self.jetpuppi_idpass[i]

            i += 1

        self.jetpuppi_size[0] = i

  
    #___________________________________________
    def processPFMissingET(self, met):
        i = 0
        for item in met:

            self.metpf_pt    [i] = item.MET
            self.metpf_phi   [i] = item.Phi
            i += 1
        self.metpf_size  [0] = i

    #___________________________________________
    def processPuppiMissingET(self, met):
        i = 0
        for item in met:

            self.metpuppi_pt    [i] = item.MET
            self.metpuppi_phi   [i] = item.Phi
            i += 1
        self.metpuppi_size  [0] = i


    def fill(self):
        self.t.Fill()

    def write(self):
        self.t.Write()

#_______________________________________________________
def dr_match(p1, p2, drmin):
    dr = p1.P4().DeltaR(p2.P4())
    return dr < drmin


#_____________________________________________________________________________________________________________
def main():

    ROOT.gSystem.Load("libDelphes")
    try:
      ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
      ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
    except:
      pass

    parser = argparse.ArgumentParser()
    parser.add_argument ('-i', '--input', help='input Delphes file',  default='delphes.root')
    parser.add_argument ('-o', '--output', help='output flat tree',  default='tree.root')
    parser.add_argument ('-n', '--nev', help='number of events', type=int, default=-1)
    parser.add_argument ('-d', '--debug', help='debug flag',  action='store_true',  default=False)

    args = parser.parse_args()

    inputFile = args.input
    outputFile = args.output
    nevents = args.nev
    debug = args.debug

    chain = ROOT.TChain("Delphes")
    chain.Add(inputFile)
    

    # Create object of class ExRootTreeReader
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

    ## for now only M for electrons, LT for muons and LT for photons are defined !!
    ## should dervie new parameterisations for other working points

    branchEvent           = treeReader.UseBranch('Event')   
    branchWeight          = treeReader.UseBranch('Weight')   

    branchVertex          = treeReader.UseBranch('Vertex')   
    branchParticle        = treeReader.UseBranch('Particle') 
    branchGenJet          = treeReader.UseBranch('GenJet')   
    branchGenMissingET    = treeReader.UseBranch('GenMissingET')   

    branchPhoton          = treeReader.UseBranch('Photon')
    branchPhotonLoose     = treeReader.UseBranch('PhotonLoose')
    branchPhotonMedium    = treeReader.UseBranch('PhotonMedium')
    branchPhotonTight     = treeReader.UseBranch('PhotonTight')

    branchElectron        = treeReader.UseBranch('Electron')
    branchElectronLoose   = treeReader.UseBranch('ElectronLoose')
    branchElectronMedium  = treeReader.UseBranch('ElectronMedium')
    branchElectronTight   = treeReader.UseBranch('ElectronTight')

    branchMuon            = treeReader.UseBranch('Muon')
    branchMuonLoose       = treeReader.UseBranch('MuonLoose')
    branchMuonMedium      = treeReader.UseBranch('MuonMedium')
    branchMuonTight       = treeReader.UseBranch('MuonTight')

    branchPuppiJet        = treeReader.UseBranch('JetPUPPI')
    branchPuppiJetLoose   = treeReader.UseBranch('JetPUPPILoose')
    branchPuppiJetTight   = treeReader.UseBranch('JetPUPPITight')
    branchFatJet          = treeReader.UseBranch('JetPUPPIAK8') 

    branchPuppiMissingET  = treeReader.UseBranch('PuppiMissingET')

    # NEED these branches to access jet constituents
    branchPuppiCandidate  = treeReader.UseBranch('ParticleFlowCandidate')
    branchRho             = treeReader.UseBranch('Rho')

    treeProducer = TreeProducer(debug)

    if nevents > 0:
        numberOfEntries = nevents

    ################ Start event loop #######################
    for entry in range(0, numberOfEntries):

        # Load selected branches with data from specified event
        treeReader.ReadEntry(entry)

        if (entry+1)%1000 == 0:
            print ' ... processed {} events ...'.format(entry+1)

        treeProducer.processEvent(entry)
        treeProducer.processWeights(branchEvent, branchWeight)
        treeProducer.processVertices(branchVertex)
        treeProducer.processElectrons(branchElectron, branchElectronLoose, branchElectronMedium, branchElectronTight)
        treeProducer.processMuons(branchMuon, branchMuonLoose, branchMuonMedium, branchMuonTight)
        treeProducer.processPuppiJets(branchPuppiJet, branchPuppiJetLoose, branchPuppiJetTight,treeProducer.elec_size,treeProducer.elec_pt,treeProducer.elec_eta,treeProducer.elec_phi,treeProducer.elec_mass,treeProducer.muon_size,treeProducer.muon_pt,treeProducer.muon_eta,treeProducer.muon_phi,treeProducer.muon_mass)
        treeProducer.processPuppiMissingET(branchPuppiMissingET)


        # At least one good vertex
        if(treeProducer.true_int <1): continue

        # At least 4 jets
        if(treeProducer.jetpuppi_size[0]<4): continue

        ## ee/em/mm events - Leading lepton pT>25
        ee, em, mm = False, False, False

        if(treeProducer.elec_size[0]==2 and treeProducer.muon_size[0]==0):
            ee=True
        if(treeProducer.elec_size[0]==1 and treeProducer.muon_size[0]==1):
            em=True
        if(treeProducer.elec_size[0]==0 and treeProducer.muon_size[0]==2):
            mm=True

        if(not(ee or em or mm)): continue
        if(treeProducer.elec_pt[0]<25. and ee): continue
        if(treeProducer.elec_pt[0]<25. and treeProducer.muon_pt[0]<25. and em): continue
        if(treeProducer.muon_pt[0]<25. and mm): continue

        # OS
        if(ee and treeProducer.elec_charge[0]*treeProducer.elec_charge[1]>0): continue
        if(em and treeProducer.elec_charge[0]*treeProducer.muon_charge[0]>0): continue
        if(mm and treeProducer.muon_charge[0]*treeProducer.muon_charge[1]>0): continue

        # mll cuts
        v1,v2,v = TLorentzVector(), TLorentzVector(), TLorentzVector()
        if(ee):
            v1.SetPtEtaPhiM(treeProducer.elec_pt[0],treeProducer.elec_eta[0],treeProducer.elec_phi[0],treeProducer.elec_mass[0])
            v2.SetPtEtaPhiM(treeProducer.elec_pt[1],treeProducer.elec_eta[1],treeProducer.elec_phi[1],treeProducer.elec_mass[1])
        if(em):
            v1.SetPtEtaPhiM(treeProducer.elec_pt[0],treeProducer.elec_eta[0],treeProducer.elec_phi[0],treeProducer.elec_mass[0])
            v2.SetPtEtaPhiM(treeProducer.muon_pt[0],treeProducer.muon_eta[0],treeProducer.muon_phi[0],treeProducer.muon_mass[0])
        if(mm):
            v1.SetPtEtaPhiM(treeProducer.muon_pt[0],treeProducer.muon_eta[0],treeProducer.muon_phi[0],treeProducer.muon_mass[0])
            v2.SetPtEtaPhiM(treeProducer.muon_pt[1],treeProducer.muon_eta[1],treeProducer.muon_phi[1],treeProducer.muon_mass[1])

        v=v1+v2
        if(v.M()<20): continue
        if((v.M()>76 and v.M()<106) and (ee or mm)): continue
        if(treeProducer.metpuppi_pt[0]<40. and (ee or mm)): continue

        ## fill tree 
        treeProducer.fill()
    

    out_root = ROOT.TFile(outputFile,"RECREATE")
    out_root.mkdir("myana")
    out_root.cd("myana")
    treeProducer.write()
 

#_______________________________________________________________________________________
if __name__ == "__main__":
    main()

