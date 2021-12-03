DelphesNtuplizer
=============

This package allows you to produce to flat Ntuples from Delphes PhaseII samples. Some information on the Delphes samples for Snowmass 2021 can be found in [this](https://twiki.cern.ch/twiki/bin/viewauth/CMS/DelphesInstructions#Delphes_instructions_for_Snowmas) TWiki.


Table of contents
=================
  * [Clone](#clone)
  * [Initialisation](#initilisation)
  * [Produce validation Delphes samples](#producing-delphes)
  * [Produce Delphes Flat trees](#producing-flatrees)
  * [Steps to produce trees for ttHbb DL analysis](#producing-flatrees-ttHbbDL)


Clone 
=====

If you do not attempt to contribute to this repository, simply clone it:
```
git clone git@github.com:recotoolsbenchmarks/DelphesNtuplizer.git
```

If you aim at contributing to the repository, you need to fork this repository (via the fork button) and then clone the forked repository:
```
git clone git@github.com:YOURGITUSERNAME/DelphesNtuplizer.git
cd DelphesNtuplizer
git remote add upstream git@github.com:recotoolsbenchmarks/DelphesNtuplizer.git
```
You can then regularly update your fork via:
```
git fetch upstream && git merge upstream/master
```

If you want to submit a new feature to ```recotoolsbenchmarks/DelphesNtuplizer``` you have to do it via pull-request (PR):
So, first commit and push your changes to ```YOURGITUSERNAME/DelphesNtuplizer``` and then make a PR via the github interface. 


Initialisation
==============

This package requires Delphes to be installed, and CMSSW for gcc, ROOT, FWLite and other dependencies:

```
cd DelphesNtuplizer
cmsrel CMSSW_10_0_5
cd CMSSW_10_0_5
cmsenv
cd ..
git clone https://github.com/delphes/delphes.git
cd delphes
./configure
sed -i -e 's/c++0x/c++1y/g' Makefile
make -j 10
cp libDelphes.so ..
```
Make a dummy test to check Delphes runs properly on GEN-SIM samples (once few events have been processed you can stop processing with CTRL+C):

```
./DelphesCMSFWLite cards/gen_card.tcl test_gensim.root /eos/cms/store/relval/CMSSW_11_2_0_pre9/RelValZMM_14/GEN-SIM-RECO/PU_112X_mcRun4_realistic_v4_2026D66PU200-v1/00000/d3b9f895-0c25-4171-9c21-d6cc9b365891.root
```

Produce validation Delphes samples 
===================================

For producing new validation samples you need to set-up the proper environment

```
cd CMSSW_10_0_5
cmsenv
cd ../delphes
```

The Delphes card stored in this repository is configured so that all needed information for validation is available in the Delphes output (gen particles, PF candidates etc ...)
As a result, the output size increases substantially compared to normal Delphes samples. 

To produce Delphes validation samples run this command (by changing the appropiate input GEN-SIM file of interest): 

```
./DelphesCMSFWLite ../cards/CMS_PhaseII_200PU_Snowmass2021_v0.tcl delphes.root /eos/cms/store/relval/CMSSW_11_2_0_pre9/RelValZMM_14/GEN-SIM-RECO/PU_112X_mcRun4_realistic_v4_2026D66PU200-v1/00000/d3b9f895-0c25-4171-9c21-d6cc9b365891.root
```

Produce Delphes flat trees
==========================

Set up the proper environment:

```
cd CMSSW_10_0_5
cmsenv
cd ..
```

The following command will produce a flat Ntuple, with 10 events.

``` 
python bin/Ntuplizer.py -i delphes/delphes.root -o flat_tree.root -n 10
```

Steps to produce trees for ttHbb DL analysis
============================================

1. Start with the **Initialisation** step

2. Fot the ttHbb DL projection studies, the `/TTJets_DiLept_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/PhaseIISummer17wmLHEGENOnly-93X_upgrade2023_realistic_v5_ext1-v3/GEN` dataset is used. Files of this dataset are
listed in **files.txt**. To set-up the proper environment and run a local test do : 

   
```
cd CMSSW_10_0_5
cmsenv
cd ../delphes
./DelphesCMSFWLite ../cards/CMS_PhaseII_200PU_Snowmass2021_v0.tcl delphes.root root://xrootd-cms.infn.it//store/mc/PhaseIISummer17wmLHEGENOnly/TTJets_DiLept_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/GEN/93X_upgrade2023_realistic_v5_ext1-v3/10000/007DE3C2-C653-E811-8EA1-0242AC1C0500.root

```


One can run for multiple files using HTCondor. List of input files and directory to save the output files as well as local paths have to be set within the script **createManyJobs.sh**. File list is stored in **files.txt**. To submit jobs, run:

```
voms-proxy-init -voms cms -rfc -out ${HOME}/.x509up_${UID} -valid 192:00
export X509_USER_PROXY=${HOME}/.x509up_${UID}
source createManyJobs.sh
```


The macro **bin/Ntuplizer_ttHbb.py** has been adapted to produce flat trees relevant to the phase-space of the ttHbb DL analysis.

Set up the proper environment:

```
cd CMSSW_10_0_5
cmsenv
cd ..
```

The following command will produce a flat Ntuple, with 10 events.
   
```
python bin/Ntuplizer_ttHbb.py -i delphes.root -o delphes/flat_tree.root -n 10
```

