#!/bin/bash
###################################
# define paths
###################################
#filelist="files_2nd-step.txt"
filelist="files_2nd-step_ttH.txt"
curDir="/afs/cern.ch/work/c/ckoraka/Delphes/DelphesNtuplizer"
outputDir="/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ckoraka/Delphes_2nd-Step_ttH"

##################################
# cleaning & creating directories
##################################

if [[ ! -d "$outputDir" ]]
then
    mkdir ${outputDir}
fi

if [[ -d "${curDir}/jobs" ]]
then
   rm -rf  ${curDir}/jobs
fi

if [[ ! -d "${curDir}/jobs" ]]
then
    mkdir ${curDir}/jobs
fi

if [[ -d "${curDir}/jobs/err" ]]
then
   rm -rf  ${curDir}/jobs/err
fi

if [[ ! -d "${curDir}/jobs/err" ]]
then
    mkdir ${curDir}/jobs/err
fi

if [[ -d "${curDir}/jobs/out" ]]
then
   rm -rf  ${curDir}/jobs/out
fi

if [[ ! -d "${curDir}/jobs/out" ]]
then
    mkdir ${curDir}/jobs/out
fi

if [[ -d "${curDir}/jobs/log" ]]
then
   rm -rf  ${curDir}/jobs/log
fi

if [[ ! -d "${curDir}/jobs/log" ]]
then
    mkdir ${curDir}/jobs/log
fi

n=1
while read line;do
# reading each line
input=$line
echo ===== Processing file $n : $input =====

################################
# create Condor job file
################################

echo ==== Creating condor job ======
jobfile=delphes_${n}.sh
condorfile=delphes_${n}.sub
echo  --- creating executable file ${outputDir}/jobs/$jobfile ..
echo  --- creating job ${outputDir}/jobs/$condorfile ..

cat > $jobfile <<@EOI
#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
scp -r /afs/cern.ch/work/c/ckoraka/Delphes/DelphesNtuplizer .
cd DelphesNtuplizer/CMSSW_10_0_5/src
source set_env.sh
cd ../..
python bin/Ntuplizer.py -i ${input} -o output.root 
mv output.root ${outputDir}/delphes_${n}.root
@EOI

echo ---created executable file---

cat > $condorfile <<@EOI
executable              = $jobfile
output                  = out/$jobfile.out
error                   = err/$jobfile.err
log                     = log/$jobfile.log
+JobFlavour             = "workday"
use_x509userproxy       = true
x509userproxy = $X509_USER_PROXY
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
queue
@EOI

############################
#submit job
############################
echo ==== Submitting condor jobs ======
mv $jobfile    ${curDir}/jobs
mv $condorfile ${curDir}/jobs
chmod ugo+rwx  ${curDir}/jobs/$condorfile
chmod ugo+rwx  ${curDir}/jobs/$jobfile
# Both executable & condor submission scripts must be in the same folder / cd to that folder every time if needed for proper submission
cd jobs
condor_submit $condorfile
cd ..
sleep 1s

n=$((n+1))

done < $filelist
