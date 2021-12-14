#!/bin/tcsh
echo
setenv joblist `condor_q -w  | grep ckoraka | awk '{print $11}'`


foreach job ($joblist)
condor_rm  ${job}
echo killing ${job} ..
end
echo
~    
