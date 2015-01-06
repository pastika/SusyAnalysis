
#!\/bin\/tcsh
Universe     = vanilla
Notification = Complete
Arguments    =   `pwd` 
Executable = CMSRUNSCRIPT
Should_Transfer_Files = YES
WhenTOTransferOutput  = ON_EXIT
Output = FILENAME.stdout
Error = FILENAME.stderr
Log =  FILENAME.condorlog
notify_user = vasurang@FNAL.GOV
Queue 1


