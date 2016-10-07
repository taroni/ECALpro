#Do not modify these
nEventsPerJob      = '-1'
outputFile         = 'EcalZNtp'           # without .root suffix
calibMapName       = 'calibMap.root'

#PATH
eosPath = '/store/user/crovelli'
#
#adding following variables to use commands like "eos ls" and "eos ls -l" commands instead of cmsLs.
#See also here for more details --> https://twiki.cern.ch/twiki/bin/view/CMSPublic/CERNStorageTools 
myeoscmd = '/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select '  #this call directly the eos command (note that eos is an alias, see link above)
myeosls = myeoscmd + 'ls '  #to avoid use of cmsLs that is deprecated since January 2016   
myeoslsl = myeosls + '-l '
myeosmkdir = myeoscmd + 'mkdir '
myeosstage = 'cmsStage -f '
# I called it myeosstage instead of myeoscp to remember that it substitutes cmsStage command
# as a convention, when adding commands like: command = myeoscmd + "some_option ", just leave a space AFTER the some_option, not before
# note that code used cmsStage -f, but eos cp doesn't support -f option
# also, code will copy *.root files from /tmp/ (where they are initially created) to eosPath, but eosPath must be preceeded by "root://eoscms/eos/cms" to have eos cp
# work as expected. So the destination will be root://eoscms/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/... . For this reason, we define here
#myPrefixToEosPath = 'root://eoscms//eos/cms'
myPrefixToEosPath = ''
# will modify calibJobHandler.py with this prefix to destination

#CRAB
isCRAB           = False               # If not is batch
CRAB_Data_Path   = '/SinglePion_FlatPt-1To15_AsymptNoPU/emanuele-SinglePion_FlatPt-1To15_AsymptNoPU-9709e5e865f17288f5a53621cf8e9935/USER'
CRAB_CopyCert    = '/afs/cern.ch/user/l/lpernie/private/x509up_u12147'
storageSite      = "T2_CH_CERN"
unitsPerJob = 10   #DBS File per Job
isOtherT2        = False
if(isCRAB):
   eosPath = '/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/' #For reason of space is better the group area
   if(isOtherT2):
       eosPath = '/pnfs/roma1.infn.it/data/cms/store/user/mciprian/piZero2016/'
       voGroup     = "itcms"
       storageSite = "T2_IT_Rome"
       outLFN      = "/store/user/mciprian/piZero2016/"

#InputList and Folder name
inputlist_n      = 'InputList/fileZskim_BtoF_purified_271036-280385_NoL1T.list'
dirname          = 'TestZ_debug'
Silent           = False                 # True->Fill modules is silent; False->Fill modules has a standard output

#TAG, QUEUE and ITERS
NameTag          = 'TestZ_'                   # Tag to the names to avoid overlap
queueForDaemon   = 'cmscaf1nw'          # Option suggested: 2nw/2nd, 1nw/1nd, cmscaf1nw/cmscaf1nd... even cmscaf2nw
queue            = 'cmscaf1nd'
nIterations      = 14
SubmitFurtherIterationsFromExisting = False
startingCalibMap = '' # used  only if SubmitFurtherIterationsFromExisting is True
if (SubmitFurtherIterationsFromExisting):  # choose path of the calibMap you want to start from
   startingCalibMap = "/store/user/crovelli/TestZ_debug/iter_0/TestZ_calibMap.root"
#N files
ijobmax          = 5                     # 5 number of files per job
nHadd            = 35                    # 35 number of files per hadd
fastHadd         = True                  # From 7_4_X we can use this faster mathod. But files have to be copied on /tmp/ to be converted in .db
if( isCRAB and isOtherT2 ):
   fastHadd      = False                 # No fastHadd on a different T2
nFit             = 2000                  # number of fits done in parallel
Barrel_or_Endcap = 'ALL_PLEASE'          # Option: 'ONLY_BARREL','ONLY_ENDCAP','ALL_PLEASE'

# electron selection
minInvMass = '70'
maxInvMass = '110'
oppositeCharge = 'False'
elePtCut = '20'
eleEtaCut = '2.5'
maxDReleSc = '0.15'
ZCalib_InvMass = 'SCTRMass'
electronSelection = '-1'
useMassInsteadOfEpsilon = 'True'

# data / json / GT
isMC               = False
json_file          = 'Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' if isMC==False else ''            #/afs/cern.ch/cms/CAF/CMSALCA/ALCA_ECALCALIB/json_ecalonly/
globaltag          = '80X_dataRun2_Prompt_v8' if isMC==False else '80X_mcRun2_asymptotic_v5' #old is GR_P_V56
