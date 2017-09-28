#
nEventsPerJob      = '-1'
#outputFile         = 'ZNtpDataFull2016bis'           # without .root suffix
outputFile         = 'ZNtpMcMoriond'                  # without .root suffix
calibMapName       = 'calibMap.root'

#PATH
#eosPath = '/eos/cms/store/group/dpg_ecal/alca_ecalcalib/zeeIc'
eosPath = '/eos/cms/store/user/taroni/alca_ecalcalib/zeeIc'
#
#See also here for more details --> https://twiki.cern.ch/twiki/bin/view/CMSPublic/CERNStorageTools 
myeoscmd = 'eos ' # from July 2017 we can use eos on lxbatch from inside scripts 
myeosls = myeoscmd + 'ls '  #to avoid use of cmsLs that is deprecated since January 2016   
myeoslsl = myeosls + '-l '
myeosmkdir = myeoscmd + 'mkdir '
myeosstage = myeoscmd + 'cp ' 
prefixSourceFile = 'root://cms-xrd-global.cern.ch/'

#CRAB
isCRAB           = False               # If not is batch
CRAB_Data_Path   = '/SinglePion_FlatPt-1To15_AsymptNoPU/emanuele-SinglePion_FlatPt-1To15_AsymptNoPU-9709e5e865f17288f5a53621cf8e9935/USER'
CRAB_CopyCert    = '/afs/cern.ch/user/t/taroni/x509up_u29820'
storageSite      = "T2_CH_CERN"
unitsPerJob = 10   #DBS File per Job
isOtherT2        = False

#InputList and Folder name
#inputlist_n      = 'InputList/fileZskim_BtoH_rereco_purifiedWithFinalJson.list'
#dirname          = 'ZDataFull2016bis'
inputlist_n      = 'InputList/fileDYToLLmoriond.list'
dirname          = 'ZMcMoriond'
Silent           = False                 # True->Fill modules is silent; False->Fill modules has a standard output

#TAG, QUEUE and ITERS
#NameTag          = 'ZDataFull2016bis_'          # Tag to the names to avoid overlap
NameTag          = 'ZMcMoriond_'             # Tag to the names to avoid overlap
queueForDaemon   = 'cmscaf1nw'               # Option suggested: 2nw/2nd, 1nw/1nd, cmscaf1nw/cmscaf1nd... even cmscaf2nw
queue            = 'cmscaf1nd'
nIterations      = 14
SubmitFurtherIterationsFromExisting = False
startingCalibMap = '' # used  only if SubmitFurtherIterationsFromExisting is True
if (SubmitFurtherIterationsFromExisting):  # choose path of the calibMap you want to start from
   startingCalibMap = "/store/user/crovelli/TestZ_miscal/miscalibMaps.root"
#N files
ijobmax          = 15                    # 2 on data; 5 on MC
nHadd            = 35                    # 35 number of files per hadd
fastHadd         = True                  # From 7_4_X we can use this faster mathod. But files have to be copied on /tmp/ to be converted in .db
nFit             = 2000                  # number of fits done in parallel
Barrel_or_Endcap = 'ALL_PLEASE'          # Option: 'ONLY_BARREL','ONLY_ENDCAP','ALL_PLEASE'

# electron selection
minInvMass = '70'
maxInvMass = '110'
oppositeCharge = 'False'
elePtCut = '25'
eleEtaCut = '2.5'
maxDReleSc = '0.15'
ZCalib_InvMass = 'SCTRMass'
puweightfile = 'pileupWeights__GoldenRereco_mcMoriond__69200'
electronSelection = '-1'
useMassInsteadOfEpsilon = 'True'

# data / json / GT
isMC               = True
json_file          = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt' if isMC==False else ''            #/afs/cern.ch/cms/CAF/CMSALCA/ALCA_ECALCALIB/json_ecalonly/
globaltag          = '80X_dataRun2_2016LegacyRepro_v3' if isMC==False else '80X_mcRun2_asymptotic_2016_TrancheIV_v7'
