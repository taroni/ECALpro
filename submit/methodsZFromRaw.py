from parametersZFromRaw import *

def printFillFromRawCfgZ1( outputfile, iteration ):
    outputfile.write('import FWCore.ParameterSet.Config as cms\n')
    outputfile.write("import os, sys, imp, re\n")
    outputfile.write('CMSSW_VERSION=os.getenv("CMSSW_VERSION")\n')
    outputfile.write('process = cms.Process("analyzerFillEpsilonForZ")\n')    
    outputfile.write('\n')  

    outputfile.write("#Import standard configurations\n")
    outputfile.write('process.load("Configuration.StandardSequences.Services_cff")\n')  
    outputfile.write('process.load("FWCore.MessageService.MessageLogger_cfi")\n\n')
    outputfile.write('process.load("Configuration.StandardSequences.GeometryRecoDB_cff")\n')  # chiara: quale?
    outputfile.write('process.load("Configuration.StandardSequences.EndOfProcess_cff")\n')  
    outputfile.write('process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")\n')   
    outputfile.write('process.load("RecoLuminosity.LumiProducer.bunchSpacingProducer_cfi")\n')
    outputfile.write('\n')  

    outputfile.write('process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")\n')  
    outputfile.write("process.GlobalTag.globaltag = '" + globaltag + "'\n")          
    outputfile.write('\n')  

    outputfile.write('### From RAW to RECO\n')
    outputfile.write('process.load("Configuration.StandardSequences.RawToDigi_Data_cff")\n') 
    outputfile.write('process.load("RecoLocalCalo.Configuration.ecalLocalRecoSequence_cff")\n') 
    outputfile.write('\n')

    outputfile.write('### Recalibration Module using previously produced map\n')
    outputfile.write('process.load("CalibCode.FillEpsilonPlot.calibRechitsFromRawForZProducer_cfi")\n')
    outputfile.write("process.ecalRecalRecHit.CurrentIteration = cms.untracked.int32(" + str(iteration) + ")\n")
    if (isCRAB):
        outputfile.write("process.ecalRecalRecHit.calibMapPath = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + NameTag + calibMapName + "')\n")
        outputfile.write("process.ecalRecalRecHit.isCRAB  = cms.untracked.bool(True)\n")
    else:
        if (SubmitFurtherIterationsFromExisting and iteration == 0):
            outputfile.write("process.ecalRecalRecHit.calibMapPath = cms.untracked.string('root://eoscms//eos/cms" + startingCalibMap + "')\n")
        else:
            outputfile.write("process.ecalRecalRecHit.calibMapPath = cms.untracked.string('root://eoscms//eos/cms" + eosPath + "/" + dirname + "/iter_" + str(iteration-1) + "/" + NameTag + calibMapName + "')\n")

    outputfile.write("process.ecalRecalRecHit.Barrel_orEndcap = cms.untracked.string('" + Barrel_or_Endcap + "')\n")
    outputfile.write('\n')

    outputfile.write('### PF rechits\n')
    outputfile.write('process.load("RecoParticleFlow.PFClusterProducer.particleFlowRecHitECAL_cfi")\n') 
    outputfile.write('process.particleFlowRecHitECAL.producers[0].src = cms.InputTag("ecalRecalRecHit","EcalRecalRecHitsEB")\n')
    outputfile.write('process.particleFlowRecHitECAL.producers[1].src = cms.InputTag("ecalRecalRecHit","EcalRecalRecHitsEE")\n')
    outputfile.write('process.load("RecoParticleFlow.PFClusterProducer.particleFlowRecHitPS_cfi")\n') 
    outputfile.write('\n')

    outputfile.write('### PF clusters\n') 
    outputfile.write('process.load("RecoParticleFlow.PFClusterProducer.particleFlowClusterECALUncorrected_cfi")\n')   
    outputfile.write('process.load("RecoParticleFlow.PFClusterProducer.particleFlowClusterECAL_cfi")\n')   
    outputfile.write('process.load("RecoParticleFlow.PFClusterProducer.particleFlowClusterPS_cfi")\n')    
    outputfile.write('process.load("RecoEcal.EgammaClusterProducers.particleFlowSuperClusterECAL_cfi")\n')   
    outputfile.write('\n') 

    outputfile.write('### PF superclusters\n')       
    outputfile.write('process.particleFlowSuperClusterECAL.PFBasicClusterCollectionBarrel = cms.string("recalibParticleFlowBasicClusterECALBarrel")\n')    
    outputfile.write('process.particleFlowSuperClusterECAL.PFSuperClusterCollectionBarrel = cms.string("recalibParticleFlowSuperClusterECALBarrel")\n')    
    outputfile.write('process.particleFlowSuperClusterECAL.PFBasicClusterCollectionEndcap = cms.string("recalibParticleFlowBasicClusterECALEndcap")\n')    
    outputfile.write('process.particleFlowSuperClusterECAL.PFSuperClusterCollectionEndcap = cms.string("recalibParticleFlowSuperClusterECALEndcap")\n')    
    outputfile.write('process.particleFlowSuperClusterECAL.PFBasicClusterCollectionPreshower = cms.string("recalibParticleFlowBasicClusterECALPreshower")\n')    
    outputfile.write('process.particleFlowSuperClusterECAL.PFSuperClusterCollectionEndcapWithPreshower = cms.string("recalibParticleFlowSuperClusterECALEndcapWithPreshower")\n')
    outputfile.write('\n') 

    outputfile.write("process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(" + nEventsPerJob +") )\n")  
    outputfile.write("process.MessageLogger.cerr.FwkReport.reportEvery = 1000\n")

    outputfile.write("process.options = cms.untracked.PSet(\n")
    outputfile.write("   wantSummary = cms.untracked.bool(True),\n")
    outputfile.write(")\n")
    outputfile.write("process.source = cms.Source('PoolSource',\n")
    outputfile.write("    fileNames = cms.untracked.vstring(\n")

def printFillFromRawCfgZ2( outputfile, pwd , iteration, outputDir, ijob ):
    outputfile.write("    )\n")
    outputfile.write(")\n")
    outputfile.write("process.analyzerFillEpsilonForZ = cms.EDAnalyzer('FillEpsilonPlotForZ')\n")    
    outputfile.write("process.analyzerFillEpsilonForZ.OutputDir = cms.untracked.string('" +  outputDir + "')\n")   
    outputfile.write("process.analyzerFillEpsilonForZ.OutputFile = cms.untracked.string('" + NameTag +  outputFile + "_" + str(ijob) + ".root')\n")            
    if (isCRAB):         
        outputfile.write("process.analyzerFillEpsilonForZ.calibMapPath = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + NameTag + calibMapName + "')\n")
        outputfile.write("process.analyzerFillEpsilonForZ.isCRAB  = cms.untracked.bool(True)\n")
    else:
        if (SubmitFurtherIterationsFromExisting and iteration == 0):        
            outputfile.write("process.analyzerFillEpsilonForZ.calibMapPath = cms.untracked.string('root://eoscms//eos/cms" + startingCalibMap + "')\n")            
        else:
            outputfile.write("process.analyzerFillEpsilonForZ.calibMapPath = cms.untracked.string('root://eoscms//eos/cms" + eosPath + "/" + dirname + "/iter_" + str(iteration-1) + "/" + NameTag + calibMapName + "')\n")       
    outputfile.write("\n")        

    outputfile.write("process.analyzerFillEpsilonForZ.CurrentIteration = cms.untracked.int32(" + str(iteration) + ")\n")
    outputfile.write("process.analyzerFillEpsilonForZ.minInvMassCut = cms.untracked.double(" + minInvMass + ")\n")   
    outputfile.write("process.analyzerFillEpsilonForZ.maxInvMassCut = cms.untracked.double(" + maxInvMass + ")\n")   
    outputfile.write("process.analyzerFillEpsilonForZ.requireOppositeCharge = cms.untracked.bool( " + oppositeCharge + ")\n")
    outputfile.write("process.analyzerFillEpsilonForZ.elePtCut = cms.untracked.double( " + elePtCut + ")\n")      
    outputfile.write("process.analyzerFillEpsilonForZ.eleEtaCut = cms.untracked.double( " + eleEtaCut + ")\n")
    outputfile.write("process.analyzerFillEpsilonForZ.maxDReleSc = cms.untracked.double( " + maxDReleSc + ")\n") 
    outputfile.write("process.analyzerFillEpsilonForZ.ZCalib_InvMass = cms.untracked.string('" + ZCalib_InvMass + "')\n")  
    outputfile.write("process.analyzerFillEpsilonForZ.electronSelection = cms.untracked.int32(" + electronSelection + ")\n") 
    outputfile.write("process.analyzerFillEpsilonForZ.Barrel_orEndcap = cms.untracked.string('" + Barrel_or_Endcap + "')\n")        
    outputfile.write("process.analyzerFillEpsilonForZ.useMassInsteadOfEpsilon = cms.untracked.bool( " + useMassInsteadOfEpsilon + ")\n")

    if(len(json_file)>0):              
        outputfile.write("process.analyzerFillEpsilonForZ.JSONfile = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + json_file + "')\n")
    if isMC:                                        
        outputfile.write("process.analyzerFillEpsilonForZ.isMC = cms.untracked.bool(True)\n")           

    outputfile.write("\n")
    outputfile.write("### Not to be changed\n")
    outputfile.write("process.analyzerFillEpsilonForZ.EBRecHitCollectionTag = cms.untracked.InputTag('ecalRecalRecHit','EcalRecalRecHitsEB')\n")
    outputfile.write("process.analyzerFillEpsilonForZ.EERecHitCollectionTag = cms.untracked.InputTag('ecalRecalRecHit','EcalRecalRecHitsEE')\n")
    outputfile.write("process.analyzerFillEpsilonForZ.EBSuperClusterCollectionTag = cms.untracked.InputTag('particleFlowSuperClusterECAL','recalibParticleFlowSuperClusterECALBarrel','analyzerFillEpsilonForZ')\n")    
    outputfile.write("process.analyzerFillEpsilonForZ.EESuperClusterCollectionTag = cms.untracked.InputTag('particleFlowSuperClusterECAL','recalibParticleFlowSuperClusterECALEndcapWithPreshower','analyzerFillEpsilonForZ')\n")    
    outputfile.write("process.analyzerFillEpsilonForZ.ElectronCollectionTag = cms.untracked.InputTag('gedGsfElectrons','','RECO')\n") 
    if isMC:   
        outputfile.write("process.analyzerFillEpsilonForZ.mcProducer = cms.untracked.string('genParticles')\n")
    else:
        outputfile.write("process.analyzerFillEpsilonForZ.mcProducer = cms.untracked.string('')\n")
    outputfile.write("\n")
    outputfile.write("process.p = cms.Path()\n")
    outputfile.write("process.p *= process.bunchSpacingProducer\n")   
    outputfile.write("process.p *= process.ecalDigis\n")   
    outputfile.write("process.p *= process.ecalPreshowerDigis\n")   
    outputfile.write("process.p *= process.ecalMultiFitUncalibRecHit\n")
    outputfile.write("process.p *= process.ecalDetIdToBeRecovered\n")
    outputfile.write("process.p *= process.ecalRecHit\n")
    outputfile.write("process.p *= process.ecalCompactTrigPrim\n")
    outputfile.write("process.p *= process.ecalTPSkim\n")
    outputfile.write("process.p *= process.ecalPreshowerRecHit\n")
    outputfile.write("process.p *= process.ecalRecalRecHit\n")    
    outputfile.write("process.p *= process.particleFlowRecHitPS\n")
    outputfile.write("process.p *= process.particleFlowClusterPS\n")
    outputfile.write("process.p *= process.particleFlowRecHitECAL\n")
    outputfile.write("process.p *= process.particleFlowClusterECALUncorrected\n")
    outputfile.write("process.p *= process.particleFlowClusterECAL\n")
    outputfile.write("process.p *= process.particleFlowSuperClusterECAL\n")
    outputfile.write("process.p *= process.analyzerFillEpsilonForZ\n")

def printFitCfgZ( outputfile, iteration, outputDir, nIn, nFin, EBorEE, nFit ):
    outputfile.write("import FWCore.ParameterSet.Config as cms\n")
    outputfile.write("process = cms.Process('FitEpsilonPlotForZ')\n")
    outputfile.write("process.load('FWCore.MessageService.MessageLogger_cfi')\n")
    outputfile.write("process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )\n")
    outputfile.write("process.source =   cms.Source('EmptySource')\n")
    outputfile.write("process.fitEpsilonForZ = cms.EDAnalyzer('FitEpsilonPlotForZ')\n")
    outputfile.write("process.fitEpsilonForZ.OutputFile = cms.untracked.string('" + NameTag + EBorEE + "_" + str(nFit) + "_" + calibMapName + "')\n")
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
        outputfile.write("process.fitEpsilonForZ.OutputDir = cms.untracked.string('$TMPDIR')\n")
    else:
        outputfile.write("process.fitEpsilonForZ.OutputDir = cms.untracked.string('" +  outputDir + "')\n")
    outputfile.write("process.fitEpsilonForZ.CurrentIteration = cms.untracked.int32(" + str(iteration) + ")\n")
    outputfile.write("process.fitEpsilonForZ.NInFit = cms.untracked.int32(" + str(nIn) + ")\n")
    outputfile.write("process.fitEpsilonForZ.NFinFit = cms.untracked.int32(" + str(nFin) + ")\n")
    outputfile.write("process.fitEpsilonForZ.EEorEB = cms.untracked.string('" + EBorEE + "')\n")

    outputfile.write("process.fitEpsilonForZ.StoreForTest = cms.untracked.bool( True )\n")
    outputfile.write("process.fitEpsilonForZ.Barrel_orEndcap = cms.untracked.string('" + Barrel_or_Endcap + "')\n")
    outputfile.write("process.fitEpsilonForZ.useMassInsteadOfEpsilon = cms.untracked.bool( " + useMassInsteadOfEpsilon + ")\n")
    if not(isCRAB): #If CRAB you have to put the correct path, and you do it on calibJobHandler.py, not on ./submitCalibration.py
        outputfile.write("process.fitEpsilonForZ.EpsilonPlotFileName = cms.untracked.string('root://eoscms//eos/cms" + eosPath + "/" + dirname + "/iter_" + str(iteration) + "/" + NameTag + "epsilonPlots.root')\n")
        if (SubmitFurtherIterationsFromExisting and iteration == 0):
            outputfile.write("process.fitEpsilonForZ.calibMapPath = cms.untracked.string('root://eoscms//eos/cms" + startingCalibMap + "')\n")
        else:
            outputfile.write("process.fitEpsilonForZ.calibMapPath = cms.untracked.string('root://eoscms//eos/cms" + eosPath + "/" + dirname + "/iter_" + str(iteration-1) + "/" + NameTag + calibMapName + "')\n")
    outputfile.write("process.p = cms.Path(process.fitEpsilonForZ)\n")

def printSubmitFitSrc(outputfile, cfgName, source, destination, pwd, logpath):
    outputfile.write("#!/bin/bash\n")
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
        outputfile.write("export SCRAM_ARCH=slc6_amd64_gcc491\n")
        outputfile.write("source $VO_CMS_SW_DIR/cmsset_default.sh\n")
        # outputfile.write("source /cvmfs/cms.cern.ch/crab3/crab.sh\n")
        outputfile.write("export X509_USER_PROXY=/localgrid/lpernie/x509up_u20580\n")
    outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    outputfile.write("echo 'cmsRun " + cfgName + " 2>&1 | awk {quote}/FIT_EPSILON:/ || /WITHOUT CONVERGENCE/ || /HAS CONVERGED/{quote}' > " + logpath  + "\n")
    outputfile.write("cmsRun " + cfgName + " 2>&1 | awk '/FIT_EPSILON:/ || /WITHOUT CONVERGENCE/ || /HAS CONVERGED/' >> " + logpath  + "\n")
    outputfile.write("echo 'ls " + source + " >> " + logpath + " 2>&1' \n" )
    outputfile.write("ls " + source + " >> " + logpath + " 2>&1 \n" )
    if not(isCRAB): #If CRAB you have to put the correct path, anbd you do it on calibJobHandler.py, not on ./submitCalibration.py
       destrooplot = destination.replace("calibMap","fitRes")
       if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
           outputfile.write("echo 'eos cp " + source + " " + destination + "' >> " + logpath  + "\n")
           outputfile.write("echo 'eos cp /tmp/Fit_Stored.root " + destrooplot + "' >> " + logpath  + "\n")
       else:
           outputfile.write("echo 'cmsStage -f " + source + " " + destination + "' >> " + logpath  + "\n")           
           outputfile.write("echo 'cmsStage -f /tmp/Fit_Stored.root " + destrooplot + "' >> " + logpath  + "\n")
       outputfile.write(myeosstage + source + " " + destination + " >> " + logpath + " 2>&1 \n")
       outputfile.write(myeosstage + " /tmp/Fit_Stored.root " + destrooplot + " >> " + logpath + " 2>&1 \n")
       outputfile.write("echo 'rm -f " + source + "' >> " + logpath + " \n")
       outputfile.write("rm -f " + source + " >> " + logpath + " 2>&1 \n")
       outputfile.write("echo 'rm -f /tmp/Fit_Stored.root' >> " + logpath + " \n")
       outputfile.write("rm -f /tmp/Fit_Stored.root >> " + logpath + " 2>&1 \n")

def printSubmitSrc(outputfile, cfgName, source, destination, pwd, logpath):
    outputfile.write("#!/bin/bash\n")
    outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    # outputfile.write("source /cvmfs/cms.cern.ch/crab3/crab.sh\n") this line produces problem when running in CMSSW_8_0_3, anyway we don't use crab
    if ( not isOtherT2 and isCRAB ):
        outputfile.write("setenv X509_USER_PROXY " + CRAB_CopyCert + "\n")
    if not(Silent):
        outputfile.write("echo 'cmsRun " + cfgName + "'\n")
        outputfile.write("cmsRun " + cfgName + "\n")
        if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
            outputfile.write("echo 'eos cp " + source + " " + destination + "'\n")
        else:
            outputfile.write("echo 'cmsStage -f " + source + " " + destination + "'\n")
        outputfile.write(myeosstage + source + " " + destination + "\n")
        outputfile.write("echo 'rm -f " + source + "'\n")
        outputfile.write("rm -f " + source + "\n")
    else:
        outputfile.write("echo 'cmsRun " + cfgName + " 2>&1 | awk {quote}/FILL_COUT:/{quote}' > " + logpath  + "\n")
        outputfile.write("cmsRun " + cfgName + " 2>&1 | awk '/FILL_COUT:/' >> " + logpath  + "\n")
        outputfile.write("echo 'ls " + source + " >> " + logpath + " 2>&1' \n" )
        outputfile.write("ls " + source + " >> " + logpath + " 2>&1 \n" )
        if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
            outputfile.write("echo 'eos cp " + source + " " + destination + "' >> " + logpath  + "\n")
        else:
            outputfile.write("echo 'cmsStage -f " + source + " " + destination + "' >> " + logpath  + "\n")
        outputfile.write(myeosstage + source + " " + destination + " >> " + logpath + " 2>&1 \n")
        outputfile.write("echo 'rm -f " + source + "' >> " + logpath + " \n")
        outputfile.write("rm -f " + source + " >> " + logpath + " 2>&1 \n")

def printParallelHadd(outputfile, outFile, list, destination, pwd):
    import os, sys, imp, re
    CMSSW_VERSION=os.getenv("CMSSW_VERSION")
    outputfile.write("#!/bin/bash\n")
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
       outputfile.write("export SCRAM_ARCH=slc6_amd64_gcc491\n")
       outputfile.write("source $VO_CMS_SW_DIR/cmsset_default.sh\n")
       # outputfile.write("source /cvmfs/cms.cern.ch/crab3/crab.sh\n")
       outputfile.write("export X509_USER_PROXY=/localgrid/lpernie/x509up_u20580\n")
    if(re.match("CMSSW_5_.*_.*",CMSSW_VERSION)):
         print "WARNING!!!! ----> I'm ging to use a harcoded path: /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src"
         print "This because you are in a release CMSSW_5_*_*, that do not allow a hadd with a @file.list."
         outputfile.write("cd /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src\n")
    else:
         outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
       outputfile.write("echo 'hadd -f -k $TMPDIR/" + outFile + " @" + list + "'\n")
       outputfile.write("hadd -f -k $TMPDIR/" + outFile + " @" + list  + "\n")
       outputfile.write("echo 'srmcp file:///$TMPDIR/" + outFile + " " + destination + "/" + outFile + "'\n")
       outputfile.write("srmcp file:///$TMPDIR/" + outFile + " " + destination + "/" + outFile + "\n")
       outputfile.write("rm -f $TMPDIR/" + outFile + "\n")
    else:
       outputfile.write("echo 'hadd -f -k /tmp/" + outFile + " @" + list + "'\n")
       outputfile.write("hadd -f -k /tmp/" + outFile + " @" + list  + "\n")
       if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:       
           outputfile.write("echo 'eos cp /tmp/" + outFile + " " + destination + "'\n")
       else:
           outputfile.write("echo 'cmsStage -f /tmp/" + outFile + " " + destination + "'\n")
       outputfile.write(myeosstage + "/tmp/" + outFile + " " + destination + "\n")
       outputfile.write("rm -f /tmp/" + outFile + "\n")

def printFinalHadd(outputfile, list, destination, pwd):
    import os, sys, imp, re
    CMSSW_VERSION=os.getenv("CMSSW_VERSION")
    outputfile.write("#!/bin/bash\n")
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
       outputfile.write("export SCRAM_ARCH=slc6_amd64_gcc491\n")
       outputfile.write("source $VO_CMS_SW_DIR/cmsset_default.sh\n")
       # outputfile.write("source /cvmfs/cms.cern.ch/crab3/crab.sh\n")
       outputfile.write("export X509_USER_PROXY=/localgrid/lpernie/x509up_u20580\n")
    if(re.match("CMSSW_5_.*_.*",CMSSW_VERSION)):
         print "WARNING!!!! ----> I'm ging to use a harcoded path: /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src"
         print "This because you are in a release CMSSW_5_*_*, that do not allow a hadd with a @file.list."
         outputfile.write("cd /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src\n")
    else:
         outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
       outputfile.write("echo 'hadd -f -k $TMPDIR/" + NameTag + "epsilonPlots.root @" + list + "'\n")
       outputfile.write("hadd -f -k $TMPDIR/" + NameTag + "epsilonPlots.root @" + list  + "\n")
       outputfile.write("echo 'srmcp file:///$TMPDIR/" + NameTag + "epsilonPlots.root " + destination + "/epsilonPlots.root" + "'\n")
       outputfile.write("srmcp file:///$TMPDIR/" + NameTag + "epsilonPlots.root " + destination + "/epsilonPlots.root" + "\n")
       outputfile.write("rm -f $TMPDIR/" + NameTag + "epsilonPlots.root\n")
    else:
       outputfile.write("echo 'hadd -f -k /tmp/" + NameTag + "epsilonPlots.root @" + list + "'\n")
       outputfile.write("hadd -f -k /tmp/" + NameTag + "epsilonPlots.root @" + list  + "\n")
       if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:       
           outputfile.write("echo 'eos cp /tmp/" + NameTag + "epsilonPlots.root " + destination + "'\n")
       else:
           outputfile.write("echo 'cmsStage -f /tmp/" + NameTag + "epsilonPlots.root " + destination + "'\n")
       outputfile.write(myeosstage + "/tmp/" + NameTag + "epsilonPlots.root " + destination + "\n")
       outputfile.write("rm -f /tmp/" + NameTag + "epsilonPlots.root\n")

def printParallelHaddFAST(outputfile, outFile, listReduced, destination, pwd, numList):
    import os, sys, imp, re
    CMSSW_VERSION=os.getenv("CMSSW_VERSION")
    outputfile.write("#!/bin/bash\n")
    if(re.match("CMSSW_5_.*_.*",CMSSW_VERSION)):
         print "WARNING!!!! ----> I'm ging to use a harcoded path: /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src"
         print "This because you are in a release CMSSW_5_*_*, that do not allow a hadd with a @file.list."
         outputfile.write("cd /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src\n")
    else:
         outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    outputfile.write("rm -rf /tmp/" + NameTag + outputFile + "_*\n")
    outputfile.write("rm -rf /tmp/" + NameTag + "FinalFile*\n")
#if we leave "cmsStage -f" to cpy file from eos to /tmp, then ok, otherwise, with "eos cp" files on eos must be preceeded by "root://eoscms//eos/cms". In the lines below $0 is a file read from listreduced, which will be of the form /store/blabla/file.root 
    if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
        outputfile.write("echo \"Copying files locally: awk '{print \"eos cp root://eoscms/eos/cms\\\"$\\0 \\\" /tmp/\"}' " + listReduced + " | bash\"\n")
        outputfile.write("awk '{print \"" + myeosstage + "root://eoscms/eos/cms\"$0 \" /tmp/\"}' " + listReduced + " | bash\n")
    elif "cmsStage -f" in myeosstage:
        outputfile.write("echo \"Copying files locally: awk '{print \"cmsStage -f \\\"$\\0 \\\"  /tmp/\"}' " + listReduced + " | bash\"\n")
        outputfile.write("awk '{print \"" + myeosstage + "\"$0 \" /tmp/\"}' " + listReduced + " | bash\n")        
    else:
        outputfile.write("echo \"Copying files locally: awk '{print \"cmsStage -f \\\"$\\0 \\\"  /tmp/\"}' " + listReduced + " | bash\"\n")
        outputfile.write("awk '{print \"cmsStage -f \"$0 \" /tmp/\"}' " + listReduced + " | bash\n")
    outputfile.write("files=`cat " + listReduced + "`\n")
    outputfile.write("filesHadd=''\n")
    outputfile.write("for file in $files;\n")
    outputfile.write("do\n")
    if( isCRAB ):
        outputfile.write("   SUBSTRING=`echo ${file} | awk -F / '{ print $14 }'`\n")
    else:
        outputfile.write("   SUBSTRING=`echo ${file} | awk -F / '{ print $9 }'`\n")  # since I added a directory on eos, must print 10th position, not 9th
    outputfile.write('   filesHadd="$filesHadd /tmp/$SUBSTRING"\n')
    outputfile.write("done\n")
    outputfile.write("echo \"hadd -k /tmp/" + NameTag + "epsilonPlots_" + str(numList) + ".root $filesHadd\"\n")
    outputfile.write("hadd -k /tmp/" + NameTag + "epsilonPlots_" + str(numList) + ".root $filesHadd\n")
    if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
        outputfile.write("echo \"eos cp /tmp/" + NameTag + "epsilonPlots_" + str(numList) + ".root " + destination + "\"\n")
    else:
        outputfile.write("echo \"cmsStage -f /tmp/" + NameTag + "epsilonPlots_" + str(numList) + ".root " + destination + "\"\n")
    outputfile.write(myeosstage + "/tmp/" + NameTag + "epsilonPlots_" + str(numList) + ".root " + destination + "\n")
    outputfile.write("rm -rf /tmp/" + NameTag + outputFile + "_*\n")
    outputfile.write("rm -rf /tmp/" + NameTag + "epsilonPlots*\n")

def printFinalHaddFAST(outputfile, listReduced, destination, pwd):
    import os, sys, imp, re
    CMSSW_VERSION=os.getenv("CMSSW_VERSION")
    outputfile.write("#!/bin/bash\n")
    if(re.match("CMSSW_5_.*_.*",CMSSW_VERSION)):
         print "WARNING!!!! ----> I'm ging to use a harcoded path: /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src"
         print "This because you are in a release CMSSW_5_*_*, that do not allow a hadd with a @file.list."
         outputfile.write("cd /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src\n")
    else:
         outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    outputfile.write("rm -rf /tmp/" + NameTag + "epsilonPlots*\n")
    outputfile.write("rm -rf /tmp/" + NameTag + "FinalFile*\n")
#if we leave "cmsStage -f" to cpy file from eos to /tmp, then ok, otherwise, with "eos cp" files on eos must be preceeded by "root://eoscms//eos/cms". In the lines below $0 is a file read from listreduced, which will be of the form /store/blabla/file.root 
    if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
        outputfile.write("echo \"Copying files locally: awk '{print \"eos cp root://eoscms/eos/cms\\\"$\\0 \\\" /tmp/\"}' " + listReduced + " | bash\"\n")
        outputfile.write("awk '{print \"" + myeosstage + "root://eoscms/eos/cms\"$0 \" /tmp/\"}' " + listReduced + " | bash\n")
    elif "cmsStage -f" in myeosstage:
        outputfile.write("echo \"Copying files locally: awk '{print \"cmsStage -f \\\"$\\0 \\\" /tmp/\"}' " + listReduced + " | bash\"\n")
        outputfile.write("awk '{print \"" + myeosstage + "\"$0 \" /tmp/\"}' " + listReduced + " | bash\n")        
    else:
        outputfile.write("echo \"Copying files locally: awk '{print \"cmsStage -f \\\"$\\0 \\\" /tmp/\"}' " + listReduced + " | bash\"\n")
        outputfile.write("awk '{print \"cmsStage -f \"$0 \" /tmp/\"}' " + listReduced + " | bash\n")
    outputfile.write("files=`cat " + listReduced + "`\n")
    outputfile.write("filesHadd=''\n")
    outputfile.write("for file in $files;\n")
    outputfile.write("do\n")
    if( isCRAB ):
        outputfile.write("   SUBSTRING=`echo ${file} | awk -F / '{ print $14 }'`\n")
    else:
        outputfile.write("   SUBSTRING=`echo ${file} | awk -F / '{ print $7 }'`\n")   # since I added a directory on eos, must print 10th position, not 9th
    outputfile.write('   filesHadd="$filesHadd /tmp/$SUBSTRING"\n')
    outputfile.write("done\n")
    outputfile.write("echo \"hadd -k /tmp/" + NameTag + "epsilonPlots.root $filesHadd\"\n")
    outputfile.write("hadd -k /tmp/" + NameTag + "epsilonPlots.root $filesHadd\n")
    if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
        outputfile.write("echo \"eos cp /tmp/" + NameTag + "epsilonPlots.root " + destination + "\"\n")
    else:
        outputfile.write("echo \"cmsStage -f /tmp/" + NameTag + "epsilonPlots.root " + destination + "\"\n")
    outputfile.write(myeosstage + "/tmp/" + NameTag + "epsilonPlots.root " + destination + "\n")
    outputfile.write("rm -rf /tmp/" + NameTag + "epsilonPlots*\n")



def printFinalHaddRegroup(outputfile, listReduced, destination, pwd, grouping=10):
    import os, sys, imp, re, ntpath
    CMSSW_VERSION=os.getenv("CMSSW_VERSION")
    outputfile.write("#!/bin/bash\n")
    outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    outputfile.write("rm -rf /tmp/" + NameTag + "epsilonPlots*\n")
    outputfile.write("rm -rf /tmp/" + NameTag + "FinalFile*\n")
    
    fileWithList = open(listReduced,"r")
    files = fileWithList.readlines()
    idx=0
    grouped_files = []
    while len(files)>0:
        filesToMerge = files[:grouping]
        mergedfile_n = "/tmp/hadded_epsilon_"+str(idx)+".root"
        outputfile.write("echo Copying files locally\n")
        strippedFiles = []
        for f in filesToMerge:
            f = f.strip()
            if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
                outputfile.write(myeosstage + "root://eoscms/eos/cms" + f + " /tmp \n")
            elif "cmsStage -f" in myeosstage:
                outputfile.write(myeosstage + f + " /tmp \n")        
            else:
                outputfile.write("cmsStage -f " + f + " /tmp \n")
            strippedFiles.append(ntpath.basename(f))
        outputfile.write("filesHadd=\"/tmp/" + " /tmp/".join(strippedFiles) + "\"\n")
        outputfile.write("echo \"hadd -k " + mergedfile_n + " $filesHadd\"\n")
        outputfile.write("hadd -k " + mergedfile_n + " $filesHadd\n")
        outputfile.write("rm -rf /tmp/" + NameTag + "epsilonPlots*\n\n")

        grouped_files.append(mergedfile_n)
        idx += 1
        files = files[grouping:]

    outputfile.write("echo now hadding the intermediate hadded files: " + " ".join(grouped_files) + "\n")
    outputfile.write("hadd -k /tmp/" +  NameTag + "epsilonPlots.root " + " ".join(grouped_files) + "\n")

    if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
        outputfile.write("echo \"eos cp /tmp/" + NameTag + "epsilonPlots.root " + destination + "\"\n")
    else:
        outputfile.write("echo \"cmsStage -f /tmp/" + NameTag + "epsilonPlots.root " + destination + "\"\n")
    outputfile.write(myeosstage + "/tmp/" + NameTag + "epsilonPlots.root " + destination + "\n")
    outputfile.write("rm -rf /tmp/" + NameTag + "epsilonPlots*\n")
    outputfile.write("rm -rf /tmp/" + NameTag + "hadded_epsilon*\n")

