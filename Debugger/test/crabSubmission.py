import sys
from sys import argv
import os
import subprocess
from shutil import copyfile

samples = {#"WJetsToLNu":"/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16DR80-FlatPU20to70HcalNZSRAW_withHLT_80X_mcRun2_asymptotic_v14-v1/GEN-SIM-RAW",
	"WJetsToLNu":"/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16DR80-FlatPU20to70HcalNZSRAW_withHLT_80X_mcRun2_asymptotic_v14-v1/AODSIM",
	   #"Data_ZMu":"/SingleMuon/Run2016G-ZMu-PromptReco-v1/RAW-RECO",
	   #"DYJetsToLNu":"/DYToLL_M_1_TuneCUETP8M1_13TeV_pythia8/RunIISpring16DR80-FlatPU20to70HcalNZSRAW_withHLT_80X_mcRun2_asymptotic_v14_ext1-v1/GEN-SIM-RAW",
	"DYJetsToLNu":"/DYToLL_M_1_TuneCUETP8M1_13TeV_pythia8/RunIISpring16DR80-FlatPU20to70HcalNZSRAW_withHLT_80X_mcRun2_asymptotic_v14_ext1-v1/AODSIM",
	   #"TTbar":"/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16DR80-FlatPU20to70HcalNZSRAW_withHLT_80X_mcRun2_asymptotic_v14_ext3-v1/GEN-SIM-RAW",
	"TTbar":"/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16DR80-FlatPU20to70HcalNZSRAW_withHLT_80X_mcRun2_asymptotic_v14_ext3-v1/AODSIM",
	#"JPsi": "/JPsiToMuMu_Pt20to100-pythia8-gun/RunIISpring16DR80-PU2016_Classic_withHLT_80X_mcRun2_asymptotic_v14-v1/GEN-SIM-RAW",
	"JPsi": "/JPsiToMuMu_Pt20to100-pythia8-gun/RunIISpring16DR80-PU2016_Classic_withHLT_80X_mcRun2_asymptotic_v14-v1/AODSIM",
	
#          "Muminus_Pt10":"/Muminus_Pt10-gun/RunIISpring16DR80-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/RAWAODSIM",
#	   "Muminus_Pt100":"/Muminus_Pt100-gun/RunIISpring16DR80-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/RAWAODSIM",
#	   "Muminus_Pt1000":"/Muminus_Pt1000-gun/RunIISpring16DR80-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/RAWAODSIM",
#          "Muplus_Pt10":"/Muplus_Pt10-gun/RunIISpring16DR80-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/RAWAODSIM",
#	   "Muplus_Pt100":"/Muplus_Pt100-gun/RunIISpring16DR80-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/RAWAODSIM",
#	   "Muplus_Pt1000":"/Muplus_Pt1000-gun/RunIISpring16DR80-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/RAWAODSIM"
	   }

user = "folguera"

lumiMasks = {"Data_ZMu":"https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt"}

mcConfig = '''
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 5
config.Data.useParent   = True
'''

dataConfig = '''
config.Data.splitting = "LumiBased"
config.Data.unitsPerJob = 50
'''

template = '''
from CRABClient.UserUtilities import config
config = config()
	
config.General.requestName = "%s"
config.General.workArea = "%s"	
config.General.transferOutputs = True
config.General.transferLogs = False
	
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "%s"
config.JobType.priority = 1

config.Data.inputDataset = "%s"
config.Data.inputDBS = "global"
%s

config.Data.publication = False
config.Data.ignoreLocality = True
config.Data.outLFNDirBase = "/store/user/folguera/MuonHLT/Effs_161121/%s/%s/"


config.Site.storageSite = "T2_CH_CERN"


'''


def main():
	configuration = argv[1]
	l3config  = argv[2]
	for sname in samples: 
		sample = samples[sname]
		print "Preparing to run over "+sample
		name=sname

		if "Data" in name:
			splittingConfig = dataConfig
			lumiMask = lumiMasks[name]
		else:
			splittingConfig = mcConfig
			lumiMask = ""

		workArea = os.getcwd()
		workDir = os.getcwd()+"/crab_effs/"+name+"_"+l3config.replace('.py','')
		if not os.path.exists(workDir):
			os.makedirs(workDir)	
		
		if "/" in configuration:
			copyfile(configuration, workDir+"/"+configuration.split("/")[-1])
			configuration = configuration.split("/")[-1]
		else:
			copyfile(configuration,workDir+"/"+configuration)
			copyfile(configuration,workDir+"/"+l3config)
			

		config = template%(str(name),str(workDir),str(workArea+"/"+configuration),str(sample),splittingConfig,l3config.replace('.py',''),name)	

		configFile = open(workDir+"/crabConfig_%s.py"%name, "w")
		configFile.write(config)
		configFile.close()
			
main()
