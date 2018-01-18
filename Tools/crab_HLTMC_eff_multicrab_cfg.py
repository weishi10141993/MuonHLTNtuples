from CRABClient.UserUtilities import config

config = config()

config.JobType.pluginName   = 'Analysis'
config.JobType.outputFiles  = ['muonNtuple.root']#, 'DQMIO.root']

config.Data.unitsPerJob     = 10000
config.Data.totalUnits      = 2000000
config.Data.splitting       = 'EventAwareLumiBased'

config.Data.useParent       = True #!!!!

config.Site.storageSite     = 'T3_US_TAMU'
config.JobType.numCores     = 4

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    tag = 'efficiency_iterL3PlusL1_v21'

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea   = tag
    config.Data.outLFNDirBase = '/store/user/wshi/' + tag

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    datasets = {}
    
    datasets['DY']        = ('/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_90X_upgrade2017_realistic_v6_C1-v1/AODSIM','hltConfig.py')
    #datasets['DY_Zpt']    = ('/DYJetsToLL_M-50_Zpt-150toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_90X_upgrade2017_realistic_v6_C1-v1/AODSIM','hltConfig.py')    
    datasets['Displaced'] = ('/DisplacedSUSY_StopToBL_M-400_CTau-10_TuneCUETP8M1_13TeV_pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_90X_upgrade2017_realistic_v6_C1-v1/AODSIM','hltConfig.py')
    #datasets['JPsi']      = ('/JpsiToMuMu_JpsiPt8_TuneCUEP8M1_13TeV-pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_90X_upgrade2017_realistic_v6_C1-v1/AODSIM', 'hltConfig.py')

    for k, v in datasets.iteritems():
        print v[0]
        config.JobType.psetName    = v[1]
        config.General.requestName = k
        config.Data.inputDataset   = v[0]
        config.Data.outputDatasetTag   = 'efficiency_iterL3PlusL1_v21_'+k
        print 'submitting config:'
        print config
        submit(config)


