import ROOT
from ROOT import TF1
import math
from functools import partial
import CombineHarvester.CombineTools.plotting as plot
import json
import argparse
import os
import numpy as np
from glob import glob


def prepare_scenario0(dc_input, lumi_values, pol1_up, pol1_noup, pol2_up, pol2_noup, prev_lumi = 0, freeze_theo = False):
   for index, lumi_value in enumerate(lumi_values):
        datacard = dc_input
        dc_S0 = open(datacard, 'r')
        lines = dc_S0.readlines()
        name = 'S0_dc_lumi' + str(int(lumi_value)) 
        dc_S0_new = open(name + '.txt', "w+")
        for line in lines:
                dc_S0_new.write(line)
	eff = 1. - (pol1_noup.Eval(lumi_values[0])-pol1_noup.Eval(int(lumi_value)))
        dc_S0_new.write("lumiscale_" + str(lumi_value) + " rateParam * * {}\n".format(((lumi_value/59.4)-(prev_lumi/59.4))*eff))
        dc_S0_new.write("nuisance edit freeze lumiscale_" + str(lumi_value) + "\n")
        dc_S0_new.close()
        print('Create the new datacard: {}\n'.format(name))
        
	if index != 0:
		os.system("combineCards.py S0_dc_lumi*.txt > combined_" + name + ".txt")
	else:
		os.system("cp " + name + ".txt combined_" + name + ".txt")
        
	scale_factor = (pol2_noup.Eval(int(lumi_value))/pol2_noup.Eval(lumi_values[0]))
	jes = (scale_factor-1)*0.3 +1
        os.system("text2workspace.py combined_" + name + ".txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO 'map=.*/WWewk:r_vbs[1,-10,10]' --X-rescale-nuisance 'CMS_scale_e_2018' " + str(scale_factor) + " --X-rescale-nuisance 'CMS_scale_met_2018' " + str(scale_factor) + " --X-rescale-nuisance 'CMS_scale_JESRelative*' " + str(jes) + " -v 0 -o combined_" + name + ".root")

        os.system("combine combined_" + name + ".root -M MultiDimFit -t -1 --setParameters r_vbs=1 --saveWorkspace -n S0_postfit_lumi" + str(int(lumi_value)))

        os.system("combineTool.py higgsCombineS0_postfit_lumi" + str(int(lumi_value)) + ".MultiDimFit.mH120.root -M MultiDimFit -t -1 --setParameters r_vbs=1 -n S0_lumi" + str(int(lumi_value)) + """ --algo grid --points 150 --split-points 20 --snapshotName MultiDimFit --setParameterRanges r_vbs=0,2 --job-mode condor --sub-opts='+JobFlavour="workday"' --task-name S0_lumi""" + str(int(lumi_value)))	
	prev_lumi = lumi_value
   return

def prepare_scenario0plus(dc_input, lumi_values, pol1_up, pol1_noup, pol2_up, pol2_noup, prev_lumi = 0):
   for index, lumi_value in enumerate(lumi_values):
        datacard = dc_input
        dc_S0plus = open(datacard, 'r')
        lines = dc_S0plus.readlines()
        name = 'S0plus_dc_lumi' + str(int(lumi_value)) 
        dc_S0plus_new = open(name + '.txt', "w+")
        for line in lines:
                dc_S0plus_new.write(line)
	eff = 1. - (pol1_noup.Eval(lumi_values[0])-pol1_up.Eval(int(lumi_value)))
        dc_S0plus_new.write("lumiscale_" + str(lumi_value) + " rateParam * * {}\n".format(((lumi_value/59.4)-(prev_lumi/59.4))*eff))
        dc_S0plus_new.write("nuisance edit freeze lumiscale_" + str(lumi_value) + "\n")
        dc_S0plus_new.close()
        print('Create the new datacard: {}\n'.format(name))
        
	if index != 0:
		os.system("combineCards.py S0plus_dc_lumi*.txt > combined_" + name + ".txt")
	else:
		os.system("cp " + name + ".txt combined_" + name + ".txt")
        
	scale_factor = pol2_up.Eval(int(lumi_value))/pol2_noup.Eval(lumi_values[0])
	if scale_factor < 1: scale_factor = 1
	jes = (scale_factor-1)*0.3 +1

        os.system("text2workspace.py combined_" + name + ".txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO 'map=.*/WWewk:r_vbs[1,-10,10]' --X-rescale-nuisance 'CMS_scale_e_2018' " + str(scale_factor) + " --X-rescale-nuisance 'CMS_scale_met_2018' " + str(scale_factor) + " --X-rescale-nuisance 'CMS_scale_JESRelative*' " + str(jes) + " --X-rescale-nuisance 'CMS_eff_m_2018' 0.7 -v 0 -o combined_" + name + ".root")

        os.system("combine combined_" + name + ".root -M MultiDimFit -t -1 --setParameters r_vbs=1 --saveWorkspace -n S0plus_postfit_lumi" + str(int(lumi_value)))

        os.system("combineTool.py higgsCombineS0plus_postfit_lumi" + str(int(lumi_value)) + ".MultiDimFit.mH120.root -M MultiDimFit -t -1 --setParameters r_vbs=1 -n S0plus_lumi" + str(int(lumi_value)) + """ --algo grid --points 150 --split-points 20 --snapshotName MultiDimFit --setParameterRanges r_vbs=0,2 --job-mode condor --sub-opts='+JobFlavour="workday"' --task-name S0plus_lumi""" + str(int(lumi_value)))

	prev_lumi = lumi_value
   return


def prepare_scenario1(dc_input, lumi_values, pol1_up, pol1_noup, pol2_up, pol2_noup, prev_lumi = 0):
   for index, lumi_value in enumerate(lumi_values):
        datacard = dc_input
        dc_S1 = open(datacard, 'r')
        lines = dc_S1.readlines()
        name = 'S1_dc_lumi' + str(int(lumi_value)) 
        dc_S1_new = open(name + '.txt', "w+")
        for line in lines:
                dc_S1_new.write(line)
	eff = 1. - (pol1_noup.Eval(lumi_values[0])-pol1_noup.Eval(int(lumi_value)))
        dc_S1_new.write("lumiscale_" + str(lumi_value) + " rateParam * * {}\n".format(((lumi_value/59.4)-(prev_lumi/59.4))*eff))
        dc_S1_new.write("nuisance edit freeze lumiscale_" + str(lumi_value) + "\n")
        dc_S1_new.close()
        print('Create the new datacard: {}\n'.format(name))
        
	if index != 0:
		os.system("combineCards.py S1_dc_lumi*.txt > combined_" + name + ".txt")
	else:
		os.system("cp " + name + ".txt combined_" + name + ".txt")
        
	scale_factor = (pol2_noup.Eval(int(lumi_value))/pol2_noup.Eval(lumi_values[0]))
	jes = (scale_factor-1)*0.3 +1
        os.system("text2workspace.py combined_" + name + ".txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO 'map=.*/WWewk:r_vbs[1,-10,10]' --X-rescale-nuisance 'CMS_scale_e_2018' " + str(scale_factor) + " --X-rescale-nuisance 'CMS_scale_met_2018' " + str(scale_factor) + " --X-rescale-nuisance 'CMS_scale_JESRelative*' " + str(jes) + " --X-rescale-nuisance 'QCDscale*' 0.5 --X-rescale-nuisance 'pdf_*' 0.5 -v 0 -o combined_" + name + ".root")

        os.system("combine combined_" + name + ".root -M MultiDimFit -t -1 --setParameters r_vbs=1 --saveWorkspace -n S1_postfit_lumi" + str(int(lumi_value)))

        os.system("combineTool.py higgsCombineS1_postfit_lumi" + str(int(lumi_value)) + ".MultiDimFit.mH120.root -M MultiDimFit -t -1 --setParameters r_vbs=1 -n S1_lumi" + str(int(lumi_value)) + """ --algo grid --points 150 --split-points 20 --snapshotName MultiDimFit --setParameterRanges r_vbs=0,2 --job-mode condor --sub-opts='+JobFlavour="workday"' --task-name S1_lumi""" + str(int(lumi_value)))

	prev_lumi = lumi_value
   return

def prepare_scenario1plus(dc_input, lumi_values, pol1_up, pol1_noup, pol2_up, pol2_noup, prev_lumi = 0):
   for index, lumi_value in enumerate(lumi_values):
        datacard = dc_input
        dc_S1plus = open(datacard, 'r')
        lines = dc_S1plus.readlines()
        name = 'S1plus_dc_lumi' + str(int(lumi_value)) 
        dc_S1plus_new = open(name + '.txt', "w+")
        for line in lines:
                dc_S1plus_new.write(line)
	eff = 1. - (pol1_noup.Eval(lumi_values[0])-pol1_up.Eval(int(lumi_value)))
        dc_S1plus_new.write("lumiscale_" + str(lumi_value) + " rateParam * * {}\n".format(((lumi_value/59.4)-(prev_lumi/59.4))*eff))
        dc_S1plus_new.write("nuisance edit freeze lumiscale_" + str(lumi_value) + "\n")
        dc_S1plus_new.close()
        print('Create the new datacard: {}\n'.format(name))
        
	if index != 0:
		os.system("combineCards.py S1plus_dc_lumi*.txt > combined_" + name + ".txt")
	else:
		os.system("cp " + name + ".txt combined_" + name + ".txt")
        
	scale_factor = pol2_up.Eval(int(lumi_value))/pol2_noup.Eval(lumi_values[0])
	jes = (scale_factor-1)*0.3 +1
        os.system("text2workspace.py combined_" + name + ".txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO 'map=.*/WWewk:r_vbs[1,-10,10]' --X-rescale-nuisance 'CMS_scale_e_2018' " + str(scale_factor) + " --X-rescale-nuisance 'CMS_scale_met_2018' " + str(scale_factor) + " --X-rescale-nuisance 'CMS_scale_JESRelative*' " + str(jes) + " --X-rescale-nuisance 'QCDscale*' 0.5 --X-rescale-nuisance 'pdf_*' 0.5 --X-rescale-nuisance 'CMS_eff_m_2018' 0.7 -v 0 -o combined_" + name + ".root")

        os.system("combine combined_" + name + ".root -M MultiDimFit -t -1 --setParameters r_vbs=1 --saveWorkspace -n S1plus_postfit_lumi" + str(int(lumi_value)))

        os.system("combineTool.py higgsCombineS1plus_postfit_lumi" + str(int(lumi_value)) + ".MultiDimFit.mH120.root -M MultiDimFit -t -1 --setParameters r_vbs=1 -n S1plus_lumi" + str(int(lumi_value)) + """ --algo grid --points 150 --split-points 20 --snapshotName MultiDimFit --setParameterRanges r_vbs=0,2 --job-mode condor --sub-opts='+JobFlavour="workday"' --task-name S1plus_lumi""" + str(int(lumi_value)))

	prev_lumi = lumi_value
   return


def prepare_scenario2(dc_input, lumi_values, pol1_up, pol1_noup, pol2_up, pol2_noup, prev_lumi = 0):
   for index, lumi_value in enumerate(lumi_values):
        datacard = dc_input
        dc_S2 = open(datacard, 'r')
        lines = dc_S2.readlines()
        name = 'S2_dc_lumi' + str(int(lumi_value)) 
        dc_S2_new = open(name + '.txt', "w+")
        for line in lines:
                dc_S2_new.write(line)
	eff = 1. - (pol1_noup.Eval(lumi_values[0])-pol1_noup.Eval(int(lumi_value)))
        dc_S2_new.write("lumiscale_" + str(lumi_value) + " rateParam * * {}\n".format(((lumi_value/59.4)-(prev_lumi/59.4))*eff))
        dc_S2_new.write("nuisance edit freeze lumiscale_" + str(lumi_value) + "\n")
        dc_S2_new.close()
        print('Create the new datacard: {}\n'.format(name))
        
	if index != 0:
		os.system("combineCards.py S2_dc_lumi*.txt > combined_" + name + ".txt")
	else:
		os.system("cp " + name + ".txt combined_" + name + ".txt")
        
	scale_factor = round((pol2_noup.Eval(int(lumi_value))/pol2_noup.Eval(lumi_values[0])), 1)
	print(scale_factor)
	jes = (scale_factor-1)*0.3 +1
        os.system("text2workspace.py combined_" + name + ".txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO 'map=.*/WWewk:r_vbs[1,-10,10]' --X-rescale-nuisance 'CMS_scale_e_2018' " + str(scale_factor) + " --X-rescale-nuisance 'CMS_scale_met_2018' " + str(scale_factor) + " --X-rescale-nuisance 'CMS_scale_JESRelative*' " + str(jes) + " -v -3 -o combined_" + name + ".root")

        os.system("combine combined_" + name + ".root -M MultiDimFit -t -1 --setParameters r_vbs=1 --saveWorkspace -n S2_postfit_lumi" + str(int(lumi_value)))

        os.system("combineTool.py higgsCombineS2_postfit_lumi" + str(int(lumi_value)) + ".MultiDimFit.mH120.root -M MultiDimFit -t -1 --setParameters r_vbs=1 -n S2_lumi" + str(int(lumi_value)) + """ --algo grid --points 150 --split-points 20 --snapshotName MultiDimFit --freezeNuisanceGroups theory --setParameterRanges r_vbs=0,2 --job-mode condor --sub-opts='+JobFlavour="workday"' --task-name S2_lumi""" + str(int(lumi_value)))
	
	prev_lumi = lumi_value
   return

def prepare_scenario2plus(dc_input, lumi_values, pol1_up, pol1_noup, pol2_up, pol2_noup, prev_lumi = 0):
   for index, lumi_value in enumerate(lumi_values):
        datacard = dc_input
        dc_S2plus = open(datacard, 'r')
        lines = dc_S2plus.readlines()
        name = 'S2plus_dc_lumi' + str(int(lumi_value)) 
        dc_S2plus_new = open(name + '.txt', "w+")
        for line in lines:
                dc_S2plus_new.write(line)
	eff = 1. - (pol1_noup.Eval(lumi_values[0])-pol1_up.Eval(int(lumi_value)))
        dc_S2plus_new.write("lumiscale_" + str(lumi_value) + " rateParam * * {}\n".format(((lumi_value/59.4)-(prev_lumi/59.4))*eff))
        dc_S2plus_new.write("nuisance edit freeze lumiscale_" + str(lumi_value)+ "\n")
        dc_S2plus_new.close()
        print('Create the new datacard: {}\n'.format(name))
        
	if index != 0:
		os.system("combineCards.py S2plus_dc_lumi*.txt > combined_" + name + ".txt")
	else:
		os.system("cp " + name + ".txt combined_" + name + ".txt")
        
	scale_factor = round((pol2_up.Eval(int(lumi_value))/pol2_noup.Eval(lumi_values[0])), 1)
	print(scale_factor)
	jes = (scale_factor-1)*0.3 +1
        os.system("text2workspace.py combined_" + name + ".txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO 'map=.*/WWewk:r_vbs[1,-10,10]' --X-rescale-nuisance 'CMS_scale_e_2018' " + str(scale_factor) + " --X-rescale-nuisance 'CMS_scale_met_2018' " + str(scale_factor) + " --X-rescale-nuisance 'CMS_scale_JESRelative*' " + str(jes) + " --X-rescale-nuisance 'CMS_eff_m_2018' 0.6 -v 0 -o combined_" + name + ".root")

        os.system("combine combined_" + name + ".root -M MultiDimFit -t -1 --setParameters r_vbs=1 --saveWorkspace -n S2plus_postfit_lumi" + str(int(lumi_value)))

        os.system("combineTool.py higgsCombineS2plus_postfit_lumi" + str(int(lumi_value)) + ".MultiDimFit.mH120.root -M MultiDimFit -t -1 --setParameters r_vbs=1 -n S2plus_lumi" + str(int(lumi_value)) + """ --algo grid --points 150 --split-points 20 --snapshotName MultiDimFit --freezeNuisanceGroups theory --setParameterRanges r_vbs=0,2 --job-mode condor --sub-opts='+JobFlavour="workday"' --task-name S2plus_lumi""" + str(int(lumi_value)))

	prev_lumi = lumi_value
   return

def prepare_scenario_freezeall(dc_input, lumi_values, prev_lumi = 0):
   for index, lumi_value in enumerate(lumi_values):
        datacard = dc_input
        dc_freeze_all = open(datacard, 'r')
        lines = dc_freeze_all.readlines()
        name = 'freeze_all_dc_lumi' + str(int(lumi_value)) 
        dc_freeze_all_new = open(name + '.txt', "w+")
        for line in lines:
                dc_freeze_all_new.write(line)
        if index == 0:
		dc_freeze_all_new.write("lumiscale_" + str(lumi_value) + " rateParam * * {}\n".format(str((lumi_values[index])/59.4)))	
	else:
		dc_freeze_all_new.write("lumiscale_" + str(lumi_value) + " rateParam * * {}\n".format(str((lumi_values[index]-lumi_values[index-1])/59.4)))
        dc_freeze_all_new.write("nuisance edit freeze lumiscale_" + str(lumi_value) + "\n")
        dc_freeze_all_new.close()
        print('Create the new datacard: {}\n'.format(name))
        
	if index != 0:
		os.system("combineCards.py freeze_all_dc_lumi*.txt > combined_" + name + ".txt")
	else:
		os.system("cp " + name + ".txt combined_" + name + ".txt")

        os.system("text2workspace.py combined_" + name + ".txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO 'map=.*/WWewk:r_vbs[1,-10,10]' -v 0 -o combined_" + name + ".root")

        os.system("combine combined_" + name + ".root -M MultiDimFit -t -1 --setParameters r_vbs=1 --saveWorkspace -n freeze_all_postfit_lumi" + str(int(lumi_value)))

        os.system("combineTool.py higgsCombinefreeze_all_postfit_lumi" + str(int(lumi_value)) + ".MultiDimFit.mH120.root -M MultiDimFit -t -1 --setParameters r_vbs=1 -n freeze_all_lumi" + str(int(lumi_value)) + """ --algo grid --points 150 --split-points 20 --snapshotName MultiDimFit --freezeParameters allConstrainedNuisances --setParameterRanges r_vbs=0,2 --job-mode condor --sub-opts='+JobFlavour="workday"' --task-name freeze_all_lumi""" + str(int(lumi_value)))
	
	old_name = name
	prev_lumi = lumi_value
   return

