#!/usr/bin/env python
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

from scenarios_functions import *
from plot_functions import *


#ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

plot.ModTDRStyle(width=700, l=0.13)
ROOT.gStyle.SetNdivisions(510, "XYZ")
ROOT.gStyle.SetMarkerSize(1)

NAMECOUNTER = 0
os.system("cmsenv")

parser = argparse.ArgumentParser()

parser.add_argument('main', help='Main input file for the scan')
parser.add_argument('--y-max', type=float, default=5, help='y-axis maximum')
parser.add_argument('--output', '-o', help='output name without file extension', default='energy_resolution')
parser.add_argument('--main-label', default='Original Datacard', type=str, help='legend label for the main scan')
parser.add_argument('--scenario', default='0', type=str, help='choose which scenario to extrapolate')
parser.add_argument('--folder', default='', type=str, help='choose which will be the working folder. Put / at the end of the name')
parser.add_argument('--main-color', default=1, type=int, help='line and marker color for main scan')
parser.add_argument('--test', default=0, type=int, help='choose if test the program or not')
parser.add_argument('--lumi-step', default=1000, type=int, help='step di campionamento luminosita')
parser.add_argument('--logo', default='CMS')
parser.add_argument('--logo-sub', default='Internal')
args = parser.parse_args()

print '--------------------------------------'
print  args.output
print '--------------------------------------'

#metto in un vettore quali scenari si volgiono estrapolare
scenarios = []
if args.scenario is not None:
    scenarios = (args.scenario).split(':')
print('Requested scenarios: {}'.format(scenarios))


# energy resolution plot
y_upgrade    =       [0.6, 0.8, 1.5, 2.2] # pag 153 TDR
y_no_upgrade = [0.80, 0.9, 1.2, 2.3, 3.3] # tabella 9.3 drive
x_lumi       = [150.,  300.,  1000., 3000., 4500.]

pol2_up, pol2_noup = energy_resolution(x_lumi, y_upgrade, y_no_upgrade, args.y_max)

y_no_upgrade_eff    = [0.94, 0.875, 0.86, 0.80, 0.78] #, 0.74]
y_upgrade_eff =       [0.94, 0.935, 0.92, 0.90, 0.88] #, 0,87] 

pol1_up, pol1_noup = electron_efficiency(x_lumi, y_upgrade_eff, y_no_upgrade_eff, 1.1)



#####################################
	#costruzione scenari 
#####################################

lumi_value = [150, 200, 300, 450]
lumi_value.extend(np.arange(600, 3100, args.lumi_step))

print('Requested extrapolation @ luminosity {}'.format(lumi_value))

if '0' in scenarios:
  prepare_scenario0(args.main, lumi_value, pol1_up, pol1_noup, pol2_up, pol2_noup)
  prepare_scenario0plus(args.main, lumi_value, pol1_up, pol1_noup, pol2_up, pol2_noup)

if '1' in scenarios:
  prepare_scenario1(args.main, lumi_value, pol1_up, pol1_noup, pol2_up, pol2_noup)
  prepare_scenario1plus(args.main, lumi_value, pol1_up, pol1_noup, pol2_up, pol2_noup)

if '2' in scenarios:
  prepare_scenario2(args.main, lumi_value, pol1_up, pol1_noup, pol2_up, pol2_noup)
  prepare_scenario2plus(args.main, lumi_value, pol1_up, pol1_noup, pol2_up, pol2_noup)

prepare_scenario_freezeall(args.main, lumi_value)

if args.test:
    scenario_test()

