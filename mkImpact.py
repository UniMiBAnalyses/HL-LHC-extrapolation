#!/usr/bin/env python
import json
import argparse
import os

#python mkImpact.py -if combined_S2_dc_lumi3150.root --start 1
#python mkImpact.py -if combined_S2_dc_lumi3150.root -o freeze-theo_met --start 2

parser = argparse.ArgumentParser()

#parser.add_argument('main', help='Main input file for the scan')
parser.add_argument('--output', '-o', help='output name without file extension', default='3000/fb-met')
parser.add_argument('--input_file', '-if', help='input file name (root-file!!!)')
parser.add_argument('--POI', help='use this parameter of interest', default='r_vbs')
parser.add_argument('--start', type=int, help='1 if condor is not already launched, 2 if it is', default='1')
parser.add_argument('--web_dir_path', type=str, help='absolute path of the web directory', default='/eos/home-m/mchiusi/www/OSWW_2018')
parser.add_argument('--options', '-opt', type=str, help='add combine options, e.g. --freezeNuisanceGroups theory', default='')


args = parser.parse_args()

print '--------------------------------------'
print 'Making Impact plot ' + args.output
print '--------------------------------------'

if args.start == 1:
    os.system('combineTool.py -M Impacts -d ' + args.input_file + ' ' + args.options + ' --points 1500 -m 125 -t -1 --setParameters r_vbs=1 -n impact --redefineSignalPOIs=r_vbs --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --doInitialFit')
    os.system("""combineTool.py -M Impacts -d """ + args.input_file + ' ' + args.options + """ -m 125 -t -1 --points 1500 --setParameters r_vbs=1 -n impact --redefineSignalPOIs=r_vbs --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --doFits --job-mode condor --sub-opts='+JobFlavour="workday"'""")

if args.start == 2:
    os.system("combineTool.py -M Impacts -d " + args.input_file + ' ' + args.options + " -m 125 -t -1 --setParameters r_vbs=1 -n impact --redefineSignalPOIs=r_vbs --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 -o " + args.output + ".json")
    os.system("plotImpacts.py -i " + args.output + ".json -o impact_" + args.output)
    os.system("cp  impact_" + args.output + ".pdf " + args.web_dir_path)
