#!/usr/bin/env python
import ROOT
import math
from functools import partial
import CombineHarvester.CombineTools.plotting as plot
import json
import argparse
import os
import numpy as np
from glob import glob
from datetime import datetime

os.system("cmsenv")
#ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

plot.ModTDRStyle(width=700, l=0.13)
ROOT.gStyle.SetNdivisions(510, "XYZ")
ROOT.gStyle.SetMarkerSize(0.7)
ROOT.gStyle.SetPalette(107)

NAMECOUNTER = 0

def read(scan, param, files, ycut):
    goodfiles = [f for f in files if plot.TFileIsGood(f)]
    limit = plot.MakeTChain(goodfiles, 'limit')
    graph = plot.TGraphFromTree(limit, param, '2*deltaNLL', 'quantileExpected > -1.5')
    graph.SetName(scan)
    graph.Sort()
    plot.RemoveGraphXDuplicates(graph)
    plot.RemoveGraphYAbove(graph, ycut)
    # graph.Print()
    return graph


def Eval(obj, x, params):
    return obj.Eval(x[0])


def BuildScan(scan, param, files, color, yvals, ycut):
    graph = read(scan, param, files, ycut)
    if graph.GetN() <= 1:
        graph.Print()
        raise RuntimeError('Attempting to build %s scan from TGraph with zero or one point (see above)' % files)
    bestfit = None
    for i in xrange(graph.GetN()):
        if graph.GetY()[i] == 0.:
            bestfit = graph.GetX()[i]
    graph.SetMarkerColor(color)
    spline = ROOT.TSpline3("spline3", graph)
    global NAMECOUNTER
    func = ROOT.TF1('splinefn'+str(NAMECOUNTER), partial(Eval, spline), graph.GetX()[0], graph.GetX()[graph.GetN() - 1], 1)
    NAMECOUNTER += 1
    func.SetLineColor(color)
    func.SetLineWidth(3)
    assert(bestfit is not None)
    crossings = {}
    cross_1sig = None
    cross_2sig = None
    other_1sig = []
    other_2sig = []
    val = None
    val_2sig = None
    for yval in yvals:
        crossings[yval] = plot.FindCrossingsWithSpline(graph, func, yval)
        for cr in crossings[yval]:
            cr["contains_bf"] = cr["lo"] <= bestfit and cr["hi"] >= bestfit
    for cr in crossings[yvals[0]]:
        if cr['contains_bf']:
            val = (bestfit, cr['hi'] - bestfit, cr['lo'] - bestfit)
            cross_1sig = cr
        else:
            other_1sig.append(cr)
    if len(yvals) > 1:
        for cr in crossings[yvals[1]]:
            if cr['contains_bf']:
                val_2sig = (bestfit, cr['hi'] - bestfit, cr['lo'] - bestfit)
                cross_2sig = cr
            else:
                other_2sig.append(cr)
    else:
        val_2sig = (0., 0., 0.)
        cross_2sig = cross_1sig
    return val[1], val[2]

parser = argparse.ArgumentParser()

#parser.add_argument('main', help='Main input file for the scan')
parser.add_argument('--y-cut', type=float, default=555., help='Remove points with y > y-cut')
parser.add_argument('--y-max', type=float, default=1.45, help='y-axis maximum')
parser.add_argument('--output', '-o', help='output name without file extension', default='dc_extrapolation_S0-S1')
parser.add_argument('--POI', help='use this parameter of interest', default='r_vbs')
parser.add_argument('--main-label', default='Original Datacard', type=str, help='legend label for the main scan')
parser.add_argument('--main-color', default=1, type=int, help='line and marker color for main scan')
parser.add_argument('--compare', default=0, type=int, help='1 if you want to activate compare section')
parser.add_argument('--others', default=[], nargs='*', help='add secondary scans processed as main: FILE:LABEL:COLOR')
parser.add_argument('--write-legend', nargs='*', help='insert legend for the all datacards')
parser.add_argument('--logo', default='CMS')
parser.add_argument('--logo-sub', default='Internal')
args = parser.parse_args()

print '--------------------------------------'
print  args.output
print '--------------------------------------'

fixed_name = args.POI

#definizione multi canvas e ambiente grafico
canv = ROOT.TCanvas(args.output, args.output)
pads = plot.OnePad()
mg = ROOT.TMultiGraph()

legend_l = 0.835
#legend = ROOT.TLegend(0.5, legend_l, 0.8, 0.925, '', 'NBNDC')
legend = ROOT.TLegend(0.46, 0.74, 0.95, 0.90)
legend.SetTextSize(0.033)

yvals = [1., 4.]
fit_results = {}

while True:
    root_files_S0 = glob("higgsCombineS0_lumi*.POINTS.*MultiDimFit.mH120.root")
    if len(root_files_S0) == 0:
	break
    lumi = root_files_S0[0].rpartition('.POINTS')[0]
    lumi = lumi.partition('_lumi')[2]

    os.system("hadd higgsCombineS0_lumi" + lumi + ".MultiDimFit.mH120.root higgsCombineS0_lumi" + lumi + ".POINTS.*MultiDimFit.mH120.root")
    os.system("mkdir points_lumi" + lumi)
    os.system("mv higgsCombineS0_lumi" + lumi + ".POINTS.*MultiDimFit.mH120.root points_lumi" + lumi)


root_files_S0 = glob("higgsCombineS0_lumi*")
if root_files_S0:
    print '--------------------------------------'
    print 'scenario 0'
    print '--------------------------------------'
    
    #estrazione dei valori di crossing @ 68%
    global x_points 
    x_points = [ ]
    y_points_S0 = [ ]

    for root_file_S0 in root_files_S0:
	lumi = root_file_S0.rpartition('.MultiDimFit')[0]
	lumi = lumi.partition('_lumi')[2]

	print([root_file_S0])
	x_points.append(float(lumi))
	val1, val2 = BuildScan(args.output, args.POI, [root_file_S0], args.main_color+7, yvals, args.y_cut)
	y_points_S0.append([val1, abs(val2)])

    #crea dei vettori contenenti i valori delle incertezze su mu
    S0_plus = [row[0] for row in y_points_S0]
    S0_minus = [row[1] for row in y_points_S0]
    
    #ordinamento
    zipped_lists = zip(x_points, S0_plus)
    sorted_zipped_lists = sorted(zipped_lists)
    sorted_S0_plus = [element for _, element in sorted_zipped_lists]
    
    zipped_lists = zip(x_points, S0_minus)
    sorted_zipped_lists = sorted(zipped_lists)
    sorted_S0_minus = [element for _, element in sorted_zipped_lists]
    
    x_points = sorted(x_points)

    fit_results['S0'] = sorted_S0_plus
    
    #grafico
    g_S0_plus = ROOT.TGraph(len(x_points), np.array(x_points), 1.+np.array(sorted_S0_plus))
    g_S0_minus = ROOT.TGraph(len(x_points), np.array(x_points), 1.-np.array(sorted_S0_minus))
    
    g_S0_plus.SetMarkerStyle(8)
    g_S0_minus.SetMarkerStyle(8)
    g_S0_plus.SetMarkerColor(8)
    g_S0_minus.SetMarkerColor(8)   
 
    legend.AddEntry(g_S0_plus, 'S0', 'LP')
    #scrivo sulla canvas
    mg.Add(g_S0_plus)
    mg.Add(g_S0_minus)



while True:
    root_files_S0plus = glob("higgsCombineS0plus_lumi*.POINTS.*MultiDimFit.mH120.root")
    if len(root_files_S0plus) == 0:
	break
    lumi = root_files_S0plus[0].rpartition('.POINTS')[0]
    lumi = lumi.partition('_lumi')[2]

    os.system("hadd higgsCombineS0plus_lumi" + lumi + ".MultiDimFit.mH120.root higgsCombineS0plus_lumi" + lumi + ".POINTS.*MultiDimFit.mH120.root")
    os.system("mv higgsCombineS0plus_lumi" + lumi + ".POINTS.*MultiDimFit.mH120.root points_lumi" + lumi)

root_files_S0plus = glob("higgsCombineS0plus_lumi*")
if root_files_S0plus:
    print '--------------------------------------'
    print 'scenario 0+'
    print '--------------------------------------'
    
    #estrazione dei valori di crossing @ 68%
    x_points = [ ]
    y_points_S0plus = [ ]

    for root_file_S0plus in root_files_S0plus:
	lumi = root_file_S0plus.rpartition('.MultiDimFit')[0]
	lumi = lumi.partition('_lumi')[2]
	
	x_points.append(float(lumi))
	val1, val2 = BuildScan(args.output, args.POI, [root_file_S0plus], args.main_color+2, yvals, args.y_cut)
	y_points_S0plus.append([val1, abs(val2)])

    #crea dei vettori contenenti i valori delle incertezze su mu
    S0plus_plus = [row[0] for row in y_points_S0plus]
    S0plus_minus = [row[1] for row in y_points_S0plus]
    
    #ordinamento
    zipped_lists = zip(x_points, S0plus_plus)
    sorted_zipped_lists = sorted(zipped_lists)
    sorted_S0plus_plus = [element for _, element in sorted_zipped_lists]
    
    zipped_lists = zip(x_points, S0plus_minus)
    sorted_zipped_lists = sorted(zipped_lists)
    sorted_S0plus_minus = [element for _, element in sorted_zipped_lists]
    
    x_points = sorted(x_points)
    
    fit_results['S0+'] = sorted_S0plus_plus
    
    #grafico
    g_S0plus_plus = ROOT.TGraph(len(x_points), np.array(x_points), 1.+np.array(sorted_S0plus_plus))
    g_S0plus_minus = ROOT.TGraph(len(x_points), np.array(x_points), 1.-np.array(sorted_S0plus_minus))
    
    g_S0plus_plus.SetMarkerStyle(8)
    g_S0plus_minus.SetMarkerStyle(8)
    g_S0plus_plus.SetMarkerColor(32)
    g_S0plus_minus.SetMarkerColor(32)   
    
    #scrivo sulla canvas
    mg.Add(g_S0plus_plus)
    mg.Add(g_S0plus_minus)

    legend.AddEntry(g_S0plus_plus, 'S0+', 'LP')


while True:
    root_files_S1 = glob("higgsCombineS1_lumi*.POINTS.*MultiDimFit.mH120.root")
    if len(root_files_S1) == 0:
	break
    lumi = root_files_S1[0].rpartition('.POINTS')[0]
    lumi = lumi.partition('_lumi')[2]

    os.system("hadd higgsCombineS1_lumi" + lumi + ".MultiDimFit.mH120.root higgsCombineS1_lumi" + lumi + ".POINTS.*MultiDimFit.mH120.root")
    os.system("mv higgsCombineS1_lumi" + lumi + ".POINTS.*MultiDimFit.mH120.root points_lumi" + lumi)


root_files_S1 = glob("higgsCombineS1_lumi*")
if root_files_S1:
    print '--------------------------------------'
    print 'scenario 1'
    print '--------------------------------------'
    
    #estrazione dei valori di crossing @ 68%
    x_points = [ ]
    y_points_S1 = [ ]

    for root_file_S1 in root_files_S1:
	lumi = root_file_S1.rpartition('.MultiDimFit')[0]
	lumi = lumi.partition('_lumi')[2]

	x_points.append(float(lumi))
	val1, val2 = BuildScan(args.output, args.POI, [root_file_S1], args.main_color+45, yvals, args.y_cut)
	y_points_S1.append([val1, abs(val2)])

    #crea dei vettori contenenti i valori delle incertezze su mu
    S1_plus = [row[0] for row in y_points_S1]
    S1_minus = [row[1] for row in y_points_S1]
    
    #ordinamento
    zipped_lists = zip(x_points, S1_plus)
    sorted_zipped_lists = sorted(zipped_lists)
    sorted_S1_plus = [element for _, element in sorted_zipped_lists]
    
    zipped_lists = zip(x_points, S1_minus)
    sorted_zipped_lists = sorted(zipped_lists)
    sorted_S1_minus = [element for _, element in sorted_zipped_lists]
    
    x_points = sorted(x_points)
    
    fit_results['S1'] = sorted_S1_plus
    
    #grafico
    g_S1_plus = ROOT.TGraph(len(x_points), np.array(x_points), 1.+np.array(sorted_S1_plus))
    g_S1_minus = ROOT.TGraph(len(x_points), np.array(x_points), 1.-np.array(sorted_S1_minus))
    
    g_S1_plus.SetMarkerStyle(8)
    g_S1_minus.SetMarkerStyle(8)
    g_S1_plus.SetMarkerColor(42)
    g_S1_minus.SetMarkerColor(42)   

    legend.SetNColumns(2) 
    legend.AddEntry(g_S1_plus, 'S1', 'LP')
    #scrivo sulla canvas
    mg.Add(g_S1_plus)
    mg.Add(g_S1_minus)


while True:
    root_files_S1plus = glob("higgsCombineS1plus_lumi*.POINTS.*MultiDimFit.mH120.root")
    if len(root_files_S1plus) == 0:
	break
    lumi = root_files_S1plus[0].rpartition('.POINTS')[0]
    lumi = lumi.partition('_lumi')[2]

    os.system("hadd higgsCombineS1plus_lumi" + lumi + ".MultiDimFit.mH120.root higgsCombineS1plus_lumi" + lumi + ".POINTS.*MultiDimFit.mH120.root")
    os.system("mv higgsCombineS1plus_lumi" + lumi + ".POINTS.*MultiDimFit.mH120.root points_lumi" + lumi)

root_files_S1plus = glob("higgsCombineS1plus_lumi*")
if root_files_S1plus:
    print '--------------------------------------'
    print 'scenario 1+'
    print '--------------------------------------'
    
    #estrazione dei valori di crossing @ 68%
    x_points = [ ]
    y_points_S1plus = [ ]

    for root_file_S1plus in root_files_S1plus:
	lumi = root_file_S1plus.rpartition('.MultiDimFit')[0]
	lumi = lumi.partition('_lumi')[2]
	
	x_points.append(float(lumi))
	val1, val2 = BuildScan(args.output, args.POI, [root_file_S1plus], args.main_color+1, yvals, args.y_cut)
	y_points_S1plus.append([val1, abs(val2)])

    #crea dei vettori contenenti i valori delle incertezze su mu
    S1plus_plus = [row[0] for row in y_points_S1plus]
    S1plus_minus = [row[1] for row in y_points_S1plus]
    
    #ordinamento
    zipped_lists = zip(x_points, S1plus_plus)
    sorted_zipped_lists = sorted(zipped_lists)
    sorted_S1plus_plus = [element for _, element in sorted_zipped_lists]
    
    zipped_lists = zip(x_points, S1plus_minus)
    sorted_zipped_lists = sorted(zipped_lists)
    sorted_S1plus_minus = [element for _, element in sorted_zipped_lists]
    
    x_points = sorted(x_points)
    
    fit_results['S1+'] = sorted_S1plus_plus
   
    #grafico
    g_S1plus_plus = ROOT.TGraph(len(x_points), np.array(x_points), 1.+np.array(sorted_S1plus_plus))
    g_S1plus_minus = ROOT.TGraph(len(x_points), np.array(x_points), 1.-np.array(sorted_S1plus_minus))
    
    g_S1plus_plus.SetMarkerStyle(8)
    g_S1plus_minus.SetMarkerStyle(8)
    g_S1plus_plus.SetMarkerColor(46)
    g_S1plus_minus.SetMarkerColor(46)   
    
    #scrivo sulla canvas
    mg.Add(g_S1plus_plus)
    mg.Add(g_S1plus_minus)

    legend.AddEntry(g_S1plus_plus, 'S1+', 'LP')
    

while True:
    root_files_S2 = glob("higgsCombineS2_lumi*.POINTS.*MultiDimFit.mH120.root")
    if len(root_files_S2) == 0:
        break
    lumi = root_files_S2[0].rpartition('.POINTS')[0]
    lumi = lumi.partition('_lumi')[2]

    os.system("hadd higgsCombineS2_lumi" + lumi + ".MultiDimFit.mH120.root higgsCombineS2_lumi" + lumi + ".POINTS.*MultiDimFit.mH120.root")
    os.system("mkdir points_lumi" + lumi)
    os.system("mv higgsCombineS2_lumi" + lumi + ".POINTS.*MultiDimFit.mH120.root points_lumi" + lumi)


root_files_S2 = glob("higgsCombineS2_lumi*")
if root_files_S2:
    print '--------------------------------------'
    print 'scenario 2'
    print '--------------------------------------'

    #estrazione dei valori di crossing @ 68%
    x_points = [ ]
    y_points_S2 = [ ]

    for root_file_S2 in root_files_S2:
        lumi = root_file_S2.rpartition('.MultiDimFit')[0]
        lumi = lumi.partition('_lumi')[2]

        print([root_file_S2])
        x_points.append(float(lumi))
        val1, val2 = BuildScan(args.output, args.POI, [root_file_S2], args.main_color+7, yvals, args.y_cut)
        y_points_S2.append([val1, abs(val2)])

    #crea dei vettori contenenti i valori delle incertezze su mu
    S2_plus = [row[0] for row in y_points_S2]
    S2_minus = [row[1] for row in y_points_S2]

    #ordinamento
    zipped_lists = zip(x_points, S2_plus)
    sorted_zipped_lists = sorted(zipped_lists)
    sorted_S2_plus = [element for _, element in sorted_zipped_lists]

    zipped_lists = zip(x_points, S2_minus)
    sorted_zipped_lists = sorted(zipped_lists)
    sorted_S2_minus = [element for _, element in sorted_zipped_lists]

    x_points = sorted(x_points)
    
    fit_results['S2'] = sorted_S2_plus

    #grafico
    g_S2_plus = ROOT.TGraph(len(x_points), np.array(x_points), 1.+np.array(sorted_S2_plus))
    g_S2_minus = ROOT.TGraph(len(x_points), np.array(x_points), 1.-np.array(sorted_S2_minus))

    g_S2_plus.SetMarkerStyle(8)
    g_S2_minus.SetMarkerStyle(8)
    g_S2_plus.SetMarkerColor(7)
    g_S2_minus.SetMarkerColor(7)

    legend.AddEntry(g_S2_plus, 'S2', 'LP')
    #scrivo sulla canvas
    mg.Add(g_S2_plus)
    mg.Add(g_S2_minus)


while True:
    root_files_S2plus = glob("higgsCombineS2plus_lumi*.POINTS.*MultiDimFit.mH120.root")
    if len(root_files_S2plus) == 0:
	break
    lumi = root_files_S2plus[0].rpartition('.POINTS')[0]
    lumi = lumi.partition('_lumi')[2]

    os.system("hadd higgsCombineS2plus_lumi" + lumi + ".MultiDimFit.mH120.root higgsCombineS2plus_lumi" + lumi + ".POINTS.*MultiDimFit.mH120.root")
    os.system("mv higgsCombineS2plus_lumi" + lumi + ".POINTS.*MultiDimFit.mH120.root points_lumi" + lumi)

root_files_S2plus = glob("higgsCombineS2plus_lumi*")
if root_files_S2plus:
    print '--------------------------------------'
    print 'scenario 2+'
    print '--------------------------------------'
    
    #estrazione dei valori di crossing @ 68%
    x_points = [ ]
    y_points_S2plus = [ ]

    for root_file_S2plus in root_files_S2plus:
	lumi = root_file_S2plus.rpartition('.MultiDimFit')[0]
	lumi = lumi.partition('_lumi')[2]
	
	x_points.append(float(lumi))
	val1, val2 = BuildScan(args.output, args.POI, [root_file_S2plus], args.main_color+2, yvals, args.y_cut)
	y_points_S2plus.append([val1, abs(val2)])

    #crea dei vettori contenenti i valori delle incertezze su mu
    S2plus_plus = [row[0] for row in y_points_S2plus]
    S2plus_minus = [row[1] for row in y_points_S2plus]
    
    #ordinamento
    zipped_lists = zip(x_points, S2plus_plus)
    sorted_zipped_lists = sorted(zipped_lists)
    sorted_S2plus_plus = [element for _, element in sorted_zipped_lists]
    
    zipped_lists = zip(x_points, S2plus_minus)
    sorted_zipped_lists = sorted(zipped_lists)
    sorted_S2plus_minus = [element for _, element in sorted_zipped_lists]
    
    x_points = sorted(x_points)
    
    fit_results['S2+'] = sorted_S2plus_plus
    
    #grafico
    g_S2plus_plus = ROOT.TGraph(len(x_points), np.array(x_points), 1.+np.array(sorted_S2plus_plus))
    g_S2plus_minus = ROOT.TGraph(len(x_points), np.array(x_points), 1.-np.array(sorted_S2plus_minus))
    
    g_S2plus_plus.SetMarkerStyle(8)
    g_S2plus_minus.SetMarkerStyle(8)
    g_S2plus_plus.SetMarkerColor(9)
    g_S2plus_minus.SetMarkerColor(9)   
    
    #scrivo sulla canvas
    mg.Add(g_S2plus_plus)
    mg.Add(g_S2plus_minus)

    legend.AddEntry(g_S2plus_plus, 'S2+', 'LP')
    

while True:
    root_files_freeze_all = glob("higgsCombinefreeze_all_lumi*.POINTS.*MultiDimFit.mH120.root")
    if len(root_files_freeze_all) == 0:
	break
    lumi = root_files_freeze_all[0].rpartition('.POINTS')[0]
    lumi = lumi.partition('_lumi')[2]

    os.system("hadd higgsCombinefreeze_all_lumi" + lumi + ".MultiDimFit.mH120.root higgsCombinefreeze_all_lumi" + lumi + ".POINTS.*MultiDimFit.mH120.root")
    os.system("mkdir points_lumi" + lumi)
    os.system("mv higgsCombinefreeze_all_lumi" + lumi + ".POINTS.*MultiDimFit.mH120.root points_lumi" + lumi)

root_files_freeze_all = glob("higgsCombinefreeze_all_lumi*")
if root_files_freeze_all:
    print '--------------------------------------'
    print 'freeze all nuisances'
    print '--------------------------------------'
    
    #estrazione dei valori di crossing @ 68%
    x_points = [ ]
    y_points_freeze_all = [ ]

    for root_file_freeze_all in root_files_freeze_all:
	lumi = root_file_freeze_all.rpartition('.MultiDimFit')[0]
	lumi = lumi.partition('_lumi')[2]
	
	x_points.append(float(lumi))
	val1, val2 = BuildScan(args.output, args.POI, [root_file_freeze_all], args.main_color+2, yvals, args.y_cut)
	y_points_freeze_all.append([val1, abs(val2)])

    #crea dei vettori contenenti i valori delle incertezze su mu
    freeze_all_plus = [row[0] for row in y_points_freeze_all]
    freeze_all_minus = [row[1] for row in y_points_freeze_all]
    
    #ordinamento
    zipped_lists = zip(x_points, freeze_all_plus)
    sorted_zipped_lists = sorted(zipped_lists)
    sorted_freeze_all_plus = [element for _, element in sorted_zipped_lists]
    
    zipped_lists = zip(x_points, freeze_all_minus)
    sorted_zipped_lists = sorted(zipped_lists)
    sorted_freeze_all_minus = [element for _, element in sorted_zipped_lists]
    
    x_points = sorted(x_points)
    
    fit_results['freeze_all'] = sorted_freeze_all_plus
    
    #grafico
    g_freeze_all_plus = ROOT.TGraph(len(x_points), np.array(x_points), 1.+np.array(sorted_freeze_all_plus))
    g_freeze_all_minus = ROOT.TGraph(len(x_points), np.array(x_points), 1.-np.array(sorted_freeze_all_minus))
    
    g_freeze_all_plus.SetMarkerStyle(8)
    g_freeze_all_minus.SetMarkerStyle(8)
    g_freeze_all_plus.SetMarkerColor(16)
    g_freeze_all_minus.SetMarkerColor(16)   
    
    #scrivo sulla canvas
    mg.Add(g_freeze_all_plus)
    mg.Add(g_freeze_all_minus)

    legend.AddEntry(g_freeze_all_plus, 'only stat', 'LP')
    
    #definisco la griglia di abbellimento del plot :)
    g_shadow = ROOT.TGraph(2*len(x_points))
    for i in range(len(x_points)):
             g_shadow.SetPoint(i,np.array(x_points)[i],1.+np.array(sorted_freeze_all_plus)[i])
             g_shadow.SetPoint(len(x_points)+i,np.array(x_points)[len(x_points)-i-1],1.-np.array(sorted_freeze_all_plus)[len(x_points)-i-1])

    g_shadow.SetFillStyle(3013)
    g_shadow.SetFillColor(16)
    g_shadow.Draw('F')
    
    #draw all
    mg.Draw('ALP PLC PMC')
    g_shadow.Draw('F')
    g_freeze_all_plus.Draw('l')
    g_freeze_all_minus.Draw('l')

axishist = plot.GetAxisHist(pads[0])

axishist.SetMaximum(float(args.y_max))
axishist.GetYaxis().SetTitle("#mu error")
axishist.GetXaxis().SetTitle("Luminosity [fb^{-1}]")

#traccio la linea su 1
line = ROOT.TLine()
line.SetLineColor(2)
line.SetLineStyle(7)
line.SetLineWidth(2)
plot.DrawHorizontalLine(pads[0], line, 1)

plot.DrawCMSLogo(pads[0], args.logo, args.logo_sub, 11, 0.045, 0.035, 1.2,  cmsTextSize = 1.)
legend.Draw()

#canv.Print('.pdf')
canv.Print('.png')


if args.compare:
    canv_compare = ROOT.TCanvas(args.output + 'compare', args.output + 'comparison')
    pads_compare = plot.OnePad()
    mg_compare = ROOT.TMultiGraph()

    legend_compare = ROOT.TLegend(0.62, 0.50, 0.92, 0.65)
    legend_compare.SetTextSize(0.033)

    valuesS0 = fit_results['S0']
    
    valuesS0plus = fit_results['S0+']
    valuesS2 = fit_results['S2']
    valuesS2plus = fit_results['S2+']

    g_compare_S0S0plus = ROOT.TGraph(len(x_points), np.array(x_points), np.divide(np.array(valuesS0), np.array(valuesS0plus)))
    g_compare_S2S2plus = ROOT.TGraph(len(x_points), np.array(x_points), np.divide(np.array(valuesS2), np.array(valuesS2plus)))
    g_compare_S0S2 = ROOT.TGraph(len(x_points), np.array(x_points), np.divide(np.array(valuesS0), np.array(valuesS2)))

    g_compare_S0S0plus.SetMarkerStyle(8)
    g_compare_S2S2plus.SetMarkerStyle(8)
    g_compare_S0S2.SetMarkerStyle(8)
    g_compare_S0S0plus.SetMarkerColor(2)
    g_compare_S0S2.SetMarkerColor(3)   
    g_compare_S2S2plus.SetMarkerColor(4)
    
    mg_compare.Add(g_compare_S0S0plus)
    mg_compare.Add(g_compare_S2S2plus)
    mg_compare.Add(g_compare_S0S2)

    legend_compare.AddEntry(g_compare_S0S0plus, 'comparison S0-S0+', 'LP')
    legend_compare.AddEntry(g_compare_S2S2plus, 'comparison S2-S2+', 'LP')
    legend_compare.AddEntry(g_compare_S0S2, 'comparison S0-S2', 'LP')
    mg_compare.Draw('ALP PLC PMC')

    axi = plot.GetAxisHist(pads_compare[0])
    #axi.SetMaximum(float())
    axi.GetYaxis().SetTitle("#Delta#mu error")
    axi.GetXaxis().SetTitle("Luminosity [fb^{-1}]")
 
    plot.DrawCMSLogo(pads_compare[0], args.logo, args.logo_sub, 11, 0.045, 0.035, 1.2,  cmsTextSize = 1.)
    legend_compare.Draw()
    
    #canv_compare.cd().SetLogy()   
    canv_compare.Print('.png')
    canv_compare.Print('.root')


# mk folder
folder_name = datetime.now().strftime('%d-%m-%Y') + "--h" + str(datetime.now().hour)
os.mkdir(folder_name)
os.system('mv S* ' + folder_name)
os.system('mv higgsCombineS* ' + folder_name)
os.system('mv condor_S* ' + folder_name)
os.system('mv combined_S* ' + folder_name)
os.system('mv *freeze_all* ' + folder_name)
os.system('mv points_lumi*/ ' + folder_name)
