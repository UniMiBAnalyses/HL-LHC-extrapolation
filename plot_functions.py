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


def energy_resolution(x_lumi, y_upgrade, y_no_upgrade, y_max, logo='CMS'):

    canv = ROOT.TCanvas('energy_resolution', 'energy_resolution')
    canv.SetGrid()
    ROOT.gStyle.SetOptStat(6)
    ROOT.gStyle.SetOptFit(0000)
    pads = plot.OnePad()
    mg = ROOT.TMultiGraph()
    
    legend_l = 0.7
    legend = ROOT.TLegend(0.3, legend_l, 0.83, 0.85)#, '', 'NBNDC')
    legend.SetTextSize(0.033)
    
    #grafico
    g_upgrade = ROOT.TGraph(len(x_lumi)-1, np.array(x_lumi)[1:], np.array(y_upgrade))
    g_no_upgrade = ROOT.TGraph(len(x_lumi)-1, np.array(x_lumi)[1:], np.array(y_no_upgrade)[1:])
    g_base = ROOT.TGraph(2, np.array(x_lumi)[:1], np.array(y_no_upgrade)[:1])
    
    g_upgrade.SetMarkerStyle(8)
    g_upgrade.SetMarkerColor(8)
    g_upgrade.SetFillColor(8)
    g_no_upgrade.SetMarkerStyle(8)
    g_no_upgrade.SetMarkerColor(2)
    g_no_upgrade.SetFillColor(2)
    g_base.SetMarkerStyle(8)
    g_base.SetMarkerColor(4)
    
    pol2_up = ROOT.TF1("f_up", "[1] * x + [0]", 300, 4500)
    pol2_up.SetLineColor(8)
    pol2_up.SetLineWidth(2)
    
    pol2_noup = ROOT.TF1("f_noup", "[1] * x + [0]", 300, 4500)
    pol2_noup.SetLineColor(2)
    pol2_noup.SetLineWidth(2)
    
    print ("Fitting...")
    g_upgrade.Fit('f_up', 'RO')
    
    g_no_upgrade.Fit('f_noup', 'RO')
    
    pol2_noup.SetLineColor(2)
    pol2_noup.SetLineWidth(1)
    
    #scrivo sulla canvas
    mg.Add(g_base)
    mg.Add(g_upgrade)
    mg.Add(g_no_upgrade)
    
    legend.AddEntry(g_no_upgrade, 'ECAL w/o upgrade', 'lp')
    legend.AddEntry(g_upgrade, 'ECAL w/ upgrade', 'lp')
    legend.AddEntry(g_base, 'ECAL Phase 1', 'lp')
    
    #draw all
    mg.Draw('AP0')
    
    axishist = plot.GetAxisHist(pads[0])
    axishist.SetMaximum(y_max)
    axishist.GetYaxis().SetTitle("Energy resolution (%)")
    axishist.GetXaxis().SetTitle("Luminosity [fb^{-1}]")
    
    #plot.DrawCMSLogo(pads[0], logo, '', 11, 0.045, 0.035, 1.2,  cmsTextSize = 1.)
    legend.Draw()
    
    pads[0].SetGrid()
    canv.Print('.pdf')
    canv.Print('.png')

    return pol2_up, pol2_noup


def electron_efficiency(x_lumi, y_upgrade, y_no_upgrade, y_max, logo='CMS'):

    canv = ROOT.TCanvas('electron_eff', 'electron_eff')
    canv.SetGrid()
    ROOT.gStyle.SetOptStat(6)
    ROOT.gStyle.SetOptFit(0000)
    pads = plot.OnePad()
    mg = ROOT.TMultiGraph()
    
    legend_l = 0.7
    legend = ROOT.TLegend(0.3, legend_l, 0.83, 0.85, '', 'NBNDC')
    legend.SetTextSize(0.033)
    
    #grafico
    g_upgrade = ROOT.TGraph(len(x_lumi), np.array(x_lumi), np.array(y_upgrade))
    g_no_upgrade = ROOT.TGraph(len(x_lumi), np.array(x_lumi), np.array(y_no_upgrade))
    g_color = ROOT.TGraph(2, np.array(x_lumi)[:1], np.array([0.94]))# np.array(y_no_upgrade)[:1])
    
    g_upgrade.SetMarkerStyle(8)
    g_upgrade.SetMarkerColor(8)
    g_no_upgrade.SetMarkerStyle(8)
    g_no_upgrade.SetMarkerColor(2)
    g_color.SetMarkerStyle(8)
    g_color.SetMarkerColor(4)
    
    pol1_up = ROOT.TF1("f_up", "[1] * x + [0]", 300, 4500)
    pol1_up.SetLineColor(8)
    pol1_up.SetLineWidth(2)
    
    pol1_base = ROOT.TF1("f_base", "([1] * x + [0])", 150, 300)
    pol1_base.SetLineColor(4)
    pol1_base.SetLineWidth(2)

    pol1_noup = ROOT.TF1("f_noup", "([0] * x + [1])", 300, 4500)
    pol1_noup.SetLineColor(2)
    pol1_noup.SetLineWidth(2)
    
    print ("Fitting...")
    g_upgrade.Fit('f_up', 'RO')
    
    g_no_upgrade.Fit('f_noup', 'RO')
    #g_no_upgrade.Fit('f_base', 'RO')
    
    #scrivo sulla canvas
    mg.Add(g_upgrade)
    mg.Add(g_no_upgrade)
    mg.Add(g_color)
    
    legend.AddEntry(g_no_upgrade, 'ECAL w/o upgrade', 'lp')
    legend.AddEntry(g_upgrade, 'ECAL w/ upgrade', 'lp')
    legend.AddEntry(g_color, 'ECAL Phase 1', 'lp')
    
    #draw all
    mg.Draw('AP0')
    
    axishist = plot.GetAxisHist(pads[0])
    axishist.SetMaximum(y_max)
    axishist.SetMinimum(0.7)
    axishist.GetYaxis().SetTitle("Efficiency ")
    axishist.GetXaxis().SetTitle("Luminosity [fb^{-1}]")
    
    #plot.DrawCMSLogo(pads[0], logo, '', 11, 0.045, 0.035, 1.2,  cmsTextSize = 1.)
    legend.Draw()
   
    pads[0].SetGrid() 
    canv.Print('.pdf')
    canv.Print('.png')

    return pol1_up, pol1_noup
