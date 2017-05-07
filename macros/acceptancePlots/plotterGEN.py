import tdrstyle 
import ROOT,os, sys
from ROOT import TFile, TH1,THStack,TCanvas,TLegend
from ROOT import *
#http://www.programcreek.com/python/example/55410/ROOT.TH1
# https://root-forum.cern.ch/t/iteration-over-ame = 'leadingMuonPt' -directory/13615/2
#https://root-forum.cern.ch/t/drawing-histograms-with-for-loop-with-python/23837/4

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.gStyle.SetOptStat(0)



#print " Starting "
ROOT.gROOT.Reset()
canv1 = ROOT.TCanvas("c1","c1",200,10,1050,750);
canv1.cd()
legend = ROOT.TLegend(0.90-.38,0.7,0.9,0.9);
draw_opts = "histSAME e1"
i = 0
#hName = 'leadingMuonPt'
#hName = 'trailingMuonPt'
#hName = 'leadingMuonPtMinSel'
#hName = 'trailingMuonPtMinSel'
hName = 'gammaPt'
DName = ["acceptance_5_5_12", "acceptance_6_2_12","acceptance_8_8_12","acceptance_12_8_12","acceptance_20_12_12","acceptance_20_2_22","acceptance_20_4_22","acceptance_20_5_22"]
#DName = ["acceptance_2_2_23"]
yname = "Nr Events"
xname = "pT [GeV/c]"


f = TFile.Open("root://cmsxrootd.fnal.gov///store/user/edacosta/ZToJPsiGamma-TuneCUETP8M1_13TeV-pythia8/crab_ZToJPsiGamma_GENANA/170422_231855/0000/acceptance_3.root", "READ")
dirlist = f.GetListOfKeys()
#print "Dir list ",dirlist
iter = dirlist.MakeIterator()
key = iter.Next()
#print "key list ", key.GetClassName() 
dirs = {}
td = None

while key:
    if key.GetClassName() == 'TDirectoryFile':
        td = key.ReadObj()
        dirName = td.GetName()
        if dirName in DName :
               i += 1
               #print "found directory", dirName
               dirs[dirName] = td
               hists = key.ReadObj().GetListOfKeys()
               for hist in hists :
                    histogram = hist.ReadObj()
                    histogram.SetDirectory(0) # <----- HERE
                    if hist.GetName() == hName :
                       #print hist.GetName(), i
                       #print histogram
                       histogram.SetMarkerColor(i+1)
                       histogram.SetLineColor(i+1)  
                       histogram.SetMarkerStyle(i+1)
                       #histogram.SetLineWidth(2)
                       histogram.SetLineStyle(i+1)
                       histogram.GetXaxis().SetTitle(xname)
                       histogram.GetYaxis().SetTitle(yname)
                       #histogram.Draw(draw_opts)
                       #legend.AddEntry(histogram,dirName,"lp")  
#                       if i != 1:    
#                            draw_opts += "SAME"
                       histogram.Draw(draw_opts)
                       legend.AddEntry(histogram,dirName,"lp")
    key = iter.Next()
#print "Drawing"
legend.Draw()
canv1.Update()
canv1.SetLogy()
#print "Saving"

#os.system("mkdir genPlots")
canv1.SaveAs('genPlots/'+'%s_file3.root' % hName)
canv1.SaveAs('genPlots/'+'%s_file3.C' % hName) 
canv1.SaveAs('genPlots/'+'%s_file3.png' % hName)
canv1.SaveAs('genPlots/'+'%s_file3.pdf' % hName) 
