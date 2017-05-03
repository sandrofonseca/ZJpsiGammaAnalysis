import tdrstyle 
import ROOT
from ROOT import TFile, TH1,THStack,TCanvas,TLegend
#http://www.programcreek.com/python/example/55410/ROOT.TH1
# https://root-forum.cern.ch/t/iteration-over-a-directory/13615/2
#https://root-forum.cern.ch/t/drawing-histograms-with-for-loop-with-python/23837/4

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.gStyle.SetOptStat(0)



print " Starting "
canv1 = ROOT.TCanvas('can')
canv1.cd()
stack = ROOT.THStack()
draw_opts = "hist e1"
legend = TLegend(.44,.72,.74,.95)
number = 0
hName = 'leadingMuonPt'
DName = ["acceptance_2_2_23",'acceptance_20_20_40']

f = TFile.Open("root://cmsxrootd.fnal.gov///store/user/edacosta/ZToJPsiGamma-TuneCUETP8M1_13TeV-pythia8/crab_ZToJPsiGamma_GENANA/170422_231855/0000/acceptance_1.root", "READ")
dirlist = f.GetListOfKeys()
print "Dir list ",dirlist
iter = dirlist.MakeIterator()
key = iter.Next()
print "key list ", key.GetClassName() 
dirs = {}
td = None
h ={}
while key:
    if key.GetClassName() == 'TDirectoryFile':
	td = key.ReadObj()
	dirName = td.GetName()
	if dirName in DName :
	   number += 1
	   print "found directory", dirName
	   dirs[dirName] = td
	   hists = key.ReadObj().GetListOfKeys()
	   for hist in hists :
		histogram = hist.ReadObj()
	        histogram.SetDirectory(0) # <----- HERE
	        h[hist.GetName()] = histogram
	        if hist.GetName() == hName :
		   print hist.GetName(), number
		   histogram.SetMarkerColor(number)
		   stack.Add(histogram)
		   legend.AddEntry(stack,dirName,"lp")       
    key = iter.Next()
print "Drawing"
stack.Draw("nostack" + draw_opts)
legend.Draw()
canv1.Update()
canv1.SaveAs('%s.png' % hName)


