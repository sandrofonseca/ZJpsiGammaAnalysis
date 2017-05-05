import tdrstyle 
import ROOT
from ROOT import TFile, TH1,THStack,TCanvas,TLegend
from ROOT import *
#http://www.programcreek.com/python/example/55410/ROOT.TH1
# https://root-forum.cern.ch/t/iteration-over-a-directory/13615/2
#https://root-forum.cern.ch/t/drawing-histograms-with-for-loop-with-python/23837/4

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.gStyle.SetOptStat(0)



print " Starting "
ROOT.gROOT.Reset()
canv1 = ROOT.TCanvas("c1","c1",200,10,1050,750);
canv1.cd()
legend = ROOT.TLegend(0.90-.38,0.7,0.9,0.9);
draw_opts = "hist e1"
i = 0
hName = 'leadingMuonPt'
DName = ["acceptance_2_2_23",'acceptance_20_20_40']
#DName = ["acceptance_2_2_23"]

f = TFile.Open("root://cmsxrootd.fnal.gov///store/user/edacosta/ZToJPsiGamma-TuneCUETP8M1_13TeV-pythia8/crab_ZToJPsiGamma_GENANA/170422_231855/0000/acceptance_1.root", "READ")
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
	   print "found directory", dirName
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
		   histogram.Draw(draw_opts)
		   legend.AddEntry(histogram,dirName,"lp") 	
	           if i != 1:
			draw_opts += " same"			 
    key = iter.Next()
print "Drawing"
legend.Draw()
canv1.Update()
print "Saving"
canv1.SaveAs('%s.png' % hName)
canv1.SaveAs('%s.pdf' % hName)

