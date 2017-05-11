import tdrstyle 
import ROOT,os, sys,itertools
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
hTotal = 'h_nEvts_ZJpsiGamma 2017 MC sample'
hGEN = 'h_nEvtsGEN_ZJpsiGamma 2017 MC sample'
#hName = 'leadingMuonPt'
#hName = 'trailingMuonPt'
#hName = 'leadingMuonPtMinSel'
#hName = 'trailingMuonPtMinSel'
hName = 'hefficiency'
#DName = ["acceptance_5_5_12", "acceptance_6_2_12","acceptance_8_8_12","acceptance_12_8_12","acceptance_20_12_12","acceptance_20_2_22","acceptance_20_4_22","acceptance_20_5_22"]
DName = ["acceptance_5_5_12"]
yname = "Eff"
#xname = "pT [GeV/c]"
xname = ""

f = TFile.Open("root://cmseos.fnal.gov//store/user/sfonseca/crab_ZToJPsiGamma_GENANA_ElizaEOS/acceptance_1_2_3.root", "READ")
dirlist = f.GetListOfKeys()
#print "Dir list ",dirlist
iter = dirlist.MakeIterator()
key = iter.Next()
#print "key list ", key.GetClassName() 
dirs = {}
td = None
h_num_den = ROOT.TH1F('h_num_den' , 'h_num_den' ,  1, 0., 1.)
h_num_tot = ROOT.TH1F('h_num_tot' , 'h_num_tot' ,  1, 0., 1.)

while key:
    if key.GetClassName() == 'TDirectoryFile':
        td = key.ReadObj()
        dirName = td.GetName()
        if dirName in DName :
               i += 1
               print "found directory", dirName
               dirs[dirName] = td
               htotkey = key.ReadObj().GetListOfKeys()
	       hgenkey = key.ReadObj().GetListOfKeys()	
               for htotal,hgen in itertools.product(htotkey,hgenkey) :
                    histogramT = htotal.ReadObj()
                    histogramG = hgen.ReadObj()
	            histogramT.SetDirectory(0) # <----- HERE
		    histogramG.SetDirectory(0) # <----- HERE
		    #print histogram	
                    if (htotal.GetName() == hTotal and hgen.GetName() == hGEN)  :
                       #print hist.GetName(), i
                       print "Total Evts: ", histogramT, "Gen Evts:", histogramG
	               h_num_tot = histogramT.Clone() 
		       h_num_den = histogramG.Clone()	
		       print "Total:" , h_num_tot.GetBinContent(1) 
		       print "Gen: ", h_num_den.GetBinContent(1)       
		       eff = ROOT.TEfficiency(h_num_tot, h_num_den)	
                       eff.SetMarkerColor(i+1)
                       eff.SetLineColor(i+1)  
                       eff.SetMarkerStyle(i+1)
                       #histogram.SetLineWidth(2)
                       eff.SetLineStyle(i+1)
                       #eff.GetXaxis().SetTitle(xname)
                       #eff.GetYaxis().SetTitle(yname)
                       #histogram.Draw(draw_opts)
                       #legend.AddEntry(histogram,dirName,"lp")  
#                       if i != 1:    
#                            draw_opts += "SAME"
                     #  eff.Draw(draw_opts)
                       legend.AddEntry(eff,dirName,"lp")
		        	
                       			
    key = iter.Next()
#print "Drawing"
legend.Draw()
canv1.Update()
#canv1.SetLogy()
#print "Saving"

#os.system("mkdir genPlots")
canv1.SaveAs('genPlots/'+'%s_file.root' % hName)
canv1.SaveAs('genPlots/'+'%s_file.C' % hName) 
canv1.SaveAs('genPlots/'+'%s_file.png' % hName)
canv1.SaveAs('genPlots/'+'%s_file.pdf' % hName) 
