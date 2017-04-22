
def replaceStringAll(pathFragment):
  index = 1
  returnString = ""
  for eg in range(12,41):
  # for eg in range(1):
    for lmu in range(21):
      for tmu in range(21):
        if (lmu >= tmu >= 2) and ((4 <= (lmu - tmu) <= 10) or lmu == tmu): 
          lmuStr = str(lmu)
          tmuStr = str(tmu)
          indexStr = str(index)
          egStr = str(eg)
          # print "INDEX: "+indexStr
          # print "L_mu: "+lmuStr
          # print "T_mu: "+str(tmu)
          # print "EG: "+egStr
          # print "###########################"
          index += 1 
          newPathFragment = pathFragment
          newPathFragment = newPathFragment.replace("@@LMU@@", str(lmu))
          newPathFragment = newPathFragment.replace("@@TMU@@", str(tmu))
          newPathFragment = newPathFragment.replace("@@EG@@", str(eg))
          returnString += newPathFragment 
  return returnString
  # return index





## process
processString = """
process.acceptance_@@LMU@@_@@TMU@@_@@EG@@ = cms.EDAnalyzer('AcceptanceStudies',
    verbose = cms.bool(True),
    configName = cms.string("ZJpsiGamma 2017 MC sample"), 
    minMuPt = cms.double(0.0),# in GeV
    maxMuEta = cms.double(2.4),
    minMuonLeadPt = cms.double(@@LMU@@.0),# in GeV
    minMuonTrailPt = cms.double(@@TMU@@.0), # in GeV
    GammaMinPtCut = cms.double(@@EG@@.0),# in GeV
    DeltaRLeadMuPhotonSel = cms.double(1.0),# deltaR>DeltaRLeadMuPhotonSel
    DeltaRTrailPhotonSel  = cms.double(1.0),# deltaR>DeltaRTrailPhotonSel 
    minJPsiMass = cms.double(2.95),# in GeV
    maxJPsiMass = cms.double(3.25)# in GeV
)\n\n
"""


print replaceStringAll(processString)



pathString = " process.acceptance_@@LMU@@_@@TMU@@_@@EG@@ +"
print "process.p = cms.Path(" + replaceStringAll(pathString) + ")"

