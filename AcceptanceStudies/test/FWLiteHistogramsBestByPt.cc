#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"


int main(int argc, char* argv[]) 
{
// define what muon you are using; this is necessary as FWLite is not 
// capable of reading edm::Views
  using reco::GenParticle;
  using pat::Muon;
  using pat::Photon;  

// load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  FWLiteEnabler::enable();

// initialize command line parser
  optutl::CommandLineParser parser ("Analyze FWLite Histograms");

// set defaults
  parser.integerValue ("maxEvents"  ) = 1000;
  parser.integerValue ("outputEvery") =   10;
  parser.stringValue  ("outputFile" ) = "analyzeFWLiteHistograms.root";

// parse arguments
  parser.parseArguments (argc, argv);
  int maxEvents_ = parser.integerValue("maxEvents");
  unsigned int outputEvery_ = parser.integerValue("outputEvery");
  std::string outputFile_ = parser.stringValue("outputFile");
  std::vector<std::string> inputFiles_ = parser.stringVector("inputFiles");

// book a set of histograms
  fwlite::TFileService fs = fwlite::TFileService(outputFile_.c_str());
  //std::cout << "Opening file: " << ievt << std::endl;
  TFileDirectory dir = fs.mkdir("plots");

//genMuons
  TH1F* leadingMuonPt_  = dir.make<TH1F>("leadingMuonPt"  , "mu1_pt"  ,   100,   0., 200.);
  TH1F* leadingMuonEta_ = dir.make<TH1F>("leadingMuonEta" , "mu1_eta" ,   100,  -5.,   5.);
  //TH2F(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup)
  //_histRvEtaLow = dir.make<TH2F>("histRvEtaLow","R v. #eta",32,0.9,2.5,26,0,6.5);
  TH2F* leadingMuonPtVsleadingMuonEta_ = dir.make<TH2F>("leadingMuonPtVsleadingMuonEta", "Pt Vs Ets",100,0.,100.,100,-5., 5.);
  TH1F* leadingMuonPhi_ = dir.make<TH1F>("leadingMuonPhi" , "mu1_phi" ,   100,  -4.,   4.);  
  TH1F* secondMuonPt_  = dir.make<TH1F>("secondMuonPt"  , "mu2_pt"  ,   100,   0., 200.);
  TH1F* secondMuonEta_ = dir.make<TH1F>("secondMuonEta" , "mu2_eta" ,   100,  -5.,   5.);
  TH1F* secondMuonPhi_ = dir.make<TH1F>("secondMuonPhi" , "mu2_phi" ,   100,  -4.,   4.);  
  TH1F* mumuMass_= dir.make<TH1F>("mumuMass", "Jpsi_mass",    100,  2.8,  3.2);

  //MinSel
  TH1F* leadingMuonPtMinSel_  = dir.make<TH1F>("leadingMuonPtMinSel"  , "mu1_pt"  ,   100,   0., 200.);
  TH1F* leadingMuonEtaMinSel_ = dir.make<TH1F>("leadingMuonEtaMinSel" , "mu1_eta" ,   100,  -5.,   5.);
  TH1F* leadingMuonPhiMinSel_ = dir.make<TH1F>("leadingMuonPhiMinSel" , "mu1_phi" ,   100,  -4.,   4.);
  //TH2F(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup) 
  TH2F* leadingMuonPtVsleadingMuonEtaMinSel_ = dir.make<TH2F>("leadingMuonPtVsleadingMuonEtaMinSel", "mu1_pt x mu1_eta" ,100,0.,100.,100,-5., 5.); 
  TH1F* secondMuonPtMinSel_  = dir.make<TH1F>("secondMuonPtMinSel"  , "mu2_pt"  ,   100,   0., 200.);
  TH1F* secondMuonEtaMinSel_ = dir.make<TH1F>("secondMuonEtaMinSel" , "mu2_eta" ,   100,  -5.,   5.);
  TH1F* secondMuonPhiMinSel_ = dir.make<TH1F>("secondMuonPhiMinSel" , "mu2_phi" ,   100,  -4.,   4.);  
  TH1F* mumuMassMinSel_= dir.make<TH1F>("mumuMassMinSel", "Jpsi_mass",    100,  2.8,  3.2);

//genGammas
  TH1F* gammaPt_  = dir.make<TH1F>("gammaPt"  , "pt"  ,   100,   0., 300.);
  TH1F* gammaEta_ = dir.make<TH1F>("gammaEta" , "eta" ,   100,  -5.,   5.);
  TH1F* gammaPhi_ = dir.make<TH1F>("gammaPhi" , "phi" ,   100,  -4.,   4.);
  TH1F* mumugammaMass_= dir.make<TH1F>("zMass", "mass",    100,  65.,  115.);

//MinSel
  TH1F* gammaPtMinSel_  = dir.make<TH1F>("gammaPtMinSel"  , "pt"  ,   100,   0., 300.);
  TH1F* gammaEtaMinSel_ = dir.make<TH1F>("gammaEtaMinSel" , "eta" ,   100,  -5.,   5.);
  TH1F* gammaPhiMinSel_ = dir.make<TH1F>("gammaPhiMinSel" , "phi" ,   100,  -4.,   4.);
  TH1F* mumugammaMassMinSel_= dir.make<TH1F>("zMassMinSel", "mass",    100,  65.,  115.);
/*  //slimmedMuons
  TH1F* slimmedMuonPt_  = dir.make<TH1F>("slimmedMuonPt"  , "pt"  ,   100,   0., 300.);
  TH1F* slimmedMuonEta_ = dir.make<TH1F>("slimmedMuonEta" , "eta" ,   100,  -3.,   3.);
  TH1F* slimmedMuonPhi_ = dir.make<TH1F>("slimmedMuonPhi" , "phi" ,   100,  -5.,   5.);  
  TH1F* slimmedmumuMass_= dir.make<TH1F>("slimmedmumuMass", "mass",    90,  2.,  4.);
  //slimmedGammas
  TH1F* slimmedgammaPt_  = dir.make<TH1F>("slimmedgammaPt"  , "pt"  ,   100,   0., 300.);
  TH1F* slimmedgammaEta_ = dir.make<TH1F>("slimmedgammaEta" , "eta" ,   100,  -3.,   3.);
  TH1F* slimmedgammaPhi_ = dir.make<TH1F>("slimmedgammaPhi" , "phi" ,   100,  -5.,   5.);
  TH1F* slimmedmumugammaMass_= dir.make<TH1F>("slimmedmumugammaMass", "mass",    100,  85.,  95.);*/



  // loop the events
  int ievt=0; 
  int nevt_tot=0; 
  int nevt_mumu=0;
  int nevt_mumugamma=0;
  int nevt_minsel=0;
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
// open input file (can be located on castor)
    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
    std::cout << "Opening file: " << inputFiles_[iFile].c_str() << std::endl;
    if( inFile ){
      fwlite::Event ev(inFile);
      for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
        edm::EventBase const & event = ev;
// break loop if maximal number of events is reached 
        if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
// simple event counter
        if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false) 
          std::cout << "  processing event: " << ievt << std::endl;

    ++nevt_tot; 
//################################ (begin) genMuons ################################
// Handle to the muon collection
        edm::Handle<std::vector<GenParticle> > muons;
        event.getByLabel(std::string("genParticles"), muons);
        //event.getByLabel(std::string("prunedGenParticles"), muons);

        float LeadingMuonPt = 0.0;
        reco::GenParticle bestMu1;
        reco::GenParticle bestMu2;
// loop muon collection and fill histograms
        for(std::vector<GenParticle>::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1){
//muonPt_ ->Fill( mu1->pt () );
//muonEta_->Fill( mu1->eta() );
//muonPhi_->Fill( mu1->phi() );
          if(abs(mu1->pdgId())== 13 && mu1->status()== 1){
//for(std::vector<GenParticle>::const_iterator mu2=muons->begin(); mu2!=muons->end(); ++mu2){
            if(mu1->pt() >= LeadingMuonPt){
              bestMu1 = *mu1->clone();
              LeadingMuonPt = mu1->pt();
            }
          }
        }

        float SecondMuonPt = 0.0;
        for(std::vector<GenParticle>::const_iterator mu2=muons->begin(); mu2!=muons->end(); ++mu2){
//muonPt_ ->Fill( mu1->pt () );
//muonEta_->Fill( mu1->eta() );
//muonPhi_->Fill( mu1->phi() );
          if(abs(mu2->pdgId())== 13 && mu2->status()== 1){
//for(std::vector<GenParticle>::const_iterator mu2=muons->begin(); mu2!=muons->end(); ++mu2){
            if(mu2->pt() >= SecondMuonPt && mu2->pt() < LeadingMuonPt){
              bestMu2 = *mu2->clone();
              SecondMuonPt = mu2->pt();
            }
          }
        }        

//std::cout << "  oi2 " << std::endl;
        ++nevt_mumu;
        mumuMass_->Fill( (bestMu1.p4()+bestMu2.p4()).mass() );
        leadingMuonPt_ ->Fill( bestMu1.pt () );
        leadingMuonEta_->Fill( bestMu1.eta() );
        leadingMuonPhi_->Fill( bestMu1.phi() );
        secondMuonPt_ ->Fill( bestMu2.pt () );
        secondMuonEta_->Fill( bestMu2.eta() );
        secondMuonPhi_->Fill( bestMu2.phi() );
        leadingMuonPtVsleadingMuonEta_ ->Fill( bestMu1.pt (),bestMu1.eta() );
//################################ (end) genMuons ################################



//################################ (begin) genGammas ################################  
// Handle to the gamma collection
        edm::Handle<std::vector<GenParticle> > photon;
        event.getByLabel(std::string("genParticles"), photon);
        //event.getByLabel(std::string("prunedGenParticles"), photon);
        

        reco::GenParticle bestGamma;
// loop gamma collection and fill histograms
        float LeadingGammaPt = 0.0;
        for(std::vector<GenParticle>::const_iterator gamma=photon->begin(); gamma!=photon->end(); ++gamma){
//muonPt_ ->Fill( mu1->pt () );
//muonEta_->Fill( mu1->eta() );
//muonPhi_->Fill( mu1->phi() );
          if(abs(gamma->pdgId())== 22 && gamma->status()== 1){
            if(gamma->pt() >= LeadingGammaPt){
              bestGamma = *gamma->clone();
              LeadingGammaPt = gamma->pt();
//std::cout << zMassDif << std::endl;
            }
          }
        }

//std::cout << "  oi2 " << std::endl;
        mumugammaMass_->Fill( (bestMu1.p4()+bestMu2.p4()+bestGamma.p4()).mass() );
        gammaPt_ ->Fill( bestGamma.pt () );
        gammaEta_->Fill( bestGamma.eta() );
        gammaPhi_->Fill( bestGamma.phi() );
        ++nevt_mumugamma;
//Some cuts asked by Andrey
////pT_mu1>20 Gev, pT_gamma > 20 GeV
       if(bestMu1.pt ()>20.0 && bestGamma.pt ()> 20.0 ){
        ++nevt_minsel;
        mumugammaMassMinSel_->Fill( (bestMu1.p4()+bestMu2.p4()+bestGamma.p4()).mass() );
        gammaPtMinSel_ ->Fill( bestGamma.pt () );
        gammaEtaMinSel_->Fill( bestGamma.eta() );
        gammaPhiMinSel_->Fill( bestGamma.phi() );
        mumuMassMinSel_->Fill( (bestMu1.p4()+bestMu2.p4()).mass() );
        leadingMuonPtMinSel_ ->Fill( bestMu1.pt () );
        leadingMuonPtVsleadingMuonEtaMinSel_ ->Fill( bestMu1.pt (),bestMu1.eta() );
        leadingMuonEtaMinSel_->Fill( bestMu1.eta() );
        leadingMuonPhiMinSel_->Fill( bestMu1.phi() );
        secondMuonPtMinSel_ ->Fill( bestMu2.pt () );
        secondMuonEtaMinSel_->Fill( bestMu2.eta() );
        secondMuonPhiMinSel_->Fill( bestMu2.phi() );
}
//################################ (end) genGammas ################################


// //################################ (begin) slimmedMuons ################################
//   // Handle to the muon collection
//         edm::Handle<std::vector<Muon> > slimmedMuons;
//         event.getByLabel(std::string("slimmedMuons"), slimmedMuons);


//         pat::Muon bestslimmedMu1;
//         pat::Muon bestslimmedMu2;

//   jpsiMassDif = 99999999999.;
//   // loop muon collection and fill histograms
//         for(std::vector<Muon>::const_iterator mu1=slimmedMuons->begin(); mu1!=slimmedMuons->end(); ++mu1){
//           //muonPt_ ->Fill( mu1->pt () );
//           //muonEta_->Fill( mu1->eta() );
//           //muonPhi_->Fill( mu1->phi() );
//             for(std::vector<Muon>::const_iterator mu2=slimmedMuons->begin(); mu2!=slimmedMuons->end(); ++mu2){
//   if(mu2>mu1){ // prevent double conting
//   if( mu1->charge()*mu2->charge()<0 ){ // check only muon pairs of unequal charge 
//     if(mu1->pt() > 10.0 && mu2->pt() > 10.0 ){
//       if(abs((mu1->p4()+mu2->p4()).mass() - jpsiMass) < jpsiMassDif){
//         bestslimmedMu1 = *mu1->clone();
//         bestslimmedMu2 = *mu2->clone();
//         jpsiMassDif = abs((mu1->p4()+mu2->p4()).mass() - jpsiMass);
//         //std::cout << "  oi1 " << std::endl;
//       }
//     }    
//       //mumuMass_->Fill( (mu1->p4()+mu2->p4()).mass() );

// /*  for(std::vector<GenParticle>::const_iterator mu3=muons->begin(); mu3!=muons->end(); ++mu3){
//     if(abs(mu3->pdgId())== 22 && mu3->status()== 1 && mu3->mother()->pdgId() == 25){
//       mumugammaMass_->Fill( (mu1->p4()+mu2->p4()+mu3->p4()).mass() );
//     }
//   }*/
// }
// }
// }


// }

//     //std::cout << "  oi2 " << std::endl;
//     slimmedmumuMass_->Fill( (bestslimmedMu1.p4()+bestslimmedMu2.p4()).mass() );
//     slimmedMuonPt_ ->Fill( bestslimmedMu1.pt () );
//     slimmedMuonEta_->Fill( bestslimmedMu1.eta() );
//     slimmedMuonPhi_->Fill( bestslimmedMu1.phi() );
//     slimmedMuonPt_ ->Fill( bestslimmedMu2.pt () );
//     slimmedMuonEta_->Fill( bestslimmedMu2.eta() );
//     slimmedMuonPhi_->Fill( bestslimmedMu2.phi() );

// //################################ (end) slimmedMuons ################################        





// //################################ (begin) slimmedGammas ################################  
//   // Handle to the gamma collection
//   edm::Handle<std::vector<Photon> > slimmedPhotons;
//   event.getByLabel(std::string("slimmedPhotons"), slimmedPhotons);

// pat::Photon bestslimmedGamma;
// zMassDif = 99999999999.;
//   // loop gamma collection and fill histograms
// for(std::vector<Photon>::const_iterator gamma=slimmedPhotons->begin(); gamma!=slimmedPhotons->end(); ++gamma){
//           //muonPt_ ->Fill( mu1->pt () );
//           //muonEta_->Fill( mu1->eta() );
//           //muonPhi_->Fill( mu1->phi() );
//     if(gamma->pt() > 10.0){
//       if(abs((bestslimmedMu1.p4()+bestslimmedMu2.p4()+gamma->p4()).mass() - zMass) < zMassDif){
//         bestslimmedGamma = *gamma->clone();
//         zMassDif = abs((bestslimmedMu1.p4()+bestslimmedMu2.p4()+gamma->p4()).mass() - zMass);
//         std::cout << (bestslimmedMu1.p4()+bestslimmedMu2.p4()+gamma->p4()).mass() << std::endl;
//       }
//     }


// }

// //std::cout << "  oi2 " << std::endl;
// slimmedmumugammaMass_->Fill( (bestslimmedMu1.p4()+bestslimmedMu2.p4()+bestslimmedGamma.p4()).mass() );
// slimmedgammaPt_ ->Fill( bestslimmedGamma.pt () );
// slimmedgammaEta_->Fill( bestslimmedGamma.eta() );
// slimmedgammaPhi_->Fill( bestslimmedGamma.phi() );
// //################################ (end) slimmedGammas ################################


      }

std::cout<<"Nr total of events: "<<nevt_tot <<std::endl;
std::cout<<"Nr of events after mumu sel: "<<nevt_mumu <<std::endl;
std::cout<<"Nr of events after mumugamma sel: "<<nevt_mumugamma <<std::endl;
std::cout<<"Nr of events after min sel: "<<nevt_minsel <<std::endl;



      // close input file
      inFile->Close();
    }
// break loop if maximal number of events is reached:
// this has to be done twice to stop the file loop as well
    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
  }
  return 0;
}
