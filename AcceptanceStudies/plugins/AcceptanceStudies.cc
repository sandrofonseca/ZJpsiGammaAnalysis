// -*- C++ -*-
//
// Package:    ZJpsiGammaAnalysis/AcceptanceStudies
// Class:      AcceptanceStudies
// 
/**\class AcceptanceStudies AcceptanceStudies.cc ZJpsiGammaAnalysis/AcceptanceStudies/plugins/AcceptanceStudies.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Sandro Fonseca de Souza
//         Created:  Tue, 18 Apr 2017 17:47:36 GMT
//
//


// system include files
#include <memory>
#include <cmath>

#include <TLorentzVector.h>
#include <TVector3.h>

// user include files

// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "Math/GenVector/VectorUtil.h"





class AcceptanceStudies : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
	explicit AcceptanceStudies(const edm::ParameterSet&);
	~AcceptanceStudies();


private:
	virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;
	bool PATGenFilter(edm::Handle <pat::PackedGenParticleCollection> genColl, const edm::Event &iEvent);
		///      bool PATRecoFilter(edm::Handle<vector<pat::Muon>> muons,edm::Handle<edm::View<pat::Photon> > photonHandle,const edm::Event &iEvent);


		// ----------member data ---------------------------

		//edm::EDGetTokenT<pat::MuonCollection> muonCollToken;
	edm::EDGetTokenT<pat::PackedGenParticleCollection> genCollToken;
		//  edm::EDGetTokenT<edm::View<pat::Photon> > photonCollection_;




		//genMuons
	TH1F* leadingMuonPt_ ;
	TH1F* leadingMuonEta_ ;
		//TH2F* leadingMuonPtVsleadingMuonEta_ ;
	TH1F* leadingMuonPhi_;
	TH1F* trailingMuonPt_ ;
	TH1F* trailingMuonEta_ ;
	TH1F* trailingMuonPhi_ ;
	TH1F* mumuMass_;

		//MinSel
	TH1F* leadingMuonPtMinSel_ ;
	TH1F* leadingMuonEtaMinSel_ ;
	TH1F* leadingMuonPhiMinSel_ ;
		//  TH2F* leadingMuonPtVsleadingMuonEtaMinSel_ ;
	TH1F* trailingMuonPtMinSel_  ;
	TH1F* trailingMuonEtaMinSel_ ;
	TH1F* trailingMuonPhiMinSel_ ;
	TH1F* mumuMassMinSel_;

		//genGammas
	TH1F* gammaPt_ ;
	TH1F* gammaEta_ ;
	TH1F* gammaPhi_ ;
	TH1F* mumugammaMass_;

		//MinSel
	TH1F* gammaPtMinSel_  ;
	TH1F* gammaEtaMinSel_ ;
	TH1F* gammaPhiMinSel_ ;
	TH1F* mumugammaMassMinSel_; 
	TH1F* Jpsi_Mass_;

		// Evts Counters
	int nEvts;
	int nEvtsGEN;
	int nEvtsRECO;
	TH1D* TnEvts;
	TH1D* TnEvtsGEN;

  // Histos map
	std::map<std::string, TH1D*> nEvtsHistosMap;

		// GEN configs
	bool verbose_;
	std::string configName_;
	double minMuPt_;
	double maxMuEta_;
	double muonLeadPt_, muonTrailPt_;
	double minJPsiMass_ ;
	double maxJPsiMass_,GammaMinPtCut_,drLeadMuPhotonSel_,drTrailPhotonSel_; 


};

AcceptanceStudies::AcceptanceStudies(const edm::ParameterSet& iConfig):
verbose_ (iConfig.getParameter< bool > ("verbose")),
configName_ (iConfig.getParameter< std::string > ("configName")),
minMuPt_ (iConfig.getParameter<double>("minMuPt")),
maxMuEta_ (iConfig.getParameter<double>("maxMuEta")), 
muonLeadPt_ (iConfig.getParameter<double>("minMuonLeadPt")),
muonTrailPt_ (iConfig.getParameter<double>("minMuonTrailPt")),
minJPsiMass_ (iConfig.getParameter<double>("minJPsiMass")),
maxJPsiMass_ (iConfig.getParameter<double>("maxJPsiMass")),
GammaMinPtCut_ (iConfig.getParameter<double>("GammaMinPtCut")),
drLeadMuPhotonSel_ (iConfig.getParameter<double>("DeltaRLeadMuPhotonSel")),
drTrailPhotonSel_ (iConfig.getParameter<double>("DeltaRTrailPhotonSel"))




{
	// edm::InputTag theMuonLabel("slimmedMuons");//Reco-PAT collections
	edm::InputTag theGenMuonLabel("packedGenParticles");// Gen particles

	//muonCollToken = consumes<pat::MuonCollection>(theMuonLabel);
	genCollToken = consumes<pat::PackedGenParticleCollection>(theGenMuonLabel);

	edm::Service<TFileService> fs;
	// Define Evts count
	nEvts = 0;
	nEvtsGEN = 0;
	nEvtsRECO = 0;

	// Books evts counters
	TnEvts = fs->make<TH1D>( ("h_nEvts_"+configName_).c_str() , ("h_nEvts_"+configName_+";  x; NEvts").c_str(), 1, 0., 1.);
	TnEvtsGEN = fs->make<TH1D>( ("h_nEvtsGEN_"+configName_).c_str() , ("h_nEvtsGEN_"+configName_+";  x; nEvtsGEN").c_str(), 1, 0., 1.);


	//TFileDirectory fs->= fs->mkfs->"plots");

	//genMuons
	TH1F* leadingMuonPt_  = fs->make<TH1F>("leadingMuonPt"  , "mu1_pt"  ,   100,   0., 200.);
	TH1F* leadingMuonEta_ = fs->make<TH1F>("leadingMuonEta" , "mu1_eta" ,   100,  -5.,   5.);
	// TH2F* leadingMuonPtVsleadingMuonEta_ = fs->make<TH2F>("leadingMuonPtVsleadingMuonEta", "Pt Vs Ets",100,0.,100.,100,-5., 5.);
	TH1F* leadingMuonPhi_ = fs->make<TH1F>("leadingMuonPhi" , "mu1_phi" ,   100,  -4.,   4.);
	TH1F* trailingMuonPt_  = fs->make<TH1F>("trailingMuonPt"  , "mu2_pt"  ,   100,   0., 200.);
	TH1F* trailingMuonEta_ = fs->make<TH1F>("trailingMuonEta" , "mu2_eta" ,   100,  -5.,   5.);
	TH1F* trailingMuonPhi_ = fs->make<TH1F>("trailingMuonPhi" , "mu2_phi" ,   100,  -4.,   4.);
	TH1F* mumuMass_= fs->make<TH1F>("mumuMass", "Dimuon_Mass",    100,  2.8,  10.2);

	//MinSel
	TH1F* leadingMuonPtMinSel_  = fs->make<TH1F>("leadingMuonPtMinSel"  , "mu1_pt"  ,   100,   0., 200.);
	TH1F* leadingMuonEtaMinSel_ = fs->make<TH1F>("leadingMuonEtaMinSel" , "mu1_eta" ,   100,  -5.,   5.);
	TH1F* leadingMuonPhiMinSel_ = fs->make<TH1F>("leadingMuonPhiMinSel" , "mu1_phi" ,   100,  -4.,   4.);
	//TH2F* leadingMuonPtVsleadingMuonEtaMinSel_ = fs->make<TH2F>("leadingMuonPtVsleadingMuonEtaMinSel", "mu1_pt x mu1_eta" ,100,0.,100.,100,-5., 5.);
	TH1F* trailingMuonPtMinSel_  = fs->make<TH1F>("trailingMuonPtMinSel"  , "mu2_pt"  ,   100,   0., 200.);
	TH1F* trailingMuonEtaMinSel_ = fs->make<TH1F>("trailingMuonEtaMinSel" , "mu2_eta" ,   100,  -5.,   5.);
	TH1F* trailingMuonPhiMinSel_ = fs->make<TH1F>("trailingMuonPhiMinSel" , "mu2_phi" ,   100,  -4.,   4.);
	TH1F* mumuMassMinSel_= fs->make<TH1F>("mumuMassMinSel", "Dimuon_Mass",    100,  2.8,  10.2);
	TH1F* Jpsi_Mass_= fs->make<TH1F>("Jpsi_Mass", "Jpsi_Mass",    100,  2.8,  3.2);

	//genGammas
	TH1F* gammaPt_  = fs->make<TH1F>("gammaPt"  , "pt"  ,   100,   0., 300.);
	TH1F* gammaEta_ = fs->make<TH1F>("gammaEta" , "eta" ,   100,  -5.,   5.);
	TH1F* gammaPhi_ = fs->make<TH1F>("gammaPhi" , "phi" ,   100,  -4.,   4.);
	TH1F* mumugammaMass_= fs->make<TH1F>("zMass", "mass",    100,  65.,  115.);

	//MinSel
	TH1F* gammaPtMinSel_  = fs->make<TH1F>("gammaPtMinSel"  , "pt"  ,   100,   0., 300.);
	TH1F* gammaEtaMinSel_ = fs->make<TH1F>("gammaEtaMinSel" , "eta" ,   100,  -5.,   5.);
	TH1F* gammaPhiMinSel_ = fs->make<TH1F>("gammaPhiMinSel" , "phi" ,   100,  -4.,   4.);
	TH1F* mumugammaMassMinSel_= fs->make<TH1F>("zMassMinSel", "mass",    100,  65.,  115.);


}


AcceptanceStudies::~AcceptanceStudies()
{

}


//
// member functions
//

// ------------ method called for each event  ------------
void
AcceptanceStudies::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace std;
	using namespace reco;
	using namespace pat;

	// GEN Muons
	//  
	edm::Handle <pat::PackedGenParticleCollection> genColl;
	iEvent.getByToken(genCollToken, genColl);
	//
	// RECO Muons
	//
	//edm::Handle<vector<pat::Muon>> muons;
	//iEvent.getByToken(muonCollToken, muons);

	//Reco Photon

	//edm::Handle<edm::View<pat::Photon> > photonHandle;
	//iEvent.getByToken(photonCollection_, photonHandle);

	bool genTest = PATGenFilter(genColl,iEvent);

	//bool recoTest = PATRecoFilter(muons, photonHandle, iEvent);

	nEvts++;
	if (genTest == true) nEvtsGEN++;
	// if (recoTest == true) nEvtsRECO++;


}// end function
///////-----------------
bool AcceptanceStudies::PATGenFilter(edm::Handle <pat::PackedGenParticleCollection> genColl, const edm::Event &iEvent)
{
	int nDimuon = 0, nJpsi = 0, nPhoton = 0; 
	using namespace std;
	vector<pat::PackedGenParticle> myLeptons;
	vector<pat::PackedGenParticle> myPhotons; 	
	// Gen matching
	for (auto genMuon = genColl->cbegin(); genMuon != genColl->cend(); ++genMuon) {
		const pat::PackedGenParticle& mcMuon = (*genMuon);
		if ( not (abs(mcMuon.pdgId()) == 13 ) ) continue; // make sure it is a muon
		if ( mcMuon.pt() > minMuPt_  && fabs(mcMuon.eta()) <= maxMuEta_ && mcMuon.status() == 1 ) {;
			if(verbose_) cout<< "mcMuon.pt() "  << mcMuon.pt() << endl;
			myLeptons.push_back(*genMuon);
		}
	}

	std::sort(myLeptons.begin(),myLeptons.end(), [](const pat::PackedGenParticle &a, const pat::PackedGenParticle &b){
		return a.pt() > b.pt();
	});


	if(verbose_) std::cout<<"PAT myLeptons.size() all  " << myLeptons.size() << std::endl;


	pat::PackedGenParticle leadingMuon;
	pat::PackedGenParticle trailingMuon;
	if (myLeptons.size() >= 2) {
		nDimuon++;
		if(verbose_) std::cout<<"PAT  Muons Multiplicity:  " << myLeptons.size() << std::endl; 
		// if(verbose_) std::cout<<"PAT Dimuons Multiplicity:  " << nDimuon << std::endl;
		leadingMuon = myLeptons[0];
		trailingMuon = myLeptons[1];
		//Dimuons  selection
		// if ((leadingMuon.charge() != trailingMuon.charge())) {
		// } else {return false;}

		if(verbose_) std::cout<< "Leading Muon pt, eta, phi, charge = " << leadingMuon.pt() << " "<< leadingMuon.eta() << " "<< leadingMuon.phi() << " " << leadingMuon.charge() << std::endl;
		if(verbose_) std::cout<< "Trailing Muon  pt, eta, phi,charge = " << trailingMuon.pt() << " " << trailingMuon.eta() << " " << trailingMuon.phi() << " " << trailingMuon.charge()<< std::endl;

		//Invariant mass of dimuons
		double Mll = (leadingMuon.p4() + trailingMuon.p4()).mass();
		double MllpT = (leadingMuon.p4() + trailingMuon.p4()).pt();
		double Mlleta = (leadingMuon.p4() + trailingMuon.p4()).eta();
		double Mllphi = (leadingMuon.p4() + trailingMuon.p4()).phi();

        //leadingMuonPt_->Fill(leadingMuon.pt());
        //leadingMuonEta_->Fill(leadingMuon.eta());
        //leadingMuonPhi_->Fill(leadingMuon.phi());
		if(verbose_) std::cout<< "Dimuons Invariant Mass Mll: " << Mll << std::endl; 
        //trailingMuonPt_->Fill(trailingMuon.pt());
        //trailingMuonEta_->Fill(trailingMuon.eta()); 
        //trailingMuonPhi_->Fill(trailingMuon.phi()); 
        //mumuMass_->Fill(Mll);

		if(verbose_) std::cout<< "Dimuons Invariant Mass Mll, pT, eta, phi: " << Mll << " " << MllpT << " " << Mlleta << " " << Mllphi << std::endl;
		if (leadingMuon.pt() >= muonLeadPt_ && trailingMuon.pt() >= muonTrailPt_ && fabs(leadingMuon.eta()) <= 2.4 && fabs(trailingMuon.eta()) <= 2.4 ) {
        /*    leadingMuonPtMinSel_->Fill(leadingMuon.pt());
            leadingMuonEtaMinSel_->Fill(leadingMuon.eta());
            leadingMuonPhiMinSel_->Fill(leadingMuon.phi());

            trailingMuonPtMinSel_->Fill(trailingMuon.pt());
            trailingMuonEtaMinSel_->Fill(trailingMuon.eta()); 
            trailingMuonPhiMinSel_->Fill(trailingMuon.phi());
            mumuMassMinSel_->Fill(Mll);
	*/	
		// ***
			//   // jpsi peak
			//     // ***
			//
			// if (Mll > minJPsiMass_ && Mll < maxJPsiMass_){
			if (true) {
				nJpsi++;                           
				if(verbose_) std::cout<<" Invariant Mass in JPsi peak, pT, eta, phi " << Mll << " " << MllpT << " " << Mlleta << " " << Mllphi << std::endl;
          //         Jpsi_Mass_->Fill(Mll);
				// if(verbose_) std::cout<<" Jpsi Multiplicity:  " <<  nJpsi << std::endl;
			} else {return false;}// jpsi selection
		} else {return false;}//lead and trail muon pT cut



	} else {return false;}//dimuons selection

	// Photon Loop

	for (auto genPhoton = genColl->cbegin(); genPhoton != genColl->cend(); ++genPhoton) {
		const pat::PackedGenParticle& mcPhoton = (*genPhoton);
			if ( not (mcPhoton.pdgId() == 22) ) continue; // make sure it is a photon
			//if(verbose_) cout<<" mcPhoton.pdgId()" << mcPhoton.pdgId() << "mcPhoton.pt() "  << mcPhoton.pt() << " mcPhoton.isPromptFinalState() " <<  mcPhoton.status() << endl;
			if(mcPhoton.pt() > GammaMinPtCut_ && mcPhoton.status() == 1 && fabs(mcPhoton.eta()) <= 2.5 ) {	
				if(verbose_) cout<<" mcPhoton.pdgId() " << mcPhoton.pdgId() << " mcPhoton.pt() "  << mcPhoton.pt() << " mcPhoton.State() " <<  mcPhoton.status() << endl;

				myPhotons.push_back(*genPhoton);
			}
		}//end loop photon



		std::sort(myPhotons.begin(),myPhotons.end(), [](const pat::PackedGenParticle &a, const pat::PackedGenParticle &b){
			return a.pt() > b.pt();
		});

		if (  myPhotons.size() >= 1 ) {
			nPhoton++;
			if(verbose_) std::cout<<" Photon Multiplicity:  " <<  nPhoton << std::endl;
			const pat::PackedGenParticle Gamma = myPhotons[0];         
			DeltaR<pat::PackedGenParticle, pat::PackedGenParticle> deltaR;
			double drLeadMuPhoton = deltaR(leadingMuon,Gamma);
			double drTrailPhoton = deltaR(trailingMuon,Gamma);
			if(verbose_) std::cout << " photon: pT, eta, phi " << Gamma.pt() << " "<< Gamma.eta() << " " << Gamma.phi() <<std::endl;                         
			if(verbose_) std::cout<< " DeltaR(LeadMu,Photon) " << drLeadMuPhoton << " DeltaR(TrailMu,Photon) " << drTrailPhoton <<std::endl;
/*
            gammaPt_->Fill(Gamma.pt());  
            gammaEta_->Fill(Gamma.eta());
            gammaPhi_->Fill(Gamma.phi());
*/
			double Mllg = (leadingMuon.p4() + trailingMuon.p4() + Gamma.p4()).mass();
			double MllgpT = (leadingMuon.p4() + trailingMuon.p4() + Gamma.p4()).pt();
			double Mllgeta = (leadingMuon.p4() + trailingMuon.p4() + Gamma.p4()).eta();    
			double Mllgphi = (leadingMuon.p4() + trailingMuon.p4() + Gamma.p4()).phi();
  //          mumugammaMass_->Fill(Mllg);

			// if (drLeadMuPhoton > drLeadMuPhotonSel_ && drTrailPhoton > drTrailPhotonSel_){
			if (true){

				//double Mllg = (leadingMuon.p4() + trailingMuon.p4() + Gamma.p4()).mass();
				//double MllgpT = (leadingMuon.p4() + trailingMuon.p4() + Gamma.p4()).pt();
				//double Mllgeta = (leadingMuon.p4() + trailingMuon.p4() + Gamma.p4()).eta();
				//double Mllgphi = (leadingMuon.p4() + trailingMuon.p4() + Gamma.p4()).phi();    
    //            gammaPtMinSel_->Fill(Gamma.pt());  
      //          gammaEtaMinSel_->Fill(Gamma.eta()); 
        //        gammaPhiMinSel_->Fill(Gamma.phi()); 
          //      mumugammaMassMinSel_->Fill(Mllg);
				if(verbose_) std::cout<< "Invariant Mass Mllg, pT, eta, phi: " << Mllg << " " << MllgpT << " " << Mllgeta << " " << Mllgphi << std::endl;  
			} else {return false;}// deltaR cuts 
		} else {return false;}//photon selection

		return true;
	}

// ------------ method called once each job just before starting event loop  ------------
	void 
	AcceptanceStudies::beginJob()
	{
	}

// ------------ method called once each job just after ending the event loop  ------------
	void 
	AcceptanceStudies::endJob() 
	{

 // set evts counters
		TnEvts->SetBinContent(1, nEvts);
		TnEvtsGEN->SetBinContent(1, nEvtsGEN);

	}


//define this as a plug-in
	DEFINE_FWK_MODULE(AcceptanceStudies);
