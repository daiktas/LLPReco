//Description: Produces and fill in muon LLPDNNX features

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "LLPReco/DataFormats/interface/XTagInfo.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "TVector3.h"


class MuonXTagInfoProducer : public edm::stream::EDProducer<> {
public:
    explicit MuonXTagInfoProducer(const edm::ParameterSet&);
    ~MuonXTagInfoProducer();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    struct CandidateHash
    {
        long operator() (const reco::CandidatePtr& cand) const 
        {
            return cand.id().id() * 100000 + cand.key();
        }
    };
    
    private:
        virtual void beginStream(edm::StreamID) override;
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        virtual void endStream() override;


        edm::EDGetTokenT<edm::View<pat::Jet>> jet_token_;
        edm::EDGetTokenT<reco::VertexCollection> vtx_token_;
        edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> sv_token_;
        edm::EDGetTokenT<std::vector<llpdnnx::ChargedCandidateFeatures>>  cpf_token_;
        //edm::EDGetTokenT<pat::MuonCollection> _muons_MiniAODToken_;
};

MuonXTagInfoProducer::MuonXTagInfoProducer(const edm::ParameterSet& iConfig) :
    jet_token_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
    vtx_token_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    cpf_token_(consumes<std::vector<llpdnnx::ChargedCandidateFeatures>>(iConfig.getParameter<edm::InputTag>("cpf_features")))
    //muon_token_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
{
    produces<llpdnnx::MuonCandidateFeatures>();
}


MuonXTagInfoProducer::~MuonXTagInfoProducer(){ }
void MuonXTagInfoProducer::beginStream(edm::StreamID) { }
// ------------ method called to produce the data  ------------
    void
MuonXTagInfoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    auto output_tag_infos = std::make_unique<std::vector<llpdnnx::MuonCandidateFeatures>>();
    edm::Handle<edm::View<pat::Jet>> jets;
    iEvent.getByToken(jet_token_, jets);

    edm::Handle<reco::VertexCollection> vtxs;
    iEvent.getByToken(vtx_token_, vtxs);

    if (vtxs->empty()) {
        // produce empty TagInfos in case no primary vertex
        iEvent.put(std::move(output_tag_infos));
        return;  // exit event
    }

    edm::Handle<std::vector<llpdnnx::ChargedCandidateFeatures>> cpf_features;
    iEvent.getByToken(cpf_token_, cpf_features);

    std::unordered_map<reco::CandidatePtr, const pat::Muon*, CandidateHash> muonMap;

    std::vector<pat::Jet> centralJets;

    for (std::size_t ijet = 0; ijet < jets->size(); ijet++) 
    {
        const pat::Jet& jet = jets->at(ijet);

        // Cut on eta
        if (std::abs(jet.eta()) > 2.4) {
        	continue;
        }

        else{
        	centralJets.emplace_back(jet);
        }
    }


    if (cpf_features->size() != centralJets.size()){
    	throw cms::Exception("Mismatch between cpf and jets collections!");
    }

    for (std::size_t ijet = 0; ijet < centralJets.size(); ijet++) 

        llpdnnx::ChargedCandidateFeatures cpf = cpf_features->at(ijet);

    
        
        // Start with global jet features
        //float uncorrectedPt = jet.correctedP4("Uncorrected").pt();
        const pat::Muon& muon = cpf.matchedMuonPtr->get();
        llpdnnx::MuonCandidateFeatures mu_features;
		t::Muon & muon = *findMuon->second;

                if (not muon.isGlobalMuon() || reco::deltaR(muon ,jet) > 0.4) continue ;
                cpf_features.cpf_matchedMuon = 1;
                mu_features.mu_isGlobal = muon.isGlobalMuon();                                   
                mu_features.mu_isTight = muon.isTightMuon(pv);                                     
                mu_features.mu_isMedium = muon.isMediumMuon();       
                mu_features.mu_isLoose = muon.isLooseMuon();
                mu_features.mu_isStandAlone = muon.isStandAloneMuon();

                mu_features.mu_ptrel = muon.pt()/uncorrectedPt; 
                mu_features.mu_deta = std::fabs(muon.eta()-jet.eta());                                                 
                mu_features.mu_dphi = std::fabs(reco::deltaPhi(muon.phi(),jet.phi()));                                                 
                mu_features.mu_charge = muon.charge();        
                mu_features.mu_energy = muon.energy()/muon.pt();                                           
                mu_features.mu_et = muon.et();   
                mu_features.mu_jetDeltaR = reco::deltaR(muon ,jet); 
                mu_features.mu_numberOfMatchedStations = muon.numberOfMatchedStations();

                mu_features.mu_2dIP = muon.dB() ; 
                mu_features.mu_2dIPSig = muon.dB()/muon.edB(); 
                mu_features.mu_3dIP = muon.dB(pat::Muon::PV3D); 
                mu_features.mu_3dIPSig = muon.dB(pat::Muon::PV3D)/muon.edB(pat::Muon::PV3D);


                reco::Candidate::Vector muonMom = muon.bestTrack()->momentum();

                mu_features.mu_EtaRel =reco::btau::etaRel(jetDir, muonMom);
                mu_features.mu_dxy = muon.bestTrack()->dxy(pv.position());
                mu_features.mu_dxyError = muon.bestTrack()->dxyError();
                mu_features.mu_dxySig = muon.bestTrack()->dxy(pv.position())/muon.bestTrack()->dxyError(); 
                mu_features.mu_dz = muon.bestTrack()->dz(pv.position());
                mu_features.mu_dzError = muon.bestTrack()->dzError();
                mu_features.mu_numberOfValidPixelHits = muon.bestTrack()->hitPattern().numberOfValidPixelHits();
                mu_features.mu_numberOfpixelLayersWithMeasurement = muon.bestTrack()->hitPattern().pixelLayersWithMeasurement();
                mu_features.mu_numberOfstripLayersWithMeasurement = muon.bestTrack()->hitPattern().stripLayersWithMeasurement();
	

                mu_features.mu_chi2 = muon.bestTrack()->chi2();
                mu_features.mu_ndof = muon.bestTrack()->ndof();

                mu_features.mu_caloIso =  muon.caloIso()/muon.pt();
                mu_features.mu_ecalIso =  muon.ecalIso()/muon.pt(); 
                mu_features.mu_hcalIso =  muon.hcalIso()/muon.pt();     


                mu_features.mu_sumPfChHadronPt  = muon.pfIsolationR04().sumChargedHadronPt/muon.pt();
                mu_features.mu_sumPfNeuHadronEt  = muon.pfIsolationR04().sumNeutralHadronEt/muon.pt();
                mu_features.mu_Pfpileup  = muon.pfIsolationR04().sumPUPt/muon.pt();
                mu_features.mu_sumPfPhotonEt = muon.pfIsolationR04().sumPhotonEt/muon.pt();


                mu_features.mu_sumPfChHadronPt03  = muon.pfIsolationR03().sumChargedHadronPt/muon.pt();
                mu_features.mu_sumPfNeuHadronEt03  = muon.pfIsolationR03().sumNeutralHadronEt/muon.pt();
                mu_features.mu_Pfpileup03  = muon.pfIsolationR03().sumPUPt/muon.pt();
                mu_features.mu_sumPfPhotonEt03 = muon.pfIsolationR03().sumPhotonEt/muon.pt();       


                mu_features.mu_timeAtIpInOut = muon.time().timeAtIpInOut;
                mu_features.mu_timeAtIpInOutErr = muon.time().timeAtIpInOutErr;
                mu_features.mu_timeAtIpOutIn = muon.time().timeAtIpOutIn; 

                features.mu_features.emplace_back(mu_features);
            }

            std::stable_sort(features.mu_features.begin(),features.mu_features.end(),[](const auto& d1, const auto& d2)
            {
                if (d1.mu_2dIPSig>0 and d2.mu_2dIPSig>0)
                {
                    if (std::fabs(d1.mu_2dIPSig-d2.mu_2dIPSig)>std::numeric_limits<float>::epsilon())
                    {
                        return std::fabs(d1.mu_2dIPSig)>std::fabs(d2.mu_2dIPSig); //sort decreasing
                    }
                }
                return d1.mu_ptrel>d2.mu_ptrel; //sort decreasing
            });

        });
        */
        //output_tag_infos->emplace_back(features, jet_ref);
    }

    iEvent.put(std::move(output_tag_infos));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuonXTagInfoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("jets", edm::InputTag("ak4PFJetsCHS"));
    desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("inclusiveCandidateSecondaryVertices"));
    desc.add<edm::InputTag>("shallow_tag_infos", edm::InputTag("pfDeepCSVTagInfos"));
}
void MuonXTagInfoProducer::endStream() {};

//define this as a plug-in
DEFINE_FWK_MODULE(MuonXTagInfoProducer);
