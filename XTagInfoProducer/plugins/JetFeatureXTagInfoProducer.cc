//Description: Produces and fill in LLPDNNX global jet features

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
#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"

#include "LLPReco/DataFormats/interface/XTagInfo.h"
#include "LLPReco/XTagInfoProducer/interface/JetSubstructure.h"


class JetFeatureXTagInfoProducer : public edm::stream::EDProducer<> {
public:
    explicit JetFeatureXTagInfoProducer(const edm::ParameterSet&);
    ~JetFeatureXTagInfoProducer();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    private:
        virtual void beginStream(edm::StreamID) override;
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        virtual void endStream() override;

        edm::EDGetTokenT<edm::View<pat::Jet>> jet_token_;
        edm::EDGetTokenT<reco::VertexCollection> vtx_token_;
        edm::EDGetTokenT<edm::View<reco::ShallowTagInfo>> shallow_tag_info_token_;
};

JetFeatureXTagInfoProducer::JetFeatureXTagInfoProducer(const edm::ParameterSet& iConfig) :
    jet_token_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
    vtx_token_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    shallow_tag_info_token_(consumes<edm::View<reco::ShallowTagInfo>>(iConfig.getParameter<edm::InputTag>("shallow_tag_infos")))
{
    produces<std::vector<llpdnnx::JetFeatures>>();
}


JetFeatureXTagInfoProducer::~JetFeatureXTagInfoProducer(){ }
void JetFeatureXTagInfoProducer::beginStream(edm::StreamID) { }
// ------------ method called to produce the data  ------------
    void
JetFeatureXTagInfoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    auto output_tag_infos = std::make_unique<std::vector<llpdnnx::JetFeatures>>();
    edm::Handle<edm::View<pat::Jet>> jets;
    iEvent.getByToken(jet_token_, jets);

    edm::Handle<reco::VertexCollection> vtxs;
    iEvent.getByToken(vtx_token_, vtxs);

    if (vtxs->empty()) {
        // produce empty TagInfos in case no primary vertex
        iEvent.put(std::move(output_tag_infos));
        return;  // exit event
    }
    edm::Handle<edm::View<reco::ShallowTagInfo>> shallow_tag_infos;
    iEvent.getByToken(shallow_tag_info_token_, shallow_tag_infos);

    for (std::size_t ijet = 0; ijet < jets->size(); ijet++) 
    {
        const pat::Jet& jet = jets->at(ijet);
      
        // Cut on eta
        if (std::abs(jet.eta()) > 2.4) continue;
    
        // create data containing structure
        llpdnnx::JetFeatures features;
        
        // Start with global jet features
        float uncorrectedPt = jet.correctedP4("Uncorrected").pt();
        
        features.pt = uncorrectedPt;  // uncorrected
        features.eta = jet.eta();

        features.phi = jet.phi();
        features.mass = jet.mass();
        features.energy = jet.energy();
        features.area = jet.jetArea();
        
        features.n60 = jet.n60();
        features.n90 = jet.n90();
        
        features.chargedEmEnergyFraction = jet.chargedEmEnergyFraction();
        features.chargedHadronEnergyFraction = jet.chargedHadronEnergyFraction();
        features.chargedMuEnergyFraction = jet.chargedMuEnergyFraction();
        features.electronEnergyFraction = jet.electronEnergyFraction();
        
        llpdnnx::JetSubstructure jetSubstructure(jet);

        features.tau1 = jetSubstructure.nSubjettiness(1);
        features.tau2 = jetSubstructure.nSubjettiness(2);
        features.tau3 = jetSubstructure.nSubjettiness(3);

        features.relMassDropMassAK = jetSubstructure.massDropMass(llpdnnx::JetSubstructure::ClusterType::AK)/features.mass;
        features.relMassDropMassCA = jetSubstructure.massDropMass(llpdnnx::JetSubstructure::ClusterType::CA)/features.mass;
        features.relSoftDropMassAK = jetSubstructure.softDropMass(llpdnnx::JetSubstructure::ClusterType::AK)/features.mass;
        features.relSoftDropMassCA = jetSubstructure.softDropMass(llpdnnx::JetSubstructure::ClusterType::CA)/features.mass;
       
        auto eventShapes = jetSubstructure.eventShapeVariables();
        features.thrust = jetSubstructure.thrust();
        features.sphericity = eventShapes.sphericity();
        features.circularity = eventShapes.circularity();
        features.isotropy = eventShapes.isotropy();
        features.eventShapeC = eventShapes.C();
        features.eventShapeD = eventShapes.D();

        // Add CSV variables
        const edm::View<reco::ShallowTagInfo>& taginfos = *shallow_tag_infos;
        edm::Ptr<reco::ShallowTagInfo> match;

        for (auto it = taginfos.begin(); it != taginfos.end(); ++it) {
            float dR = reco::deltaR(it->jet()->p4(),jet.p4());
            if (dR<0.01) {
                match = taginfos.ptrAt(it - taginfos.begin());
                break;
            }
        }
        reco::ShallowTagInfo tag_info;
        if (match.isNonnull()) {
            tag_info = *match;
        }  // will be default values otherwise

        reco::TaggingVariableList vars = tag_info.taggingVariables();
        features.csv_trackSumJetEtRatio = vars.get(reco::btau::trackSumJetEtRatio, -1);
        features.csv_trackSumJetDeltaR = vars.get(reco::btau::trackSumJetDeltaR, -1);
        features.csv_vertexCategory = vars.get(reco::btau::vertexCategory, -1);
        features.csv_trackSip2dValAboveCharm = vars.get(reco::btau::trackSip2dValAboveCharm, -1);
        features.csv_trackSip2dSigAboveCharm = vars.get(reco::btau::trackSip2dSigAboveCharm, -1);
        features.csv_trackSip3dValAboveCharm = vars.get(reco::btau::trackSip3dValAboveCharm, -1);
        features.csv_trackSip3dSigAboveCharm = vars.get(reco::btau::trackSip3dSigAboveCharm, -1);
        features.csv_jetNTracksEtaRel = vars.get(reco::btau::jetNTracksEtaRel, -1);
        features.csv_jetNSelectedTracks = vars.get(reco::btau::jetNSelectedTracks, -1);


        output_tag_infos->emplace_back(features);
    }

    iEvent.put(std::move(output_tag_infos));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void JetFeatureXTagInfoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
    desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
    desc.add<edm::InputTag>("shallow_tag_infos", edm::InputTag("pfDeepCSVTagInfos"));
}
void JetFeatureXTagInfoProducer::endStream() {};

//define this as a plug-in
DEFINE_FWK_MODULE(JetFeatureXTagInfoProducer);
