//Description: Produces and fill in npf LLPDNNX features


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

#include "LLPReco/DataFormats/interface/XTagInfo.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "TVector3.h"


class NeutralCandidateXTagInfoProducer : public edm::stream::EDProducer<> {
public:
    explicit NeutralCandidateXTagInfoProducer(const edm::ParameterSet&);
    ~NeutralCandidateXTagInfoProducer();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    private:
        virtual void beginStream(edm::StreamID) override;
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        virtual void endStream() override;

        edm::EDGetTokenT<edm::View<pat::Jet>> jet_token_;
        edm::EDGetTokenT<reco::VertexCollection> vtx_token_;
        edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> sv_token_;

};

NeutralCandidateXTagInfoProducer::NeutralCandidateXTagInfoProducer(const edm::ParameterSet& iConfig) :
    jet_token_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
    vtx_token_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    sv_token_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices")))
{
    produces<std::vector<std::vector<llpdnnx::NeutralCandidateFeatures>>>();
}


NeutralCandidateXTagInfoProducer::~NeutralCandidateXTagInfoProducer(){ }
void NeutralCandidateXTagInfoProducer::beginStream(edm::StreamID) { }
// ------------ method called to produce the data  ------------
    void
NeutralCandidateXTagInfoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    auto output_tag_infos = std::make_unique<std::vector<std::vector<llpdnnx::NeutralCandidateFeatures>>>();
    edm::Handle<edm::View<pat::Jet>> jets;
    iEvent.getByToken(jet_token_, jets);

    edm::Handle<reco::VertexCollection> vtxs;
    iEvent.getByToken(vtx_token_, vtxs);

    edm::Handle<reco::VertexCompositePtrCandidateCollection> svs;
    iEvent.getByToken(sv_token_, svs);

    if (vtxs->empty()) {
        // produce empty TagInfos in case no primary vertex
        iEvent.put(std::move(output_tag_infos));
        return;  // exit event
    }

    for (std::size_t ijet = 0; ijet < jets->size(); ijet++) 
    {
        const pat::Jet& jet = jets->at(ijet);

        // Cut on eta
        if (std::abs(jet.eta()) > 2.4) continue;
            
        // Start with global jet features
        float uncorrectedPt = jet.correctedP4("Uncorrected").pt();

        // Fill neutral hadron info
        std::vector<llpdnnx::NeutralCandidateFeatures> jet_tag_infos;

        for (unsigned int idaughter = 0; idaughter < jet.numberOfDaughters(); ++idaughter)
        {
            const pat::PackedCandidate* constituent = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(idaughter));
            if ((not constituent) or constituent->charge()!=0)
            {
                continue;
            }
            llpdnnx::NeutralCandidateFeatures npf_features;

            npf_features.npf_ptrel = constituent->pt()/uncorrectedPt;
            npf_features.npf_deta = std::fabs(constituent->eta()-jet.eta());
            npf_features.npf_dphi = std::fabs(reco::deltaPhi(constituent->phi(),jet.phi()));
            npf_features.npf_puppi_weight = constituent->puppiWeight();
            npf_features.npf_deltaR = reco::deltaR(*constituent,jet);
            npf_features.npf_isGamma = abs(constituent->pdgId())==22;
            npf_features.npf_hcal_fraction = constituent->hcalFraction();

            npf_features.npf_drminsv = 0.4;
            for (const auto& sv: *svs.product())
            {
                float dR = reco::deltaR(sv,*constituent);
                npf_features.npf_drminsv = std::min(npf_features.npf_drminsv,dR);
            }

            if (jet.mass()<1e-10) 
            {
                npf_features.npf_relmassdrop = -1;
            }
            else
            {
                npf_features.npf_relmassdrop = (jet.p4()- constituent->p4()).mass()/jet.mass();
            }
            jet_tag_infos.emplace_back(npf_features);
            
        }
        std::stable_sort(jet_tag_infos.begin(),jet_tag_infos.end(),[](llpdnnx::NeutralCandidateFeatures d1, llpdnnx::NeutralCandidateFeatures d2)
        {

            if (std::fabs(d1.npf_drminsv-d2.npf_drminsv)>std::numeric_limits<float>::epsilon())
            {
                return d1.npf_drminsv<d2.npf_drminsv; //sort increasing
            }
            else
            {
                return d1.npf_ptrel>d2.npf_ptrel; //sort decreasing
            }
            return false;
        });
    output_tag_infos->emplace_back(jet_tag_infos);

    }

    iEvent.put(std::move(output_tag_infos));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void NeutralCandidateXTagInfoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
    desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
    desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("inclusiveCandidateSecondaryVertices"));
}
void NeutralCandidateXTagInfoProducer::endStream() {};

//define this as a plug-in
DEFINE_FWK_MODULE(NeutralCandidateXTagInfoProducer);
