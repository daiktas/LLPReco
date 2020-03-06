//Description: Produces and fill in sv LLPDNNX features

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

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "LLPReco/DataFormats/interface/XTagInfo.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "TVector3.h"


class SecondaryVertexXTagInfoProducer : public edm::stream::EDProducer<> {
public:
    explicit SecondaryVertexXTagInfoProducer(const edm::ParameterSet&);
    ~SecondaryVertexXTagInfoProducer();
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
};

SecondaryVertexXTagInfoProducer::SecondaryVertexXTagInfoProducer(const edm::ParameterSet& iConfig) :
    jet_token_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
    vtx_token_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    sv_token_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices")))
{
    produces<std::vector<std::vector<llpdnnx::SecondaryVertexFeatures>>>();
}


SecondaryVertexXTagInfoProducer::~SecondaryVertexXTagInfoProducer(){ }
void SecondaryVertexXTagInfoProducer::beginStream(edm::StreamID) { }
// ------------ method called to produce the data  ------------
    void
SecondaryVertexXTagInfoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    auto output_tag_infos = std::make_unique<std::vector<std::vector<llpdnnx::SecondaryVertexFeatures>>>();
    edm::Handle<edm::View<pat::Jet>> jets;
    iEvent.getByToken(jet_token_, jets);

    edm::Handle<reco::VertexCollection> vtxs;
    iEvent.getByToken(vtx_token_, vtxs);

    if (vtxs->empty()) {
        // produce empty TagInfos in case no primary vertex
        iEvent.put(std::move(output_tag_infos));
        return;  // exit event
    }
    const auto& pv = vtxs->at(0);
    edm::ESHandle<TransientTrackBuilder> builder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

    edm::Handle<reco::VertexCompositePtrCandidateCollection> svs;
    iEvent.getByToken(sv_token_, svs);

    for (std::size_t ijet = 0; ijet < jets->size(); ijet++) 
    {
        const pat::Jet& jet = jets->at(ijet);

        float uncorrectedPt = jet.correctedP4("Uncorrected").pt();

        std::unordered_set<reco::CandidatePtr, CandidateHash> jetConstituentSet;
        for (unsigned int idaughter = 0; idaughter < jet.numberOfDaughters(); ++idaughter)
        {
            jetConstituentSet.insert(jet.daughterPtr(idaughter));
        }

        std::unordered_set<reco::CandidatePtr, CandidateHash> candidatesMatchedToSV;
        // fill features from secondary vertices

        std::vector<llpdnnx::SecondaryVertexFeatures> jet_tag_infos;

        for (unsigned int isv = 0; isv < svs->size(); ++isv)
        {
            const reco::VertexCompositePtrCandidate& sv = svs->at(isv);
            
            if (reco::deltaR(sv,jet)>0.4)
            {
                continue;
            }
            bool matchingTrack = false;
            for (auto const& candidateFromVertex: sv.daughterPtrVector())
            {
                if (jetConstituentSet.find(candidateFromVertex)!=jetConstituentSet.end())
                {
                    candidatesMatchedToSV.insert(candidateFromVertex);
                    matchingTrack = true;
                }
            }
            if (not matchingTrack) continue;
            
            llpdnnx::SecondaryVertexFeatures sv_features;

            sv_features.sv_ptrel = sv.pt()/uncorrectedPt;
            sv_features.sv_deta = std::fabs(sv.eta()-jet.eta());
            sv_features.sv_dphi = std::fabs(reco::deltaPhi(sv.phi(),jet.phi()));
            sv_features.sv_deltaR = reco::deltaR(sv,jet);
            sv_features.sv_mass = sv.mass();
            sv_features.sv_ntracks = sv.numberOfDaughters();
            sv_features.sv_chi2 = sv.vertexChi2();
            sv_features.sv_ndof = sv.vertexNdof();


            reco::Vertex::CovarianceMatrix covsv; 
            sv.fillVertexCovariance(covsv);
            reco::Vertex svtx(sv.vertex(), covsv);

            VertexDistanceXY distXY;
            Measurement1D distanceXY = distXY.distance(svtx, pv);
            sv_features.sv_dxy = distanceXY.value();
            sv_features.sv_dxysig = distanceXY.value()/distanceXY.error();

            VertexDistance3D dist3D;
            Measurement1D distance3D = dist3D.distance(svtx, pv);
            sv_features.sv_d3d = distance3D.value();
            sv_features.sv_d3dsig = distance3D.value()/distance3D.error();

            reco::Candidate::Vector distance(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
            sv_features.sv_costhetasvpv = sv.momentum().Unit().Dot(distance.Unit());
            sv_features.sv_enratio = sv.energy()/jet.pt();


            jet_tag_infos.emplace_back(sv_features);
        }

        std::stable_sort(jet_tag_infos.begin(),jet_tag_infos.end(),[](llpdnnx::SecondaryVertexFeatures d1, llpdnnx::SecondaryVertexFeatures d2)
        {
            return d1.sv_dxysig>d2.sv_dxysig; //sort decreasing
        });

    output_tag_infos->emplace_back(jet_tag_infos);


    }

    iEvent.put(std::move(output_tag_infos));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SecondaryVertexXTagInfoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
    desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
    desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("inclusiveCandidateSecondaryVertices"));
    descriptions.addDefault(desc);
}
void SecondaryVertexXTagInfoProducer::endStream() {};

//define this as a plug-in
DEFINE_FWK_MODULE(SecondaryVertexXTagInfoProducer);
