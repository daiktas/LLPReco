//Description: Produces and fill in cpf LLPDNNX features

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
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "LLPReco/DataFormats/interface/XTagInfo.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "TVector3.h"


class ChargedCandidateXTagInfoProducer : public edm::stream::EDProducer<> {
public:
    explicit ChargedCandidateXTagInfoProducer(const edm::ParameterSet&);
    ~ChargedCandidateXTagInfoProducer();
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
        edm::EDGetTokenT<edm::View<reco::Candidate>> candidateToken_;
        edm::EDGetTokenT<pat::MuonCollection> muon_token_;
        edm::EDGetTokenT<pat::ElectronCollection> electron_token_;
};

ChargedCandidateXTagInfoProducer::ChargedCandidateXTagInfoProducer(const edm::ParameterSet& iConfig) :
    jet_token_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
    vtx_token_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    sv_token_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))),
    muon_token_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    electron_token_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons")))
{
    produces<std::unordered_map<edm::RefToBase<reco::Jet>, std::vector<llpdnnx::ChargedCandidateFeatures>>>();
}


ChargedCandidateXTagInfoProducer::~ChargedCandidateXTagInfoProducer(){ }
void ChargedCandidateXTagInfoProducer::beginStream(edm::StreamID) { }
// ------------ method called to produce the data  ------------
    void
ChargedCandidateXTagInfoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    auto output_tag_infos = std::make_unique<std::unordered_map<edm::RefToBase<reco::Jet>, std::vector<llpdnnx::ChargedCandidateFeatures>>>();
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

    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muon_token_, muons);

    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electron_token_, electrons);

    std::unordered_map<reco::CandidatePtr, const pat::Muon*, CandidateHash> muonMap;
    std::unordered_map<reco::CandidatePtr, const pat::Electron*, CandidateHash> electronMap;

    for (const pat::Muon& muon: *muons)
    {
        for (unsigned int i = 0 ; i < muon.numberOfSourceCandidatePtrs(); ++i )
        {
            muonMap[muon.sourceCandidatePtr(i)] = &muon;
        }
    }


    for (const pat::Electron& electron: *electrons)
    {
        for (unsigned int i = 0 ; i < electron.numberOfSourceCandidatePtrs(); ++i )
        {
            electronMap[electron.sourceCandidatePtr(i)] = &electron;
        }
    }

    for (std::size_t ijet = 0; ijet < jets->size(); ijet++) 
    {
        const pat::Jet& jet = jets->at(ijet);
        edm::RefToBase<reco::Jet> jet_ref(jets->refAt(ijet));


        std::unordered_set<reco::CandidatePtr, CandidateHash> jetConstituentSet;
        for (unsigned int idaughter = 0; idaughter < jet.numberOfDaughters(); ++idaughter)
        {
            jetConstituentSet.insert(jet.daughterPtr(idaughter));
        }

        // Cut on eta
        if (std::abs(jet.eta()) > 2.4) continue;

        float uncorrectedPt = jet.correctedP4("Uncorrected").pt();
        //edm::RefToBase<reco::Jet> jet_ref(jets->refAt(ijet)); //upcast

        std::unordered_set<reco::CandidatePtr, CandidateHash> candidatesMatchedToSV;
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

        }

        std::vector<llpdnnx::ChargedCandidateFeatures> jet_tag_infos;

            
        // Fill cpf info
        for (unsigned int idaughter = 0; idaughter < jet.numberOfDaughters(); ++idaughter)
        {
            const pat::PackedCandidate* constituent = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(idaughter));
            if ((not constituent) or constituent->charge()==0 or (not constituent->hasTrackDetails()))
            {
               continue;
            }

            llpdnnx::ChargedCandidateFeatures cpf_features;

            cpf_features.cpf_ptrel = constituent->pt()/uncorrectedPt;
            cpf_features.cpf_deta = std::fabs(constituent->eta()-jet.eta());
            cpf_features.cpf_dphi = std::fabs(reco::deltaPhi(constituent->phi(),jet.phi()));

            cpf_features.cpf_drminsv = 0.4;
            for (const auto& sv: *svs.product())
            {
                float dR = reco::deltaR(sv,*constituent);
                cpf_features.cpf_drminsv = std::min(cpf_features.cpf_drminsv,dR);
            }

            cpf_features.cpf_vertex_association = constituent->pvAssociationQuality();
            cpf_features.cpf_fromPV = constituent->fromPV();
            cpf_features.cpf_puppi_weight = constituent->puppiWeight();
            cpf_features.cpf_track_chi2 = constituent->pseudoTrack().chi2();
            cpf_features.cpf_track_ndof = constituent->pseudoTrack().ndof();
            cpf_features.cpf_track_quality = constituent->pseudoTrack().qualityMask();
	        cpf_features.cpf_track_numberOfValidPixelHits = constituent->pseudoTrack().hitPattern().numberOfValidPixelHits();
	        cpf_features.cpf_track_pixelLayersWithMeasurement  = constituent->pseudoTrack().hitPattern().pixelLayersWithMeasurement();
	        cpf_features.cpf_track_numberOfValidStripHits = constituent->pseudoTrack().hitPattern().numberOfValidStripHits();
	        cpf_features.cpf_track_stripLayersWithMeasurement = constituent->pseudoTrack().hitPattern().stripLayersWithMeasurement();
		

            if (jet.mass()<1e-10)
            {
                cpf_features.cpf_relmassdrop = -1;
            }
            else
            {
                cpf_features.cpf_relmassdrop = (jet.p4()-constituent->p4()).mass()/jet.mass();
            }
            
            reco::TransientTrack transientTrack = builder->build(constituent->pseudoTrack());
            reco::Candidate::Vector jetDir = jet.momentum().Unit();
            GlobalVector jetRefTrackDir(jet.px(),jet.py(),jet.pz());

            Measurement1D meas_ip2d=IPTools::signedTransverseImpactParameter(transientTrack, jetRefTrackDir, pv).second;
            Measurement1D meas_ip3d=IPTools::signedImpactParameter3D(transientTrack, jetRefTrackDir, pv).second;
            Measurement1D jetdist=IPTools::jetTrackDistance(transientTrack, jetRefTrackDir, pv).second;
            reco::Candidate::Vector trackMom = constituent->pseudoTrack().momentum();
            double trackMag = std::sqrt(trackMom.Mag2());
            TVector3 trackMom3(trackMom.x(),trackMom.y(),trackMom.z());
            TVector3 jetDir3(jetDir.x(),jetDir.y(),jetDir.z());

            cpf_features.cpf_trackEtaRel=reco::btau::etaRel(jetDir, trackMom);
            cpf_features.cpf_trackPtRel=trackMom3.Perp(jetDir3);
            cpf_features.cpf_trackPPar=jetDir.Dot(trackMom);
            cpf_features.cpf_trackDeltaR=reco::deltaR(trackMom, jetDir);
            cpf_features.cpf_trackPtRatio=cpf_features.cpf_trackPtRel/trackMag;
            cpf_features.cpf_trackPParRatio=cpf_features.cpf_trackPPar/trackMag;

            cpf_features.cpf_trackSip2dVal=meas_ip2d.value();
            cpf_features.cpf_trackSip2dSig=meas_ip2d.significance();
            cpf_features.cpf_trackSip3dVal=meas_ip3d.value();
            cpf_features.cpf_trackSip3dSig=meas_ip3d.significance();
            if (std::isnan(cpf_features.cpf_trackSip2dSig) or std::isnan(cpf_features.cpf_trackSip3dSig))
            {
                cpf_features.cpf_trackSip2dSig=-1.;
                cpf_features.cpf_trackSip3dSig=-1.;
            }

            cpf_features.cpf_trackJetDistVal = jetdist.value();
            cpf_features.cpf_trackJetDistSig = jetdist.significance();

            cpf_features.cpf_matchedMuon = 0;
            cpf_features.cpf_matchedElectron = 0;
            
            if (candidatesMatchedToSV.find(jet.daughterPtr(idaughter))!=candidatesMatchedToSV.end())
            {
                cpf_features.cpf_matchedSV = 1;
            }
            else
            {
                cpf_features.cpf_matchedSV = 0;
            }

            // reco::CandidatePtr

            
            //find matching muons
            auto findMuon = muonMap.find(jet.daughterPtr(idaughter));  

            if (findMuon != muonMap.end())
            {
                const pat::Muon& muon = *findMuon->second;
                if (not muon.isGlobalMuon() or reco::deltaR(muon, jet) > 0.4) continue;
                cpf_features.cpf_matchedMuon = 1;
                cpf_features.matchedMuonPtr = findMuon->first;

            }


            //find matching electrons
            auto findElectron = electronMap.find(jet.daughterPtr(idaughter));  
            if(findElectron != electronMap.end())
            {
                const pat::Electron& electron = *findElectron->second;
                if(reco::deltaR(electron, jet) > 0.4) continue; 
                cpf_features.cpf_matchedElectron = 1;
                cpf_features.matchedElectronPtr = findElectron->first;
            }

            jet_tag_infos.emplace_back(cpf_features);


        } //end loop over charged consistuents
        
    
        std::stable_sort(jet_tag_infos.begin(), jet_tag_infos.end(),[](llpdnnx::ChargedCandidateFeatures d1, llpdnnx::ChargedCandidateFeatures d2)
        {
            if (d1.cpf_trackSip2dSig>0 and d2.cpf_trackSip2dSig>0)
            {
                return std::fabs(d1.cpf_trackSip2dSig)>std::fabs(d2.cpf_trackSip2dSig); //sort decreasing
            }
            else if (d1.cpf_trackSip2dSig<0 and d2.cpf_trackSip2dSig>0)
            {
                return false;
            }
            else if (d1.cpf_trackSip2dSig>0 and d2.cpf_trackSip2dSig<0)
            {
                return true;
            }
            else if (std::fabs(d1.cpf_drminsv-d2.cpf_drminsv)>std::numeric_limits<float>::epsilon())
            {
                return d1.cpf_drminsv<d2.cpf_drminsv; //sort increasing
            }
            else
            {
                return d1.cpf_ptrel>d2.cpf_ptrel;  //sort decreasing
            }
            
            return false;
        });

    output_tag_infos->emplace(std::make_pair(jet_ref, jet_tag_infos));
      
    }

    iEvent.put(std::move(output_tag_infos));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ChargedCandidateXTagInfoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
    desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
    desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("inclusiveCandidateSecondaryVertices"));
    desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
    desc.add<edm::InputTag>("electrons", edm::InputTag("slimmedElectrons"));
    descriptions.addDefault(desc);
}
void ChargedCandidateXTagInfoProducer::endStream() {};

//define this as a plug-in
DEFINE_FWK_MODULE(ChargedCandidateXTagInfoProducer);
