#include "DataFormats/Common/interface/Wrapper.h"
#include "LLPReco/DataFormats/interface/XTagFeatures.h"
#include "LLPReco/DataFormats/interface/XTagInfo.h"
#include "LLPReco/DataFormats/interface/DisplacedGenVertex.h"
#include "LLPReco/DataFormats/interface/LLPLabel.h"
#include "LLPReco/DataFormats/interface/LLPGenDecayInfo.h"
#include "LLPReco/DataFormats/interface/LLPGhostFlavourInfo.h"

namespace {

    struct dictionary
    {
        llpdnnx::JetFeatures JetFeaturesProduct;
        std::vector<llpdnnx::JetFeatures> JetFeaturesVectorProduct;
        edm::Wrapper<std::vector<llpdnnx::JetFeatures>> JetFeaturesVectorWrapperProduct;

        llpdnnx::ChargedCandidateFeatures ChargedCandidateFeaturesProduct;
        std::vector<llpdnnx::ChargedCandidateFeatures> ChargedCandidateFeaturesVectorProduct;
        edm::Wrapper<std::vector<llpdnnx::ChargedCandidateFeatures>> ChargedCandidateFeaturesVectorWrapperProduct;

        llpdnnx::NeutralCandidateFeatures NeutralCandidateFeaturesProduct;
        std::vector<llpdnnx::NeutralCandidateFeatures> NeutralCandidateFeaturesVectorProduct;
        edm::Wrapper<std::vector<llpdnnx::NeutralCandidateFeatures>> NeutralCandidateFeaturesVectorWrapperProduct;

        llpdnnx::SecondaryVertexFeatures SecondaryVertexFeaturesProduct;
        std::vector<llpdnnx::SecondaryVertexFeatures> SecondaryVertexFeaturesVectorProduct;
        edm::Wrapper<std::vector<llpdnnx::SecondaryVertexFeatures>> SecondaryVertexFeaturesVectorWrapperProduct;
        
        llpdnnx::MuonCandidateFeatures MuonCandidateFeaturesProduct;
        std::vector<llpdnnx::MuonCandidateFeatures> MuonCandidateFeaturesVectorProduct;
        edm::Wrapper<std::vector<llpdnnx::MuonCandidateFeatures>> MuonCandidateFeaturesVectorWrapperProduct;
        llpdnnx::ElectronCandidateFeatures ElectronCandidateFeaturesProduct;
        std::vector<llpdnnx::ElectronCandidateFeatures> ElectronCandidateFeaturesVectorProduct;
        edm::Wrapper<std::vector<llpdnnx::ElectronCandidateFeatures>> ElectronCandidateFeaturesVectorWrapperProduct;
        
        llpdnnx::DisplacedGenVertexCollection DisplacedGenVertexCollectionProduct;
        edm::Wrapper<llpdnnx::DisplacedGenVertexCollection> DisplacedGenVertexCollectionWrapperProduct;
        
        edm::Ptr<llpdnnx::DisplacedGenVertex> DisplacedGenVertexProduct;
        edm::Wrapper<edm::Ptr<llpdnnx::DisplacedGenVertex>> DisplacedGenVertexWrapperProduct;
        
        edm::Ptr<llpdnnx::DisplacedGenVertexCollection> DisplacedGenVertexPtrProduct;
        edm::Wrapper<edm::Ptr<llpdnnx::DisplacedGenVertexCollection>> DisplacedGenVertexPtrWrapperProduct;

        edm::PtrVector<llpdnnx::DisplacedGenVertexCollection> DisplacedGenVertexPtrVectorProduct;
        edm::Wrapper<edm::PtrVector<llpdnnx::DisplacedGenVertexCollection>> DisplacedGenVertexPtrVectorWrapperProduct;

        /*
        std::vector<reco::FeaturesTagInfo<llpdnnx::XTagFeatures>> dummy0;
        edm::Wrapper<std::vector<reco::FeaturesTagInfo<llpdnnx::XTagFeatures>>> dummy1;

        reco::FeaturesTagInfo<llpdnnx::XTagFeatures> dummy2;
        edm::Wrapper<reco::FeaturesTagInfo<llpdnnx::XTagFeatures>> dummy3;
      
        edm::PtrVector<reco::GenParticle> dummy17;
        edm::Wrapper<edm::PtrVector<reco::GenParticle>> dummy18;
        
        llpdnnx::LLPLabel dummy19;
        reco::FeaturesTagInfo<llpdnnx::LLPLabel> dummy20;

        std::vector<reco::FeaturesTagInfo<llpdnnx::LLPLabel>> dummy21;
        edm::Wrapper<std::vector<reco::FeaturesTagInfo<llpdnnx::LLPLabel>>> dummy22;
        
        llpdnnx::LLPGenDecayInfo dummy23;
        std::vector<llpdnnx::LLPGenDecayInfo> dummy24;
        edm::Wrapper<std::vector<llpdnnx::LLPGenDecayInfo>> dummy25;
        
        llpdnnx::LLPGhostFlavourInfo dummy26;
        std::vector<llpdnnx::LLPGhostFlavourInfo> dummy27;
        edm::Wrapper<std::vector<llpdnnx::LLPGhostFlavourInfo>> dummy28;
        edm::ValueMap<llpdnnx::LLPGhostFlavourInfo> dummy29;
        edm::Wrapper<edm::ValueMap<llpdnnx::LLPGhostFlavourInfo>> dummy30;

        llpdnnx::ElectronCandidateFeatures dummy31;
        llpdnnx::MuonCandidateFeatures dummy32;
        */

    };

}
