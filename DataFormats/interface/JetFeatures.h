#ifndef LLPReco_DataFormats_JetFeatures_h
#define LLPReco_DataFormats_JetFeatures_h

namespace llpdnnx {

struct JetFeatures 
{

    int jetIdx;
    float pt;
    float eta;
    float phi;
    float mass;
    float energy;
    
    float area;
    
    int n60;
    int n90;
    
    float chargedEmEnergyFraction;
    float chargedHadronEnergyFraction;
    float chargedMuEnergyFraction;
    float electronEnergyFraction;

    float tau1;
    float tau2;
    float tau3;
    
    float relMassDropMassAK;
    float relMassDropMassCA;
    float relSoftDropMassAK;
    float relSoftDropMassCA;
    
    float thrust;
    float sphericity;
    float circularity;
    float isotropy;
    float eventShapeC;
    float eventShapeD;

    float csv_trackSumJetEtRatio;      // ratio of track sum transverse energy over jet energy
    float csv_trackSumJetDeltaR;       // pseudoangular distance between jet axis and track fourvector sum
    int csv_vertexCategory;          // category of secondary vertex (Reco, Pseudo, No)
    float csv_trackSip2dValAboveCharm; // track 2D signed impact parameter of first track lifting mass above charm
    float csv_trackSip2dSigAboveCharm; // track 2D signed impact parameter significance of first track lifting mass above charm
    float csv_trackSip3dValAboveCharm; // track 3D signed impact parameter of first track lifting mass above charm
    float csv_trackSip3dSigAboveCharm; // track 3D signed impact parameter significance of first track lifting mass above charm
    // track info
    int csv_jetNTracksEtaRel; // tracks associated to jet for which trackEtaRel is calculated
    int csv_jetNSelectedTracks;
};

}

#endif 
