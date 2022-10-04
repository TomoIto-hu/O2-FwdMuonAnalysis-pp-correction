#include "ReconstructionDataFormats/GlobalFwdTrack.h"
#include "DataFormatsMFT/TrackMFT.h"

#include<iostream>
#include<memory>

#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "ReconstructionDataFormats/GlobalFwdTrack.h"
#include "ReconstructionDataFormats/MatchInfoFwd.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TTree.h"
#include <TGeoGlobalMagField.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>


double getZField(double x, double y, double z){
  const auto grp = o2::parameters::GRPObject::loadFrom("o2sim_grp.root");
  std::unique_ptr<o2::parameters::GRPObject> mGRP = nullptr;
  mGRP.reset(grp);
  o2::base::Propagator::initFieldFromGRP(grp);
  auto field = static_cast<o2::field::MagneticField *>(TGeoGlobalMagField::Instance()->GetField());

  double position[3] = {x, y, z};
  auto Bz = field->getBz(position);
  printf("B field z = %f [kGauss]\n", Bz);
  return Bz;
}

double getToXYplane(double x, double y, double z, double theta, double phi){
  auto dZ = (0.-z);
  auto xp = x + dZ*tan(theta)*cos(phi); //prospagated x coordinate
  auto yp = y + dZ*tan(theta)*sin(phi); //prospagated y coordinate
  auto DCAp = pow(pow(xp,2)+pow(yp,2),0.5); //DCA
  //printf("DCA = %f\n", DCAp);
  return DCAp;
  }

void FwdTrackHandsOn_MFT(const std::string trkFile = "globalfwdtracks.root",
                     const std::string o2sim_KineFile = "o2sim_Kine.root",
                     const std::string sig_KineFile = "sgn_Kine.root",
                     const std::string MFTtrkFile = "mfttracks.root",
                     const std::string MCHtrkFile = "mchtracks.root"){

  // Set global muon track (MFT + MCH) information
  //std::unique_ptr<TFile> trkFileIn(new TFile(trkFile.c_str()));
  //std::unique_ptr<TTree> gmTrackTree((TTree*)trkFileIn->Get("GlobalFwdTracks"));
  //std::vector<o2::dataformats::GlobalFwdTrack> trackGMVec, *trackGMVecP = &trackGMVec;
  //gmTrackTree->SetBranchAddress("fwdtracks", &trackGMVecP);
  //vector<o2::MCCompLabel>* mcLabels_gm = nullptr;
  //gmTrackTree->SetBranchAddress("MCTruth", &mcLabels_gm);
  //gmTrackTree->GetEntry(0);

  // Set MFT track information
  std::unique_ptr<TFile> MFTtrkFileIn(new TFile(MFTtrkFile.c_str()));
  std::unique_ptr<TTree> MFTTrackTree((TTree*)MFTtrkFileIn->Get("o2sim"));
  std::vector<o2::mft::TrackMFT> trackMFTVec, *trackMFTVecP = &trackMFTVec;
  MFTTrackTree->SetBranchAddress("MFTTrack", &trackMFTVecP);
  std::vector<o2::MCCompLabel>* mcLabels_mft = nullptr;
  MFTTrackTree->SetBranchAddress("MFTTrackMCTruth", &mcLabels_mft);
  MFTTrackTree->GetEntry(0);

  // Set MCH track MCH information
  //std::unique_ptr<TFile> MCHtrkFileIn(new TFile(MCHtrkFile.c_str()));
  //std::unique_ptr<TTree> MCHTrackTree((TTree*)MCHtrkFileIn->Get("o2sim"));
  //std::vector<o2::mch::TrackMCH> trackMCHVec, *trackMCHVecP = &trackMCHVec;
  //MCHTrackTree->SetBranchAddress("tracks", &trackMCHVecP);
  //vector<o2::MCCompLabel>* mcLabels_mch = nullptr;
  //MCHTrackTree->SetBranchAddress("tracklabels", &mcLabels_mch);
  //MCHTrackTree->GetEntry(0);

  // Set MC track information
  std::unique_ptr<TFile> o2sim_KineFileIn(new TFile(o2sim_KineFile.c_str()));
  std::unique_ptr<TTree> o2SimKineTree((TTree*)o2sim_KineFileIn->Get("o2sim"));
  vector<o2::MCTrackT<float>>* mcTr = nullptr;
  o2SimKineTree->SetBranchAddress("MCTrack", &mcTr);
  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  o2SimKineTree->SetBranchAddress("MCEventHeader.", &eventHeader);

  //Set global run condition
  auto field_z = getZField(0, 0, -61.4); // Get field at Center of MFT

  //Prepare output objects
  TFile* output = new TFile("AnalysisResults_MFT.root","recreate");


  std::unique_ptr<TTree> TotalTreeObject_MFT(new TTree("TotalTreeObject_MFT","TotalTreeObject_MFT"));
  double P_MFT,Pt_MFT,Eta_MFT,Theta_MFT,Phi_MFT,vx_MFT,vy_MFT,vz_MFT,DfIP_MFT;
  double P_MCMFT,Pt_MCMFT,Eta_MCMFT,Theta_MCMFT,Phi_MCMFT,vx_MCMFT,vy_MCMFT,vz_MCMFT,DfIP_MCMFT,MotherId_MCMFT;
  //double P_GrandMCGM,Pt_GrandMCGM,Eta_GrandMCGM,vx_GrandMCGM,vy_GrandMCGM,vz_GrandMCGM,DfIP_GrandMCGM;
  int PDG_MCMFT,Mother_MCMFT,GrandMother_MCMFT;
  int SourceCheck;
  TotalTreeObject_MFT -> Branch("P_MFT",&P_MFT,"P_MFT/D");
  TotalTreeObject_MFT -> Branch("Pt_MFT",&Pt_MFT,"Pt_MFT/D");
  TotalTreeObject_MFT -> Branch("Eta_MFT",&Eta_MFT,"Eta_MFT/D");
  TotalTreeObject_MFT -> Branch("Theta_MFT",&Theta_MFT,"Theta_MFT/D");
  TotalTreeObject_MFT -> Branch("Phi_MFT",&Phi_MFT,"Phi_MFT/D");
  TotalTreeObject_MFT -> Branch("vx_MFT",&vx_MFT,"vx_MFT/D");
  TotalTreeObject_MFT -> Branch("vy_MFT",&vy_MFT,"vy_MFT/D");
  TotalTreeObject_MFT -> Branch("vz_MFT",&vz_MFT,"vz_MFT/D");
  TotalTreeObject_MFT -> Branch("DfIP_MFT",&DfIP_MFT,"DfIP_MFT/D");
  TotalTreeObject_MFT -> Branch("PDG_MCMFT",&PDG_MCMFT,"PDG_MCMFT/I");
  TotalTreeObject_MFT -> Branch("vx_MCMFT",&vx_MCMFT,"vx_MCMFT/D");
  TotalTreeObject_MFT -> Branch("vy_MCMFT",&vy_MCMFT,"vy_MCMFT/D");
  TotalTreeObject_MFT -> Branch("vz_MCMFT",&vz_MCMFT,"vz_MCMFT/D");
  TotalTreeObject_MFT -> Branch("P_MCMFT",&P_MCMFT,"P_MCMFT/D");
  TotalTreeObject_MFT -> Branch("Pt_MCMFT",&Pt_MCMFT,"Pt_MCMFT/D");
  TotalTreeObject_MFT -> Branch("Eta_MCMFT",&Eta_MCMFT,"Eta_MCMFT/D");
  TotalTreeObject_MFT -> Branch("Theta_MCMFT",&Theta_MCMFT,"Theta_MCMFT/D");
  TotalTreeObject_MFT -> Branch("Phi_MCMFT",&Phi_MCMFT,"Phi_MCMFT/D");
  TotalTreeObject_MFT -> Branch("DfIP_MCMFT",&DfIP_MCMFT,"DfIP_MCMFT/D");
  TotalTreeObject_MFT -> Branch("MotherId_MCMFT",&MotherId_MCMFT,"MotherId_MCMFT/I");
  TotalTreeObject_MFT -> Branch("GrandMother_MCMFT",&GrandMother_MCMFT,"GrandMother_MCMFT/I");
  TotalTreeObject_MFT -> Branch("SourceCheck",&SourceCheck,"SourceCheck/I");
  TotalTreeObject_MFT -> Branch("Mother_MCMFT",&Mother_MCMFT,"Mother_MCMFT/I");

  int nTracks = mcLabels_mft->size();
  int Count_miss=0;
  int pdgcode_MC;
  std::cout<<"n Track = "<<nTracks<<std::endl;

  for (int iTrack = 0;iTrack<nTracks; ++iTrack) {
    double P_MC,Pt_MC,vx_MC,vy_MC,vz_MC,Eta_MC,Theta_MC,Phi_MC;
    double P_GrandMC,Pt_GrandMC,vx_GrandMC,vy_GrandMC,vz_GrandMC,Eta_GrandMC,Theta_GrandMC,Phi_GrandMC;
    int SaveId,MotherId,Mother_MC,GrandMotherId,GrandMother_MC;
    SourceCheck=0;
    auto& mftTrack = trackMFTVecP->at(iTrack);
    //Get label information
    const auto& label_mft = mcLabels_mft->at(iTrack);

    //Get MC track id
    auto trkID = label_mft.getTrackID();
    //Get Event id
    auto evtID = label_mft.getEventID();

    //Load MC track information
    o2::MCTrackT<float>* thisTrack;
    if (label_mft.getSourceID() == 0) {
      SourceCheck = 1;
      o2SimKineTree->GetEntry(evtID);
      thisTrack = &(mcTr->at(trkID));

    //MC track information
    pdgcode_MC = thisTrack->GetPdgCode();//get particle pdg code
    //Propagate global muon track to primary vertex position z
    //gmTrack.propagateToZquadratic(0, field_z);
    //gmTrack.propagateToZlinear(0);
    mftTrack.propagateToZhelix(0, field_z);

    //Mother Track
    SaveId =  trkID;
    MotherId = trkID;
    Mother_MC = pdgcode_MC;
    while(fabs(pdgcode_MC) == fabs(Mother_MC)){
    	SaveId = MotherId;
      vx_MC = thisTrack->GetStartVertexCoordinatesX();//get particle production x-posiiton
    	vy_MC = thisTrack->GetStartVertexCoordinatesY();//get particle production y-posiiton
    	vz_MC = thisTrack->GetStartVertexCoordinatesZ();//get particle production z-posiiton
    	Pt_MC = thisTrack->GetPt();//get particle pt
    	P_MC = thisTrack->GetP();//get particle total momentum
    	Eta_MC = thisTrack->GetEta();
    	Theta_MC = thisTrack->GetTheta();
      Phi_MC = thisTrack->GetPhi();
    	MotherId = thisTrack->getMotherTrackId();
    	if(MotherId==-1) break;
    	thisTrack = &(mcTr->at(MotherId));
    	Mother_MC = thisTrack->GetPdgCode();
    }
    	GrandMotherId = MotherId;
    	GrandMother_MC = Mother_MC;
    	while(fabs(GrandMother_MC) == fabs(Mother_MC)){
		GrandMotherId = thisTrack->getMotherTrackId();
    		if(GrandMotherId==-1) break;
    		thisTrack = &(mcTr->at(GrandMotherId));
    		GrandMother_MC = thisTrack->GetPdgCode();
    	}
      }else{
      //gmTrack.propagateToZquadratic(0, field_z);
      //gmTrack.propagateToZlinear(0);
      mftTrack.propagateToZhelix(0, field_z);
    }
    //Mother_MC = thisTrack->GetPdgCode();
    //Global muon track information
    auto MFTEta = mftTrack.o2::track::TrackParFwd::getEta();//get pseudo rapidity
    auto MFTPt = mftTrack.o2::track::TrackParFwd::getPt();//get transverse momentum
    auto MFTP = mftTrack.o2::track::TrackParFwd::getP();//get total momentum
    auto MFTPhi = mftTrack.o2::track::TrackParFwd::getPhi();//get phi
    auto MFTTheta = mftTrack.o2::track::TrackParFwd::getTheta();
    auto MFTvx = mftTrack.o2::track::TrackParFwd::getX();//get X
    auto MFTvy = mftTrack.o2::track::TrackParFwd::getY();//get Y
    auto MFTvz = mftTrack.o2::track::TrackParFwd::getZ();//get Z


      P_MFT = MFTP;
      Pt_MFT = MFTPt;
      vx_MFT = MFTvx;
      vy_MFT = MFTvy;
      vz_MFT = MFTvz;
      Eta_MFT = MFTEta;
      Theta_MFT = MFTTheta;
      Phi_MFT = MFTPhi;
      PDG_MCMFT = pdgcode_MC;
      DfIP_MFT = pow(pow(MFTvx,2)+pow(MFTvy,2)+pow(MFTvz,2),0.5);

      P_MCMFT = P_MC;
      Pt_MCMFT = Pt_MC;
      vx_MCMFT = vx_MC;
      vy_MCMFT = vy_MC;
      vz_MCMFT = vz_MC;
      Eta_MCMFT = Eta_MC;
      Theta_MCMFT = Theta_MC;
      Phi_MCMFT = Phi_MC;
      //DfIP_MCMFT = getToXYplane(vx_MC,vy_MC,vz_MC,Theta_MC,Phi_MC);
      Mother_MCMFT = Mother_MC;
      MotherId_MCMFT = MotherId;
      //cout<<"Mother = "<<Mother_MC<<endl;
      //cout<<"Mother_MCGM = "<<Mother_MCGM<<endl;

      //DfIP_GrandMCGM = getToXYplane(vx_GrandMC,vy_GrandMC,vz_GrandMC,Theta_GrandMC,Phi_GrandMC);
      GrandMother_MCMFT = GrandMother_MC;
      TotalTreeObject_MFT->Fill();
  }
  //Save histograms

  output->WriteTObject(TotalTreeObject_MFT.get());

  std::cout<<"Analysis complete"<<std::endl;
}
