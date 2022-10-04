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
  printf("DCA = %f\n", DCAp);
  return DCAp;
  }

void GMTrackHandsOn_Comp(const std::string trkFile = "globalfwdtracks.root",
                     const std::string o2sim_KineFile = "o2sim_Kine.root",
                     const std::string sig_KineFile = "sgn_Kine.root",
                     const std::string MFTtrkFile = "AnalysisResults_MFT.root",
                     const std::string MCHtrkFile = "mchtracks.root"){

  // Set global muon track (MFT + MCH) information
  std::unique_ptr<TFile> trkFileIn(new TFile(trkFile.c_str()));
  std::unique_ptr<TTree> gmTrackTree((TTree*)trkFileIn->Get("GlobalFwdTracks"));
  std::vector<o2::dataformats::GlobalFwdTrack> trackGMVec, *trackGMVecP = &trackGMVec;
  gmTrackTree->SetBranchAddress("fwdtracks", &trackGMVecP);
  vector<o2::MCCompLabel>* mcLabels_gm = nullptr;
  gmTrackTree->SetBranchAddress("MCTruth", &mcLabels_gm);
  gmTrackTree->GetEntry(0);

  // Set MFT track information
  std::unique_ptr<TFile> MFTtrkFileIn(new TFile(MFTtrkFile.c_str()));
  std::unique_ptr<TTree> MFTTrackTree((TTree*)MFTtrkFileIn->Get("TotalTreeObject_MFT"));

  // Set MCH track MCH information
  std::unique_ptr<TFile> MCHtrkFileIn(new TFile(MCHtrkFile.c_str()));
  std::unique_ptr<TTree> MCHTrackTree((TTree*)MCHtrkFileIn->Get("o2sim"));
  std::vector<o2::mch::TrackMCH> trackMCHVec, *trackMCHVecP = &trackMCHVec;
  MCHTrackTree->SetBranchAddress("tracks", &trackMCHVecP);
  vector<o2::MCCompLabel>* mcLabels_mch = nullptr;
  MCHTrackTree->SetBranchAddress("tracklabels", &mcLabels_mch);
  MCHTrackTree->GetEntry(0);

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
  TFile* output = new TFile("AnalysisResults_FwdDet.root","recreate");
  std::unique_ptr<TH2F> fHistGMTrackEtaPt(new TH2F("fHistGMTrackEtaPt","",150,-4.0,-2.5,100,0,10));
  std::unique_ptr<TH2F> fHistGMTrackEtaPt_CorrectMatch(new TH2F("fHistGMTrackEtaPt_CorrectMatch","",150,-4.0,-2.5,100,0,10));

  std::unique_ptr<TTree> TotalTreeObject_GM(new TTree("TotalTreeObject_GM","TotalTreeObject_GM"));
  double P_GM,Pt_GM,Eta_GM,Theta_GM,Phi_GM,vx_GM,vy_GM,vz_GM,DfIP_GM;
  double Eta_p_GM,Theta_p_GM,Phi_p_GM,vx_p_GM,vy_p_GM,vz_p_GM,invTanl_p_GM,Qpt_p_GM,Hz_p_GM,k_p_GM;
  double P_MCGM,Pt_MCGM,Eta_MCGM,Theta_MCGM,Phi_MCGM,vx_MCGM,vy_MCGM,vz_MCGM,DfIP_MCGM,MotherId_MCGM;
  int PDG_MCGM,Mother_MCGM,GrandMother_MCGM;
  int SourceCheck,CorrectCheck;
  double MIDCheck,MatchingScore;
  TotalTreeObject_GM -> Branch("Eta_p_GM",&Eta_p_GM,"Eta_p_GM/D");
  TotalTreeObject_GM -> Branch("Theta_p_GM",&Theta_p_GM,"Theta_p_GM/D");
  TotalTreeObject_GM -> Branch("Phi_p_GM",&Phi_p_GM,"Phi_p_GM/D");
  TotalTreeObject_GM -> Branch("vx_p_GM",&vx_p_GM,"vx_p_GM/D");
  TotalTreeObject_GM -> Branch("vy_p_GM",&vy_p_GM,"vy_p_GM/D");
  TotalTreeObject_GM -> Branch("vz_p_GM",&vz_p_GM,"vz_p_GM/D");
  TotalTreeObject_GM -> Branch("invTanl_p_GM",&invTanl_p_GM,"invTanl_p_GM/D");
  TotalTreeObject_GM -> Branch("Qpt_p_GM",&Qpt_p_GM,"Qpt_p_GM/D");
  TotalTreeObject_GM -> Branch("k_p_GM",&k_p_GM,"k_p_GM/D");
  TotalTreeObject_GM -> Branch("Hz_p_GM",&Hz_p_GM,"Hz_p_GM/D");
  TotalTreeObject_GM -> Branch("P_GM",&P_GM,"P_GM/D");
  TotalTreeObject_GM -> Branch("Pt_GM",&Pt_GM,"Pt_GM/D");
  TotalTreeObject_GM -> Branch("Eta_GM",&Eta_GM,"Eta_GM/D");
  TotalTreeObject_GM -> Branch("Theta_GM",&Theta_GM,"Theta_GM/D");
  TotalTreeObject_GM -> Branch("Phi_GM",&Phi_GM,"Phi_GM/D");
  TotalTreeObject_GM -> Branch("vx_GM",&vx_GM,"vx_GM/D");
  TotalTreeObject_GM -> Branch("vy_GM",&vy_GM,"vy_GM/D");
  TotalTreeObject_GM -> Branch("vz_GM",&vz_GM,"vz_GM/D");
  TotalTreeObject_GM -> Branch("DfIP_GM",&DfIP_GM,"DfIP_GM/D");
  TotalTreeObject_GM -> Branch("PDG_MCGM",&PDG_MCGM,"PDG_MCGM/I");
  TotalTreeObject_GM -> Branch("vx_MCGM",&vx_MCGM,"vx_MCGM/D");
  TotalTreeObject_GM -> Branch("vy_MCGM",&vy_MCGM,"vy_MCGM/D");
  TotalTreeObject_GM -> Branch("vz_MCGM",&vz_MCGM,"vz_MCGM/D");
  TotalTreeObject_GM -> Branch("P_MCGM",&P_MCGM,"P_MCGM/D");
  TotalTreeObject_GM -> Branch("Pt_MCGM",&Pt_MCGM,"Pt_MCGM/D");
  TotalTreeObject_GM -> Branch("Eta_MCGM",&Eta_MCGM,"Eta_MCGM/D");
  TotalTreeObject_GM -> Branch("Theta_MCGM",&Theta_MCGM,"Theta_MCGM/D");
  TotalTreeObject_GM -> Branch("Phi_MCGM",&Phi_MCGM,"Phi_MCGM/D");
  TotalTreeObject_GM -> Branch("DfIP_MCGM",&DfIP_MCGM,"DfIP_MCGM/D");
  TotalTreeObject_GM -> Branch("MotherId_MCGM",&MotherId_MCGM,"MotherId_MCGM/I");
  TotalTreeObject_GM -> Branch("GrandMother_MCGM",&GrandMother_MCGM,"GrandMother_MCGM/I");
  TotalTreeObject_GM -> Branch("SourceCheck",&SourceCheck,"SourceCheck/I");
  TotalTreeObject_GM -> Branch("CorrectCheck",&CorrectCheck,"CorrectCheck/I");
  TotalTreeObject_GM -> Branch("Mother_MCGM",&Mother_MCGM,"Mother_MCGM/I");
  TotalTreeObject_GM -> Branch("MIDCheck",&MIDCheck,"MIDCheck/D");
  TotalTreeObject_GM -> Branch("MatchingScore",&MatchingScore,"MatchingScore/D");

  //Get MFT Branch
  double P_MFT,Pt_MFT,DfIP_MFT,vz_MFT,vy_MFT,vx_MFT,Eta_MFT;
  double P_MCMFT,Pt_MCMFT,DfIP_MCMFT,vz_MCMFT,vy_MCMFT,vx_MCMFT,Eta_MCMFT;
  double Theta_MFT,Theta_MCMFT,Phi_MFT,Phi_MCMFT;
  double P_GrandMCMFT,Pt_GrandMCMFT,DfIP_GrandMCMFT,vz_GrandMCMFT,vy_GrandMCMFT,vx_GrandMCMFT,Eta_GrandMCMFT;
  int MotherId_MCMFT;
  int Mother_MCMFT,PDG_MCMFT,GrandMother_MCMFT;
  MFTTrackTree ->SetBranchAddress("DfIP_MFT",&DfIP_MFT);
  MFTTrackTree ->SetBranchAddress("P_MFT",&P_MFT);
  MFTTrackTree ->SetBranchAddress("Pt_MFT",&Pt_MFT);
  MFTTrackTree ->SetBranchAddress("vx_MFT",&vx_MFT);
  MFTTrackTree ->SetBranchAddress("vy_MFT",&vy_MFT);
  MFTTrackTree ->SetBranchAddress("vz_MFT",&vz_MFT);
  MFTTrackTree ->SetBranchAddress("Eta_MFT",&Eta_MFT);
  MFTTrackTree ->SetBranchAddress("Theta_MFT",&Theta_MFT);
  MFTTrackTree ->SetBranchAddress("Phi_MFT",&Phi_MFT);
  MFTTrackTree ->SetBranchAddress("DfIP_MCMFT",&DfIP_MCMFT);
  MFTTrackTree ->SetBranchAddress("P_MCMFT",&P_MCMFT);
  MFTTrackTree ->SetBranchAddress("Pt_MCMFT",&Pt_MCMFT);
  MFTTrackTree ->SetBranchAddress("vx_MCMFT",&vx_MCMFT);
  MFTTrackTree ->SetBranchAddress("vy_MCMFT",&vy_MCMFT);
  MFTTrackTree ->SetBranchAddress("vz_MCMFT",&vz_MCMFT);
  MFTTrackTree ->SetBranchAddress("Eta_MCMFT",&Eta_MCMFT);
  MFTTrackTree ->SetBranchAddress("Theta_MCMFT",&Theta_MCMFT);
  MFTTrackTree ->SetBranchAddress("Phi_MCMFT",&Phi_MCMFT);
  MFTTrackTree ->SetBranchAddress("PDG_MCMFT",&PDG_MCMFT);
  MFTTrackTree ->SetBranchAddress("Mother_MCMFT",&Mother_MCMFT);
  MFTTrackTree ->SetBranchAddress("MotherId_MCMFT",&MotherId_MCMFT);
  MFTTrackTree ->SetBranchAddress("GrandMother_MCMFT",&GrandMother_MCMFT);
  //Set MFT Branch
  TotalTreeObject_GM -> Branch("P_MFT",&P_MFT,"P_MFT/D");
  TotalTreeObject_GM -> Branch("Pt_MFT",&Pt_MFT,"Pt_MFT/D");
  TotalTreeObject_GM -> Branch("Eta_MFT",&Eta_MFT,"Eta_MFT/D");
  TotalTreeObject_GM -> Branch("Theta_MFT",&Theta_MFT,"Theta_MFT/D");
  TotalTreeObject_GM -> Branch("Phi_MFT",&Phi_MFT,"Phi_MFT/D");
  TotalTreeObject_GM -> Branch("vx_MFT",&vx_MFT,"vx_MFT/D");
  TotalTreeObject_GM -> Branch("vy_MFT",&vy_MFT,"vy_MFT/D");
  TotalTreeObject_GM -> Branch("vz_MFT",&vz_MFT,"vz_MFT/D");
  TotalTreeObject_GM -> Branch("DfIP_MFT",&DfIP_MFT,"DfIP_MFT/D");
  TotalTreeObject_GM -> Branch("PDG_MCMFT",&PDG_MCMFT,"PDG_MCMFT/I");
  TotalTreeObject_GM -> Branch("vx_MCMFT",&vx_MCMFT,"vx_MCMFT/D");
  TotalTreeObject_GM -> Branch("vy_MCMFT",&vy_MCMFT,"vy_MCMFT/D");
  TotalTreeObject_GM -> Branch("vz_MCMFT",&vz_MCMFT,"vz_MCMFT/D");
  TotalTreeObject_GM -> Branch("P_MCMFT",&P_MCMFT,"P_MCMFT/D");
  TotalTreeObject_GM -> Branch("Pt_MCMFT",&Pt_MCMFT,"Pt_MCMFT/D");
  TotalTreeObject_GM -> Branch("Eta_MCMFT",&Eta_MCMFT,"Eta_MCMFT/D");
  TotalTreeObject_GM -> Branch("Theta_MCMFT",&Theta_MCMFT,"Theta_MCMFT/D");
  TotalTreeObject_GM -> Branch("Phi_MCMFT",&Phi_MCMFT,"Phi_MCMFT/D");
  TotalTreeObject_GM -> Branch("DfIP_MCMFT",&DfIP_MCMFT,"DfIP_MCMFT/D");
  TotalTreeObject_GM -> Branch("MotherId_MCMFT",&MotherId_MCMFT,"MotherId_MCMFT/I");
  TotalTreeObject_GM -> Branch("GrandMother_MCMFT",&GrandMother_MCMFT,"GrandMother_MCMFT/I");
  TotalTreeObject_GM -> Branch("Mother_MCMFT",&Mother_MCMFT,"Mother_MCMFT/I");

  //MCH Branch
  double P_MCH,Pt_MCH,Eta_MCH,vx_MCH,vy_MCH,vz_MCH,DfIP_MCH,Px_MCH,Py_MCH,Pz_MCH,Theta_MCH,Phi_MCH;
  //double P_MCMCH,Pt_MCMCH,Eta_MCMCH,vx_MCMCH,vy_MCMCH,vz_MCMCH,DfIP_MCMCH,MotherId_MCMCH;
  TotalTreeObject_GM -> Branch("P_MCH",&P_MCH,"P_MCH/D");
  TotalTreeObject_GM -> Branch("Pt_MCH",&Pt_MCH,"Pt_MCH/D");
  TotalTreeObject_GM -> Branch("Eta_MCH",&Eta_MCH,"Eta_MCH/D");
  TotalTreeObject_GM -> Branch("vx_MCH",&vx_MCH,"vx_MCH/D");
  TotalTreeObject_GM -> Branch("vy_MCH",&vy_MCH,"vy_MCH/D");
  TotalTreeObject_GM -> Branch("vz_MCH",&vz_MCH,"vz_MCH/D");
  TotalTreeObject_GM -> Branch("DfIP_MCH",&DfIP_MCH,"DfIP_MCH/D");
  TotalTreeObject_GM -> Branch("Px_MCH",&Px_MCH,"Px_MCH/D");
  TotalTreeObject_GM -> Branch("Py_MCH",&Py_MCH,"Py_MCH/D");
  TotalTreeObject_GM -> Branch("Pz_MCH",&Pz_MCH,"Pz_MCH/D");
  TotalTreeObject_GM -> Branch("Theta_MCH",&Theta_MCH,"Theta_MCH/D");
  TotalTreeObject_GM -> Branch("Phi_MCH",&Phi_MCH,"Phi_MCH/D");

  int nTracks = mcLabels_gm->size();
  int Count_miss=0;
  int pdgcode_MC;
  std::cout<<"n Track = "<<nTracks<<std::endl;

  std::vector<int> MatchingSaver(1);
  std::vector<double> MatchingChecker(1);
  MatchingSaver[0] = -1;
  MatchingChecker[0] = -1;
  int mchTrackId;
  double Chi2Check;
  for(int iRoop =0;iRoop<nTracks; ++iRoop){
    auto& gmTrack = trackGMVecP->at(iRoop);
    mchTrackId = gmTrack.getMCHTrackID();
    Chi2Check = gmTrack.getMFTMCHMatchingChi2();
    std::cout<<"Roop = "<<iRoop<<std::endl;
    std::cout<<"mchTrackId = "<<mchTrackId<<std::endl;
    std::cout<<"Chi2 = "<<Chi2Check<<std::endl;
    MatchingSaver.push_back(mchTrackId);
    MatchingChecker.push_back(Chi2Check);
    for(int i = 1;i<MatchingSaver.size();i++){
      std::cout<<"MatchingSaver save : "<<MatchingSaver[i]<<std::endl;
      std::cout<<"MatchingChecker save : "<<MatchingChecker[i]<<std::endl;
      if(mchTrackId == MatchingSaver[i] && i!=iRoop+1) MatchingSaver[iRoop+1]<MatchingSaver[i]? MatchingChecker[iRoop+1] = -1 : MatchingChecker[i] = -1;
    }
    std::cout<<"mchTrackId = "<<MatchingSaver[iRoop+1]<<std::endl;
    std::cout<<"Chi2 = "<<MatchingChecker[iRoop+1]<<std::endl;
  }

  for (int iTrack = 0;iTrack<nTracks; ++iTrack) {
    if(MatchingChecker[iTrack+1]==-1) continue;
    double P_MC,Pt_MC,vx_MC,vy_MC,vz_MC,Eta_MC,Theta_MC,Phi_MC;
    double P_GrandMC,Pt_GrandMC,vx_GrandMC,vy_GrandMC,vz_GrandMC,Eta_GrandMC,Theta_GrandMC,Phi_GrandMC;
    int SaveId,MotherId,Mother_MC,GrandMotherId,GrandMother_MC;
    CorrectCheck=0;
    SourceCheck=0;
    //std::cout<<"now Track"<<iTrack<<std::endl;
    auto& gmTrack = trackGMVecP->at(iTrack);
    auto& mchTrack = trackMCHVecP->at(MatchingSaver[iTrack+1]);
    //std::cout<<"Total GMTrack = "<<trackGMVecP->size()<<std::endl;
    //std::cout<<"Total MCHTrack = "<<trackMCHVecP->size()<<std::endl;
    //Get label information
    const auto& label_gm = mcLabels_gm->at(iTrack);
    const auto& label_mch = mcLabels_mch->at(MatchingSaver[iTrack+1]);
    //Get Matching information
    auto mftTrackId = gmTrack.getMFTTrackID();
    auto mchTrackId = gmTrack.getMCHTrackID();
    auto midTrackId = gmTrack.getMIDTrackID();
    MFTTrackTree->GetEntry(mftTrackId);

    //Get MC track id
    auto trkID = label_gm.getTrackID();
    //Get Event id
    auto evtID = label_gm.getEventID();


    //Load MC track information
    o2::MCTrackT<float>* thisTrack;
    if (label_gm.getSourceID() == 0) {
      SourceCheck = 1;
      o2SimKineTree->GetEntry(evtID);
      thisTrack = &(mcTr->at(trkID));

    //MC track information
    pdgcode_MC = thisTrack->GetPdgCode();//get particle pdg code

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
      }
     //Global track informmation for propagation
      auto GMEta_p = gmTrack.o2::track::TrackParFwd::getEta();//get pseudo rapidity
      auto GMPhi_p = gmTrack.o2::track::TrackParFwd::getPhi();//get phi
      auto GMTheta_p = gmTrack.o2::track::TrackParFwd::getTheta();
      auto GMvx_p = gmTrack.o2::track::TrackParFwd::getX();//get X
      auto GMvy_p = gmTrack.o2::track::TrackParFwd::getY();//get Y
      auto GMvz_p = gmTrack.o2::track::TrackParFwd::getZ();//get Z
      auto GMTanl_p = gmTrack.o2::track::TrackParFwd::getTanl();//get tan long
      auto GMinvTanl_p = 1.0/GMTanl_p;
      auto GMinvQpt_p = gmTrack.o2::track::TrackParFwd::getInvQPt();
      auto GMQpt_p = 1.0/GMinvQpt_p;
      auto GMk_p = TMath::Abs(o2::constants::math::B2C*field_z);
      auto GMHz_p = std::copysign(1,field_z);
      //Propagate global muon track to primary vertex position z
    gmTrack.propagateToZhelix(0, field_z);
    //Mother_MC = thisTrack->GetPdgCode();
    //Global muon track information
    auto GMEta = gmTrack.o2::track::TrackParFwd::getEta();//get pseudo rapidity
    auto GMPt = gmTrack.o2::track::TrackParFwd::getPt();//get transverse momentum
    auto GMP = gmTrack.o2::track::TrackParFwd::getP();//get total momentum
    auto GMPhi = gmTrack.o2::track::TrackParFwd::getPhi();//get phi
    auto GMTheta = gmTrack.o2::track::TrackParFwd::getTheta();
    auto GMvx = gmTrack.o2::track::TrackParFwd::getX();//get X
    auto GMvy = gmTrack.o2::track::TrackParFwd::getY();//get Y
    auto GMvz = gmTrack.o2::track::TrackParFwd::getZ();//get Z
    MIDCheck = gmTrack.getMIDMatchingChi2();//get MCH-MID matching Chi2
    MatchingScore = gmTrack.getMFTMCHMatchingChi2();//get matching score
    //std::cout<<"MID Chi2 = "<<MIDCheck<<std::endl;
    //std::cout<<"Matching Score = "<<MatchingScore<<std::endl;
    //Fill all pair results
    fHistGMTrackEtaPt->Fill(GMEta,GMPt);

    auto MCHPx = mchTrack.o2::mch::TrackMCH::getPx();//get py
    auto MCHPy = mchTrack.o2::mch::TrackMCH::getPy();//get Py
    auto MCHPz = mchTrack.o2::mch::TrackMCH::getPz();//get Pz
    auto MCHP = mchTrack.o2::mch::TrackMCH::getP();//get total momentum
    auto MCHvx = mchTrack.o2::mch::TrackMCH::getX();//get X
    auto MCHvy = mchTrack.o2::mch::TrackMCH::getY();//get Y
    auto MCHvz = mchTrack.o2::mch::TrackMCH::getZ();//get Z

    auto MCHPt = TMath::Sqrt(pow(MCHPx,2)+pow(MCHPy,2));
    auto MCHLamda = TMath::ATan(MCHPz/MCHPt);
    auto MCHTheta = TMath::PiOver2() - MCHLamda;
    auto MCHEta = -TMath::Log(TMath::Tan(MCHTheta/2));
    auto MCHPhi = TMath::ATan(MCHPy/MCHPx);

    //std::cout<<"MCH vertex z = "<<MCHvz<<std::endl;

    //Fill correct match results
    if(label_gm.getSourceID() == 0){
      if (label_gm.isCorrect()) {
	CorrectCheck = 1;
        fHistGMTrackEtaPt_CorrectMatch->Fill(GMEta,GMPt);
      }
    }
      //GMTrack for propagation
      vx_p_GM = GMvx_p;
      vy_p_GM = GMvy_p;
      vz_p_GM = GMvz_p;
      Eta_p_GM = GMEta_p;
      Theta_p_GM = GMTheta_p;
      Phi_p_GM = GMPhi_p;
      invTanl_p_GM = GMinvTanl_p;
      Qpt_p_GM = GMQpt_p;
      k_p_GM = GMk_p;
      Hz_p_GM = GMHz_p;
      //GMTrack
      P_GM = GMP;
      Pt_GM = GMPt;
      vx_GM = GMvx;
      vy_GM = GMvy;
      vz_GM = GMvz;
      Eta_GM = GMEta;
      Theta_GM = GMTheta;
      Phi_GM = GMPhi;
      PDG_MCGM = pdgcode_MC;
      DfIP_GM = pow(pow(GMvx,2)+pow(GMvy,2)+pow(GMvz,2),0.5);

      //MCH
      P_MCH = MCHP;
      Pt_MCH = MCHPt;
      vx_MCH = MCHvx;
      vy_MCH = MCHvy;
      vz_MCH = MCHvz;
      Eta_MCH = MCHEta;
      Px_MCH = MCHPx;
      Py_MCH = MCHPy;
      Pz_MCH = MCHPz;
      Phi_MCH = MCHPhi;
      Theta_MCH = MCHTheta;
      DfIP_MCH = getToXYplane(vx_MCH,vy_MCH,vz_MCH,Theta_MCH,Phi_MCH);

      P_MCGM = P_MC;
      Pt_MCGM = Pt_MC;
      vx_MCGM = vx_MC;
      vy_MCGM = vy_MC;
      vz_MCGM = vz_MC;
      Eta_MCGM = Eta_MC;
      Theta_MCGM = Theta_MC;
      Phi_MCGM = Phi_MC;
      DfIP_MCGM = getToXYplane(vx_MC,vy_MC,vz_MC,Theta_MC,Phi_MC);
      Mother_MCGM = Mother_MC;
      MotherId_MCGM = MotherId;
      //cout<<"Mother = "<<Mother_MC<<endl;
      //cout<<"Mother_MCGM = "<<Mother_MCGM<<endl;

      //P_GrandMCGM = P_GrandMC;
      //Pt_GrandMCGM = Pt_GrandMC;
      //vx_GrandMCGM = vx_GrandMC;
      //vy_GrandMCGM = vy_GrandMC;
      //vz_GrandMCGM = vz_GrandMC;
      //Eta_GrandMCGM = Eta_GrandMC;
      //DfIP_GrandMCGM = getToXYplane(vx_GrandMC,vy_GrandMC,vz_GrandMC,Theta_GrandMC,Phi_GrandMC);
      GrandMother_MCGM = GrandMother_MC;
      TotalTreeObject_GM->Fill();
  }
  //Save histograms
  output->WriteTObject(fHistGMTrackEtaPt.get());
  output->WriteTObject(fHistGMTrackEtaPt_CorrectMatch.get());
  output->WriteTObject(TotalTreeObject_GM.get());

  std::cout<<"Analysis complete"<<std::endl;
}
