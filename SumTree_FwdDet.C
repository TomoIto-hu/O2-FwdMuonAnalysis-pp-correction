#include<iostream>
#include<memory>

#include "TFile.h"
#include "TTree.h"
#include <TGeoGlobalMagField.h>
#include <TGraph.h>
#include <TGraph2DPainter.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TEfficiency.h>
#include <THStack.h>

struct Propagated_parametor{
  double propagated_X,propagated_Y,propagated_DCA,propagated_Phi;
} typedef Propagated_parametor;

Propagated_parametor getToXYplane_linear(double x, double y, double z, double theta, double phi){
  Propagated_parametor Prop_p_linear;
  auto dZ = (0.-z);
  //std::cout<<"dZ = "<<dZ<<std::endl;
  Prop_p_linear.propagated_X = x + dZ*tan(theta)*cos(phi); //prospagated x coordinate
  Prop_p_linear.propagated_Y = y + dZ*tan(theta)*sin(phi); //prospagated y coordinate
  Prop_p_linear.propagated_DCA = pow(pow(Prop_p_linear.propagated_X,2)+pow(Prop_p_linear.propagated_Y,2),0.5); //DCA
  Prop_p_linear.propagated_Phi = phi;
  //printf("DCA = %f\n", DCAp);
  return Prop_p_linear;
}

Propagated_parametor getToXYplane_zhelix(double x, double y, double z, double theta, double phi, double qptOverk, double hz){
  Propagated_parametor Prop_p_zhelix;
  auto dZ = (0.-z);
  auto theta_helix = -dZ*tan(theta)/qptOverk;
  Prop_p_zhelix.propagated_X = x + hz*(1-cos(theta_helix))*qptOverk*sin(phi) - qptOverk*sin(theta_helix)*cos(phi);
  Prop_p_zhelix.propagated_Y = y - hz*(1-cos(theta_helix))*qptOverk*cos(phi) - qptOverk*sin(theta_helix)*sin(phi);
  Prop_p_zhelix.propagated_DCA = pow(pow(Prop_p_zhelix.propagated_X,2)+pow(Prop_p_zhelix.propagated_Y,2),0.5); //DCA
  Prop_p_zhelix.propagated_Phi = phi + hz*theta_helix;
  //printf("DCA = %f\n", DCAp);
  return Prop_p_zhelix;
}


void SumTree_FwdDet(){
  // Open root file and get Object
  // MC Information
  std::unique_ptr<TH1F> MCTrackPDG_dis_All(new TH1F("MCTrackPDG_dis_All","",10000,0-0.05,1000-0.05));
  std::unique_ptr<TH1F> MCTrackPDG_dis_Muon(new TH1F("MCTrackPDG_dis_Muon","",10000,0-0.05,1000-0.05));
  std::unique_ptr<TH1F> MCTrackDfIP_dis_Heavy(new TH1F("MCTrackDfIP_dis_Heavy","",10000,0-0.005,100-0.005));
  std::unique_ptr<TH1F> MCTrackPt_Muon_nonPiK(new TH1F("MCTrackPt_Muon_nonPiK","",1000,0-0.05,100-0.05));

  // GM Information All
  std::unique_ptr<TH1F> GMTrackP_All(new TH1F("GMTrackP_All","",10000,0-0.05,1000-0.05));
  std::unique_ptr<TH1F> GMTrackEta_All(new TH1F("GMTrackEta_All","",20000,-10-0.0005,10-0.0005));
  std::unique_ptr<TH1F> GMTrackX_All(new TH1F("GMTrackX_All","",20000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> GMTrackY_All(new TH1F("GMTrackY_All","",20000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> GMTrackZ_All(new TH1F("GMTrackZ_All","",20000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> GMTrackZ_Correct(new TH1F("GMTrackZ_Correct","",20000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> GMTrackZ_InCorrect(new TH1F("GMTrackZ_inCorrect","",20000,-1000-0.005,1000-0.005));
  std::unique_ptr<TH1F> GMTrackDfIP_All(new TH1F("GMTrackDfIP_All","",10000,0,100));
  //std::unique_ptr<TH1F> GMTrackDfIP_All(new TH1F("GMTrackDfIP_All","",20000,0,100));
  std::unique_ptr<TH1F> GMTrackDfIP_Muon_All(new TH1F("GMTrackDfIP_Muon_All","",10000,0,100));
  //std::unique_ptr<TH1F> GMTrackPt_DCAcut_021(new TH1F("GMTrackPt_DCAcut_021","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> GMTrackPDG_All(new TH1F("GMTrackPDG_All","",80000,0-0.05,8000-0.05));
  std::unique_ptr<TH1F> GMTrackMother_All(new TH1F("GMTrackMother_All","",80000,0-0.05,8000-0.05));
  std::unique_ptr<TH1F> GMTrackMother_Correct(new TH1F("GMTrackMother_Correct","",10000,0-0.05,1000-0.05));
  std::unique_ptr<TH1F> GMTrackDfIP_InCorrect(new TH1F("GMTrackDfIP_InCorrect","",10000,0-0.005,100-0.005));
  std::unique_ptr<TH1F> GMTrackDfIP_Correct(new TH1F("GMTrackDfIP_Correct","",10000,0,100));
  //std::unique_ptr<TH1F> GMTrackDfIP_Correct(new TH1F("GMTrackDfIP_Correct","",20000,0,100));
  std::unique_ptr<TH1F> GMTrackScore_All(new TH1F("GMTrackScore_All","",2000,0,1000));
  std::unique_ptr<TH1F> GMTrackScore_Correct(new TH1F("GMTrackScore_Correct","",2000,0,1000));
  std::unique_ptr<TH1F> GMTrackScore_InCorrect(new TH1F("GMTrackScore_InCorrect","",2000,0,1000));
  // GM Information PiK
  std::unique_ptr<TH1F> GMTrackDfIP_PiK(new TH1F("GMTrackDfIP_PiK","",10000,0,100));
  //std::unique_ptr<TH1F> GMTrackDfIP_PiK(new TH1F("GMTrackDfIP_PiK","",20000,0,100));
  std::unique_ptr<TH1F> GMTrackDfIP_Correct_PiK(new TH1F("GMTrackDfIP_Correct_PiK","",10000,0,100));
  std::unique_ptr<TH1F> GMTrackDfIP_InCorrect_PiK(new TH1F("GMTrackDfIP_InCorrect_PiK","",10000,0,100));
  //std::unique_ptr<TH1F> GMTrackDfIP_Correct_PiK(new TH1F("GMTrackDfIP_Correct_PiK","",20000,0,100));
  std::unique_ptr<TH1F> GMTrackMCDfIP_PiK(new TH1F("GMTrackMCDfIP_PiK","",10000,0,100));
  std::unique_ptr<TH1F> GMTrackMCDfIP_Correct_PiK(new TH1F("GMTrackMCDfIP_Correct_PiK","",10000,0,100));
  //GMTrack propagatopn linear
  std::unique_ptr<TH1F> GMTrackX_propagated_All(new TH1F("GMTrackX_propagated_All","",20000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> GMTrackY_propagated_All(new TH1F("GMTrackY_propagated_All","",20000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> GMTrackDCA_propagated_All(new TH1F("GMTrackDCA_propagated_All","",1000,0,100));
  //GM propagated helix
  std::unique_ptr<TH1F> GMTrackX_p_All(new TH1F("GMTrackX_p_All","",20000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> GMTrackY_p_All(new TH1F("GMTrackY_p_All","",20000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> GMTrackTheta_p_All(new TH1F("GMTrackTheta_p_All","",2000,-10-0.05,10-0.05));
  std::unique_ptr<TH1F> GMTrackPhi_p_All(new TH1F("GMTrackPhi_p_All","",2000,-20-0.0025,20-0.0025));
  std::unique_ptr<TH1F> GMTrackX_phelix_All(new TH1F("GMTrackX_phelix_All","",20000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> GMTrackY_phelix_All(new TH1F("GMTrackY_phelix_All","",20000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> GMTrackDCA_phelix_All(new TH1F("GMTrackDCA_phelix_All","",1000,0,100));

  // GM Information nonPiK
  std::unique_ptr<TH1F> GMTrackDfIP_nonPiK(new TH1F("GMTrackDfIP_nonPiK","",10000,0,100));
  std::unique_ptr<TH1F> GMTrackPt_nonPiK(new TH1F("GMTrackPt_nonPiK","",1000,0-0.05,100-0.05));
  //std::unique_ptr<TH1F> GMTrackDfIP_nonPiK(new TH1F("GMTrackDfIP_nonPiK","",20000,0,100));
  std::unique_ptr<TH1F> GMTrackDfIP_Correct_nonPiK(new TH1F("GMTrackDfIP_Correct_nonPiK","",10000,0,100));
  std::unique_ptr<TH1F> GMTrackDfIP_InCorrect_nonPiK(new TH1F("GMTrackDfIP_InCorrect_nonPiK","",10000,0,100));
  //std::unique_ptr<TH1F> GMTrackDfIP_Correct_nonPiK(new TH1F("GMTrackDfIP_Correct_nonPiK","",20000,0,100));
  std::unique_ptr<TH1F> GMTrackMCDfIP_nonPiK(new TH1F("GMTrackMCDfIP_nonPiK","",10000,0,100));
  std::unique_ptr<TH1F> GMTrackMCDfIP_Correct_nonPiK(new TH1F("GMTrackMCDfIP_Correct_nonPiK","",10000,0,100));

  // GM Information Vector
  std::unique_ptr<TH1F> GMTrackDfIP_Correct_Vector(new TH1F("GMTrackDfIP_Correct_Vector","",10000,0-0.005,100-0.005));
  // GM Information Heavy
  std::unique_ptr<TH1F> GMTrackDfIP_Correct_Heavy(new TH1F("GMTrackDfIP_Correct_Heavy","",10000,0-0.005,100-0.005));
  // GM Information Other
  std::unique_ptr<TH1F> GMTrackDfIP_Correct_Other(new TH1F("GMTrackDfIP_Correct_Other","",10000,0-0.005,100-0.005));

  //Spectol
  //All
  std::unique_ptr<TH1F> GMTrackPt_All(new TH1F("GMTrackPt_All","",1000,0-0.05,100-0.05));
  //Pion Keon
  std::unique_ptr<TH1F> GMTrackPt_PiK(new TH1F("GMTrackPt_PiK","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> GMTrackPt_Pi(new TH1F("GMTrackPt_Pi","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> GMTrackPt_K(new TH1F("GMTrackPt_K","",1000,0-0.05,100-0.05));
  //Heavy Flavor
  std::unique_ptr<TH1F> GMTrackPt_HF(new TH1F("GMTrackPt_HF","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> GMTrackPt_D(new TH1F("GMTrackPt_D","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> GMTrackPt_B(new TH1F("GMTrackPt_B","",1000,0-0.05,100-0.05));
  //Low mass Light Vecor meson
  std::unique_ptr<TH1F> GMTrackPt_LVM(new TH1F("GMTrackPt_LVM","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> GMTrackPt_Rou(new TH1F("GMTrackPt_rou","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> GMTrackPt_Omega(new TH1F("GMTrackPt_omega","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> GMTrackPt_Phi(new TH1F("GMTrackPt_phi","",1000,0-0.05,100-0.05));

  std::unique_ptr<TH1F> GMTrackGrandMother_Rou(new TH1F("GMTrackGrandMother_Rou","",100000,0,10000));
  //Quarkonia
  std::unique_ptr<TH1F> GMTrackPt_Quarkonia(new TH1F("GMTrackPt_Quarkonia","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> GMTrackPt_JPsi(new TH1F("GMTrackPt_JPsi","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> GMTrackPt_Upsilon(new TH1F("GMTrackPt_Upsilon","",1000,0-0.05,100-0.05));
  //Others(Muon)
  std::unique_ptr<TH1F> GMTrackPt_Direct(new TH1F("GMTrackPt_Direct","",1000,0-0.05,100-0.05));
  //Others(Muon)
  std::unique_ptr<TH1F> GMTrackPt_OtherMuon(new TH1F("GMTrackPt_OtherMuon","",1000,0-0.05,100-0.05));
  //Others(nonMuon)
  std::unique_ptr<TH1F> GMTrackPt_OtherNonMuon(new TH1F("GMTrackPt_OtherNonMuon","",1000,0-0.05,100-0.05));

  //MC information for GM tracks
  std::unique_ptr<TH1F> MCGMTrackZ_All(new TH1F("MCGMTrackZ_All","",1000,-900-0.5,100-0.5));
  std::unique_ptr<TH1F> MCGMTrackZ_PiK(new TH1F("MCGMTrackZ_PiK","",1000,-900-0.5,100-0.5));
  std::unique_ptr<TH1F> MCGMTrackZ_nonPiK(new TH1F("MCGMTrackZ_nonPiK","",1000,-900-0.5,100-0.5));

  // MC + GM
  std::unique_ptr<TH1F> MCGMmatchDfIP_All(new TH1F("MCGMmatchDfIP_All","",20000,-100-0.005,100-0.005));
  std::unique_ptr<TH1F> MCGMmatchDfIP_PiK(new TH1F("MCGMmatchDfIP_PiK","",20000,-100-0.005,100-0.005));
  std::unique_ptr<TH1F> MCGMmatchDfIP_InCorrect_PiK(new TH1F("MCGMmatchDfIP_InCorrect_PiK","",20000,-100-0.005,100-0.005));
  std::unique_ptr<TH1F> MCGMmatchDfIP_nonPiK(new TH1F("MCGMmatchDfIP_nonPiK","",20000,-100-0.005,100-0.005));
  std::unique_ptr<TH1F> MCGMmatchDfIP_InCorrect_nonPiK(new TH1F("MCGMmatchDfIP_InCorrect_nonPiK","",20000,-100-0.005,100-0.005));

  //MCH Track information
  std::unique_ptr<TH1F> MCHTrackP_All(new TH1F("MCHTrackP_All","",10000,0-0.05,1000-0.05));
  std::unique_ptr<TH1F> MCHTrackPt_All(new TH1F("MCHTrackPt_All","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> MCHTrackPx_All(new TH1F("MCHTrackPx_All","",2000,-100-0.05,100-0.05));
  std::unique_ptr<TH1F> MCHTrackPy_All(new TH1F("MCHTrackPy_All","",2000,-100-0.05,100-0.05));
  std::unique_ptr<TH1F> MCHTrackPz_All(new TH1F("MCHTrackPz_All","",2000,-100-0.05,100-0.05));
  std::unique_ptr<TH1F> MCHTrackEta_All(new TH1F("MCHTrackEta_All","",2000,-10-0.005,10-0.005));
  std::unique_ptr<TH1F> MCHTrackTheta_All(new TH1F("MCHTrackTheta_All","",2000,-10-0.005,10-0.005));
  std::unique_ptr<TH1F> MCHTrackPhi_All(new TH1F("MCHTrackPhi_All","",2000,-20-0.0025,20-0.0025));
  std::unique_ptr<TH1F> MCHTrackX_All(new TH1F("MCHTrackX_All","",20000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> MCHTrackY_All(new TH1F("MCHTrackY_All","",20000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> MCHTrackZ_All(new TH1F("MCHTrackZ_All","",20000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> MCHTrackDCA_All(new TH1F("MCHTrackDCA_All","",10000,0,100));
  //MCH propagated linear
  std::unique_ptr<TH1F> MCHTrackX_propagated_All(new TH1F("MCHTrackX_propagated_All","",2000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> MCHTrackY_propagated_All(new TH1F("MCHTrackY_propagated_All","",2000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> MCHTrackDCA_propagated_All(new TH1F("MCHTrackDCA_propagated_All","",100,0,100));
  //MCH propagated helix
  std::unique_ptr<TH1F> MCHTrackX_phelix_All(new TH1F("MCHTrackX_phelix_All","",2000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> MCHTrackY_phelix_All(new TH1F("MCHTrackY_phelix_All","",2000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> MCHTrackDCA_phelix_All(new TH1F("MCHTrackDCA_phelix_All","",100,0,100));
  //MCH Track PiK & nonPiK
  std::unique_ptr<TH1F> MCHTrackDCA_PiK(new TH1F("MCHTrackDCA_PiK","",1000,0,100));
  std::unique_ptr<TH1F> MCHTrackDCA_nonPiK(new TH1F("MCHTrackDCA_nonPiK","",1000,0,100));
  std::unique_ptr<TH1F> MCHTrackX_phelix_PiK(new TH1F("MCHTrackX_phelix_PiK","",2000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> MCHTrackY_phelix_PiK(new TH1F("MCHTrackY_phelix_PiK","",2000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> MCHTrackDCA_phelix_PiK(new TH1F("MCHTrackDCA_phelix_PiK","",100,0,100));
  std::unique_ptr<TH1F> MCHTrackX_phelix_nonPiK(new TH1F("MCHTrackX_phelix_nonPiK","",2000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> MCHTrackY_phelix_nonPiK(new TH1F("MCHTrackY_phelix_nonPiK","",2000,-1000-0.05,1000-0.05));
  std::unique_ptr<TH1F> MCHTrackDCA_phelix_nonPiK(new TH1F("MCHTrackDCA_phelix_nonPiK","",100,0,100));

  //Get TObject from existing "SumResult_SumTree_FwdDet.root"
  TFile* CorrectionFile = TFile::Open("SumResult_SumTree_FwdDet.root");
  TGraphErrors *Eff_pT_Correction = (TGraphErrors*)CorrectionFile->Get("Eff_nonPiK");
  //Prepare output objects
  TFile* output = new TFile("SumResult_SumTree_FwdDet.root","recreate");
  TChain *SumTotalTree_MC=new TChain("TotalTreeObject_All","SumTotalTree_MC_All");
  TChain *SumTotalTree_GM=new TChain("TotalTreeObject_GM","SumTotalTree_GM_All");
  //Fill SumObject
  SumTotalTree_MC -> Add("MCTrackCompInfo*.root");
  SumTotalTree_GM -> Add("AnalysisResults*.root");

  //Get Branch Information :GM
  double P_GM,Pt_GM,DfIP_GM,vz_GM,vy_GM,vx_GM,Eta_GM,Theta_GM,Phi_GM;
  double vz_p_GM,vy_p_GM,vx_p_GM,Eta_p_GM,Theta_p_GM,Phi_p_GM,invTanl_p_GM,Qpt_p_GM,Hz_p_GM,k_p_GM;
  double P_MCGM,Pt_MCGM,DfIP_MCGM,vz_MCGM,vy_MCGM,vx_MCGM,Eta_MCGM;
  double P_GrandMCGM,Pt_GrandMCGM,DfIP_GrandMCGM,vz_GrandMCGM,vy_GrandMCGM,vx_GrandMCGM,Eta_GrandMCGM;
  double MIDCheck,MatchingScore;
  int CorrectCheck,SourceCheck,MotherId_MCGM;
  int Mother_MCGM,PDG_MCGM,GrandMother_MCGM;
  SumTotalTree_GM ->SetBranchAddress("vx_p_GM",&vx_p_GM);
  SumTotalTree_GM ->SetBranchAddress("vy_p_GM",&vy_p_GM);
  SumTotalTree_GM ->SetBranchAddress("vz_p_GM",&vz_p_GM);
  SumTotalTree_GM ->SetBranchAddress("Eta_p_GM",&Eta_p_GM);
  SumTotalTree_GM ->SetBranchAddress("Theta_p_GM",&Theta_p_GM);
  SumTotalTree_GM ->SetBranchAddress("Phi_p_GM",&Phi_p_GM);
  SumTotalTree_GM ->SetBranchAddress("invTanl_p_GM",&invTanl_p_GM);
  SumTotalTree_GM ->SetBranchAddress("Hz_p_GM",&Hz_p_GM);
  SumTotalTree_GM ->SetBranchAddress("k_p_GM",&k_p_GM);
  SumTotalTree_GM ->SetBranchAddress("Qpt_p_GM",&Qpt_p_GM);
  SumTotalTree_GM ->SetBranchAddress("DfIP_GM",&DfIP_GM);
  SumTotalTree_GM ->SetBranchAddress("P_GM",&P_GM);
  SumTotalTree_GM ->SetBranchAddress("Pt_GM",&Pt_GM);
  SumTotalTree_GM ->SetBranchAddress("vx_GM",&vx_GM);
  SumTotalTree_GM ->SetBranchAddress("vy_GM",&vy_GM);
  SumTotalTree_GM ->SetBranchAddress("vz_GM",&vz_GM);
  SumTotalTree_GM ->SetBranchAddress("Eta_GM",&Eta_GM);
  SumTotalTree_GM ->SetBranchAddress("Theta_GM",&Theta_GM);
  SumTotalTree_GM ->SetBranchAddress("Phi_GM",&Phi_GM);
  SumTotalTree_GM ->SetBranchAddress("DfIP_MCGM",&DfIP_MCGM);
  SumTotalTree_GM ->SetBranchAddress("P_MCGM",&P_MCGM);
  SumTotalTree_GM ->SetBranchAddress("Pt_MCGM",&Pt_MCGM);
  SumTotalTree_GM ->SetBranchAddress("vx_MCGM",&vx_MCGM);
  SumTotalTree_GM ->SetBranchAddress("vy_MCGM",&vy_MCGM);
  SumTotalTree_GM ->SetBranchAddress("vz_MCGM",&vz_MCGM);
  SumTotalTree_GM ->SetBranchAddress("Eta_MCGM",&Eta_MCGM);
  SumTotalTree_GM ->SetBranchAddress("PDG_MCGM",&PDG_MCGM);
  SumTotalTree_GM ->SetBranchAddress("Mother_MCGM",&Mother_MCGM);
  SumTotalTree_GM ->SetBranchAddress("MotherId_MCGM",&MotherId_MCGM);
  SumTotalTree_GM ->SetBranchAddress("GrandMother_MCGM",&GrandMother_MCGM);
  SumTotalTree_GM ->SetBranchAddress("CorrectCheck",&CorrectCheck);
  SumTotalTree_GM ->SetBranchAddress("SourceCheck",&SourceCheck);
  SumTotalTree_GM ->SetBranchAddress("MIDCheck",&MIDCheck);
  SumTotalTree_GM ->SetBranchAddress("MatchingScore",&MatchingScore);

  //Get Branch Information :MFT
  double P_MFT,Pt_MFT,DfIP_MFT,vz_MFT,vy_MFT,vx_MFT,Eta_MFT;
  double P_MCMFT,Pt_MCMFT,DfIP_MCMFT,vz_MCMFT,vy_MCMFT,vx_MCMFT,Eta_MCMFT;
  double Theta_MFT,Theta_MCMFT,Phi_MFT,Phi_MCMFT;
  double P_GrandMCMFT,Pt_GrandMCMFT,DfIP_GrandMCMFT,vz_GrandMCMFT,vy_GrandMCMFT,vx_GrandMCMFT,Eta_GrandMCMFT;
  int MotherId_MCMFT;
  int Mother_MCMFT,PDG_MCMFT,GrandMother_MCMFT;
  SumTotalTree_GM ->SetBranchAddress("DfIP_MFT",&DfIP_MFT);
  SumTotalTree_GM ->SetBranchAddress("P_MFT",&P_MFT);
  SumTotalTree_GM ->SetBranchAddress("Pt_MFT",&Pt_MFT);
  SumTotalTree_GM ->SetBranchAddress("vx_MFT",&vx_MFT);
  SumTotalTree_GM ->SetBranchAddress("vy_MFT",&vy_MFT);
  SumTotalTree_GM ->SetBranchAddress("vz_MFT",&vz_MFT);
  SumTotalTree_GM ->SetBranchAddress("Eta_MFT",&Eta_MFT);
  SumTotalTree_GM ->SetBranchAddress("Theta_MFT",&Theta_MFT);
  SumTotalTree_GM ->SetBranchAddress("Phi_MFT",&Phi_MFT);
  SumTotalTree_GM ->SetBranchAddress("DfIP_MCMFT",&DfIP_MCMFT);
  SumTotalTree_GM ->SetBranchAddress("P_MCMFT",&P_MCMFT);
  SumTotalTree_GM ->SetBranchAddress("Pt_MCMFT",&Pt_MCMFT);
  SumTotalTree_GM ->SetBranchAddress("vx_MCMFT",&vx_MCMFT);
  SumTotalTree_GM ->SetBranchAddress("vy_MCMFT",&vy_MCMFT);
  SumTotalTree_GM ->SetBranchAddress("vz_MCMFT",&vz_MCMFT);
  SumTotalTree_GM ->SetBranchAddress("Eta_MCMFT",&Eta_MCMFT);
  SumTotalTree_GM ->SetBranchAddress("Theta_MCMFT",&Theta_MCMFT);
  SumTotalTree_GM ->SetBranchAddress("Phi_MCMFT",&Phi_MCMFT);
  SumTotalTree_GM ->SetBranchAddress("PDG_MCMFT",&PDG_MCMFT);
  SumTotalTree_GM ->SetBranchAddress("Mother_MCMFT",&Mother_MCMFT);
  SumTotalTree_GM ->SetBranchAddress("MotherId_MCMFT",&MotherId_MCMFT);
  SumTotalTree_GM ->SetBranchAddress("GrandMother_MCMFT",&GrandMother_MCMFT);

  //Get Branch Information :MCH
  double P_MCH,Pt_MCH,Eta_MCH,vx_MCH,vy_MCH,vz_MCH,DfIP_MCH,Px_MCH,Py_MCH,Pz_MCH,Theta_MCH,Phi_MCH;
  SumTotalTree_GM ->SetBranchAddress("DfIP_MCH",&DfIP_MCH);
  SumTotalTree_GM ->SetBranchAddress("P_MCH",&P_MCH);
  SumTotalTree_GM ->SetBranchAddress("Pt_MCH",&Pt_MCH);
  SumTotalTree_GM ->SetBranchAddress("vx_MCH",&vx_MCH);
  SumTotalTree_GM ->SetBranchAddress("vy_MCH",&vy_MCH);
  SumTotalTree_GM ->SetBranchAddress("vz_MCH",&vz_MCH);
  SumTotalTree_GM ->SetBranchAddress("Eta_MCH",&Eta_MCH);
  SumTotalTree_GM ->SetBranchAddress("Phi_MCH",&Phi_MCH);
  SumTotalTree_GM ->SetBranchAddress("Theta_MCH",&Theta_MCH);
  SumTotalTree_GM ->SetBranchAddress("Px_MCH",&Px_MCH);
  SumTotalTree_GM ->SetBranchAddress("Py_MCH",&Py_MCH);
  SumTotalTree_GM ->SetBranchAddress("Pz_MCH",&Pz_MCH);

  //Get Branch Information :MC
  float P_All,Pt_All,DfIP_All,vz_All,Eta_All;
  int PDG_All,Mother_All,MotherId_All;
  SumTotalTree_MC ->SetBranchAddress("DfIP_All",&DfIP_All);
  SumTotalTree_MC ->SetBranchAddress("PDG_All",&PDG_All);
  SumTotalTree_MC ->SetBranchAddress("P_All",&P_All);
  SumTotalTree_MC ->SetBranchAddress("Mother_All",&Mother_All);
  SumTotalTree_MC ->SetBranchAddress("MotherId_All",&MotherId_All);
  SumTotalTree_MC ->SetBranchAddress("Pt_All",&Pt_All);
  SumTotalTree_MC ->SetBranchAddress("vz_All",&vz_All);
  SumTotalTree_MC ->SetBranchAddress("Eta_All",&Eta_All);

  int nEntry_MC = SumTotalTree_MC->GetEntries();
  vector<vector<int>> MotherData(2, vector<int>(1));
  for(int iEntry_MC=0 ;iEntry_MC < nEntry_MC;++iEntry_MC){
  	SumTotalTree_MC->GetEntry(iEntry_MC);
	  MCTrackPDG_dis_All->Fill(PDG_All);
		MCTrackPDG_dis_Muon->Fill(Mother_All);
    for(int i = 0; i<MotherData.at(0).size();i++){
      if(fabs(Mother_All)==MotherData.at(0).at(i)){
        MotherData.at(1).at(i)+=1;
        //std::cout<<PDGdata.at(0).at(i)<<" : "<<PDGdata.at(1).at(i)<<std::endl;
        break;
      }
      if(i==MotherData.at(0).size()-1){
        MotherData.at(0).push_back(fabs(Mother_All));
        MotherData.at(1).push_back(1);
        //std::cout<<MotherData.at(0).at(i+1)<<" : "<<MotherData.at(1).at(i+1)<<std::endl;
        break;
      }
    }
		if(Mother_All == 130 || Mother_All == 211 || Mother_All == 321 ) continue;
		MCTrackPt_Muon_nonPiK->Fill(Pt_All);

  }
  for(int j = 0; j<MotherData.at(0).size();j++){
    std::cout<<MotherData.at(0).at(j)<<" : "<<MotherData.at(1).at(j)<<std::endl;
  }
  cout<<"Entry = "<<nEntry_MC<<endl;

  int nEntry_GM = SumTotalTree_GM->GetEntries();

  double dZ,qpt_Over_k,Hz,theta_helix,xp,yp,DCAp,theta_mch;

  for(int iEntry_GM=0 ;iEntry_GM < nEntry_GM;++iEntry_GM){
        SumTotalTree_GM->GetEntry(iEntry_GM);
        //MCH phi redefinition
        if(Px_MCH<0 && Py_MCH>0)Phi_MCH=Phi_MCH+TMath::Pi();
        if(Px_MCH<0 && Py_MCH<0)Phi_MCH=Phi_MCH-TMath::Pi();
        //MCH theta redefinition
        Theta_MCH = TMath::Pi()-Theta_MCH;
        auto Prop_p_linear_GM = getToXYplane_linear(vx_p_GM,vy_p_GM,vz_p_GM,Theta_p_GM,Phi_p_GM);
        auto Prop_p_zhelix_GM = getToXYplane_zhelix(vx_p_GM,vy_p_GM,vz_p_GM,Theta_p_GM,Phi_p_GM,qpt_Over_k,Hz_p_GM);
        auto Prop_p_linear_MCH = getToXYplane_linear(vx_MCH,vy_MCH,vz_MCH,Theta_MCH,Phi_MCH);
        auto Prop_p_zhelix_MCH = getToXYplane_zhelix(vx_MCH,vy_MCH,vz_MCH,Theta_MCH,Phi_MCH,qpt_Over_k,Hz_p_GM);

        if(MIDCheck==-1) continue; //MID requirement switch
	if(P_GM < 4) continue;
  //if(DfIP_MFT > 0.1) continue;
  //if(Prop_p_zhelix_MCH.propagated_DCA > 6) continue;
        //if(MatchingScore>5.5) continue;
        //if(DfIP_GM>0.030) continue;
        GMTrackX_All->Fill(vx_GM);
        GMTrackY_All->Fill(vy_GM);
        GMTrackP_All->Fill(P_GM);
        GMTrackPt_All->Fill(Pt_GM);//Spectol All
        GMTrackEta_All->Fill(Eta_GM);
        GMTrackDfIP_All->Fill(DfIP_GM);
        GMTrackScore_All->Fill(MatchingScore);
        GMTrackTheta_p_All->Fill(Theta_p_GM);
        GMTrackPhi_p_All->Fill(Phi_p_GM);
        GMTrackX_p_All->Fill(vx_p_GM);
        GMTrackY_p_All->Fill(vy_p_GM);

        //MCGM informatio
        MCGMTrackZ_All->Fill(vz_MCGM);

        //MCH Track distribution
        MCHTrackP_All->Fill(P_MCH);
        MCHTrackPt_All->Fill(Pt_MCH);//Spectol All
        MCHTrackPx_All->Fill(Px_MCH);
        MCHTrackPy_All->Fill(Py_MCH);
        MCHTrackPz_All->Fill(Pz_MCH);
        MCHTrackEta_All->Fill(Eta_MCH);
        MCHTrackTheta_All->Fill(Theta_MCH);
        MCHTrackPhi_All->Fill(Phi_MCH);
        MCHTrackX_All->Fill(vx_MCH);
        MCHTrackY_All->Fill(vy_MCH);
        MCHTrackZ_All->Fill(vz_MCH);
        MCHTrackDCA_All->Fill(DfIP_MCH);

        //GM propagation linear
        GMTrackX_propagated_All->Fill(Prop_p_linear_GM.propagated_X);
        GMTrackY_propagated_All->Fill(Prop_p_linear_GM.propagated_Y);
        GMTrackDCA_propagated_All->Fill(Prop_p_linear_GM.propagated_DCA);
        //GM propagation zhelix
        qpt_Over_k = Qpt_p_GM/k_p_GM;
        GMTrackX_phelix_All->Fill(Prop_p_zhelix_GM.propagated_X);
        GMTrackY_phelix_All->Fill(Prop_p_zhelix_GM.propagated_Y);
        GMTrackDCA_phelix_All->Fill(Prop_p_zhelix_GM.propagated_DCA);

        //if(Theta_MCH<3.07 || Theta_MCH>3.08) continue;
	      //Propagate MCH track to XY plane linear.
        MCHTrackX_propagated_All->Fill(Prop_p_linear_MCH.propagated_X);
        MCHTrackY_propagated_All->Fill(Prop_p_linear_MCH.propagated_Y);
        MCHTrackDCA_propagated_All->Fill(Prop_p_linear_MCH.propagated_DCA);

        //Propagate MCH track to XY plane with zfield (helix).
        MCHTrackX_phelix_All->Fill(Prop_p_zhelix_MCH.propagated_X);
        MCHTrackY_phelix_All->Fill(Prop_p_zhelix_MCH.propagated_Y);
        MCHTrackDCA_phelix_All->Fill(Prop_p_zhelix_MCH.propagated_DCA);

        if(fabs(PDG_MCGM)==13){
          if(fabs(Mother_MCGM) == 130 || fabs(Mother_MCGM) == 211 || fabs(Mother_MCGM) == 321 ){
            MCHTrackX_phelix_PiK->Fill(Prop_p_zhelix_MCH.propagated_X);
            MCHTrackY_phelix_PiK->Fill(Prop_p_zhelix_MCH.propagated_Y);
            MCHTrackDCA_phelix_PiK->Fill(Prop_p_zhelix_MCH.propagated_DCA);
          }else{
            MCHTrackX_phelix_nonPiK->Fill(Prop_p_zhelix_MCH.propagated_X);
            MCHTrackY_phelix_nonPiK->Fill(Prop_p_zhelix_MCH.propagated_Y);
            MCHTrackDCA_phelix_nonPiK->Fill(Prop_p_zhelix_MCH.propagated_DCA);
          }
        }else{//not muon track
          MCHTrackX_phelix_PiK->Fill(Prop_p_zhelix_MCH.propagated_X);
          MCHTrackY_phelix_PiK->Fill(Prop_p_zhelix_MCH.propagated_Y);
          MCHTrackDCA_phelix_PiK->Fill(Prop_p_zhelix_MCH.propagated_DCA);
        }

	//MC information check
	if(SourceCheck==0){//GMTrack Source==false
    std::cout<<"No source data"<<std::endl;
		GMTrackPt_OtherNonMuon->Fill(Pt_GM);//Spectol NonMuon
    		MCHTrackDCA_PiK->Fill(DfIP_MCH);
	}else{//GMTrack source==true
		GMTrackZ_All->Fill(vz_MCGM);
		GMTrackPDG_All->Fill(fabs(PDG_MCGM));
		if(CorrectCheck==1){//GMTrack All and Correct match
			MCGMmatchDfIP_All->Fill(DfIP_MCGM-DfIP_GM);
			GMTrackZ_Correct->Fill(vz_MCGM);
    			GMTrackScore_Correct->Fill(MatchingScore);
		}else{//GMTrack All and Fake match
	  		GMTrackZ_InCorrect->Fill(vz_MCGM);
        GMTrackScore_InCorrect->Fill(MatchingScore);
		}//END GMTrack All
		if(fabs(PDG_MCGM)==13){//All muon track
			GMTrackMother_All->Fill(fabs(Mother_MCGM));
			if(CorrectCheck==1) {//Correct(Muon) match
				GMTrackMother_Correct->Fill(fabs(Mother_MCGM));
				GMTrackDfIP_Correct->Fill(DfIP_GM);
			}
			if(CorrectCheck==0) GMTrackDfIP_InCorrect->Fill(DfIP_GM);
			GMTrackDfIP_Muon_All->Fill(DfIP_GM);
			if(fabs(Mother_MCGM) == 130 || fabs(Mother_MCGM) == 211 || fabs(Mother_MCGM) == 321 ){//Mother is PiK
				GMTrackPt_PiK->Fill(Pt_GM);//Spectol Pion Keon
				GMTrackDfIP_PiK->Fill(DfIP_GM);
				GMTrackMCDfIP_PiK->Fill(DfIP_MCGM);
        //MCGM PiK information
        MCGMTrackZ_PiK->Fill(vz_MCGM);
        //MCHTrack PiK information
        MCHTrackDCA_PiK->Fill(DfIP_MCH);
        if(fabs(Mother_MCGM) == 211){//Mother is Pion
          GMTrackPt_Pi->Fill(Pt_GM);//Spectol Pion
        }else{
          GMTrackPt_K->Fill(Pt_GM);//Spectol Keon
        }
				if(CorrectCheck==1){//Mother is PiK and Correct match
					GMTrackDfIP_Correct_PiK->Fill(DfIP_GM);
					GMTrackMCDfIP_Correct_PiK->Fill(DfIP_MCGM);
					MCGMmatchDfIP_PiK->Fill(DfIP_MCGM-DfIP_GM);
				}else{//Mother is PiK and fake match
					GMTrackDfIP_InCorrect_PiK->Fill(DfIP_GM);
					MCGMmatchDfIP_InCorrect_PiK->Fill(DfIP_MCGM-DfIP_GM);
				}
			}else{//Mother is nonPiK
				GMTrackPt_nonPiK->Fill(Pt_GM);
				GMTrackDfIP_nonPiK->Fill(DfIP_GM);
				GMTrackMCDfIP_nonPiK->Fill(DfIP_MCGM);
        //MCGM nonPiK information
        MCGMTrackZ_PiK->Fill(vz_MCGM);

        //MCHTrack nonPiK Information
        MCHTrackDCA_nonPiK->Fill(DfIP_MCH);
        if(fabs(Mother_MCGM) == 421 || fabs(Mother_MCGM) == 411 || fabs(Mother_MCGM) == 431 || fabs(Mother_MCGM) == 511|| fabs(Mother_MCGM) == 521 || fabs(Mother_MCGM) == 531){
              //GMTrackDfIP_Correct_Heavy->Fill(DfIP_GM);
              GMTrackPt_HF->Fill(Pt_GM);//Spectol Heavy Flaver
              if(fabs(Mother_MCGM) == 521 || fabs(Mother_MCGM) == 531|| fabs(Mother_MCGM) == 511){
                    GMTrackPt_B->Fill(Pt_GM);//Spectol B-meson
              }else{
                    GMTrackPt_D->Fill(Pt_GM);//Spectol D-meson
              }
        }else if(fabs(Mother_MCGM) == 113 || fabs(Mother_MCGM) == 223 || fabs(Mother_MCGM) == 333){
          //GMTrackDfIP_Correct_Vector->Fill(DfIP_GM);
              GMTrackPt_LVM->Fill(Pt_GM);//Spetol Low mass Light Vector meson
              if(fabs(Mother_MCGM) == 333)GMTrackPt_Phi->Fill(Pt_GM);//Spectol Phi
              if(fabs(Mother_MCGM) == 113){
			GMTrackPt_Rou->Fill(Pt_GM);//Spectol Rou
			GMTrackGrandMother_Rou->Fill(fabs(GrandMother_MCGM));
	      }
              if(fabs(Mother_MCGM) == 223)GMTrackPt_Omega->Fill(Pt_GM);//Spectol Omega
        }else	if(fabs(Mother_MCGM) == 443 || fabs(Mother_MCGM) == 553 || fabs(Mother_MCGM) == 100553){
              //GMTrackDfIP_Correct_Other->Fill(DfIP_GM);
              GMTrackPt_Quarkonia->Fill(Pt_GM);//Spectol Quarkonia
              if(fabs(Mother_MCGM) == 443)GMTrackPt_JPsi->Fill(Pt_GM);//Spectol JPsi
              if(fabs(Mother_MCGM) == 553 || fabs(Mother_MCGM) == 100553){
			GMTrackPt_Upsilon->Fill(Pt_GM);//Spectol upsilon
	      		std::cout<<"Other Muon Mother = "<<Mother_MCGM<<std::endl;
		}
        }else if(fabs(Mother_MCGM == 22)){
              GMTrackPt_Direct->Fill(Pt_GM);//Spectol Other Muon
        }else{
	      std::cout<<"Other Muon Mother = "<<Mother_MCGM<<std::endl;
              GMTrackPt_OtherMuon->Fill(Pt_GM);//Spectol Other Muon
        }

				if(CorrectCheck==1){//Mother is signal and correct match
					GMTrackDfIP_Correct_nonPiK->Fill(DfIP_GM);
					GMTrackMCDfIP_Correct_nonPiK->Fill(DfIP_MCGM);
				}else{//mother is signal and fake match
					GMTrackDfIP_InCorrect_nonPiK->Fill(DfIP_GM);
					MCGMmatchDfIP_InCorrect_nonPiK->Fill(DfIP_MCGM-DfIP_GM);
				}
			}
  		}else{//nonMuon GMTrack
			GMTrackPt_OtherNonMuon->Fill(Pt_GM);
			GMTrackDfIP_PiK->Fill(DfIP_GM);
			GMTrackMCDfIP_PiK->Fill(DfIP_MCGM);
		        //MCH Track information
      			MCHTrackDCA_PiK->Fill(DfIP_MCH);
			if(CorrectCheck == 1){//NonMuon and correct match
				GMTrackDfIP_Correct_PiK->Fill(DfIP_GM);
				GMTrackDfIP_Correct->Fill(DfIP_GM);
				MCGMmatchDfIP_PiK->Fill(DfIP_MCGM-DfIP_GM);
			}else{//nonMuon and fake match
				MCGMmatchDfIP_InCorrect_PiK->Fill(DfIP_MCGM-DfIP_GM);
				GMTrackDfIP_InCorrect_PiK->Fill(DfIP_GM);
			}
		}
  }
  }
  cout<<"Entry = "<<nEntry_GM<<endl;

  //puriry efficiency DCA cut discussion
  int n = 1000;
  double Count_True_Muon[n];
  double Count_CorrectPiK[n];
  double Count_PiK_All[n];
  double DfIP_axis[n];
  double Chi2_axis[n];
  double E_Chi2_axis[n];
  double E_DfIP[n];//Error parametor

  double Count_GMTrack_Correct_PiK[n],Count_GMTrack_Correct_nonPiK[n],Count_GMTrack_Correct_All[n],Count_GMTrack_All[n];
  double Count_GMTrack_All_PiK[n],Count_GMTrack_All_nonPiK[n];
  double Efficiency[n],Purity[n],SignalPurity[n],RealEfficiency[n],RemovalPiK[n],PurityPiK[n],SNRate[n],SNRateW[n],Collection[n];
  double Purity_All[n],SignalPurity_All[n],RealEfficiency_All[n],RemovalPiK_All[n],PurityPiK_All[n],SNRate_All[n],SNRateW_All[n],Collection_All[n];
  double E_Eff[n],E_Puri[n],E_Coll[n],E_Real[n],E_Sig[n],E_Rem[n],E_PiKP[n],E_SN[n],E_SNW[n];
  double E_Eff_All[n],E_Puri_All[n],E_Coll_All[n],E_Real_All[n],E_Sig_All[n],E_Rem_All[n],E_PiKP_All[n],E_SN_All[n],E_SNW_All[n];

  double Count_MCGMTrack_nonPiK[n],Count_MCGMTrack_PiK[n],Count_MCGMTrack_All[n],Count_MCGMTrack_PiK_All[n];
  double RemovalPiK_MCGM[n],SNRate_MCGM[n],SNRateW_MCGM[n];
  double E_Rem_MCGM[n],E_SN_MCGM[n],E_SNW_MCGM[n];

  //PiK decay point Check
  double z_axis[n],E_z_axis[n];
  double Count_Total_zMCGM_PiK;
  double Count_zMCGM_All[n],Count_zMCGM_PiK[n],Count_zMCGM_nonPiK[n];
  double RateDecay_PiK[n],E_RateDecay_PiK[n];
  double RateDecay_nonPiK[n],E_RateDecay_nonPiK[n];

  //Matching parametor
  double Count_Chi2_Correct[n],Count_Chi2_InCorrect[n];
  double GMChi2_SN[n],GMChi2_SNW[n];
  double E_GMChi2_SN[n],E_GMChi2_SNW[n];

  //MCH Track information
  double MCHDCA_axis[n],E_MCHDCA[n];
  double Count_MCH_DCA_PiK[n],Count_MCH_DCA_nonPiK[n];
  double MCH_DCA_SN[n],MCH_DCA_SNW[n];
  double E_MCH_DCA_SN[n],E_MCH_DCA_SNW[n];

  //Pt
  double Pt_axis[n],E_Pt[n];
  double Count_True_Muon_Pt[n],Count_Track_Pt[n],Count_Track_PiK_Pt[n],Count_Track_nonPiK_Pt[n],Count_Track_DCAcut[n];
  double Efficiency_nonPiK[n],Efficiency_DCAcut[n],Efficiency_nonPiK_Signal[n],Efficiency_nonPiK_Other[n];
  double E_Eff_nonPiK[n],E_Eff_DCAcut[n];

  double Efficiency_Even;

  //Pt_correction
  double pT_correction_Eff[n],pT_dis_corrected[n];
  double E_pT_dis_corrected[n];

  for(int i = 0; i<n; ++i){
	std::cout<<"roop = "<<i<<std::endl;
	//Correct
	Count_CorrectPiK[i] = GMTrackDfIP_Correct_PiK->GetEntries();
	Count_PiK_All[i] = GMTrackDfIP_PiK->GetEntries();
	Count_True_Muon[i] = MCTrackPt_Muon_nonPiK->GetEntries();
  	//std::cout<<"Count True Muon = "<<Count_True_Muon[i]<<std::endl;
  	//std::cout<<"Count CorrectPiK = "<<Count_CorrectPiK[i]<<std::endl;
  	//std::cout<<"Count PiK All = "<<Count_PiK_All[i]<<std::endl;
	//DfIP
	if(i==0){//初期値を入れる
  		Count_GMTrack_Correct_PiK[i] =  GMTrackDfIP_Correct_PiK->GetBinContent(i+1);
  		Count_GMTrack_Correct_nonPiK[i] = GMTrackDfIP_Correct_nonPiK->GetBinContent(i+1);
		  Count_GMTrack_Correct_All[i] = Count_GMTrack_Correct_PiK[i] + Count_GMTrack_Correct_nonPiK[i];
		  Count_GMTrack_All[i] = GMTrackDfIP_All->GetBinContent(i+1);
      //PiK decay point Check
      Count_zMCGM_PiK[i] = MCGMTrackZ_PiK->GetBinContent(1001-i);
      Count_zMCGM_nonPiK[i] = MCGMTrackZ_nonPiK->GetBinContent(1001-i);
      Count_zMCGM_All[i] = Count_zMCGM_PiK[n]+Count_zMCGM_nonPiK[i];

      //Matching Chi2
      Count_Chi2_Correct[i] = GMTrackScore_Correct->GetBinContent(i+1);
      Count_Chi2_InCorrect[i] = GMTrackScore_InCorrect->GetBinContent(i+1);
      //MCH DCA
      Count_MCH_DCA_PiK[i] = MCHTrackDCA_phelix_PiK->GetBinContent(i+1);
      Count_MCH_DCA_nonPiK[i] = MCHTrackDCA_phelix_nonPiK->GetBinContent(i+1);
 	}
  	Count_GMTrack_Correct_PiK[i] = Count_GMTrack_Correct_PiK[i-1] + GMTrackDfIP_Correct_PiK->GetBinContent(i+1);
  	Count_GMTrack_Correct_nonPiK[i] = Count_GMTrack_Correct_nonPiK[i-1] + GMTrackDfIP_Correct_nonPiK->GetBinContent(i+1);
	Count_GMTrack_Correct_All[i] = Count_GMTrack_Correct_PiK[i-1] + Count_GMTrack_Correct_nonPiK[i-1];
	Count_GMTrack_All[i] = Count_GMTrack_All[i-1] + GMTrackDfIP_All->GetBinContent(i+1);

  //PiK decay point Check
  Count_zMCGM_PiK[i] = Count_zMCGM_PiK[i-1] + MCGMTrackZ_PiK->GetBinContent(1001-i);
  Count_zMCGM_nonPiK[i] = Count_zMCGM_nonPiK[i-1] + MCGMTrackZ_nonPiK->GetBinContent(1001-i);
  Count_zMCGM_All[i] = Count_zMCGM_PiK[n]+Count_zMCGM_nonPiK[i];
  Count_Total_zMCGM_PiK = MCGMTrackZ_PiK->GetEntries();

  	//Matching Chi2
  	Count_Chi2_Correct[i] = Count_Chi2_Correct[i-1] + GMTrackScore_Correct->GetBinContent(i+1);
  	Count_Chi2_InCorrect[i] = Count_Chi2_InCorrect[i-1] + GMTrackScore_InCorrect->GetBinContent(i+1);

	//MCH DCA
  	Count_MCH_DCA_PiK[i] = Count_MCH_DCA_PiK[i-1] + MCHTrackDCA_phelix_PiK->GetBinContent(i+1);
  	Count_MCH_DCA_nonPiK[i] = Count_MCH_DCA_nonPiK[i-1] + MCHTrackDCA_phelix_nonPiK->GetBinContent(i+1);

	//MCGM
	Count_MCGMTrack_PiK_All[i] = GMTrackMCDfIP_PiK->GetEntries();
	if(i==0){
		Count_MCGMTrack_nonPiK[i] = GMTrackMCDfIP_nonPiK->GetBinContent(i+1);
		Count_MCGMTrack_PiK[i] = GMTrackMCDfIP_PiK->GetBinContent(i+1);
	}
	Count_MCGMTrack_nonPiK[i] = Count_MCGMTrack_nonPiK[i-1] + GMTrackMCDfIP_nonPiK->GetBinContent(i+1);
	Count_MCGMTrack_PiK[i] = Count_MCGMTrack_PiK[i-1] + GMTrackMCDfIP_PiK->GetBinContent(i+1);
	Count_MCGMTrack_All[i] = Count_MCGMTrack_PiK[i] + Count_MCGMTrack_nonPiK[i];
	//All
  	Count_GMTrack_All_PiK[i] = Count_GMTrack_All_PiK[i-1] + GMTrackDfIP_PiK->GetBinContent(i+1);
  	Count_GMTrack_All_nonPiK[i] = Count_GMTrack_All_nonPiK[i-1] + GMTrackDfIP_nonPiK->GetBinContent(i+1);

	//Pt
	Count_True_Muon_Pt[i] = MCTrackPt_Muon_nonPiK->GetBinContent(i+1);
	Count_Track_nonPiK_Pt[i] = GMTrackPt_nonPiK->GetBinContent(i+1);
	Count_Track_PiK_Pt[i] = GMTrackPt_PiK->GetBinContent(i+1);
	Count_Track_Pt[i] = Count_Track_nonPiK_Pt[i]+Count_Track_PiK_Pt[i];
	//Count_Track_DCAcut[i] = GMTrackPt_DCAcut_021->GetBinContent(i+1);

	//DfIP
	//DfIP_axis[i] = (i+1)*0.005;
	//E_DfIP[i] = 0.0025;
	DfIP_axis[i] = (i+1)*0.01 - 0.005;
	E_DfIP[i] = 0.005;
  //MCGM z axis
  z_axis[i] = 100 - (i)*1;
  E_z_axis[i] = 0.5;
  //MCH DCA
  MCHDCA_axis[i] = (i+1)*1 - 0.5;
  E_MCHDCA[i] = 0.5;
  	//Chi2
  	Chi2_axis[i] = (i+1)*0.5 - 0.25;
  	E_Chi2_axis[i] = 0.25;
	//Pt
	Pt_axis[i] = i*0.1;
	E_Pt[i] = 0.05;

	if(i<300){
	   //DfIP
	   //std::cout<<"Count_GMTrack_Correct_PiK = "<<Count_GMTrack_Correct_PiK[i]<<std::endl;
	   //std::cout<<"Count_GMTrack_Correct_nonPiK = "<<Count_GMTrack_Correct_nonPiK[i]<<std::endl;
	   //std::cout<<"Count_GMTrack_Correct_All = "<<Count_GMTrack_Correct_All[i]<<std::endl;
     //std::cout<<"Count_GMTrack_All = "<<Count_GMTrack_All[i]<<std::endl;
	   //std::cout<<"Count_GMTrack_All_PiK = "<<Count_GMTrack_All_PiK[i]<<std::endl;
     //std::cout<<"Count_GMTrack_All_nonPiK = "<<Count_GMTrack_All_nonPiK[i]<<std::endl;

     //PiK decay point Check
     std::cout<<"z = "<<z_axis[i]<<std::endl;
     std::cout<<"Count_zMCGM_PiK = "<<Count_zMCGM_PiK[i]<<std::endl;
     std::cout<<"Count_Total_zMCGM_PiK = "<<Count_Total_zMCGM_PiK<<std::endl;

	   //Pt
     //std::cout<<"Pt axis = "<<Pt_axis[i]<<std::endl;
     //std::cout<<"Count_True_Muon_Pt = "<<Count_True_Muon_Pt[i]<<std::endl;
     //std::cout<<"Count_Track_Pt = "<<Count_Track_Pt[i]<<std::endl;
     //std::cout<<"Count_Track_PiK_Pt = "<<Count_Track_PiK_Pt[i]<<std::endl;
     //std::cout<<"Count_Track_nonPiK_Pt = "<<Count_Track_nonPiK_Pt[i]<<std::endl;
     //std::cout<<"Count_Track_DCAcut = "<<Count_Track_DCAcut[i]<<std::endl;
	}
	//Out put
	Efficiency[i] = Count_GMTrack_All[i]/Count_True_Muon[i];;
	Purity_All[i] = Count_GMTrack_All_nonPiK[i]/Count_GMTrack_All[i];
	if(Count_GMTrack_Correct_All[i] >= 1) Purity[i] = Count_GMTrack_Correct_nonPiK[i]/Count_GMTrack_Correct_All[i];

	//Correct
	SignalPurity[i] = Count_GMTrack_Correct_nonPiK[i]/Count_GMTrack_All[i];//Correct-Signal/All
        RealEfficiency[i] = Efficiency[i]*SignalPurity[i];//Correct-Signal/Generated signal muon
	Collection[i] = Count_GMTrack_Correct_nonPiK[i]/40000.;
	RemovalPiK[i] = (Count_CorrectPiK[i] - Count_GMTrack_Correct_PiK[i])/Count_CorrectPiK[i]; //Rate of PiK removal
	PurityPiK[i] = Count_GMTrack_Correct_PiK[i]/Count_CorrectPiK[i]; //Rate of PiK
	if(Count_GMTrack_Correct_PiK[i] >= 1){
		SNRate[i] = Count_GMTrack_Correct_nonPiK[i]/Count_GMTrack_Correct_PiK[i];
		SNRateW[i] = pow(Count_GMTrack_Correct_nonPiK[i],1.5)/Count_GMTrack_Correct_PiK[i];
	}

	//All
	SignalPurity_All[i] = Count_GMTrack_All_nonPiK[i]/Count_GMTrack_All[i];//Correct-Signal/All
        RealEfficiency_All[i] = Efficiency[i]*SignalPurity_All[i];//Correct-Signal/Generated signal muon
	RemovalPiK_All[i] = (Count_PiK_All[i]-Count_GMTrack_All_PiK[i])/Count_PiK_All[i]; //Rate of PiK removal
	Collection_All[i] = Count_GMTrack_All_nonPiK[i]/40000.;
	PurityPiK_All[i] = Count_GMTrack_All_PiK[i]/Count_PiK_All[i]; //Rate of PiK
	SNRate_All[i] = Count_GMTrack_All_nonPiK[i]/Count_GMTrack_All_PiK[i];
	SNRateW_All[i] = pow(Count_GMTrack_All_nonPiK[i],1.5)/Count_GMTrack_All_PiK[i];

  //PiK decay point check
  RateDecay_PiK[i] = Count_zMCGM_PiK[i]/Count_Total_zMCGM_PiK;

	//MCGM
	if(Count_MCGMTrack_PiK[i]>=1){
		RemovalPiK_MCGM[i] = 1. - (Count_MCGMTrack_PiK[i]/Count_MCGMTrack_PiK_All[i]);
		SNRate_MCGM[i] = Count_MCGMTrack_nonPiK[i]/Count_MCGMTrack_PiK[i];
		SNRateW_MCGM[i] = pow(Count_MCGMTrack_nonPiK[i],1.5)/Count_MCGMTrack_PiK[i];
	}

  //Matching Chi2
  if(Count_Chi2_InCorrect[i]>=1){
    GMChi2_SN[i] = Count_Chi2_Correct[i]/Count_Chi2_InCorrect[i];
    GMChi2_SNW[i] = pow(Count_Chi2_Correct[i],1.5)/Count_Chi2_InCorrect[i];
    if(i<300){
      //std::cout<<"Chi2 cut = "<<Chi2_axis[i]<<std::endl;
      //std::cout<<"Chi2 SN = "<<GMChi2_SN[i]<<std::endl;
      //std::cout<<"Chi2 SNW = "<<GMChi2_SNW[i]<<std::endl;
      //std::cout<<"Chi2 Correct = "<<Count_Chi2_Correct[i]<<std::endl;
      //std::cout<<"Chi2 InCorrect = "<<Count_Chi2_InCorrect[i]<<std::endl;
    }
  }

  //MCH Track DCA
  if(Count_MCH_DCA_PiK[i]>=1){
    MCH_DCA_SN[i] = Count_MCH_DCA_nonPiK[i]/Count_MCH_DCA_PiK[i];
    MCH_DCA_SNW[i] = pow(Count_MCH_DCA_nonPiK[i],1.5)/Count_MCH_DCA_PiK[i];
    if(i<300){
      //std::cout<<"DCA cut = "<<DfIP_axis[i]<<std::endl;
      //std::cout<<"DCA(MCH) SN = "<<MCH_DCA_SN[i]<<std::endl;
      //std::cout<<"DCA(MCH) SNW = "<<MCH_DCA_SNW[i]<<std::endl;
      //std::cout<<"DCA(MCH) nonPiK = "<<Count_MCH_DCA_nonPiK[i]<<std::endl;
      //std::cout<<"DCA(MCH) PiK = "<<Count_MCH_DCA_PiK[i]<<std::endl;
    }
  }
	//Pt
	if(Count_True_Muon_Pt[i] >= 1){
		Efficiency_nonPiK[i] = (Count_Track_Pt[i]/Count_True_Muon_Pt[i]);
		Efficiency_nonPiK_Signal[i] = Count_Track_nonPiK_Pt[i]/Count_True_Muon_Pt[i];
		Efficiency_nonPiK_Other[i] = Count_Track_PiK_Pt[i]/Count_True_Muon_Pt[i];
		//Efficiency_DCAcut[i] = Count_Track_DCAcut[i]/Count_True_Muon_Pt[i];
	}
	Efficiency_Even += Efficiency_nonPiK[i]*(Count_True_Muon_Pt[i]/Count_True_Muon[i]);

  //Eff correction
  pT_correction_Eff[i] = Eff_pT_Correction->GetPointY(i);
  //pT Correction
  if(pT_correction_Eff[i] == 0){
    pT_dis_corrected[i] = 0;
  }else{
    pT_dis_corrected[i] = Count_Track_Pt[i]/(pT_correction_Eff[i]/100);
  }

  if(i<300){
	   //std::cout<<"Efficiency_Even = "<<Efficiency_Even<<std::endl;

	   //std::cout<<"Efficiency_nonPiK = "<<Efficiency_nonPiK[i]<<std::endl;
	   //std::cout<<"Efficiency_nonPiK_Signal = "<<Efficiency_nonPiK_Signal[i]<<std::endl;
	   //std::cout<<"Efficiency_nonPiK_Other = "<<Efficiency_nonPiK_Other[i]<<std::endl;
	   //Analysis_All = PurityPiK_All/RealEfficiency_All;
	   //std::cout<<"DCA cut = "<<DfIP_axis[i]<<std::endl;
	   //std::cout<<"Efficieny = "<<Efficiency[i]<<std::endl;
	   //std::cout<<"Purity = "<<Purity[i]<<std::endl;
	   //std::cout<<"Sig = "<<SignalPurity[i]<<std::endl;
	   //std::cout<<"Real Eff = "<<RealEfficiency[i]<<std::endl;
	   //std::cout<<"PiK Purity = "<<PurityPiK[i]<<std::endl;
	   //std::cout<<"SNRate = "<<SNRate[i]<<std::endl;
	   //std::cout<<"SNRateW = "<<SNRateW[i]<<std::endl;
	   //std::cout<<"Purity_All = "<<Purity_All[i]<<std::endl;
	   //std::cout<<"Sig_All = "<<SignalPurity_All[i]<<std::endl;
	   //std::cout<<"Real Eff_All = "<<RealEfficiency_All[i]<<std::endl;
	   //std::cout<<"PiK Purity_All = "<<PurityPiK_All[i]<<std::endl;
     //std::cout<<"DCA cut  = "<<DfIP_axis[i]<<std::endl;
	   //std::cout<<"SNRate_All = "<<SNRate_All[i]<<std::endl;
     //std::cout<<"SNRateW_All = "<<SNRateW_All[i]<<std::endl;
     std::cout<<"PiK decay rate = "<<RateDecay_PiK[i]<<std::endl;
     //std::cout<<"pT_correction_Eff/100 = "<<pT_correction_Eff[i]<<std::endl;
     //std::cout<<"pT_dis_corrected = "<<pT_dis_corrected[i]<<std::endl;
     std::cout<<""<<std::endl;

	}
	//Error
	E_Eff[i] = pow(Count_GMTrack_All[i],0.5)/Count_True_Muon[i];
 	if(Count_GMTrack_Correct_PiK[i] >= 1){
		if(Count_GMTrack_Correct_All[i] >= 1)E_Puri[i] = pow(Count_GMTrack_Correct_nonPiK[i],0.5)/Count_GMTrack_Correct_All[i];
		E_SN[i] = pow(SNRate[i]*(1.+SNRate[i])/Count_GMTrack_Correct_PiK[i],0.5);
		E_SNW[i] = pow(SNRate[i]*SNRate[i]*(1.+SNRate[i]),0.5);
	}
 	E_Coll[i] = pow(Count_GMTrack_Correct_nonPiK[i],0.5)/40000.;
 	E_Sig[i] = pow(Count_GMTrack_Correct_nonPiK[i],0.5)/Count_GMTrack_All[i];
 	E_Real[i] = pow(Count_GMTrack_Correct_nonPiK[i],0.5)/Count_GMTrack_All[i];
 	E_Rem[i] = pow(Count_CorrectPiK[i] - Count_GMTrack_Correct_PiK[i],0.5)/Count_CorrectPiK[i];
 	E_PiKP[i] = pow(Count_GMTrack_Correct_PiK[i],0.5)/Count_CorrectPiK[i];
	//All
 	if(Count_GMTrack_All_PiK[i] >= 1){
		if(Count_GMTrack_Correct_All[i] >= 1)E_Puri_All[i] = pow(Count_GMTrack_All_nonPiK[i],0.5)/Count_GMTrack_All[i];
		E_SN_All[i] = pow(SNRate_All[i]*(1.+SNRate_All[i])/Count_GMTrack_All_PiK[i],0.5);
		E_SNW_All[i] = pow(SNRate_All[i]*SNRate_All[i]*(1.+SNRate_All[i]),0.5);
	}
 	E_Coll_All[i] = pow(Count_GMTrack_All_nonPiK[i],0.5)/40000.;
 	E_Sig_All[i] = pow(Count_GMTrack_All_nonPiK[i],0.5)/Count_GMTrack_All[i];
 	E_Real_All[i] = pow(Count_GMTrack_All_nonPiK[i],0.5)/Count_GMTrack_All[i];
 	E_Rem_All[i] = pow(Count_PiK_All[i] - Count_GMTrack_All_PiK[i],0.5)/Count_PiK_All[i];
 	E_PiKP_All[i] = pow(Count_GMTrack_All_PiK[i],0.5)/Count_PiK_All[i];

  //PiK decay point Check
  E_RateDecay_PiK[i] = pow(RateDecay_PiK[i]*(1.-RateDecay_PiK[i]),0.5)/Count_Total_zMCGM_PiK;

	//MCGM
	if(Count_MCGMTrack_PiK[i]>=1){
		E_Rem_MCGM[i] = pow(Count_MCGMTrack_PiK_All[i] - Count_MCGMTrack_PiK[i],0.5)/Count_MCGMTrack_PiK_All[i];
		E_SN_MCGM[i] = pow(SNRate_MCGM[i]*(1.+SNRate_MCGM[i])/Count_MCGMTrack_PiK[i],0.5);
		E_SNW_MCGM[i] = pow(SNRate_MCGM[i]*SNRate_MCGM[i]*(1.+SNRate_MCGM[i]),0.5);
	}

  //Matching Chi2
  if(Count_Chi2_InCorrect[i]>=1){
		E_GMChi2_SN[i] = pow(GMChi2_SN[i]*(1.+GMChi2_SN[i])/Count_Chi2_InCorrect[i],0.5);
		E_GMChi2_SNW[i] = pow(GMChi2_SN[i]*GMChi2_SN[i]*(1.+GMChi2_SN[i]),0.5);
	}

  //MCH DCA
  if(Count_MCH_DCA_PiK[i]>=1){
		E_MCH_DCA_SN[i] = pow(MCH_DCA_SN[i]*(1.+MCH_DCA_SN[i])/Count_MCH_DCA_PiK[i],0.5);
		E_MCH_DCA_SNW[i] = pow(MCH_DCA_SN[i]*MCH_DCA_SN[i]*(1.+MCH_DCA_SN[i]),0.5);
	}

	//Pt
	if(Count_True_Muon_Pt[i] >=1){
		E_Eff_nonPiK[i] = pow(Efficiency_nonPiK[i]*(1.+ Efficiency_nonPiK[i])/Count_True_Muon_Pt[i],0.5);
		//E_Eff_DCAcut[i] = pow(Count_Track_DCAcut[i],0.5)/Count_True_Muon_Pt[i];
	}

  //Eff Correction
  E_pT_dis_corrected[i] = 0;

	//transfer eff->%
	if(Count_True_Muon_Pt[i] >=1){
		Efficiency_nonPiK[i] = Efficiency_nonPiK[i]*100;
		E_Eff_nonPiK[i] = E_Eff_nonPiK[i]*100;
	}
  }

  TGraphErrors* DCA_Eff = new TGraphErrors(n,DfIP_axis,Efficiency,E_DfIP,E_Eff);
  DCA_Eff -> SetName("DCA_Eff");
  TGraphErrors* DCA_Puri = new TGraphErrors(n,DfIP_axis,Purity,E_DfIP,E_Puri);
  DCA_Puri -> SetName("DCA_Puri");
  TGraphErrors* DCA_Coll = new TGraphErrors(n,DfIP_axis,Collection,E_DfIP,E_Coll);
  DCA_Coll -> SetName("DCA_Coll");
  TGraphErrors* DCA_Sig = new TGraphErrors(n,DfIP_axis,SignalPurity,E_DfIP,E_Sig);
  DCA_Sig -> SetName("DCA_Sig");
  TGraphErrors* DCA_Real = new TGraphErrors(n,DfIP_axis,RealEfficiency,E_DfIP,E_Real);
  DCA_Real -> SetName("DCA_Real");
  TGraphErrors* DCA_Rem = new TGraphErrors(n,DfIP_axis,RemovalPiK,E_DfIP,E_Rem);
  DCA_Rem -> SetName("DCA_Rem");
  TGraphErrors* DCA_PiKP = new TGraphErrors(n,DfIP_axis,PurityPiK,E_DfIP,E_PiKP);
  DCA_PiKP -> SetName("DCA_PiKP");
  TGraphErrors* DCA_SN = new TGraphErrors(n,DfIP_axis,SNRate,E_DfIP,E_SN);
  DCA_SN -> SetName("DCA_SN");
  TGraphErrors* DCA_SNW = new TGraphErrors(n,DfIP_axis,SNRateW,E_DfIP,E_SNW);
  DCA_SNW -> SetName("DCA_SNW");
  //All
  TGraphErrors* DCA_Puri_All = new TGraphErrors(n,DfIP_axis,Purity_All,E_DfIP,E_Puri_All);
  DCA_Puri_All -> SetName("DCA_Puri_All");
  TGraphErrors* DCA_Coll_All = new TGraphErrors(n,DfIP_axis,Collection_All,E_DfIP,E_Coll_All);
  DCA_Coll_All -> SetName("DCA_Coll_All");
  TGraphErrors* DCA_Sig_All = new TGraphErrors(n,DfIP_axis,SignalPurity_All,E_DfIP,E_Sig_All);
  DCA_Sig_All -> SetName("DCA_Sig_All");
  TGraphErrors* DCA_Real_All = new TGraphErrors(n,DfIP_axis,RealEfficiency_All,E_DfIP,E_Real_All);
  DCA_Real_All -> SetName("DCA_Real_All");
  TGraphErrors* DCA_Rem_All = new TGraphErrors(n,DfIP_axis,RemovalPiK_All,E_DfIP,E_Rem_All);
  DCA_Rem_All -> SetName("DCA_Rem_All");
  TGraphErrors* DCA_PiKP_All = new TGraphErrors(n,DfIP_axis,PurityPiK_All,E_DfIP,E_PiKP_All);
  DCA_PiKP_All -> SetName("DCA_PiKP_All");
  TGraphErrors* DCA_SN_All = new TGraphErrors(n,DfIP_axis,SNRate_All,E_DfIP,E_SN_All);
  DCA_SN_All -> SetName("DCA_SN_All");
  TGraphErrors* DCA_SNW_All = new TGraphErrors(n,DfIP_axis,SNRateW_All,E_DfIP,E_SNW_All);
  DCA_SNW_All -> SetName("DCA_SNW_All");
  //PiK decay point Check
  TGraphErrors* PiK_Decay_Rate_z = new TGraphErrors(n,z_axis,RateDecay_PiK,E_z_axis,E_RateDecay_PiK);
  PiK_Decay_Rate_z->SetName("PiK_Decay_Rate_z");
  //MCGM
  TGraphErrors* DCA_Rem_MCGM = new TGraphErrors(n,DfIP_axis,RemovalPiK_MCGM,E_DfIP,E_Rem_MCGM);
  DCA_Rem_MCGM -> SetName("DCA_Rem_MCGM");
  TGraphErrors* DCA_SN_MCGM = new TGraphErrors(n,DfIP_axis,SNRate_MCGM,E_DfIP,E_SN_MCGM);
  DCA_SN_MCGM -> SetName("DCA_SN_MCGM");
  TGraphErrors* DCA_SNW_MCGM = new TGraphErrors(n,DfIP_axis,SNRateW_MCGM,E_DfIP,E_SNW_MCGM);
  DCA_SNW_MCGM -> SetName("DCA_SNW_MCGM");

  //GM Matching Chi2
  TGraphErrors* GM_Chi2_SN = new TGraphErrors(n,Chi2_axis,GMChi2_SN,E_Chi2_axis,E_GMChi2_SN);
  GM_Chi2_SN -> SetName("GM_Chi2_SN");
  TGraphErrors* GM_Chi2_SNW = new TGraphErrors(n,Chi2_axis,GMChi2_SNW,E_Chi2_axis,E_GMChi2_SNW);
  GM_Chi2_SNW -> SetName("GM_Chi2_SNW");

  //MCH DCA
  TGraphErrors* DCA_SN_MCH = new TGraphErrors(n,MCHDCA_axis,MCH_DCA_SN,E_MCHDCA,E_MCH_DCA_SN);
  DCA_SN_MCH -> SetName("DCA_SN_MCH");
  TGraphErrors* DCA_SNW_MCH = new TGraphErrors(n,MCHDCA_axis,MCH_DCA_SNW,E_MCHDCA,E_MCH_DCA_SNW);
  DCA_SNW_MCH -> SetName("DCA_SNW_MCH");

  //Pt
  TGraphErrors* Eff_nonPiK = new TGraphErrors(n,Pt_axis,Efficiency_nonPiK,E_Pt,E_Eff_nonPiK);
  Eff_nonPiK -> SetName("Eff_nonPiK");
  //TGraphErrors* Eff_DCAcut = new TGraphErrors(n,Pt_axis,Efficiency_DCAcut,E_Pt,E_Eff_DCAcut);
  //Eff_DCAcut -> SetName("Eff_DCAcut");
  TGraphErrors* pT_corrected_All = new TGraphErrors(n,Pt_axis,pT_dis_corrected,E_Pt,E_pT_dis_corrected);
  pT_corrected_All -> SetName("pT_corrected_All");

  //Axis setting
  DCA_Eff->GetXaxis()->SetLimits(0,10);
  DCA_Puri->GetXaxis()->SetLimits(0,10);
  //Puri_Eff->GetXaxis()->SetLimits(0,1);
  DCA_Coll->GetXaxis()->SetLimits(0,10);
  DCA_Sig->GetXaxis()->SetLimits(0,10);
  DCA_Real->GetXaxis()->SetLimits(0,10);
  DCA_Rem->GetXaxis()->SetLimits(0,10);

  //DCA_Ana->SetMinimum(0.);
  DCA_SN_All->SetMinimum(0.);
  DCA_SNW_All->SetMinimum(0.);

  //Save histograms
  // MC Information
  output->WriteTObject(MCTrackPDG_dis_All.get());
  output->WriteTObject(MCTrackPDG_dis_Muon.get());
  output->WriteTObject(MCTrackPt_Muon_nonPiK.get());
  //GM Information All
  output->WriteTObject(GMTrackP_All.get());
  output->WriteTObject(GMTrackPt_All.get());
  //output->WriteTObject(GMTrackPt_DCAcut_021.get());
  output->WriteTObject(GMTrackEta_All.get());
  output->WriteTObject(GMTrackDfIP_All.get());
  output->WriteTObject(GMTrackDfIP_Muon_All.get());
  output->WriteTObject(GMTrackX_All.get());
  output->WriteTObject(GMTrackY_All.get());
  output->WriteTObject(GMTrackZ_All.get());
  output->WriteTObject(GMTrackZ_Correct.get());
  output->WriteTObject(GMTrackZ_InCorrect.get());
  output->WriteTObject(GMTrackPDG_All.get());
  output->WriteTObject(GMTrackMother_All.get());
  output->WriteTObject(GMTrackScore_All.get());
  output->WriteTObject(GMTrackScore_Correct.get());
  output->WriteTObject(GMTrackScore_InCorrect.get());
  //GM propagation linear
  output->WriteTObject(GMTrackX_propagated_All.get());
  output->WriteTObject(GMTrackY_propagated_All.get());
  output->WriteTObject(GMTrackDCA_propagated_All.get());
  //MCH propagated zhelix
  output->WriteTObject(GMTrackX_p_All.get());
  output->WriteTObject(GMTrackY_p_All.get());
  output->WriteTObject(GMTrackTheta_p_All.get());
  output->WriteTObject(GMTrackPhi_p_All.get());
  output->WriteTObject(GMTrackX_phelix_All.get());
  output->WriteTObject(GMTrackY_phelix_All.get());
  output->WriteTObject(GMTrackDCA_phelix_All.get());
  //GM Information PiK
  output->WriteTObject(GMTrackDfIP_PiK.get());
  output->WriteTObject(GMTrackPt_PiK.get());
  output->WriteTObject(GMTrackMCDfIP_PiK.get());
  //GM Information nonPiK
  output->WriteTObject(GMTrackDfIP_nonPiK.get());
  output->WriteTObject(GMTrackPt_nonPiK.get());
  output->WriteTObject(GMTrackMCDfIP_nonPiK.get());
  //GM Information Correct
  //All
  output->WriteTObject(GMTrackMother_Correct.get());
  output->WriteTObject(GMTrackDfIP_InCorrect.get());
  output->WriteTObject(GMTrackDfIP_Correct.get());
  //PiK nonPiK
  output->WriteTObject(GMTrackMCDfIP_Correct_PiK.get());
  output->WriteTObject(GMTrackDfIP_InCorrect_PiK.get());
  output->WriteTObject(GMTrackMCDfIP_Correct_nonPiK.get());
  output->WriteTObject(GMTrackDfIP_InCorrect_nonPiK.get());
  output->WriteTObject(GMTrackDfIP_Correct_PiK.get());
  output->WriteTObject(GMTrackDfIP_Correct_nonPiK.get());
  //Vector
  output->WriteTObject(GMTrackDfIP_Correct_Vector.get());
  //Heavy
  output->WriteTObject(GMTrackDfIP_Correct_Heavy.get());
  //Other
  output->WriteTObject(GMTrackDfIP_Correct_Other.get());
  //All
  output->WriteTObject(GMTrackPt_All.get());
  //Pion Keon
  output->WriteTObject(GMTrackPt_PiK.get());
  output->WriteTObject(GMTrackPt_Pi.get());
  output->WriteTObject(GMTrackPt_K.get());
  //Heavy Flavor
  output->WriteTObject(GMTrackPt_HF.get());
  output->WriteTObject(GMTrackPt_D.get());
  output->WriteTObject(GMTrackPt_B.get());
  //Low mass Vector meson
  output->WriteTObject(GMTrackPt_LVM.get());
  output->WriteTObject(GMTrackPt_Rou.get());
  output->WriteTObject(GMTrackPt_Omega.get());
  output->WriteTObject(GMTrackPt_Phi.get());

  output->WriteTObject(GMTrackGrandMother_Rou.get());
  //Quarkonia
  output->WriteTObject(GMTrackPt_Quarkonia.get());
  output->WriteTObject(GMTrackPt_JPsi.get());
  output->WriteTObject(GMTrackPt_Upsilon.get());
  //Direct
  output->WriteTObject(GMTrackPt_Direct.get());
  //OtehrMuon
  output->WriteTObject(GMTrackPt_OtherMuon.get());
  //NonMuon
  output->WriteTObject(GMTrackPt_OtherNonMuon.get());

  //MC information for GM tracks
  output->WriteTObject(MCGMTrackZ_All.get());
  output->WriteTObject(MCGMTrackZ_PiK.get());
  output->WriteTObject(MCGMTrackZ_nonPiK.get());
  //MC + GM
  output->WriteTObject(MCGMmatchDfIP_All.get());
  output->WriteTObject(MCGMmatchDfIP_PiK.get());
  output->WriteTObject(MCGMmatchDfIP_nonPiK.get());
  output->WriteTObject(MCGMmatchDfIP_InCorrect_PiK.get());
  output->WriteTObject(MCGMmatchDfIP_InCorrect_nonPiK.get());
  //MCH Track
  output->WriteTObject(MCHTrackP_All.get());
  output->WriteTObject(MCHTrackPt_All.get());
  output->WriteTObject(MCHTrackPx_All.get());
  output->WriteTObject(MCHTrackPy_All.get());
  output->WriteTObject(MCHTrackPz_All.get());
  output->WriteTObject(MCHTrackX_All.get());
  output->WriteTObject(MCHTrackY_All.get());
  output->WriteTObject(MCHTrackZ_All.get());
  output->WriteTObject(MCHTrackDCA_All.get());
  output->WriteTObject(MCHTrackEta_All.get());
  output->WriteTObject(MCHTrackPhi_All.get());
  output->WriteTObject(MCHTrackTheta_All.get());
  //MCH propagated linear
  output->WriteTObject(MCHTrackX_propagated_All.get());
  output->WriteTObject(MCHTrackY_propagated_All.get());
  output->WriteTObject(MCHTrackDCA_propagated_All.get());
  //MCH propagated zhelix
  output->WriteTObject(MCHTrackX_phelix_All.get());
  output->WriteTObject(MCHTrackY_phelix_All.get());
  output->WriteTObject(MCHTrackDCA_phelix_All.get());
  //MCH PiK & nonPiK
  output->WriteTObject(MCHTrackDCA_PiK.get());
  output->WriteTObject(MCHTrackDCA_nonPiK.get());
  output->WriteTObject(MCHTrackX_phelix_PiK.get());
  output->WriteTObject(MCHTrackY_phelix_PiK.get());
  output->WriteTObject(MCHTrackDCA_phelix_PiK.get());
  output->WriteTObject(MCHTrackX_phelix_nonPiK.get());
  output->WriteTObject(MCHTrackY_phelix_nonPiK.get());
  output->WriteTObject(MCHTrackDCA_phelix_nonPiK.get());
//TGraph
  //Correct
  DCA_Eff->Write();
  DCA_Puri->Write();
  //Puri_Eff->Write();
  DCA_Coll->Write();
  DCA_Sig->Write();
  DCA_Real->Write();
  DCA_Rem->Write();
  DCA_PiKP->Write();
  DCA_SN->Write();
  DCA_SNW->Write();
  //All
  DCA_Puri_All->Write();
  DCA_Coll_All->Write();
  DCA_Sig_All->Write();
  DCA_Real_All->Write();
  DCA_Rem_All->Write();
  DCA_PiKP_All->Write();
  DCA_SN_All->Write();
  DCA_SNW_All->Write();
  //PiK decay point Check
  PiK_Decay_Rate_z->Write();
  //MCGM
  DCA_Rem_MCGM->Write();
  DCA_SN_MCGM->Write();
  DCA_SNW_MCGM->Write();
  //Matching Chi2
  GM_Chi2_SN->Write();
  GM_Chi2_SNW->Write();
  //MCH DCA
  DCA_SN_MCH->Write();
  DCA_SNW_MCH->Write();
  //Pt
  Eff_nonPiK->Write();
  //Eff_DCAcut->Write();
  //pT dis corrected
  pT_corrected_All->Write();
}
