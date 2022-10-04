#include "ReconstructionDataFormats/GlobalFwdTrack.h"
#include "DataFormatsMFT/TrackMFT.h"

#include<iostream>
#include<memory>

#include "ReconstructionDataFormats/TrackFwd.h"
#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "SimulationDataFormat/TrackReference.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/BaseHits.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TTree.h"
#include <TGeoGlobalMagField.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TEfficiency.h>

double getToXYplane(double x, double y, double z, double theta, double phi){
  if(cos(theta == 0)) return pow(pow(x,2)+pow(y,2),0.5);
  auto Dv = pow(pow(x,2)+pow(y,2)+pow(z,2),0.5);
  auto Dpv = fabs(z/cos(theta));
  auto Dp = pow(Dpv*Dpv + Dv*Dv -2*Dpv*(x*sin(theta)*cos(phi) + y*sin(theta)*sin(phi) + z*cos(theta)),0.5);
  printf("Distance from Vertex to IP = %f\n", Dp);
  return Dp;
}

void MCTrackChecker_Comp( const std::string o2sim_KineFile = "o2sim_Kine.root",
                     const std::string sig_KineFile = "sgn_Kine.root"){

  // Set MC track information
  std::unique_ptr<TFile> o2sim_KineFileIn(new TFile(o2sim_KineFile.c_str()));
  std::unique_ptr<TTree> o2SimKineTree((TTree*)o2sim_KineFileIn->Get("o2sim"));
  vector<o2::MCTrackT<float>>mcTrVec, *mcTr = &mcTrVec;
  o2SimKineTree->SetBranchAddress("MCTrack", &mcTr);
  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  o2SimKineTree->SetBranchAddress("MCEventHeader.", &eventHeader);
  vector<o2::TrackReference>* trackRefs = nullptr;
  o2SimKineTree->SetBranchAddress("TrackRefs", &trackRefs);
  o2SimKineTree->GetEntry(0);

  //Parepare output objects
  TFile* output = new TFile("MCTrackCompInfo.root","recreate");
  std::unique_ptr<TH1F> MCTrackdDeltaTheta_muon_pion(new TH1F("MCTrackdDeltaTheta_muon_pion","",20000,-10,10));
  //Tree Object
  //Information_All
  float P_All,Eta_All,DfIP_All,vx_All,vy_All,vz_All,vt_All,Mother_MC,Pt_All,vdeltatheta_All;
  int PDG_All,Mother_All,SaveId_All,MotherId_All;
  std::unique_ptr<TTree> TotalTreeObject_Comp(new TTree("TotalTreeObject_All","TotalTreeObject_All"));
  TotalTreeObject_Comp->Branch("P_All",&P_All,"P_All/F");
  TotalTreeObject_Comp->Branch("Pt_All",&Pt_All,"Pt_All/F");
  TotalTreeObject_Comp->Branch("Eta_All",&Eta_All,"Eta_All/F");
  TotalTreeObject_Comp->Branch("DfIP_All",&DfIP_All,"DfIP_All/F");
  TotalTreeObject_Comp->Branch("vx_All",&vx_All,"vx_All/F");
  TotalTreeObject_Comp->Branch("vy_All",&vy_All,"vy_All/F");
  TotalTreeObject_Comp->Branch("vz_All",&vz_All,"vz_All/F");
  TotalTreeObject_Comp->Branch("vt_All",&vt_All,"vt_All/F");
  TotalTreeObject_Comp->Branch("vdeltatheta_All",&vdeltatheta_All,"vdeltatheta_All/F");
  TotalTreeObject_Comp->Branch("PDG_All",&PDG_All,"PDG_All/I");
  TotalTreeObject_Comp->Branch("Mother_All",&Mother_All,"Mother_All/I");
  TotalTreeObject_Comp->Branch("MotherId_All",&MotherId_All,"MotherId_All/I");
  TotalTreeObject_Comp->Branch("SaveId_All",&SaveId_All,"SaveId_All/I");

  float vx_MC,vy_MC,vz_MC,P_MC,Pt_MC,vt_MC,vr_MC,vtheta_MC,vphi_MC,DistanceIP3D_MC,DistanceIP2D_MC;
  int pdgcode_MC,MotherId,FirstDaughterId,SaveId;
  int C_myu,C_e,C_p,C_n,C_pie,C_gamma,C_K,C_D,C_B,C_Jpsi,C_rou,C_omega,C_phi;
  //Check muon_theta - pion_theta
  float vtheta_delta;
  int nEvent = o2SimKineTree->GetEntries();
  for (int iEvent=0; iEvent < nEvent; iEvent++){
  	//Get tree information of the "iTrack"th line
  	//int iEvent = 0;
  	o2SimKineTree->GetEntry(iEvent);
	std::vector<int> DuplicationChecker(1);
        DuplicationChecker[0] = -1;
  	int nMCTracks = (*mcTr).size();
  	for (Int_t iMCTrack=0; iMCTrack<(*mcTr).size();++iMCTrack) {
	       	o2::MCTrackT<float> *thisTrack = &(*mcTr).at(iMCTrack);
  	     	SaveId = iMCTrack;
		MotherId = iMCTrack;
      		pdgcode_MC = thisTrack->GetPdgCode();//get particle pdg code
      		int Mother_MC = pdgcode_MC;
		while(fabs(pdgcode_MC) == fabs(Mother_MC)){
			SaveId = MotherId;
			vx_MC = thisTrack->GetStartVertexCoordinatesX();
                	vy_MC = thisTrack->GetStartVertexCoordinatesY();
                	vz_MC = thisTrack->GetStartVertexCoordinatesZ();
                	vt_MC = thisTrack->GetStartVertexCoordinatesT();
                	vr_MC = thisTrack->GetEta();
                  vtheta_MC = thisTrack->GetTheta();
                	P_MC = thisTrack->GetP();
                	Pt_MC = thisTrack->GetPt();
			MotherId = thisTrack->getMotherTrackId();
			if(MotherId==-1) break;
			thisTrack = &(mcTr->at(MotherId));
			Mother_MC = thisTrack->GetPdgCode();
		}
		//thisTrack = &(mcTr->at(SaveId));

		//std::cout<<"size = "<<DuplicationChecker.size()<<std::endl;
		//std::cout<<"SaveId = "<<SaveId<<std::endl;
		int DuplicationId = 0;
                for(int i = 0;i<DuplicationChecker.size();i++){
			//std::cout<<"[i] = "<<DuplicationChecker[i]<<std::endl;
                	if(DuplicationChecker[i] == SaveId){
			//std::cout<<""<<std::endl;
                        	DuplicationId = i;
                        	break;
                	}
                }

                if(DuplicationChecker[DuplicationId] == SaveId) continue;
                DuplicationChecker.push_back(SaveId);

      DistanceIP3D_MC = pow(pow(vx_MC,2)+pow(vy_MC,2)+pow(vz_MC,2),0.5);
		  DistanceIP2D_MC = getToXYplane(vx_MC,vy_MC,vz_MC,vtheta_MC,vphi_MC);

      if(fabs(pdgcode_MC)!=13) continue;
      if(vz_MC<-46) continue;
      if(P_MC<4) continue;
      if(vr_MC<-3.6 || vr_MC>-2.45) continue;

      //Check effect of lorenz boost
      if(fabs(pdgcode_MC)==13){
        if(fabs(Mother_MC==211) || fabs(Mother_MC==321) || fabs(Mother_MC==130)){
          std::cout<<"muon p = "<<P_MC<<std::endl;
          std::cout<<"muon theta = "<<vtheta_MC<<std::endl;
          vtheta_delta = vtheta_MC;
          P_MC = thisTrack->GetP();
          vtheta_MC = thisTrack->GetTheta();
          std::cout<<"pion keon p = "<<P_MC<<std::endl;
          std::cout<<"pion keon theta = "<<vtheta_MC<<std::endl;
          vtheta_delta = vtheta_MC - vtheta_delta;
          std::cout<<"delta_theta = "<<vtheta_delta<<std::endl;
          std::cout<<""<<std::endl;
          vtheta_delta = vtheta_delta*360/(2*TMath::Pi());
          MCTrackdDeltaTheta_muon_pion->Fill(vtheta_delta);
        }
      }
		  //Fill Branch
		  P_All = P_MC;
		  vz_All = vz_MC;
		  vy_All = vy_MC;
		  vx_All = vx_MC;
		  vt_All = vt_MC;
		  Pt_All = Pt_MC;
      Eta_All = vr_MC;
      DfIP_All = DistanceIP2D_MC;
      PDG_All = fabs(pdgcode_MC);
      Mother_All = fabs(Mother_MC);
      MotherId_All = MotherId;
      vdeltatheta_All = vtheta_delta;
      SaveId_All = SaveId;
      TotalTreeObject_Comp->Fill();
    }
  }
  //Save Objects
  output->WriteTObject(TotalTreeObject_Comp.get());
  output->WriteTObject(MCTrackdDeltaTheta_muon_pion.get());
  std::cout<<"Analysis complete"<<std::endl;
}
