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

void SumTree_MFT(){
  // Open root file and get Object
  // MC Information
  std::unique_ptr<TH1F> MCTrackPDG_dis_All(new TH1F("MCTrackPDG_dis_All","",10000,0-0.05,1000-0.05));
  std::unique_ptr<TH1F> MCTrackPDG_dis_Muon(new TH1F("MCTrackPDG_dis_Muon","",10000,0-0.05,1000-0.05));
  std::unique_ptr<TH1F> MCTrackDfIP_dis_Heavy(new TH1F("MCTrackDfIP_dis_Heavy","",10000,0-0.005,100-0.005));
  std::unique_ptr<TH1F> MCTrackPt_Muon_SigMuon(new TH1F("MCTrackPt_Muon_SigMuon","",1000,0-0.05,100-0.05));

  // MFT Information All
  std::unique_ptr<TH1F> MFTTrackP_All(new TH1F("MFTTrackP_All","",10000,0-0.05,1000-0.05));
  std::unique_ptr<TH1F> MFTTrackEta_All(new TH1F("MFTTrackEta_All","",20000,-10-0.0005,10-0.0005));
  std::unique_ptr<TH1F> MFTTrackX_All(new TH1F("MFTTrackX_All","",20000,-1000-0.005,1000-0.005));
  std::unique_ptr<TH1F> MFTTrackY_All(new TH1F("MFTTrackY_All","",20000,-1000-0.005,1000-0.005));
  std::unique_ptr<TH1F> MFTTrackZ_All(new TH1F("MFTTrackZ_All","",20000,-1000-0.005,1000-0.005));
  std::unique_ptr<TH1F> MFTTrackDfIP_All(new TH1F("MFTTrackDfIP_All","",10000,0,100));
  //std::unique_ptr<TH1F> MFTTrackDfIP_All(new TH1F("MFTTrackDfIP_All","",20000,0,100));
  std::unique_ptr<TH1F> MFTTrackDfIP_Muon_All(new TH1F("MFTTrackDfIP_Muon_All","",10000,0,100));
  //std::unique_ptr<TH1F> MFTTrackPt_DCAcut_021(new TH1F("MFTTrackPt_DCAcut_021","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> MFTTrackPDG_All(new TH1F("MFTTrackPDG_All","",80000,0-0.05,8000-0.05));
  std::unique_ptr<TH1F> MFTTrackMother_All(new TH1F("MFTTrackMother_All","",80000,0-0.05,8000-0.05));

  // MFT Information nonSigMuon
  std::unique_ptr<TH1F> MFTTrackDfIP_nonSigMuon(new TH1F("MFTTrackDfIP_nonSigMuon","",10000,0,100));
  //std::unique_ptr<TH1F> MFTTrackDfIP_nonSigMuon(new TH1F("MFTTrackDfIP_nonSigMuon","",20000,0,100));
  std::unique_ptr<TH1F> MFTTrackMCDfIP_nonSigMuon(new TH1F("MFTTrackMCDfIP_nonSigMuon","",10000,0,100));

  // MFT Information SigMuon
  std::unique_ptr<TH1F> MFTTrackDfIP_SigMuon(new TH1F("MFTTrackDfIP_SigMuon","",10000,0,100));
  std::unique_ptr<TH1F> MFTTrackPt_SigMuon(new TH1F("MFTTrackPt_SigMuon","",1000,0-0.05,100-0.05));
  //std::unique_ptr<TH1F> MFTTrackDfIP_SigMuon(new TH1F("MFTTrackDfIP_SigMuon","",20000,0,100));
  std::unique_ptr<TH1F> MFTTrackMCDfIP_SigMuon(new TH1F("MFTTrackMCDfIP_SigMuon","",10000,0,100));

  //Spectol
  //All
  std::unique_ptr<TH1F> MFTTrackPt_All(new TH1F("MFTTrackPt_All","",1000,0-0.05,100-0.05));
  //Pion Keon
  std::unique_ptr<TH1F> MFTTrackPt_nonSigMuon(new TH1F("MFTTrackPt_nonSigMuon","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> MFTTrackPt_Pi(new TH1F("MFTTrackPt_Pi","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> MFTTrackPt_K(new TH1F("MFTTrackPt_K","",1000,0-0.05,100-0.05));
  //Heavy Flavor
  std::unique_ptr<TH1F> MFTTrackPt_HF(new TH1F("MFTTrackPt_HF","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> MFTTrackPt_D(new TH1F("MFTTrackPt_D","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> MFTTrackPt_B(new TH1F("MFTTrackPt_B","",1000,0-0.05,100-0.05));
  //Low mass Light Vecor meson
  std::unique_ptr<TH1F> MFTTrackPt_LVM(new TH1F("MFTTrackPt_LVM","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> MFTTrackPt_Rou(new TH1F("MFTTrackPt_rou","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> MFTTrackPt_Omega(new TH1F("MFTTrackPt_omega","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> MFTTrackPt_Phi(new TH1F("MFTTrackPt_phi","",1000,0-0.05,100-0.05));

  std::unique_ptr<TH1F> MFTTrackGrandMother_Rou(new TH1F("MFTTrackGrandMother_Rou","",100000,0,10000));
  //Quarkonia
  std::unique_ptr<TH1F> MFTTrackPt_Quarkonia(new TH1F("MFTTrackPt_Quarkonia","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> MFTTrackPt_JPsi(new TH1F("MFTTrackPt_JPsi","",1000,0-0.05,100-0.05));
  std::unique_ptr<TH1F> MFTTrackPt_Upsilon(new TH1F("MFTTrackPt_Upsilon","",1000,0-0.05,100-0.05));
  //Others(Muon)
  std::unique_ptr<TH1F> MFTTrackPt_Direct(new TH1F("MFTTrackPt_Direct","",1000,0-0.05,100-0.05));
  //Others(Muon)
  std::unique_ptr<TH1F> MFTTrackPt_OtherMuon(new TH1F("MFTTrackPt_OtherMuon","",1000,0-0.05,100-0.05));
  //Others(nonMuon)
  std::unique_ptr<TH1F> MFTTrackPt_OtherNonMuon(new TH1F("MFTTrackPt_OtherNonMuon","",1000,0-0.05,100-0.05));

  // MC + MFT
  std::unique_ptr<TH1F> MCMFTmatchDfIP_All(new TH1F("MCMFTmatchDfIP_All","",20000,-100-0.005,100-0.005));
  std::unique_ptr<TH1F> MCMFTmatchDfIP_nonSigMuon(new TH1F("MCMFTmatchDfIP_nonSigMuon","",20000,-100-0.005,100-0.005));
  std::unique_ptr<TH1F> MCMFTmatchDfIP_SigMuon(new TH1F("MCMFTmatchDfIP_SigMuon","",20000,-100-0.005,100-0.005));




  //Prepare output objects
  TFile* output = new TFile("SumResult_SumTree_MFT.root","recreate");
  TChain *SumTotalTree_MC=new TChain("TotalTreeObject_All","SumTotalTree_MC_All");
  TChain *SumTotalTree_MFT=new TChain("TotalTreeObject_MFT","SumTotalTree_MFT_All");
  //Fill SumObject
  SumTotalTree_MC -> Add("MCTrackCompInfo*.root");
  SumTotalTree_MFT -> Add("AnalysisResults*.root");

  //Get Branch Information :MFT
  double P_MFT,Pt_MFT,DfIP_MFT,vz_MFT,vy_MFT,vx_MFT,Eta_MFT;
  double P_MCMFT,Pt_MCMFT,DfIP_MCMFT,vz_MCMFT,vy_MCMFT,vx_MCMFT,Eta_MCMFT;
  double P_GrandMCMFT,Pt_GrandMCMFT,DfIP_GrandMCMFT,vz_GrandMCMFT,vy_GrandMCMFT,vx_GrandMCMFT,Eta_GrandMCMFT;
  int SourceCheck,MotherId_MCMFT;
  int Mother_MCMFT,PDG_MCMFT,GrandMother_MCMFT;
  SumTotalTree_MFT ->SetBranchAddress("DfIP_MFT",&DfIP_MFT);
  SumTotalTree_MFT ->SetBranchAddress("P_MFT",&P_MFT);
  SumTotalTree_MFT ->SetBranchAddress("Pt_MFT",&Pt_MFT);
  SumTotalTree_MFT ->SetBranchAddress("vx_MFT",&vx_MFT);
  SumTotalTree_MFT ->SetBranchAddress("vy_MFT",&vy_MFT);
  SumTotalTree_MFT ->SetBranchAddress("vz_MFT",&vz_MFT);
  SumTotalTree_MFT ->SetBranchAddress("Eta_MFT",&Eta_MFT);
  SumTotalTree_MFT ->SetBranchAddress("DfIP_MCMFT",&DfIP_MCMFT);
  SumTotalTree_MFT ->SetBranchAddress("P_MCMFT",&P_MCMFT);
  SumTotalTree_MFT ->SetBranchAddress("Pt_MCMFT",&Pt_MCMFT);
  SumTotalTree_MFT ->SetBranchAddress("vx_MCMFT",&vx_MCMFT);
  SumTotalTree_MFT ->SetBranchAddress("vy_MCMFT",&vy_MCMFT);
  SumTotalTree_MFT ->SetBranchAddress("vz_MCMFT",&vz_MCMFT);
  SumTotalTree_MFT ->SetBranchAddress("Eta_MCMFT",&Eta_MCMFT);
  SumTotalTree_MFT ->SetBranchAddress("PDG_MCMFT",&PDG_MCMFT);
  SumTotalTree_MFT ->SetBranchAddress("Mother_MCMFT",&Mother_MCMFT);
  SumTotalTree_MFT ->SetBranchAddress("MotherId_MCMFT",&MotherId_MCMFT);
  //SumTotalTree_MFT ->SetBranchAddress("DfIP_GrandMCMFT",&DfIP_GrandMCMFT);
  //SumTotalTree_MFT ->SetBranchAddress("P_GrandMCMFT",&P_GrandMCMFT);
  //SumTotalTree_MFT ->SetBranchAddress("Pt_GrandMCMFT",&Pt_GrandMCMFT);
  //SumTotalTree_MFT ->SetBranchAddress("vx_GrandMCMFT",&vx_GrandMCMFT);
  //SumTotalTree_MFT ->SetBranchAddress("vy_GrandMCMFT",&vy_GrandMCMFT);
  //SumTotalTree_MFT ->SetBranchAddress("vz_GrandMCMFT",&vz_GrandMCMFT);
  //SumTotalTree_MFT ->SetBranchAddress("Eta_GrandMCMFT",&Eta_GrandMCMFT);
  SumTotalTree_MFT ->SetBranchAddress("GrandMother_MCMFT",&GrandMother_MCMFT);
  SumTotalTree_MFT ->SetBranchAddress("SourceCheck",&SourceCheck);


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

  for(int iEntry_MC=0 ;iEntry_MC < nEntry_MC;++iEntry_MC){
  	SumTotalTree_MC->GetEntry(iEntry_MC);
	  MCTrackPDG_dis_All->Fill(PDG_All);
		MCTrackPDG_dis_Muon->Fill(Mother_All);
		if(Mother_All == 130 || Mother_All == 211 || Mother_All == 321 ) continue;
		MCTrackPt_Muon_SigMuon->Fill(Pt_All);
  }
  cout<<"Entry = "<<nEntry_MC<<endl;

  int nEntry_MFT = SumTotalTree_MFT->GetEntries();

  int q;//charge

  for(int iEntry_MFT=0 ;iEntry_MFT < nEntry_MFT;++iEntry_MFT){
        SumTotalTree_MFT->GetEntry(iEntry_MFT);
	//if(P_MFT < 4) continue;
  //if(DfIP_MFT>0.1) continue;
        //if(MatchingScore>5.25) continue;
        //if(DfIP_MFT>0.030) continue;
        MFTTrackX_All->Fill(vx_MFT);
        MFTTrackY_All->Fill(vy_MFT);
        MFTTrackP_All->Fill(P_MFT);
        MFTTrackPt_All->Fill(Pt_MFT);//Spectol All
        MFTTrackEta_All->Fill(Eta_MFT);
        MFTTrackDfIP_All->Fill(DfIP_MFT);

	//MC information check
	if(SourceCheck==0){//MFTTrack Source==false
    std::cout<<"No source data"<<std::endl;
		MFTTrackPt_OtherNonMuon->Fill(Pt_MFT);//Spectol NonMuon
	}else{//MFTTrack source==true
		MFTTrackZ_All->Fill(vz_MCMFT);
		MFTTrackPDG_All->Fill(fabs(PDG_MCMFT));

		if(fabs(PDG_MCMFT)==13){//All muon track
			MFTTrackMother_All->Fill(fabs(Mother_MCMFT));
			MFTTrackDfIP_Muon_All->Fill(DfIP_MFT);
			if(fabs(Mother_MCMFT) == 130 || fabs(Mother_MCMFT) == 211 || fabs(Mother_MCMFT) == 321 ){//Mother is nonSigMuon
				MFTTrackPt_nonSigMuon->Fill(Pt_MFT);//Spectol Pion Keon
				MFTTrackDfIP_nonSigMuon->Fill(DfIP_MFT);
				MFTTrackMCDfIP_nonSigMuon->Fill(DfIP_MCMFT);
        if(fabs(Mother_MCMFT) == 211){//Mother is Pion
          MFTTrackPt_Pi->Fill(Pt_MFT);//Spectol Pion
        }else{
          MFTTrackPt_K->Fill(Pt_MFT);//Spectol Keon
        }
			}else{//Mother is SigMuon
				MFTTrackPt_SigMuon->Fill(Pt_MFT);
				MFTTrackDfIP_SigMuon->Fill(DfIP_MFT);
				MFTTrackMCDfIP_SigMuon->Fill(DfIP_MCMFT);
        if(fabs(Mother_MCMFT) == 421 || fabs(Mother_MCMFT) == 411 || fabs(Mother_MCMFT) == 431 || fabs(Mother_MCMFT) == 511|| fabs(Mother_MCMFT) == 521 || fabs(Mother_MCMFT) == 531){
              MFTTrackPt_HF->Fill(Pt_MFT);//Spectol Heavy Flaver
              if(fabs(Mother_MCMFT) == 521 || fabs(Mother_MCMFT) == 531|| fabs(Mother_MCMFT) == 511){
                    MFTTrackPt_B->Fill(Pt_MFT);//Spectol B-meson
              }else{
                    MFTTrackPt_D->Fill(Pt_MFT);//Spectol D-meson
              }
        }else if(fabs(Mother_MCMFT) == 113 || fabs(Mother_MCMFT) == 223 || fabs(Mother_MCMFT) == 333){
              MFTTrackPt_LVM->Fill(Pt_MFT);//Spetol Low mass Light Vector meson
              if(fabs(Mother_MCMFT) == 333)MFTTrackPt_Phi->Fill(Pt_MFT);//Spectol Phi
              if(fabs(Mother_MCMFT) == 113){
			MFTTrackPt_Rou->Fill(Pt_MFT);//Spectol Rou
			MFTTrackGrandMother_Rou->Fill(fabs(GrandMother_MCMFT));
	      }
              if(fabs(Mother_MCMFT) == 223)MFTTrackPt_Omega->Fill(Pt_MFT);//Spectol Omega
        }else	if(fabs(Mother_MCMFT) == 443 || fabs(Mother_MCMFT) == 553 || fabs(Mother_MCMFT) == 100553){
              MFTTrackPt_Quarkonia->Fill(Pt_MFT);//Spectol Quarkonia
              if(fabs(Mother_MCMFT) == 443)MFTTrackPt_JPsi->Fill(Pt_MFT);//Spectol JPsi
              if(fabs(Mother_MCMFT) == 553 || fabs(Mother_MCMFT) == 100553){
			MFTTrackPt_Upsilon->Fill(Pt_MFT);//Spectol upsilon
	      		std::cout<<"Other Muon Mother = "<<Mother_MCMFT<<std::endl;
		}
        }else if(fabs(Mother_MCMFT == 22)){
              MFTTrackPt_Direct->Fill(Pt_MFT);//Spectol Other Muon
        }else{
	      std::cout<<"Other Muon Mother = "<<Mother_MCMFT<<std::endl;
              MFTTrackPt_OtherMuon->Fill(Pt_MFT);//Spectol Other Muon
        }
			}
  		}else{//nonMuon MFTTrack
			MFTTrackPt_OtherNonMuon->Fill(Pt_MFT);
      MFTTrackPt_nonSigMuon->Fill(Pt_MFT);//Spectol Pion Keon
			//MFTTrackDfIP_nonSigMuon->Fill(DfIP_MFT);
			MFTTrackMCDfIP_nonSigMuon->Fill(DfIP_MCMFT);
		}
  }
  }
  cout<<"Entry = "<<nEntry_MFT<<endl;

  //puriry efficiency DCA cut discussion
  int n = 1000;
  double Count_True_Muon[n];
  double Count_nonSigMuon_All[n];
  double DfIP_axis[n];
  double E_Chi2_axis[n];
  double E_DfIP[n];//Error parametor

  double Count_MFTTrack_All[n];
  double Count_MFTTrack_All_nonSigMuon[n],Count_MFTTrack_All_SigMuon[n];
  double Efficiency[n],Purity[n],SignalPurity[n],RealEfficiency[n],RemovalnonSigMuon[n],PuritynonSigMuon[n],SNRate[n],SNRateW[n],Collection[n];
  double Purity_All[n],SignalPurity_All[n],RealEfficiency_All[n],RemovalnonSigMuon_All[n],PuritynonSigMuon_All[n],SNRate_All[n],SNRateW_All[n],Collection_All[n];
  double E_Eff[n],E_Puri[n],E_Coll[n],E_Real[n],E_Sig[n],E_Rem[n],E_nonSigMuonP[n],E_SN[n],E_SNW[n];
  double E_Eff_All[n],E_Puri_All[n],E_Coll_All[n],E_Real_All[n],E_Sig_All[n],E_Rem_All[n],E_nonSigMuonP_All[n],E_SN_All[n],E_SNW_All[n];

  double Count_MCMFTTrack_SigMuon[n],Count_MCMFTTrack_nonSigMuon[n],Count_MCMFTTrack_All[n],Count_MCMFTTrack_nonSigMuon_All[n];
  double RemovalnonSigMuon_MCMFT[n],SNRate_MCMFT[n],SNRateW_MCMFT[n];
  double E_Rem_MCMFT[n],E_SN_MCMFT[n],E_SNW_MCMFT[n];


  //Pt
  double Pt_axis[n],E_Pt[n];
  double Count_True_Muon_Pt[n],Count_Track_Pt[n],Count_Track_nonSigMuon_Pt[n],Count_Track_SigMuon_Pt[n],Count_Track_DCAcut[n];
  double Efficiency_SigMuon[n],Efficiency_DCAcut[n],Efficiency_SigMuon_Signal[n],Efficiency_SigMuon_Other[n];
  double E_Eff_SigMuon[n],E_Eff_DCAcut[n];

  double Efficiency_Even;

  for(int i = 0; i<n; ++i){
	std::cout<<"roop = "<<i<<std::endl;
	Count_nonSigMuon_All[i] = MFTTrackDfIP_nonSigMuon->GetEntries();
	Count_True_Muon[i] = MCTrackPt_Muon_SigMuon->GetEntries();
  	//std::cout<<"Count True Muon = "<<Count_True_Muon[i]<<std::endl;
  	//std::cout<<"Count nonSigMuon All = "<<Count_nonSigMuon_All[i]<<std::endl;
	//DfIP
	if(i==0){//初期値を入れる
		  Count_MFTTrack_All[i] = MFTTrackDfIP_All->GetBinContent(i+1);
 	}
	Count_MFTTrack_All[i] = Count_MFTTrack_All[i-1] + MFTTrackDfIP_All->GetBinContent(i+1);


	//MCMFT
	Count_MCMFTTrack_nonSigMuon_All[i] = MFTTrackMCDfIP_nonSigMuon->GetEntries();
	if(i==0){
		Count_MCMFTTrack_SigMuon[i] = MFTTrackMCDfIP_SigMuon->GetBinContent(i+1);
		Count_MCMFTTrack_nonSigMuon[i] = MFTTrackMCDfIP_nonSigMuon->GetBinContent(i+1);
	}
	Count_MCMFTTrack_SigMuon[i] = Count_MCMFTTrack_SigMuon[i-1] + MFTTrackMCDfIP_SigMuon->GetBinContent(i+1);
	Count_MCMFTTrack_nonSigMuon[i] = Count_MCMFTTrack_nonSigMuon[i-1] + MFTTrackMCDfIP_nonSigMuon->GetBinContent(i+1);
	Count_MCMFTTrack_All[i] = Count_MCMFTTrack_nonSigMuon[i] + Count_MCMFTTrack_SigMuon[i];
	//All
  	Count_MFTTrack_All_nonSigMuon[i] = Count_MFTTrack_All_nonSigMuon[i-1] + MFTTrackDfIP_nonSigMuon->GetBinContent(i+1);
  	Count_MFTTrack_All_SigMuon[i] = Count_MFTTrack_All_SigMuon[i-1] + MFTTrackDfIP_SigMuon->GetBinContent(i+1);

	//Pt
	Count_True_Muon_Pt[i] = MCTrackPt_Muon_SigMuon->GetBinContent(i+1);
	Count_Track_SigMuon_Pt[i] = MFTTrackPt_SigMuon->GetBinContent(i+1);
	Count_Track_nonSigMuon_Pt[i] = MFTTrackPt_nonSigMuon->GetBinContent(i+1);
	Count_Track_Pt[i] = Count_Track_SigMuon_Pt[i]+Count_Track_nonSigMuon_Pt[i];
	//Count_Track_DCAcut[i] = MFTTrackPt_DCAcut_021->GetBinContent(i+1);

	//DfIP
	//DfIP_axis[i] = (i+1)*0.005;
	//E_DfIP[i] = 0.0025;
	DfIP_axis[i] = (i+1)*0.01;
	E_DfIP[i] = 0;

	//Pt
	Pt_axis[i] = i*0.1;
	E_Pt[i] = 0.05;
	if(i<300){
	   //DfIP
     std::cout<<"Count_MFTTrack_All = "<<Count_MFTTrack_All[i]<<std::endl;
	   std::cout<<"Count_MFTTrack_All_nonSigMuon = "<<Count_MFTTrack_All_nonSigMuon[i]<<std::endl;
     std::cout<<"Count_MFTTrack_All_SigMuon = "<<Count_MFTTrack_All_SigMuon[i]<<std::endl;

	   //Pt
     //std::cout<<"Pt axis = "<<Pt_axis[i]<<std::endl;
     //std::cout<<"Count_True_Muon_Pt = "<<Count_True_Muon_Pt[i]<<std::endl;
     //std::cout<<"Count_Track_Pt = "<<Count_Track_Pt[i]<<std::endl;
     //std::cout<<"Count_Track_nonSigMuon_Pt = "<<Count_Track_nonSigMuon_Pt[i]<<std::endl;
     //std::cout<<"Count_Track_SigMuon_Pt = "<<Count_Track_SigMuon_Pt[i]<<std::endl;
     //std::cout<<"Count_Track_DCAcut = "<<Count_Track_DCAcut[i]<<std::endl;
	}
	//Out put
	Efficiency[i] = Count_MFTTrack_All[i]/Count_True_Muon[i];;
	Purity_All[i] = Count_MFTTrack_All_SigMuon[i]/Count_MFTTrack_All[i];



	//All
	PuritynonSigMuon_All[i] = Count_MFTTrack_All_nonSigMuon[i]/Count_nonSigMuon_All[i]; //Rate of nonSigMuon
	SNRate_All[i] = Count_MFTTrack_All_SigMuon[i]/Count_MFTTrack_All_nonSigMuon[i];
	SNRateW_All[i] = pow(Count_MFTTrack_All_SigMuon[i],1.5)/Count_MFTTrack_All_nonSigMuon[i];

	//MCMFT
	if(Count_MCMFTTrack_nonSigMuon[i]>=1){
		RemovalnonSigMuon_MCMFT[i] = 1. - (Count_MCMFTTrack_nonSigMuon[i]/Count_MCMFTTrack_nonSigMuon_All[i]);
		SNRate_MCMFT[i] = Count_MCMFTTrack_SigMuon[i]/Count_MCMFTTrack_nonSigMuon[i];
		SNRateW_MCMFT[i] = pow(Count_MCMFTTrack_SigMuon[i],1.5)/Count_MCMFTTrack_nonSigMuon[i];
	}

	//Pt
	if(Count_True_Muon_Pt[i] >= 1){
		Efficiency_SigMuon[i] = (Count_Track_Pt[i]/Count_True_Muon_Pt[i]);
		Efficiency_SigMuon_Signal[i] = Count_Track_SigMuon_Pt[i]/Count_True_Muon_Pt[i];
		Efficiency_SigMuon_Other[i] = Count_Track_nonSigMuon_Pt[i]/Count_True_Muon_Pt[i];
		//Efficiency_DCAcut[i] = Count_Track_DCAcut[i]/Count_True_Muon_Pt[i];
	}
	Efficiency_Even += Efficiency_SigMuon[i]*(Count_True_Muon_Pt[i]/Count_True_Muon[i]);

  if(i<300){
	   //std::cout<<"Efficiency_Even = "<<Efficiency_Even<<std::endl;

	   //std::cout<<"Efficiency_SigMuon = "<<Efficiency_SigMuon[i]<<std::endl;
	   //std::cout<<"Efficiency_SigMuon_Signal = "<<Efficiency_SigMuon_Signal[i]<<std::endl;
	   //std::cout<<"Efficiency_SigMuon_Other = "<<Efficiency_SigMuon_Other[i]<<std::endl;
	   //Analysis_All = PuritynonSigMuon_All/RealEfficiency_All;
	   //std::cout<<"DCA cut = "<<DfIP_axis[i]<<std::endl;
	   //std::cout<<"Efficieny = "<<Efficiency[i]<<std::endl;
	   //std::cout<<"Purity = "<<Purity[i]<<std::endl;
	   //std::cout<<"Sig = "<<SignalPurity[i]<<std::endl;
	   //std::cout<<"Real Eff = "<<RealEfficiency[i]<<std::endl;
	   //std::cout<<"nonSigMuon Purity = "<<PuritynonSigMuon[i]<<std::endl;
	   //std::cout<<"SNRate = "<<SNRate[i]<<std::endl;
	   //std::cout<<"SNRateW = "<<SNRateW[i]<<std::endl;
	   //std::cout<<"Purity_All = "<<Purity_All[i]<<std::endl;
	   //std::cout<<"Sig_All = "<<SignalPurity_All[i]<<std::endl;
	   //std::cout<<"Real Eff_All = "<<RealEfficiency_All[i]<<std::endl;
	   //std::cout<<"nonSigMuon Purity_All = "<<PuritynonSigMuon_All[i]<<std::endl;
     std::cout<<"DCA cut  = "<<DfIP_axis[i]<<std::endl;
	   std::cout<<"SNRate_All = "<<SNRate_All[i]<<std::endl;
     std::cout<<"SNRateW_All = "<<SNRateW_All[i]<<std::endl;
     std::cout<<""<<std::endl;
	}
	//Error
	E_Eff[i] = pow(Count_MFTTrack_All[i],0.5)/Count_True_Muon[i];
		E_SNW[i] = pow(SNRate[i]*SNRate[i]*(1.+SNRate[i]),0.5);

	//All
 	if(Count_MFTTrack_All_nonSigMuon[i] >= 1){
		E_SN_All[i] = pow(SNRate_All[i]*(1.+SNRate_All[i])/Count_MFTTrack_All_nonSigMuon[i],0.5);
		E_SNW_All[i] = pow(SNRate_All[i]*SNRate_All[i]*(1.+SNRate_All[i]),0.5);
	}
 	E_Coll_All[i] = pow(Count_MFTTrack_All_SigMuon[i],0.5)/40000.;
 	E_Sig_All[i] = pow(Count_MFTTrack_All_SigMuon[i],0.5)/Count_MFTTrack_All[i];
 	E_Real_All[i] = pow(Count_MFTTrack_All_SigMuon[i],0.5)/Count_MFTTrack_All[i];
 	E_Rem_All[i] = pow(Count_nonSigMuon_All[i] - Count_MFTTrack_All_nonSigMuon[i],0.5)/Count_nonSigMuon_All[i];
 	E_nonSigMuonP_All[i] = pow(Count_MFTTrack_All_nonSigMuon[i],0.5)/Count_nonSigMuon_All[i];

	//MCMFT
	if(Count_MCMFTTrack_nonSigMuon[i]>=1){
		E_Rem_MCMFT[i] = pow(Count_MCMFTTrack_nonSigMuon_All[i] - Count_MCMFTTrack_nonSigMuon[i],0.5)/Count_MCMFTTrack_nonSigMuon_All[i];
		E_SN_MCMFT[i] = pow(SNRate_MCMFT[i]*(1.+SNRate_MCMFT[i])/Count_MCMFTTrack_nonSigMuon[i],0.5);
		E_SNW_MCMFT[i] = pow(SNRate_MCMFT[i]*SNRate_MCMFT[i]*(1.+SNRate_MCMFT[i]),0.5);
	}

	//Pt
	if(Count_True_Muon_Pt[i] >=1){
		E_Eff_SigMuon[i] = pow(Efficiency_SigMuon[i]*(1.+ Efficiency_SigMuon[i])/Count_True_Muon_Pt[i],0.5);
		//E_Eff_DCAcut[i] = pow(Count_Track_DCAcut[i],0.5)/Count_True_Muon_Pt[i];
	}

	//transfer eff->%
	if(Count_True_Muon_Pt[i] >=1){
		Efficiency_SigMuon[i] = Efficiency_SigMuon[i]*100;
		E_Eff_SigMuon[i] = E_Eff_SigMuon[i]*100;
	}
  }

  //All
  TGraphErrors* DCA_Puri_All = new TGraphErrors(n,DfIP_axis,Purity_All,E_DfIP,E_Puri_All);
  DCA_Puri_All -> SetName("DCA_Puri_All");
  TGraphErrors* DCA_Coll_All = new TGraphErrors(n,DfIP_axis,Collection_All,E_DfIP,E_Coll_All);
  DCA_Coll_All -> SetName("DCA_Coll_All");
  TGraphErrors* DCA_Sig_All = new TGraphErrors(n,DfIP_axis,SignalPurity_All,E_DfIP,E_Sig_All);
  DCA_Sig_All -> SetName("DCA_Sig_All");
  TGraphErrors* DCA_Real_All = new TGraphErrors(n,DfIP_axis,RealEfficiency_All,E_DfIP,E_Real_All);
  DCA_Real_All -> SetName("DCA_Real_All");
  TGraphErrors* DCA_Rem_All = new TGraphErrors(n,DfIP_axis,RemovalnonSigMuon_All,E_DfIP,E_Rem_All);
  DCA_Rem_All -> SetName("DCA_Rem_All");
  TGraphErrors* DCA_nonSigMuonP_All = new TGraphErrors(n,DfIP_axis,PuritynonSigMuon_All,E_DfIP,E_nonSigMuonP_All);
  DCA_nonSigMuonP_All -> SetName("DCA_nonSigMuonP_All");
  TGraphErrors* DCA_SN_All = new TGraphErrors(n,DfIP_axis,SNRate_All,E_DfIP,E_SN_All);
  DCA_SN_All -> SetName("DCA_SN_All");
  TGraphErrors* DCA_SNW_All = new TGraphErrors(n,DfIP_axis,SNRateW_All,E_DfIP,E_SNW_All);
  DCA_SNW_All -> SetName("DCA_SNW_All");

  //MCMFT
  TGraphErrors* DCA_Rem_MCMFT = new TGraphErrors(n,DfIP_axis,RemovalnonSigMuon_MCMFT,E_DfIP,E_Rem_MCMFT);
  DCA_Rem_MCMFT -> SetName("DCA_Rem_MCMFT");
  TGraphErrors* DCA_SN_MCMFT = new TGraphErrors(n,DfIP_axis,SNRate_MCMFT,E_DfIP,E_SN_MCMFT);
  DCA_SN_MCMFT -> SetName("DCA_SN_MCMFT");
  TGraphErrors* DCA_SNW_MCMFT = new TGraphErrors(n,DfIP_axis,SNRateW_MCMFT,E_DfIP,E_SNW_MCMFT);
  DCA_SNW_MCMFT -> SetName("DCA_SNW_MCMFT");

  //Pt
  TGraphErrors* Eff_SigMuon = new TGraphErrors(n,Pt_axis,Efficiency_SigMuon,E_Pt,E_Eff_SigMuon);
  Eff_SigMuon -> SetName("Eff_SigMuon");
  //TGraphErrors* Eff_DCAcut = new TGraphErrors(n,Pt_axis,Efficiency_DCAcut,E_Pt,E_Eff_DCAcut);
  //Eff_DCAcut -> SetName("Eff_DCAcut");

  //Save histograms
  // MC Information
  output->WriteTObject(MCTrackPDG_dis_All.get());
  output->WriteTObject(MCTrackPDG_dis_Muon.get());
  output->WriteTObject(MCTrackPt_Muon_SigMuon.get());
  //MFT Information All
  output->WriteTObject(MFTTrackP_All.get());
  output->WriteTObject(MFTTrackPt_All.get());
  //output->WriteTObject(MFTTrackPt_DCAcut_021.get());
  output->WriteTObject(MFTTrackEta_All.get());
  output->WriteTObject(MFTTrackDfIP_All.get());
  output->WriteTObject(MFTTrackDfIP_Muon_All.get());
  output->WriteTObject(MFTTrackX_All.get());
  output->WriteTObject(MFTTrackY_All.get());
  output->WriteTObject(MFTTrackZ_All.get());
  output->WriteTObject(MFTTrackPDG_All.get());
  output->WriteTObject(MFTTrackMother_All.get());
  //MFT Information nonSigMuon
  output->WriteTObject(MFTTrackDfIP_nonSigMuon.get());
  output->WriteTObject(MFTTrackPt_nonSigMuon.get());
  output->WriteTObject(MFTTrackMCDfIP_nonSigMuon.get());
  //MFT Information SigMuon
  output->WriteTObject(MFTTrackDfIP_SigMuon.get());
  output->WriteTObject(MFTTrackPt_SigMuon.get());
  output->WriteTObject(MFTTrackMCDfIP_SigMuon.get());
  //All
  output->WriteTObject(MFTTrackPt_All.get());
  //Pion Keon
  output->WriteTObject(MFTTrackPt_nonSigMuon.get());
  output->WriteTObject(MFTTrackPt_Pi.get());
  output->WriteTObject(MFTTrackPt_K.get());
  //Heavy Flavor
  output->WriteTObject(MFTTrackPt_HF.get());
  output->WriteTObject(MFTTrackPt_D.get());
  output->WriteTObject(MFTTrackPt_B.get());
  //Low mass Vector meson
  output->WriteTObject(MFTTrackPt_LVM.get());
  output->WriteTObject(MFTTrackPt_Rou.get());
  output->WriteTObject(MFTTrackPt_Omega.get());
  output->WriteTObject(MFTTrackPt_Phi.get());

  output->WriteTObject(MFTTrackGrandMother_Rou.get());
  //Quarkonia
  output->WriteTObject(MFTTrackPt_Quarkonia.get());
  output->WriteTObject(MFTTrackPt_JPsi.get());
  output->WriteTObject(MFTTrackPt_Upsilon.get());
  //Direct
  output->WriteTObject(MFTTrackPt_Direct.get());
  //OtehrMuon
  output->WriteTObject(MFTTrackPt_OtherMuon.get());
  //NonMuon
  output->WriteTObject(MFTTrackPt_OtherNonMuon.get());
  //MC + MFT
  output->WriteTObject(MCMFTmatchDfIP_All.get());
  output->WriteTObject(MCMFTmatchDfIP_nonSigMuon.get());
  output->WriteTObject(MCMFTmatchDfIP_SigMuon.get());
//TGraph
  //All
  DCA_Puri_All->Write();
  DCA_Coll_All->Write();
  DCA_Sig_All->Write();
  DCA_Real_All->Write();
  DCA_Rem_All->Write();
  DCA_nonSigMuonP_All->Write();
  DCA_SN_All->Write();
  DCA_SNW_All->Write();
  //MCMFT
  DCA_Rem_MCMFT->Write();
  DCA_SN_MCMFT->Write();
  DCA_SNW_MCMFT->Write();
  //Pt
  Eff_SigMuon->Write();
  //Eff_DCAcut->Write();
}
