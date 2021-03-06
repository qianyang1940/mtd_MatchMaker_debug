//This *.C is to read mtdEVENT Tree, and do the analysis the deltaR and dletaY
//In which a small mtd relative Tree will also be stored
//Qian Yang
//2014_3_30
//1.0
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <sys/types.h>
#include <utility>
#include <getopt.h>
#include "math.h"
#include "string.h"

//#include "cut_corrz.h"
#include "cut.h"
#include "mtdEvent.h"
#include "mtdEventAnaly.h"
#ifndef __CINT__
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#endif
using namespace std;
using std::cout;
using std::endl;


bool passEvent(mtdEvent* event);
bool validTrack(StMtdTrack& track);
void clearTrack(StMtdTrack& track);
bool validHits(StMtdHits& hit);
void clearHits(StMtdHits& hit);
void destroyTrack(StMtdTrack& track);
void destroyHit(StMtdHits& hit);
void Hitsdist(mtdEvent* event,mtdCellHitVector& trackVec);
void bookHistograms();
bool loadCorrPars();
void bookmtdRelatTree(char* outFile);
bool TriggerRead(mtdEvent* event);
bool EventRead(mtdEvent* event,Int_t i);
void writeHistograms(char* outFile);
bool pairCosmic(mtdEvent* event,Int_t i);
void match(mtdEvent* event,mtdCellHitVector& matchHitCellsVec);
//*************** mtd correction parameters ***************
const Int_t maxnBins = 100;
const int BKLGNUM=30;
const int MODUNUM=30;
Float_t t0Corr[30][5][12];
Int_t   nBins[30][5];
Float_t tot[30][5][maxnBins];
Float_t tCorr[30][5][maxnBins];
Float_t zCenCorr[30][5];
Float_t tofOffsetCorr[30][5];
float TacdiffWestUp[BKLGNUM][MODUNUM];
float TacdiffWestbottom[BKLGNUM][MODUNUM];
float TacdiffEastUp[BKLGNUM][MODUNUM];
float TacdiffEastbottom[BKLGNUM][MODUNUM];
bool nNeighbors = kTRUE;
//*********************************************************

Int_t main(int argc, char **argv)
{
	if(argc!=3 && argc!=1)return 0;
        const Char_t *inFile="0.list";
	char *outFile = "huhu";

	TString *string;
        if(argc==1){
                string = new TString("0");
        }
        if(argc==3){
                inFile = argv[1];
		outFile = argv[2];
	}
	if(!loadCorrPars())return -1;
	if(argc!=1 && argc!=3)return -1;
	TChain *chain = new TChain("mtdEvent");

	Int_t ifile =0;
	char filename[512];
	ifstream *inputStream = new ifstream;
	inputStream->open(inFile);
	if(!(inputStream))
		{
			printf("can not open list file\n");
			return 0;
		}
	if(nNeighborsFlag==0)nNeighbors = kFALSE;
	else nNeighbors = kTRUE;
	for(;inputStream->good();)
		{
			inputStream->getline(filename,512);
			if(inputStream->good())
				{
					TFile *ftmp = new TFile(filename);
					if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys()))
						{
							cout<<"something is wrong"<<endl;
						}
					else
						{
							cout<<"read in"<<ifile<<"th file:"<<filename<<endl;
							chain->Add(filename);
							ifile++;
						}
					delete ftmp;
				}
		}
	delete inputStream;

	bookHistograms();
	bookmtdRelatTree(outFile);
	mtdEvent* eventtree = new mtdEvent(chain);
	Int_t entries = chain->GetEntries();
	cout<<entries<<"event"<<endl;
	Int_t entries_MB=0;
	for(Int_t i=0;i<entries;i++)
		{
		//	cout<<"step 0 begin ================================"<<endl;
			if(i%100==0)cout<<"begin"<<i<<"th entry..."<<endl;
			eventtree->GetEntry(i);
			if(!passEvent(eventtree))continue;
		//	cout<<"step 1 passEvent ended++++++++++++++++ is ok"<<endl;
			mtdCellHitVector trackVec;
			Hitsdist(eventtree,trackVec);
			mtdCellHitVector matchHitCellsVec;
			entries_MB++;
			match(eventtree,matchHitCellsVec);
		//	cout<<"step 2 match tracks ended-------------------------is ok"<<endl;
		}
		hHitEtavsPhi->Scale(1./entries_MB);
		hTrackEtavsPhi->Scale(1./entries_MB);
		hTrackEtavsPhi_1->Scale(1./entries_MB);
		hTrackEtavsPhi_Or->Scale(1./entries_MB);
		hMatchEtavsPhi->Scale(1./entries_MB);

		writeHistograms(outFile);
}
void match(mtdEvent *event,mtdCellHitVector& matchHitCellsVec)
{
	StructCellHit cellHit;
	int nMtdHits = event->nMtdHits;
	int ngTrks = event->ngTrks;	
	Double_t tStart = event->tStart;	
	int hits=0;
	int hits_cut=0;
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   matching the track with hits 
//   before matching t0 and swing effect are also need correct(marked as from gaygay)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	StMtdTrack track_loop;	
	StMtdHits hit;
	for(int i=0;i<nMtdHits;i++)
		{
			hit.bkleg = event->backleg[i];
			hit.module = event->module[i];
			hit.cell = event->cell[i];
			hit.leTimeWest = event->leTimeWest[i];
			hit.leTimeEast = event->leTimeEast[i];
			hit.totWest = event->totWest[i];
			hit.totEast = event->totEast[i];
			hit.TdiffWest = event->TdiffWest[i];
			hit.TdiffEast = event->TdiffEast[i];
			hits_cut++;

			char backleg = event->backleg[i];
			char module = event->module[i];
			char cell = event->cell[i];
			double TdiffWest = event->TdiffWest[i];
			double TdiffEast = event->TdiffEast[i];
			double zdaq = (hit.leTimeEast -hit.leTimeWest)/2./mvDrift*1e3;
		 	hits++;	
			float leTimediff = (TdiffWest-TdiffEast)/2.;
			if(!validHits(hit))
			{	
				clearHits(hit);
				continue;
			}

//			++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++from gaygay
	//T0 correction
//	Float_t tof2MTD = (leTimeWest+leTimeEast)/2.-tStart-t0Corr[backleg-1][module-1][cell-1];
			Float_t tof2MTD = (hit.leTimeWest+hit.leTimeEast)/2.-tStart;
	
	//Slewing correction
	/*
	Int_t totBin = -1;
	Float_t muTot = sqrt(totWest*totEast);
	for(Int_t j=0; j<nBins[backleg-1][module-1]; j++){
		if(muTot>=tot[backleg-1][module-1][j]
				&& muTot<tot[backleg-1][module-1][j+1]){
			totBin = j;
			break;
		}
	}
	Float_t slewCorr = 0;
	if(totBin>=0 && totBin<nBins[backleg-1][module-1]){
		Float_t x1 = tot[backleg-1][module-1][totBin];
		Float_t x2 = tot[backleg-1][module-1][totBin+1];
		Float_t y1 = tCorr[backleg-1][module-1][totBin]; 
		Float_t y2 = tCorr[backleg-1][module-1][totBin+1]; 
		if(x2-x1!=0) slewCorr = y1+(muTot-x1)*(y2-y1)/(x2-x1);
	}
	tof2MTD -= slewCorr;
	if(tof2MTD>25) tof2MTD -=25; //related tof issue during run12 calibration.
	*/
//	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++end from gaygay
			for(Int_t j=0;j<ngTrks;j++)
				{
					Int_t idVec = event->gTrkMatchNum[j];
					track_loop.mtdpt = event->gpt[j];
					track_loop.mtdeta = event->geta[j];
					track_loop.mtdphi = event->gphi[j];
					track_loop.mtdnFtPts = event->gnFtPts[j];
					track_loop.mtdnDedxPts = event->gnDedxPts[j];
					track_loop.mtddca = event->gdca[j];
					track_loop.mtdnSigmaPi = event->gnSigmaPi[j];
					for(Int_t k=0;k<idVec;k++)
					{
						track_loop.mtdBL.push_back(event->gprojMtdBackLeg[j][k]);
						track_loop.mtdmod.push_back(event->gprojMtdModule[j][k]);
						track_loop.mtdcell.push_back(event->gprojMtdCell[j][k]);
						track_loop.mtdProjphi.push_back(event->gprojMtdPhi[j][k]);
						track_loop.mtdProjZ.push_back(event->gprojMtdZ[j][k]);
						track_loop.mtdProjLength.push_back(event->gprojMtdLength[j][k]);
						track_loop.mtdtof2Mtd.push_back(event->gtof2Mtd[j][k]);
					}
					if(!(validTrack(track_loop)))
					{	
						clearTrack(track_loop);
						continue;
					}
					float pt = event->gpt[j];
					float eta = event->geta[j];
					float phi = event->gphi[j];
					float gnSigmaPi = event->gnSigmaPi[j];
					int gtrackindex = event->gtrackindex[j];
					for(Int_t k=0;k<idVec;k++)
					{
						int gprojMtdBackLeg = event->gprojMtdBackLeg[j][k];
						int gprojMtdModule = event->gprojMtdModule[j][k];
						float projMtdPhi = event->gprojMtdPhi[j][k];
						float projMtdZ = event->gprojMtdZ[j][k];
						double gtof2Mtd = event->gtof2Mtd[j][k];
						bool isMatch = kFALSE;	
						bool isMatchNebi = kFALSE;
						
						if((backleg==gprojMtdBackLeg)&&(module==gprojMtdModule))
						{
							isMatch = kTRUE;
						}
						if(nNeighbors)
						{
							if(backleg==gprojMtdBackLeg && module!=gprojMtdModule){
								if(zdaq<0){
									if((module-1>0)&&((module-1)==gprojMtdModule)) isMatchNebi = kTRUE;
								}else if(zdaq==0){
									if(abs(module-gprojMtdModule)<=1) isMatchNebi = kTRUE;
								}else{
									if((module+1<6)&&((module+1)==gprojMtdModule)) isMatchNebi = kTRUE;
								}
							}
						}
						/////////////////////////////////////////////////////
							if(isMatch||isMatchNebi)
							{
								//	cout<<"+++++++++++++++++++++++++++++++++++++++++"<<endl;
								//	cout<<"isMatch="<<isMatch<<endl;
								//	cout<<"isMatchNebi="<<isMatchNebi<<endl;
								//	cout<<"-----------------------------------------"<<endl;
								cellHit.pt = pt;
								cellHit.tof2MTD = tof2MTD;
								cellHit.leTimediff = leTimediff;
								cellHit.gtof2Mtd = gtof2Mtd;
								cellHit.nSigmaPi = gnSigmaPi;
								cellHit.backleg = backleg;
								cellHit.module = module;
								cellHit.cell = cell;
								cellHit.eta = eta;
								cellHit.projMtdZ = projMtdZ;
								cellHit.projMtdPhi = projMtdPhi;
								cellHit.phi = phi;
								cellHit.trackIdVec.push_back(gtrackindex);
								cellHit.MatchC=isMatch;
								cellHit.MatchN=isMatchNebi;
								matchHitCellsVec.push_back(cellHit);
							}	
					}//endl of the idvec					
				}//end of the track loop
		}//endl of the mathcing

	mtdCellHitVector tempVec = matchHitCellsVec;
	mtdCellHitVector erasedVec = tempVec;
//	cout<<"matchedhit="<<tempVec.size()<<endl;
//	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
// loop all matching, find out the matching: one hits for one track(may one track two hits); one hit two tracks; 
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	Int_t mFlag = 0;
	//Int_t Numberofhits=0;
	 while (tempVec.size() != 0) 
		{
			Int_t nTracks = 0;
			Int_t nTracks_sort = 0;
			idVector trackIdVec;
			const int arrynumber = 20;
			Int_t bckleg[arrynumber] ={0};
			Float_t pt[arrynumber]={0};
			Int_t modl[arrynumber] ={0};
			Int_t cells[arrynumber]={0};
			Float_t projmtdz[arrynumber] ={0.};
			Float_t leTimediff[arrynumber] ={0.};
			Float_t projmtdphi[arrynumber] ={0.};
			Float_t eta[arrynumber] ={0.};
			Float_t phi[arrynumber] ={0.};
			Int_t inode[arrynumber] ={0};
			Float_t tof2MTD[arrynumber] = {0};
			Float_t gtof2Mtd[arrynumber] ={0};
			bool MatchCflag[arrynumber]={kFALSE};
			bool MatchNflag[arrynumber]={kFALSE};

			mtdCellHitVectorIter tempIter=tempVec.begin();
			mtdCellHitVectorIter erasedIter=erasedVec.begin();
			while(erasedIter!= erasedVec.end()) {
				if(tempIter->backleg == erasedIter->backleg &&
						tempIter->module == erasedIter->module &&
						tempIter->cell == erasedIter->cell) {

					bckleg[nTracks]=erasedIter->backleg;
					modl[nTracks]=erasedIter->module;
					cells[nTracks]=erasedIter->cell;
					pt[nTracks]=erasedIter->pt;
					projmtdz[nTracks]=erasedIter->projMtdZ;
					projmtdphi[nTracks]=erasedIter->projMtdPhi;
					eta[nTracks]=erasedIter->eta;
					phi[nTracks]=erasedIter->phi;
					inode[nTracks]=erasedIter->trackIdVec.back();
					tof2MTD[nTracks]=erasedIter->tof2MTD;
					leTimediff[nTracks]=erasedIter->leTimediff;
					gtof2Mtd[nTracks]=erasedIter->gtof2Mtd;
					MatchCflag[nTracks]=erasedIter->MatchC;
					MatchNflag[nTracks]=erasedIter->MatchN;
					nTracks++;
					trackIdVec.push_back(erasedIter->trackIdVec.back());  // merge
					erasedVec.erase(erasedIter);
					erasedIter--;
				}
				erasedIter++;
			}	
			
			idVector tmpIdVec = trackIdVec;
			sort(tmpIdVec.begin(),tmpIdVec.end());
			tmpIdVec.erase(unique(tmpIdVec.begin(),tmpIdVec.end()),tmpIdVec.end());
			nTracks_sort = tmpIdVec.size();
			if(nTracks>0)
			{		
				hTracksPerHit->Fill(nTracks);
				hTracksPerHit_sort->Fill(nTracks_sort);
			}
			if(nTracks_sort==1)
			{
				double hitZ = (modl[0]-3)*87.-leTimediff[0]/mvDrift*1e3;
				double projZ=projmtdz[0];
				double deltaZ = hitZ-projZ;
				int index =999;
				index = bckleg[0]*5+modl[0];
				// deltaZ = deltaZ+corrdeltaz[index];
				double deltaR = deltaZ;
				hDeltazvsStrip->Fill(cells[0]+12*(modl[0]-1)+60*(bckleg[0]-1),deltaR);	 
				double pT = pt[0];
				hdeltaR_singlevsPt->Fill(deltaR,pT);
				hTimediffvsPt->Fill(tempIter->tof2MTD-tempIter->gtof2Mtd,pT);
			//	Int_t moduleNum = 0;
			//	if((tempIter->tof2MTD-tempIter->gtof2Mtd)>9||(tempIter->tof2MTD-tempIter->gtof2Mtd)<11)
			//	{
			//		moduleNum= (tempIter->backleg)*5+tempIter->module;
			//		hModule->Fill(moduleNum);
			//	}	
				if(nTracks==1)
				{	
				//	cout<<"matchCflag[0]="<<MatchCflag[0]<<endl;
					if(MatchCflag[0]==kTRUE)
					{
						hMatchSource->Fill(2);
						hDeltazvsStrip_ISMATCH->Fill(cells[0]+12*(modl[0]-1)+60*(bckleg[0]-1),deltaR);	 
							
					}
					else
					{
						hMatchSource->Fill(3);
						hDeltazvsStrip_ISMATCHNebi->Fill(cells[0]+12*(modl[0]-1)+60*(bckleg[0]-1),deltaR);	 
				//		cout<<"is nebi::::matchCflag[0]:matchNflag[0]="<<MatchCflag[0]<<" "<<MatchNflag[0]<<endl;
					}
				}
				if(nTracks==2)
				{
					for(int kk=0;kk<nTracks;kk++)
					{
						double hitZZ = (modl[kk]-3)*87.-leTimediff[kk]/mvDrift*1e3;
						double projZZ=projmtdz[kk];
						double deltaZZ = hitZZ-projZZ;
						if(MatchCflag[kk]==kTRUE)
						{
							hMatchSource->Fill(4);
							hDeltazvsStrip_ISMATCH_TWO->Fill(cells[kk]+12*(modl[kk]-1)+60*(bckleg[kk]-1),deltaZZ);	 
						}
						else
						{
							hMatchSource->Fill(5);
							hDeltazvsStrip_ISMATCHNebi_TWO->Fill(cells[kk]+12*(modl[kk]-1)+60*(bckleg[kk]-1),deltaZZ);	 
						}
					}// track double hit loop end
				}
			}

			if(nTracks_sort==2) 
			{
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				if(nTracks==2)
				{
					mFlag = 1;
					Double_t hitZ = (modl[0]-3)*87.-leTimediff[0]/mvDrift*1e3;
					//	Double_t projZ=((((projmtdz[0]+2.5*87.)/87.+1)-(int)((projmtdz[0]+2.5*87.)/87.+1))-0.5)*87.;
					Double_t projZ = projmtdz[0];
					Double_t deltaZ = hitZ-projZ;
					int index=999;
					index = bckleg[0]*5+modl[0];
					//	deltaZ = deltaZ + corrdeltaz[index];
					double deltaR1 = deltaZ;
					hitZ = (modl[1]-3)*87.-leTimediff[1]/mvDrift*1e3;
					projZ = projmtdz[1];
					deltaZ = hitZ-projZ;
					index =999;
					index = bckleg[1]*5+modl[1];
					//	deltaZ = deltaZ + corrdeltaz[index];
					double deltaR2 = deltaZ;
					float tdiff1 = TMath::Abs(tof2MTD[0]-gtof2Mtd[0]+0.8751);
					float tdiff2 = TMath::Abs(tof2MTD[1]-gtof2Mtd[1]+0.8751);
					float tdiff3 = tof2MTD[0]-gtof2Mtd[0];
					float tdiff4 = tof2MTD[1]-gtof2Mtd[1];

					if(MatchCflag[0]==kTRUE&&MatchCflag[1]==kTRUE)
					{
						hMatchSource->Fill(7);
					}//two track in the same hit module
					if((MatchCflag[0]+MatchCflag[1])==1)
					{
						hMatchSource->Fill(8);
						if(MatchCflag[0]==kTRUE)
						{

							hTimediff_Rcut_l->Fill(tdiff3,pt[0]);
							hTimediff_Rcut_s->Fill(tdiff4,pt[1]);
							hdeltaR_svsPt->Fill(deltaR2,pt[1]);
							hdeltaR_lvsPt->Fill(deltaR1,pt[0]);
						}
						else
						{
							hdeltaR_lvsPt->Fill(deltaR2,pt[1]);
							hdeltaR_svsPt->Fill(deltaR1,pt[0]);
							hTimediff_Rcut_l->Fill(tdiff4,pt[1]);
							hTimediff_Rcut_s->Fill(tdiff3,pt[0]);


						}
					}//one track in and an anther adject track 
					if((MatchCflag[0]+MatchCflag[1])==0)
					{
						hMatchSource->Fill(9);
						if(abs(deltaR1)<abs(deltaR2))
						{
							hdeltaR_svsPt_Two->Fill(deltaR1,pt[0]);
							hdeltaR_lvsPt_Two->Fill(deltaR2,pt[1]);
						}
						else
						{
							hdeltaR_lvsPt_Two->Fill(deltaR1,pt[0]);
							hdeltaR_svsPt_Two->Fill(deltaR2,pt[1]);
						}
					}//two adjec track
				}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
			}
			Int_t FiredStripNumber = (tempIter->cell)+(tempIter->module-1)*12+(tempIter->backleg-1)*60;//test Qian Yang
			double dphi = 12./180.*TMath::Pi();
			Float_t projMtdZ1 = tempIter->projMtdZ;
			Float_t projMtdPhi1 = tempIter->projMtdPhi;
			//Int_t projmodule1 = (int)((projMtdZ1+2.5*87.)/87.+1);
			int projMtdBackLeg1 = -1;
			projMtdBackLeg1 = (int)(projMtdPhi1/dphi);
			projMtdBackLeg1 += 24;
			if(projMtdBackLeg1>30) projMtdBackLeg1 -= 30;
			if(nTracks>0)
			{
				hMatchEtavsPhi->Fill(tempIter->module,tempIter->backleg);
				hNumberofMatchHits->Fill(FiredStripNumber,nTracks_sort);
			}

			tempVec = erasedVec;
		//	cout<<"------------------------------------------------------end"<<endl;
		}
	cout<<"step 1.9"<<endl;
}
bool loadCorrPars(){

	ifstream indata;
	indata.open("/star/u/syang/data01/mtd/run12/cuau/pid/pid_Hist/calibDatabase/t0_4DB.dat");
	if(!indata.is_open()){
		cout<<"Failed to open the T0 correction file !!!"<<endl;
		return kFALSE;
	}
	cout<<"************* Load T0 correction parameters start ************"<<endl;
	memset(t0Corr,0,sizeof(t0Corr));
	Int_t ib,it,ic;
	while(indata>>ib>>it>>ic){
		indata>>t0Corr[ib-1][it-1][ic-1];
		cout<<"BL:"<<ib<<"  \tTray:"<<it<<"  \tCell:"<<ic<<"  \tT0corr:"<<t0Corr[ib-1][it-1][ic-1]<<endl;
	};
	indata.close();
	cout<<"************* Load T0 correction parameters done ************"<<endl;
	cout<<endl;
	cout<<endl;


	indata.open("/star/u/syang/data01/mtd/run12/cuau/pid/pid_Hist/calibDatabase/totCali_4DB.dat");
	cout<<"************* Load Slewing correction parameters start ************"<<endl;
	if(!indata.is_open()){
		cout<<"Failed to open the Slewing correction file !!!"<<endl;
		return kFALSE;
	}
	memset(nBins,0,sizeof(nBins));
	memset(tot,0,sizeof(tot));
	memset(tCorr,0,sizeof(tCorr));
	while(indata>>ib>>it){
		indata>>nBins[ib-1][it-1];
		for(Int_t j=0;j<=nBins[ib-1][it-1];j++){
			indata>>tot[ib-1][it-1][j];
		}
		for(Int_t j=0;j<=nBins[ib-1][it-1];j++){
			indata>>tCorr[ib-1][it-1][j];
		}
		cout<<"BL:"<<ib<<"  \tTray:"<<it<<"  \ttotBins:"<<nBins[ib-1][it-1]<<endl;
		cout<<"ToT Bins:  \t";
		for(Int_t j=0;j<=nBins[ib-1][it-1];j++){
			cout<<tot[ib-1][it-1][j]<<"  \t";
		}
		cout<<endl;
		cout<<"tCorr Bins:  \t";
		for(Int_t j=0;j<=nBins[ib-1][it-1];j++){
			cout<<tCorr[ib-1][it-1][j]<<"  \t";
		}
		cout<<endl;
	}
	cout<<"************* Load Slewing correction parameters done ************"<<endl;
	cout<<endl;
	cout<<endl;
	indata.close();


	indata.open("/star/u/syang/data01/mtd/run12/cuau/pid/pid_Hist/calibDatabase/zCenter.dat");
	cout<<"************* Load deltaZ center correction parameters start ************"<<endl;
	if(!indata.is_open()){
		cout<<"Failed to open the deltaZ center correction file !!!"<<endl;
		return kFALSE;
	}
	memset(zCenCorr,0,sizeof(zCenCorr));
	while(indata>>ib>>it){
		indata>>zCenCorr[ib-1][it-1];
		cout<<"BL:"<<ib<<"  \tTray:"<<it<<"  \tzCenCorr:"<<zCenCorr[ib-1][it-1]<<endl;
	};
	indata.close();
	cout<<"************* Load deltaZ center correction parameters done ************"<<endl;
	cout<<endl;
	cout<<endl;


	indata.open("/star/u/syang/data01/mtd/run12/cuau/pid/pid_Hist/calibDatabase/tofOffset.dat");
	cout<<"************* Load tofOffset correction parameters start ************"<<endl;
	if(!indata.is_open()){
		cout<<"Failed to open the tofOffset correction file !!!"<<endl;
		return kFALSE;
	}
	memset(tofOffsetCorr,0,sizeof(tofOffsetCorr));
	while(indata>>ib>>it){
		indata>>tofOffsetCorr[ib-1][it-1];
		cout<<"BL:"<<ib<<"  \tTray:"<<it<<"  \ttofOffsetCorr:"<<tofOffsetCorr[ib-1][it-1]<<endl;
	};
	indata.close();
	cout<<"************* Load tofOffset correction parameters done ************"<<endl;
	cout<<endl;
	cout<<endl;

	indata.open("/star/u/tc88qy/data01/mtd/out/ana_1/Timediff/timediff_cp.dat");
	cout<<"************* Load tofOffset correction parameters start ************"<<endl;
	if(!indata.is_open()){
		cout<<"Failed to open the tofOffset correction file !!!"<<endl;
		return kFALSE;
	}
	memset(TacdiffWestUp,0,sizeof(TacdiffWestUp));
	memset(TacdiffWestbottom,0,sizeof(TacdiffWestbottom));
	memset(TacdiffEastUp,0,sizeof(TacdiffEastUp));
	memset(TacdiffEastbottom,0,sizeof(TacdiffEastbottom));
	while(indata>>ib>>it){
		float mean;
		float sigma;
		indata>>mean;
		indata>>sigma;
		TacdiffWestUp[ib-1][it-1]=mean+3*sigma;
		TacdiffWestbottom[ib-1][it-1]=mean-3*sigma;
		cout<<"WestUp:Westbottom="<<mean+3*sigma<<":"<<mean-3*sigma<<endl;
		indata>>mean;
		indata>>sigma;
		TacdiffEastUp[ib-1][it-1]=mean+3*sigma;
		TacdiffEastbottom[ib-1][it-1]=mean-3*sigma;
		cout<<"EastUp:Eastbottom="<<mean+3*sigma<<":"<<mean-3*sigma<<endl;
	};
	indata.close();
	cout<<"************* Load tofOffset correction parameters done ************"<<endl;
	cout<<endl;
	cout<<endl;
	
	cout<<"Load Parameters DONE !!!"<<endl;

	return kTRUE;
};
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool passEvent(mtdEvent* event )
{

	Int_t nTrgIds = event->nTrgIds;
	Int_t trgId[10];
	for(Int_t i=0;i<nTrgIds;i++)
		{
			trgId[i] = event->trgId[i];
		}
	Float_t vx = event->vertexX;
	Float_t vy = event->vertexY;
	Float_t vz = event->vertexZ;
	Float_t vzVpd = event->vpdVz;
	Float_t vzDiff = vz - vzVpd;
		
	bool TriggerFlag  = kFALSE;
	for(Int_t i =0 ; i<nTrgIds;i++)
		{
			if(trgId[i]==450635||trgId[i]==450645)
				{		
					TriggerFlag = kTRUE;
				}
		}
	if(TMath::Abs(vx)<1e-5&&TMath::Abs(vy)<1e-5&&TMath::Abs(vz)<1e-5)return kFALSE;

	hVertexXY->Fill(vx,vy);
	hVzDiff->Fill(vzDiff);
	hvzVpdvsvz->Fill(vzVpd,vz);
	Int_t refMult = event->refMult;
	hRefmult->Fill(refMult);
	if(!(TriggerFlag))return kFALSE;
	if(vz==-999)return kFALSE;

	if(vzDiff<-3||vzDiff>3)return kFALSE;


	return kTRUE;	
	
}
void Hitsdist(mtdEvent* event,mtdCellHitVector& trackVec)
{

	//cout<<"step 1.0"<<endl;
	StructCellHit track;
	StMtdTrack track_loop;	
	StMtdHits hit;
	MtdHits hitsVec; 
	MtdTrack trkVec; 
	Int_t nMtdHits = event->nMtdHits;
	Int_t ngTrks = event->ngTrks;	
	Int_t refMult = event->refMult;	
	Int_t hits=0;
	Int_t hits_cut=0;
	Int_t _bkleg[200]={0};
	Int_t _module[200]={0};
	Int_t _cell[200]={0};
	for(Int_t i=0;i<nMtdHits;i++)
		{
			hit.bkleg = event->backleg[i];
			hit.module = event->module[i];
			hit.cell = event->cell[i];
			hit.leTimeWest = event->leTimeWest[i];
			hit.leTimeEast = event->leTimeEast[i];
			hit.totWest = event->totWest[i];
			hit.totEast = event->totEast[i];
			hit.TdiffWest = event->TdiffWest[i];
			hit.TdiffEast = event->TdiffEast[i];
		 	hits++;	
			char backleg = event->backleg[i];
			char module = event->module[i];
			char cell = event->cell[i];
			double TdiffWest = event->TdiffWest[i];
			double TdiffEast = event->TdiffEast[i];
			double zdaq = (hit.leTimeWest -hit.leTimeEast)/2./mvDrift*1e3;
			int Strip;
			if(module<3){
				cell = cell+1;
				Strip  = (backleg-1)*60+(module-1)*12+cell;
			}
			else{
				Strip = (backleg-1)*60+(module-1)*12+(12-cell);
			};
			Float_t leTimediff = (TdiffWest-TdiffEast)/2.;
			hTdiffWestvsStrip->Fill(Strip,TdiffWest);
			hTdiffEastvsStrip->Fill(Strip,TdiffEast);
			if(!validHits(hit))
			{
				clearHits(hit);
				continue;
			}
			hTdiffWestvsStrip_Tac->Fill(Strip,TdiffWest);
			hTdiffEastvsStrip_Tac->Fill(Strip,TdiffEast);
			hHitEtavsPhi->Fill(module,backleg);
		
			hitsVec.push_back(hit);

			_bkleg[hits_cut]=hit.bkleg;
			_module[hits_cut]=hit.module;
			_cell[hits_cut]=hit.cell;
			clearHits(hit);
			hits_cut++;
		}
    hHits->Fill(hits);
	hHits_cut->Fill(hits_cut);
	//cout<<"step 2.0"<<endl;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//   for those events have two mtd hits 
	//   filling histo to get the hits pros information
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	for(int i=0;i<hits_cut;i++)
	{
		for(int j =i+1;j<hits_cut;j++)
		{
			int hitsdiff = 7;
			//int stripdiff = 0;
			int strip_loop1 =999;
			int strip_loop2 =999;
			if(_module[i]<3)
			{
				strip_loop1 =(_bkleg[i]-1)*60+(_module[i]-1)*12+(_cell[i]+1);
			}
			else
			{
				strip_loop1 =(_bkleg[i]-1)*60+(_module[i]-1)*12+(12-_cell[i]);
			};
			if(_module[j]<3)
			{
				strip_loop2 =(_bkleg[j]-1)*60+(_module[j]-1)*12+(_cell[j]+1);
			}
			else
			{
				strip_loop2 =(_bkleg[j]-1)*60+(_module[j]-1)*12+(12-_cell[j]);
			};
			if(_bkleg[i]==_bkleg[j])
			{
				hitsdiff= 2;
				if(TMath::Abs(_module[i]-_module[j])==1)
				{
					hitsdiff=1;
				}
				if(_module[i]==_module[j])
				{
					hitsdiff=0;
				}
			}
			if(TMath::Abs(_bkleg[i]-_bkleg[j])==1) 
			{
				hitsdiff =4;
				if(_module[i]==_module[j])hitsdiff=3;
			}
			if(TMath::Abs(_bkleg[i]-_bkleg[j])==2)
			{
				hitsdiff=5;
			}
			hHitDis->Fill(hitsdiff,strip_loop1-strip_loop2);
		}
	}
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//   loop all tracks and sort in the vector
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//cout<<"step 3.0"<<endl;
	Int_t numberoftrack=0;
//	cout<<"ngTrks="<<ngTrks<<endl;
	for(Int_t j=0;j<ngTrks;j++)
	{
		Int_t idVec = event->gTrkMatchNum[j];

		track_loop.mtdpt = event->gpt[j];
		track_loop.mtdeta = event->geta[j];
		track_loop.mtdphi = event->gphi[j];
		track_loop.mtdnFtPts = event->gnFtPts[j];
		track_loop.mtdnDedxPts = event->gnDedxPts[j];
		track_loop.mtddca = event->gdca[j];
		track_loop.mtdnSigmaPi = event->gnSigmaPi[j];
		track_loop.trackId = event->gtrackindex[j];
		hTrackEtavsPhi_Or->Fill(event->geta[j],event->gphi[j]);
		for(int k=0;k<idVec;k++)
		{
			track_loop.mtdBL.push_back(event->gprojMtdBackLeg[j][k]);
			track_loop.mtdmod.push_back(event->gprojMtdModule[j][k]);
			track_loop.mtdcell.push_back(event->gprojMtdCell[j][k]);
			track_loop.mtdProjphi.push_back(event->gprojMtdPhi[j][k]);
			track_loop.mtdProjZ.push_back(event->gprojMtdZ[j][k]);
			track_loop.mtdProjLength.push_back(event->gprojMtdLength[j][k]);
			track_loop.mtdtof2Mtd.push_back(event->gtof2Mtd[j][k]);
		}
		if(!(validTrack(track_loop)))
		{
			clearTrack(track_loop);
			continue;
		}
//		cout<<"idVec="<<idVec<<endl;
		trkVec.push_back(track_loop);

		numberoftrack++;
		hTracksMatchNumber->Fill(idVec);
		hTrackEtavsPhi_1->Fill(event->geta[j],event->gphi[j]);

		float projtemz = -999;
		track.pt = track_loop.mtdpt;
		track.eta = track_loop.mtdeta;
		for(int k=0;k<idVec;k++)
		{
			track.backleg = event->gprojMtdBackLeg[j][k];
			track.module = event->gprojMtdModule[j][k];
			track.cell = event->gprojMtdCell[j][k];
			track.projMtdZ = event->gprojMtdZ[j][k];
			track.projMtdPhi = event->gprojMtdPhi[j][k];
			track.trackIdVec.push_back((int)event->gtrackindex[j]);
			trackVec.push_back(track);
			float z = event->gprojMtdZ[j][k];
			int cell = event->gprojMtdCell[j][k];
			hTrackprojZ->Fill(z);
			hTrackprojCell->Fill(cell);
			if(idVec==2)
			{
				if(k==0)projtemz=event->gprojMtdZ[j][k];
				if(k==1)
				{
					hTrackprojZdiffvsEta->Fill(event->geta[j],event->gprojMtdZ[j][k]-projtemz);
					hTrackprojZdiffvsPt->Fill(event->gpt[j],event->gprojMtdZ[j][k]-projtemz);
				}
			}
			hTrackEtavsPhi->Fill(event->gprojMtdModule[j][k],event->gprojMtdBackLeg[j][k]);
		}
		clearTrack(track_loop);
	}

	hTracks->Fill(numberoftrack);
	htracksvsHits->Fill(numberoftrack,hits_cut);
	//cout<<"step 1.2"<<endl;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//   loop all track and hit get the distribution of the hits vs track. 
	//   one track one hits(in three module) two tracks and ont hits ........
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	MtdHits tempHitsVec = 	hitsVec;
	MtdHits erasetempHitsVec = tempHitsVec;
	MtdTrack temptrkVec = trkVec;	
	MtdTrack erasetemptrkVec = temptrkVec;
	MtdHits SingalhitsVec; 
	StMtdHits hits_temp;
	while(tempHitsVec.size()!=0)
	{
		MtdHitsVectorIter tempHitsVecIter = tempHitsVec.begin();

		int hitsbkgIter = tempHitsVecIter->bkleg;
		int hitsmoduleIter = tempHitsVecIter->module;
		int hitscellIter = tempHitsVecIter->cell;
		bool HitsFlag = kFALSE;
		bool HitsNeboileftFlag = kFALSE;
		bool HitsNeboirightFlag = kFALSE;
		bool HitsNeboiFlag = kFALSE;
		int hitsflag = 0;
		int hitsleftflag = 0;
		int hitsrightflag = 0;
		MtdHitsVectorIter erasetempHitsVecIter = erasetempHitsVec.begin();
		while(erasetempHitsVecIter!=erasetempHitsVec.end())
		{
			int erasehitsbkgIter = erasetempHitsVecIter->bkleg;
			int erasehitsmoduleIter = erasetempHitsVecIter->module;
			int erasehitscellIter = erasetempHitsVecIter->cell;
			if(hitsbkgIter==erasehitsbkgIter)
			{
				if(hitsmoduleIter==erasehitsmoduleIter)
				{
					hitsflag++;
					if(hitsflag==1)
					{
						erasetempHitsVec.erase(erasetempHitsVecIter);
						erasetempHitsVecIter--;
					}
				}
				if(hitsmoduleIter==(erasehitsmoduleIter-1))
				{
					hitsleftflag++;
				}
				if(hitsmoduleIter==(erasehitsmoduleIter+1))
				{
					hitsrightflag++;
				}
			}
			erasetempHitsVecIter++;
		}
	    hHitsMulti->Fill(2,hitsflag);
	    hHitsMulti->Fill(3,hitsleftflag);
	    hHitsMulti->Fill(4,hitsrightflag);
		if(hitsflag==1)
		{
			if((hitsrightflag+hitsleftflag)==0)HitsFlag=kTRUE;
			if(hitsrightflag>0&&hitsleftflag==0)HitsNeboirightFlag=kTRUE;
			if(hitsrightflag==0&&hitsleftflag>0)HitsNeboileftFlag=kTRUE;
			if(hitsrightflag>0&&hitsleftflag>0)HitsNeboiFlag=kTRUE;

			hits_temp.bkleg = hitsbkgIter;
			hits_temp.module = hitsmoduleIter;
			hits_temp.cell = hitsmoduleIter;
			hits_temp.leTimeWest = tempHitsVecIter->leTimeWest;
			hits_temp.leTimeEast = tempHitsVecIter->leTimeEast;
			hits_temp.totWest = tempHitsVecIter->totWest;
			hits_temp.totEast = tempHitsVecIter->totEast;
			hits_temp.TdiffWest = tempHitsVecIter->TdiffWest;
			hits_temp.TdiffEast = tempHitsVecIter->TdiffEast;
			hits_temp.HitsFlag = HitsFlag;
			hits_temp.HitsNeboirightFlag = HitsNeboirightFlag;
			hits_temp.HitsNeboileftFlag = HitsNeboileftFlag;
			hits_temp.HitsNeboiFlag = HitsNeboiFlag;
			SingalhitsVec.push_back(hits_temp);
		}
		tempHitsVec = erasetempHitsVec;
	}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//			add track neib flag
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		MtdTrack temptrkVec_temp = trkVec;
		MtdTrack erasetemptrkVec_temp = temptrkVec_temp;
		MtdTrack SingaltrackVec; 
		StMtdTrack track_temp;
		while(temptrkVec_temp.size()!=0)
		{
			MtdTrackVectorIter temptrkVecIter = temptrkVec_temp.begin();
			idVector modVect = temptrkVecIter->mtdmod;
			//cout<<"projtrk hits in moudle :"<<(temptrkVecIter->mtdmod).size()<<endl;	
			if(modVect.size()==0)continue;
			MtdTrackVectorIter erasetemptrkVecIter = erasetemptrkVec_temp.begin();
			bool trackFlag = kFALSE;
			bool trackNeboileftFlag = kFALSE;
			bool trackNeboirightFlag = kFALSE;
			bool trackNeboiFlag = kFALSE;
			int projbklgIter = (int)((temptrkVecIter->mtdBL).back());
			int projmoduleIter = (int)((temptrkVecIter->mtdmod).back());
			int projcellIter = (int)((temptrkVecIter->mtdcell).back());
			int protrackId = (int)(temptrkVecIter->trackId);
			float projphi = (float)((temptrkVecIter->mtdProjphi).back());
			float projz = (float)((temptrkVecIter->mtdProjZ).back());
			float projlength = (float)((temptrkVecIter->mtdProjLength).back());
			float tof2Mtd = (float)((temptrkVecIter->mtdtof2Mtd).back());
		//	cout<<"proj:bklegIter:moduleIter:cellIter="<<projbklgIter<<":"<<projmoduleIter<<":"<<projcellIter<<endl;
			idVector trackIdC_sort;
			idVector trackIdL_sort;
			idVector trackIdR_sort;
			int trackflag=0;
		   	int trackleftflag=0;
			int trackrightflag=0;	
			while(erasetemptrkVecIter!=erasetemptrkVec_temp.end())
			{
				int erabklgIter = (int)((erasetemptrkVecIter->mtdBL).back());
				int eramodulegIter = (int)((erasetemptrkVecIter->mtdmod).back());
				int eracellIter = (int)((erasetemptrkVecIter->mtdcell).back());
				int eratrackId = (int)(erasetemptrkVecIter->trackId);
		//	cout<<"era:bklegIter:moduleIter:cellIter="<<erabklgIter<<":"<<eramodulegIter<<":"<<eracellIter<<endl;
				if(projbklgIter==erabklgIter)
				{
					if(projmoduleIter==eramodulegIter)
					{
						trackflag++;
						trackIdC_sort.push_back(eratrackId);
						if(trackflag==1)
						{
							erasetemptrkVec_temp.erase(erasetemptrkVecIter);
							erasetemptrkVecIter--;
						}
					}
					if(projmoduleIter==(eramodulegIter-1))
					{
						trackleftflag++;
						trackIdL_sort.push_back(eratrackId);
						if(trackleftflag==1)
						{
							trackIdL_sort.push_back(protrackId);
						}
					}
					if(projmoduleIter==(eramodulegIter+1))
					{
						trackrightflag++;
						trackIdL_sort.push_back(eratrackId);
						if(trackrightflag==1)
						{
							trackIdR_sort.push_back(protrackId);
						}
					}

				}
				erasetemptrkVecIter++;
			}//endl of the erase loop	
	//		cout<<"trackIdC_sort.size()="<<trackIdC_sort.size()<<endl;
	//		cout<<"trackIdR_sort.size()="<<trackIdR_sort.size()<<endl;
	//		cout<<"trackIdL_sort.size()="<<trackIdL_sort.size()<<endl;
			sort(trackIdC_sort.begin(),trackIdC_sort.end());
			sort(trackIdL_sort.begin(),trackIdL_sort.end());
			sort(trackIdR_sort.begin(),trackIdR_sort.end());
			trackIdC_sort.erase(unique(trackIdC_sort.begin(),trackIdC_sort.end()),trackIdC_sort.end());
			trackIdL_sort.erase(unique(trackIdL_sort.begin(),trackIdL_sort.end()),trackIdL_sort.end());
			trackIdR_sort.erase(unique(trackIdR_sort.begin(),trackIdR_sort.end()),trackIdR_sort.end());
			int trackIdC = (int)trackIdC_sort.size();
			int trackIdL = (int)trackIdL_sort.size();
			int trackIdR = (int)trackIdR_sort.size();
	    	hTracksSur->Fill(2,trackIdC);
	    	hTracksSur->Fill(3,trackIdL);
	    	hTracksSur->Fill(4,trackIdR);
			if(trackflag==1||(trackflag==2&&trackIdC==1))
			{
				if((trackleftflag+trackrightflag)==0)trackFlag=kTRUE;
				if((trackleftflag>0)&&(trackIdL==1)&&(trackrightflag==0))trackFlag=kTRUE;
				if((trackrightflag>0)&&(trackIdR==1)&&(trackleftflag==0))trackFlag=kTRUE;
				if(((trackleftflag>0)&&(trackIdL>1))&&(trackrightflag==0||trackIdR==1))trackNeboileftFlag=kTRUE;
				if(((trackrightflag>0)&&(trackIdR>1))&&(trackleftflag==0||trackIdL==1))trackNeboirightFlag=kTRUE;
				if(((trackrightflag>0)&&(trackIdR>1))&&(trackleftflag>0||trackIdL>1))trackNeboiFlag=kTRUE;
				track_temp.mtdpt=temptrkVecIter->mtdpt;
				track_temp.mtdeta=temptrkVecIter->mtdeta;
				track_temp.mtdphi=temptrkVecIter->mtdphi;
				track_temp.mtdnFtPts= temptrkVecIter->mtdnFtPts;
				track_temp.mtdnDedxPts=temptrkVecIter->mtdnDedxPts;
				track_temp.mtddca=temptrkVecIter->mtddca;
				track_temp.mtdnSigmaPi=temptrkVecIter->mtdnSigmaPi;	
				track_temp.trackId=temptrkVecIter->trackId;
				track_temp.mtdBL.push_back(projbklgIter);
				track_temp.mtdmod.push_back(projmoduleIter);
				track_temp.mtdProjphi.push_back(projphi);
				track_temp.mtdProjZ.push_back(projz);
				track_temp.mtdProjLength.push_back(projlength);
				track_temp.mtdtof2Mtd.push_back(tof2Mtd);
				track_temp.trackFlag=trackFlag;
				track_temp.trackNeboileftFlag=trackNeboileftFlag;
				track_temp.trackNeboirightFlag=trackNeboirightFlag;
				track_temp.trackNeboiFlag=trackNeboiFlag;
				SingaltrackVec.push_back(track_temp);
			}
			trackIdC_sort.clear();
			trackIdR_sort.clear();
			trackIdL_sort.clear();
			temptrkVec_temp=erasetemptrkVec_temp;
		}//endl of the track loop

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//   loop over alll tracks; sort the information of two tracks extrapolations to the same module
	//   phi difference(phi - cellcenter)
	//   z difference 
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	mtdCellHitVector tempVec = trackVec;
	mtdCellHitVector erasedVec = tempVec;
	 while (tempVec.size() != 0) 
		{
			int nTracks = 0;
			
			idVector trackIdVec;
			mtdCellHitVectorIter tempIter=tempVec.begin();
			
			float projMtdZ1 = tempIter->projMtdZ;
			float projMtdPhi1 = tempIter->projMtdPhi;
			int projmodule1 = tempIter->module;
			int projMtdBackLeg1 = tempIter->backleg;
			if((projMtdBackLeg1>30||projMtdBackLeg1<1)||(projmodule1>5||projmodule1<1)) 
			{
				cout<<"wrong ::  projection backleg or module out of range!!!!!"<<endl;
				continue;
			}
			mtdCellHitVectorIter erasedIter=erasedVec.begin();
			while(erasedIter!= erasedVec.end()) {
				float projMtdZ2 = erasedIter->projMtdZ;
				float projMtdPhi2 = erasedIter->projMtdPhi;
				int projmodule2 = erasedIter->module;
				int projMtdBackLeg2 = erasedIter->backleg;
				if((projMtdBackLeg2>30||projMtdBackLeg2<1)||(projmodule2>5||projmodule2<1)) 
				{
					cout<<"wrong ::  projection backleg or module out of range!!!!!"<<endl;
					continue;
				}
				bool isMatch = kFALSE;
				if(projMtdBackLeg1 == projMtdBackLeg2 &&projmodule1 == projmodule2)isMatch = kTRUE; 
				if(nNeighbors)
				{
					if(projMtdBackLeg1==projMtdBackLeg2 && projmodule1!=projmodule2){
							if((projmodule1-1>0)&&((projmodule1-1)==projmodule2)) isMatch = true;
							if(abs(projmodule1-projmodule2)<=1) isMatch = true;
							if((projmodule1+1<6)&&((projmodule1+1)==projmodule2)) isMatch = true;
						}
				}
				if(isMatch)
				{
					nTracks++;
					if(nTracks==2)
					{
						Double_t phicenter = 0.5*TMath::Pi()+(projMtdBackLeg1-1)*(12.*TMath::Pi()/180.);
						if(phicenter>2.*(TMath::Pi()))phicenter -= 2.*(TMath::Pi());
						if(phicenter<0) phicenter += 2.*(TMath::Pi());

						hPhidis->Fill(1,phicenter-projMtdPhi1);
						hPhidis->Fill(2,phicenter-projMtdPhi2);
						hZdiff->Fill(projMtdZ1-projMtdZ2);
						hPhidiff->Fill(projMtdPhi1-projMtdPhi2);
					}
					erasedVec.erase(erasedIter);
					erasedIter--;
				}
				erasedIter++;
			}
			hTracksMulti->Fill(nTracks,refMult);
			tempVec = erasedVec;
		}
	 //cout<<"step 1.3"<<endl;
}

void bookHistograms()
{
	hVertexXY = new TH2F("hVertexXY","vertexXY;Vx (cm);Vy (cm);Counts",600,-3,3.,600,-3,3);
	hVzDiff = new TH1F("hVzDiff","Vz Diff.;Vz(TPC)-Vz(VPD) (cm);Counts",1000,-25,25);
	hvzVpdvsvz = new TH2D("hvzVpdvsvz","vzVpd vs vzTpc;vzVpd;vzTpc",400,-200,200,400,-200,200);
	hNumberofMatchHits = new TH2D("hNumberofMatchHits","Number of tracks for hits;#strip;#tracks",1900,0,1900,20,0,20);
	hTdiffWestvsStrip = new TH2D("hTdiffWestvsStrip","hitleadingedgeTime-triggerTime;# Strip;T_hitsleandingedge - T_trigger",2000,0,2000,5140,0,5140);
	hTdiffEastvsStrip = new TH2D("hTdiffEastvsStrip","hitleadingedgeTime-triggerTime;# Strip;T_hitsleandingedge - T_trigger",2000,0,2000,5140,0,5140);
	hTdiffWestvsStrip_Tac = new TH2D("hTdiffWestvsStrip_Tac","hitleadingedgeTime-triggerTime;# Strip;T_hitsleandingedge - T_trigger",2000,0,2000,5140,0,5140);
	hTdiffEastvsStrip_Tac = new TH2D("hTdiffEastvsStrip_Tac","hitleadingedgeTime-triggerTime;# Strip;T_hitsleandingedge - T_trigger",2000,0,2000,5140,0,5140);
	//double_t etarange[8]={-1,-0.5,-0.37,-0.11,0.11,0.37,0.5,1};
//	hHitEtavsPhi = new TH2D("HitsEtavsPhi","hits eta vs phi;#eta;#phi",7,etarange,30,-TMath::Pi(),TMath::Pi());
	hHitEtavsPhi = new TH2D("HitsEtavsPhi","hit denstiy;module;bkleg",5,1,6,30,1,31);
	hTrackEtavsPhi = new TH2D("TrackEtavsPhi","track density;modlue;bkleg",5,1,6,30,1,31);
	//hTrackEtavsPhi_1 = new TH2D("TrackEtavsPhi_1","track density;#eta;#phi",200,-1,1,360,-TMath::Pi(),TMath::Pi());
	hTrackEtavsPhi_1 = new TH2D("TrackEtavsPhi_1","track density;#eta;#phi",200,-1,1,360,0.,2.*TMath::Pi());
	hTrackEtavsPhi_Or = new TH2D("TrackEtavsPhi_Or","track density;#eta;#phi",200,-1,1,360,0.,2.*TMath::Pi());
	hMatchEtavsPhi = new TH2D("MatchEtavsPhi","hhhh;module;bkleg",5,1,6,30,1,31);
	hTrackMatchHits = new TH1D("TrackMatchHits","one track mathc with # hits",10,1,11);	
	hRefmultvsTracks = new TH2D("RefmultvsTracks","# match vs refmult;refmult;# of tracks matched one hit",100,0,400,10,0,10);
	hRefmult = new TH1D("Refmult","refmult;refmult;",400,0,400);
	hTracksPerHit = new TH1D("TracksPerHit","# of tracks matched with hit;# tracks matched with the same hit",50,0,50);
	hTracksPerHit_sort = new TH1D("TracksPerHit_sort","# of tracks matched with hit;# tracks matched with the same hit",50,0,50);
	hHits = new TH1D("Hits","# of hits per event",100,0,100);
	hHits_cut = new TH1D("Hits_cut","# of hits per event after time diff cut",100,0,100);
	hTracks = new TH1D("Tracks","# of tracs event after cut",700,0,700);
	hTracksMatchNumber = new TH1D("TracksMatchNumber","# of tracs event after cut",100,0,100);
	hHitDis = new TH2D("hHitDis","Hits distribution;;bkleg_diff*60+moudle_diff*12+cell_diff",10,0,10,360,-180,180);
	hHitDis_cut = new TH2D("hHitDis_cut","hitdis after matching",10,0,10,360,-180,180);
	htracksvsHits = new TH2D("htracksvsHits","# of tracks vs total hits;# of tracks;hits",50,0,50,50,0,50);
	hTracksMulti = new TH2D("hTracksMulti","# of tracks in the same module",10,0,10,1000,0,1000);
	hZdiff = new TH1D("hZdiff","track z Diff in the same modlue;z diff",200,-100,100);
	hPhidiff = new TH1D("hPhidiff","track phi diff in the same modlue; phi diff",360,-TMath::Pi(),TMath::Pi());
	hPhidis = new TH2D("hPhidis","phi distribution in the same modlue;track Id;modlue center - track phi",5,0,5,720,-TMath::Pi(),TMath::Pi());
	hTimediffvsPt = new TH2D("hTimediffvsPt","mearsured time - projection time;measred time - proj time;p_{T} (GeV/c)",800,-20,20,300,0,30);
	hTimediff_2tracks = new TH2D("hTimediff_2tracks","mearsured time - projection time;measred time - proj time;p_{T} (GeV/c)",1600,-40,40,300,0,30);
	hTimediff_Timecut_l = new TH2D("hTimediff_Timecut_l","mearsured time - projection time;measred time - proj time;p_{T} (GeV/c)",1600,-40,40,300,0,30);
	hTimediff_Timecut_s = new TH2D("hTimediff_Timecut_s","mearsured time - projection time;measred time - proj time;p_{T} (GeV/c)",1600,-40,40,300,0,30);
	hTimediff_Ptcut_l = new TH2D("hTimediff_Ptcut_l","mearsured time - projection time;measred time - proj time;p_{T} (GeV/c)",1600,-40,40,300,0,30);
	hTimediff_Ptcut_s = new TH2D("hTimediff_Ptcut_s","mearsured time - projection time;measred time - proj time;p_{T} (GeV/c)",1600,-40,40,300,0,30);
	hTimediff_Rcut_l = new TH2D("hTimediff_Rcut_l","mearsured time - projection time;measred time - proj time;p_{T} (GeV/c)",1600,-40,40,300,0,30);
	hTimediff_Rcut_s = new TH2D("hTimediff_Rcut_s","mearsured time - projection time;measred time - proj time;p_{T} (GeV/c)",1600,-40,40,300,0,30);
	hdeltaR_singlevsPt = new TH2D("hdeltaR_singlevsPt","one track one hit;#DeltaZ;p_{T} (GeV/c)",3000,-150,150,300,0,30);
	hdeltaR_both = new TH2D("hdeltaR_both","two track one hit;#DeltaZ;p_{T} (GeV/c)",3000,-150,150,300,0,30);
	hdeltaR_Timecut_l = new TH2D("hdeltaR_Timecut_l","two track one hit;#DeltaZ;p_{T} (GeV/c)",3000,-150,150,300,0,30);
	hdeltaR_Timecut_s = new TH2D("hdeltaR_Timecut_s","two track one hit;#DeltaZ;p_{T} (GeV/c)",3000,-150,150,300,0,30);
	hdeltaR_Ptcut_l = new TH2D("hdeltaR_Ptcut_l","two track one hit;#DeltaZ;p_{T} (GeV/c)",3000,-150,150,300,0,30);
	hdeltaR_Ptcut_s = new TH2D("hdeltaR_Ptcut_s","two track one hit;#DeltaZ;p_{T} (GeV/c)",3000,-150,150,300,0,30);
	hdeltaR_lvsPt = new TH2D("hdeltaR_lvsPt","two track one hit;#DeltaZ;p_{T} (GeV/c)",3000,-150,150,300,0,30);
	hdeltaR_svsPt = new TH2D("hdeltaR_svsPt","two track one hit;#DeltaZ;p_{T} (GeV/c)",3000,-150,150,300,0,30);
	hdeltaR_lvsPt_Two = new TH2D("hdeltaR_lvsPt_Two","two track one hit;#DeltaZ;p_{T} (GeV/c)",3000,-150,150,300,0,30);
	hdeltaR_svsPt_Two = new TH2D("hdeltaR_svsPt_Two","two track one hit;#DeltaZ;p_{T} (GeV/c)",3000,-150,150,300,0,30);
	hModule = new TH1D("hModule","module",100,0,100);
	hDeltazvsStrip = new TH2D("hDeltazvsStrip","deltaz vs stripId",1800,0,1800,300,-150,150);
	hDeltazvsStrip_ISMATCH = new TH2D("hDeltazvsStrip_ISMATCH","deltaz vs stripId",1800,0,1800,300,-150,150);
	hDeltazvsStrip_ISMATCHNebi = new TH2D("hDeltazvsStrip_ISMATCHNebi","deltaz vs stripId",1800,0,1800,300,-150,150);
	hDeltazvsStrip_ISMATCH_TWO = new TH2D("hDeltazvsStrip_ISMATCH_TWO","deltaz vs stripId",1800,0,1800,300,-150,150);
	hDeltazvsStrip_ISMATCHNebi_TWO = new TH2D("hDeltazvsStrip_ISMATCHNebi_TWO","deltaz vs stripId",1800,0,1800,300,-150,150);
	hTrackprojZ = new TH1D("hTrackprojZ","projz distribution",8000,-400,400);
	hTrackprojCell = new TH1D("hTrackprojCell","projCell distribution",30,-10,20);
	hTrackprojZdiffvsEta = new TH2D("hTrackprojZdiffvsEta","zdiff vs eta;eta,zdiff",20,-1,1,3000,-150,150);
	hTrackprojZdiffvsPt = new TH2D("hTrackprojZdiffvsPt","zdiff vs p_{T};p_{T} GeV/c,zdiff",300,0,30,3000,-150,150);
	hHitsMulti = new TH2D("hHitsMulti","hits surranding hits distrubution;;hits",10,0,10,100,0,100);
	hTracksSur = new TH2D("hTracksSur","tracks surranding trackss distrubution;;# 0f tracks",10,0,10,400,0,400);
	hMatchSource = new TH1D("hMatchSource","hMatchSource",20,0,20);
}
void writeHistograms(char* outFile)
{
	char buf[1024];
	sprintf(buf,"%s.histo.root",outFile);
	TFile f(buf,"recreate");
	hMatchSource->Write();
	hVertexXY->Write();
	hVzDiff->Write();
	hNumberofMatchHits->Write();
	hTdiffWestvsStrip_Tac->Write();
	hTdiffEastvsStrip_Tac->Write();
	hTdiffWestvsStrip->Write();
	hTdiffEastvsStrip->Write();
	hvzVpdvsvz->Write();
	hHitEtavsPhi->Write();
	hTrackEtavsPhi->Write();
	hTrackEtavsPhi_1->Write();
	hTrackEtavsPhi_Or->Write();
	hMatchEtavsPhi->Write();
	hTrackMatchHits->Write();
	hTracksPerHit->Write();
	hTracksPerHit_sort->Write();
	hRefmultvsTracks->Write();
	hRefmult->Write();
	hHits->Write();
	hHits_cut->Write();
	hTracks->Write();
	hTracksMatchNumber->Write();
	hHitDis->Write();
	hHitDis_cut->Write();
	htracksvsHits->Write();
	hTracksMulti->Write();
	hZdiff->Write();
	hPhidiff->Write();
	hPhidis->Write();
	hTimediffvsPt->Write();
	hTimediff_2tracks->Write();
	hTimediff_Timecut_l->Write();
	hTimediff_Timecut_s->Write();
	hTimediff_Ptcut_l->Write();
	hTimediff_Ptcut_s->Write();
	hTimediff_Rcut_l->Write();
	hTimediff_Rcut_s->Write();
	hdeltaR_singlevsPt->Write();
	hdeltaR_both->Write();
	hdeltaR_Timecut_l->Write();
	hdeltaR_Timecut_s->Write();
	hdeltaR_Ptcut_l->Write();
	hdeltaR_Ptcut_s->Write();
	hdeltaR_lvsPt->Write();
	hdeltaR_svsPt->Write();
	hdeltaR_lvsPt_Two->Write();
	hdeltaR_svsPt_Two->Write();
	hDeltazvsStrip->Write();
	hDeltazvsStrip_ISMATCH->Write();
	hDeltazvsStrip_ISMATCHNebi->Write();
	hDeltazvsStrip_ISMATCH_TWO->Write();
	hDeltazvsStrip_ISMATCHNebi_TWO->Write();
	hModule->Write();
	hTrackprojZ->Write();
	hTrackprojCell->Write();
	hTrackprojZdiffvsEta->Write();
	hTrackprojZdiffvsPt->Write();
	hHitsMulti->Write();
	hTracksSur->Write();
} 
void bookmtdRelatTree(char* outFile)
{
//  return 1;
}
bool validHits(StMtdHits& hit)
{
	int bkleg = (hit.bkleg)-1;
	int module = (hit.module)-1;
	float leTimeWest = hit.leTimeWest;
	float leTimeEast = hit.leTimeEast;
	float totWest = hit.totWest;
	float totEast = hit.totEast;
	float TdiffWest = hit.TdiffWest;
	float TdiffEast = hit.TdiffEast;
//	cout<<bkleg<<" "<<module<<endl;
	if(TdiffWest>(TacdiffWestUp[bkleg][module])||TdiffWest<(TacdiffWestbottom[bkleg][module])||(TdiffEast>TacdiffEastUp[bkleg][module])||TdiffWest<(TacdiffEastbottom[bkleg][module]))return kFALSE;
	return kTRUE;
}
void clearHits(StMtdHits& hit)
{
	hit.bkleg = -999;
	hit.module = -999;
	hit.cell = -999;
	hit.leTimeWest = -999;
	hit.leTimeEast = -999;
	hit.totEast = -999;
	hit.totWest = -999;
	hit.TdiffWest = -999;
	hit.TdiffEast = -999;
	
}
bool validTrack(StMtdTrack& track)
{
	float pt = track.mtdpt;
	float eta = track.mtdeta;
	float phi = track.mtdphi;
	float nFtPts = track.mtdnFtPts;
	float nDedxPts = track.mtdnDedxPts;
	float dca = track.mtddca;
   	float nSigmaPi = track.mtdnSigmaPi;
	idVector bklegVec = (idVector)(track.mtdBL);
	idVector moduleVec = (idVector)(track.mtdmod);
	idVector cellVec = (idVector)(track.mtdcell);
	if(bklegVec.size()==0)return kFALSE;
	int bkleg = (int)(bklegVec.back());
	int module = (int)(moduleVec.back());
	int cell = (int)(cellVec.back());
	if(pt<PT)return kFALSE;
	if(eta<ETA[0]||eta>ETA[1])return kFALSE;
	if(nFtPts<NFTPTS)return kFALSE;
	if(nDedxPts<NDEDXPTS)return kFALSE;
	if(nSigmaPi<NSIGMAPI[0]||nSigmaPi>NSIGMAPI[1])return kFALSE;
	if(bkleg<1||bkleg>30)return kFALSE;	
	if(module<1||module>5)return kFALSE;	
	return kTRUE;
}
void clearTrack(StMtdTrack& track)
{
	track.mtdpt = -999;
	track.mtdeta = -999;
	track.mtdphi = -999;
	track.mtdnFtPts = -999;
	track.mtdnDedxPts = -999;
	track.mtddca = -999;
	track.mtdnSigmaPi = -999;
	track.mtdBL.clear();
	track.mtdmod.clear();
	track.mtdcell.clear();
	track.mtdProjphi.clear();
	track.mtdProjZ.clear();
	track.mtdProjLength.clear();
	track.mtdtof2Mtd.clear();
}
