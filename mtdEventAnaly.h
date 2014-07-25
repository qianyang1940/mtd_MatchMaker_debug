
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
#include "mtdEvent.h"
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

#ifndef ST_NO_TEMPLATE_DEF_ARGS
typedef vector<Int_t> idVector;
#else
typedef vector<Int_t,allocator<Int_t>> idVector;
#endif
typedef idVector::iterator idVectorIter;

struct StructCellHit{
	int backleg;
	Float_t tof2MTD;
	Float_t leTimediff;
	Float_t gtof2Mtd;
	Float_t nSigmaPi;
	int module;
	int cell;
	Float_t pt;
	Float_t eta;
	Float_t phi;
	TVector3 hitPosition;
	idVector trackIdVec;
	Int_t matchFlag;
	Float_t zhit;
	float projMtdZ;
	float projMtdPhi;
	Float_t yhit;
	pair<Double_t,Double_t> tot;
	pair<Double_t,Double_t> leadingEdgeTime;
	Int_t index2MtdHit;
	Double_t theta;
	float pathLength;
};
#ifndef ST_NO_TEMPLATE_DEF_ARGSA
typedef vector<StructCellHit> mtdCellHitVector;
#else
typedef vector<StructCellHit,allocator<StructCellHit>> mtdCellHitVector;
#endif
typedef vector<StructCellHit>::iterator mtdCellHitVectorIter;
typedef pair<Double_t,Double_t> pairD;

struct StMtdTrack{
	float mtdpt;
	float mtdeta;
	float mtdphi;
	float mtdnFtPts;
	float mtdnDedxPts;
	float mtddca;
	float mtdnSigmaPi;
	int trackId;
	vector<int> mtdBL;
	vector<int> mtdmod;
	vector<int> mtdcell;
	vector<float> mtdProjphi;
	vector<float> mtdProjZ;
	vector<float> mtdProjLength;
	vector<float> mtdtof2Mtd;
};
#ifndef ST_NO_TEMPLATE_DEF_ARGSA
typedef vector<StMtdTrack> MtdTrack;
#else
typedef vector<StMtdTrack,allocator<StMtdTrack>> MtdTrack;
#endif
typedef vector<StMtdTrack>::iterator MtdTrackVectorIter;

struct StMtdHits{
	int bkleg;
	int module;
	int cell;
	float leTimeWest;
	float leTimeEast;
	float totWest;
	float totEast;
	float TdiffWest;
	float TdiffEast;
};
#ifndef ST_NO_TEMPLATE_DEF_ARGSA
typedef vector<StMtdHits> MtdHits;
#else
typedef vector<StMtdHits,allocator<StMtdHits>> MtdHits;
#endif
typedef vector<StMtdHits>::iterator MtdHitsVectorIter;

TH2F *hVertexXY;
TH1F *hVzDiff;
TH2D *hNumberofMatchHits;
TH2D *hMatchHitsvsPt;
TH2D *hvzVpdvsvz;
TH2D *hTdiffWestvsStrip_Tac;
TH2D *hTdiffEastvsStrip_Tac;
TH2D *hTdiffWestvsStrip;
TH2D *hTdiffEastvsStrip;
TH2D *hHitEtavsPhi;
TH2D *hTrackEtavsPhi;
TH2D *hTrackEtavsPhi_1;
TH2D *hMatchEtavsPhi;
TH1D *hTrackMatchHits;
TH2D *hRefmultvsTracks;
TH1D *hRefmult;
TH1D *hTracksPerHit;
TH1D *hHits;
TH1D *hHits_cut;
TH1D *hTracks;
TH1D *hTracksMatchNumber;
TH2D *hHitDis;
TH2D *hHitDis_cut;
TH2D *htracksvsHits;
TH2D *hTracksMulti;
TH1D *hZdiff;
TH1D *hPhidiff;
TH2D *hPhidis;
TH2D *hTimediffvsPt;


TH2D *hTimediff_2tracks;
TH2D *hTimediff_Timecut_l;
TH2D *hTimediff_Timecut_s;
TH2D *hTimediff_Ptcut_l;
TH2D *hTimediff_Ptcut_s;
TH2D *hTimediff_Rcut_l;
TH2D *hTimediff_Rcut_s;
TH2D *hdeltaR_singlevsPt;
TH2D *hdeltaR_both;
TH2D *hdeltaR_Timecut_l;
TH2D *hdeltaR_Timecut_s;
TH2D *hdeltaR_Ptcut_l;
TH2D *hdeltaR_Ptcut_s;
TH2D *hdeltaR_lvsPt;
TH2D *hdeltaR_svsPt;
TH1D *hModule;
TH2D *hDeltazvsStrip;
TH1D *hTrackprojZ;
TH2D *hTrackprojZdiffvsPt;
TH2D *hTrackprojZdiffvsEta;

const int moudleNum = 13;
const double corrdeltaz[moudleNum]={-5.68597,-2.6286,4.60387,-6.072,-5.58134,-2.52121,4.77746,4.98287,-5.89513,-5.43165,-2.37928,5.2931,6.35884};
