//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 14 22:27:17 2014 by ROOT version 5.34/09
// from TTree mtdEvent/mtdEvent
// found on file: ../mtd_test.root
//////////////////////////////////////////////////////////

#ifndef mtdEvent_h
#define mtdEvent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class mtdEvent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           evt;
   Int_t           nTrgIds;
   Int_t           trgId[5];   //[nTrgIds]
   Float_t         bField;
   Float_t         vertexX;
   Float_t         vertexY;
   Float_t         vertexZ;
   Double_t        triggerTime[2];
   Float_t         vpdVz;
   Float_t         tStart;
   Int_t           refMult;
   Int_t           nMtdRawHits;
   Char_t          flag[1000];   //[nMtdRawHits]
   UChar_t         backlegRaw[1000];   //[nMtdRawHits]
   UChar_t         chn[1000];   //[nMtdRawHits]
   Double_t        tdc[1000];   //[nMtdRawHits]
   Int_t           nMtdHits;
   UChar_t         backleg[200];   //[nMtdHits]
   UChar_t         module[200];   //[nMtdHits]
   UChar_t         cell[200];   //[nMtdHits]
   Double_t        leTimeWest[200];   //[nMtdHits]
   Double_t        leTimeEast[200];   //[nMtdHits]
   Double_t        totWest[200];   //[nMtdHits]
   Double_t        totEast[200];   //[nMtdHits]
   Double_t        TdiffWest[200];   //[nMtdHits]
   Double_t        TdiffEast[200];   //[nMtdHits]
   Int_t           ngTrks;
   Float_t         gpt[1000];   //[ngTrks]
   Float_t         geta[1000];   //[ngTrks]
   Float_t         gphi[1000];   //[ngTrks]
   Float_t         ppt[1000];   //[ngTrks]
   Float_t         peta[1000];   //[ngTrks]
   Float_t         pphi[1000];   //[ngTrks]
   Float_t         ghelixpx[1000];   //[ngTrks]
   Float_t         ghelixpy[1000];   //[ngTrks]
   Float_t         ghelixpz[1000];   //[ngTrks]
   Float_t         ghelixox[1000];   //[ngTrks]
   Float_t         ghelixoy[1000];   //[ngTrks]
   Float_t         ghelixoz[1000];   //[ngTrks]
   Float_t         gdedx[1000];   //[ngTrks]
   Float_t         gnSigmaPi[1000];   //[ngTrks]
   Float_t         gnSigmaK[1000];   //[ngTrks]
   Float_t         gnSigmaP[1000];   //[ngTrks]
   Float_t         gnSigmaE[1000];   //[ngTrks]
   Char_t          gq[1000];   //[ngTrks]
   Int_t           gtrackId[1000];   //[ngTrks]
   Int_t           gIndex2Primary[1000];   //[ngTrks]
   Char_t          gnFtPts[1000];   //[ngTrks]
   Char_t          gnDedxPts[1000];   //[ngTrks]
   Int_t           gchannel[1000];   //[ngTrks]
   Float_t         gyLocal[1000];   //[ngTrks]
   Float_t         gzLocal[1000];   //[ngTrks]
   Float_t         gtdc[1000];   //[ngTrks]
   Float_t         gtot[1000];   //[ngTrks]
   Float_t         gtof[1000];   //[ngTrks]
   Float_t         gpathLength[1000];   //[ngTrks]
   Float_t         gbeta[1000];   //[ngTrks]
   Float_t         gtdiff[1000];   //[ngTrks]
   Float_t         gdca[1000];   //[ngTrks]
   Int_t           gtrackindex[1000];   //[ngTrks]
   Int_t           gTrkMatchNum[1000];   //[ngTrks]
   UChar_t         gprojMtdBackLeg[1000][2];   //[ngTrks]
   UChar_t         gprojMtdModule[1000][2];   //[ngTrks]
   UChar_t         gprojMtdCell[1000][2];   //[ngTrks]
   Float_t         gprojMtdPhi[1000][2];   //[ngTrks]
   Float_t         gprojMtdZ[1000][2];   //[ngTrks]
   Float_t         gprojMtdLength[1000][2];   //[ngTrks]
   Float_t         gtof2Mtd[1000][2];   //[ngTrks]
   Int_t           gnMatchMtdHits[1000];   //[ngTrks]
   Int_t           gmMtdHitIndex[1000];   //[ngTrks]
   UChar_t         gmBackLeg[1000];   //[ngTrks]
   UChar_t         gmModule[1000];   //[ngTrks]
   UChar_t         gmCell[1000];   //[ngTrks]
   Float_t         gmLeTimeWest[1000];   //[ngTrks]
   Float_t         gmTotWest[1000];   //[ngTrks]
   Float_t         gmLeTimeEast[1000];   //[ngTrks]
   Float_t         gmTotEast[1000];   //[ngTrks]
   Float_t         gmLocalZ[1000];   //[ngTrks]
   Float_t         gmLocalY[1000];   //[ngTrks]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_nTrgIds;   //!
   TBranch        *b_trgId;   //!
   TBranch        *b_bField;   //!
   TBranch        *b_vertexX;   //!
   TBranch        *b_vertexY;   //!
   TBranch        *b_vertexZ;   //!
   TBranch        *b_triggerTime;   //!
   TBranch        *b_vpdVz;   //!
   TBranch        *b_tStart;   //!
   TBranch        *b_refMult;   //!
   TBranch        *b_nMtdRawHits;   //!
   TBranch        *b_flag;   //!
   TBranch        *b_backlegRaw;   //!
   TBranch        *b_chn;   //!
   TBranch        *b_tdc;   //!
   TBranch        *b_nMtdHits;   //!
   TBranch        *b_backleg;   //!
   TBranch        *b_module;   //!
   TBranch        *b_cell;   //!
   TBranch        *b_leTimeWest;   //!
   TBranch        *b_leTimeEast;   //!
   TBranch        *b_totWest;   //!
   TBranch        *b_totEast;   //!
   TBranch        *b_TdiffWest;   //!
   TBranch        *b_TdiffEast;   //!
   TBranch        *b_ngTrks;   //!
   TBranch        *b_gpt;   //!
   TBranch        *b_geta;   //!
   TBranch        *b_gphi;   //!
   TBranch        *b_ppt;   //!
   TBranch        *b_peta;   //!
   TBranch        *b_pphi;   //!
   TBranch        *b_ghelixpx;   //!
   TBranch        *b_ghelixpy;   //!
   TBranch        *b_ghelixpz;   //!
   TBranch        *b_ghelixox;   //!
   TBranch        *b_ghelixoy;   //!
   TBranch        *b_ghelixoz;   //!
   TBranch        *b_gdedx;   //!
   TBranch        *b_gnSigmaPi;   //!
   TBranch        *b_gnSigmaK;   //!
   TBranch        *b_gnSigmaP;   //!
   TBranch        *b_gnSigmaE;   //!
   TBranch        *b_gq;   //!
   TBranch        *b_gtrackId;   //!
   TBranch        *b_gIndex2Primary;   //!
   TBranch        *b_gnFtPts;   //!
   TBranch        *b_gnDedxPts;   //!
   TBranch        *b_gchannel;   //!
   TBranch        *b_gyLocal;   //!
   TBranch        *b_gzLocal;   //!
   TBranch        *b_gtdc;   //!
   TBranch        *b_gtot;   //!
   TBranch        *b_gtof;   //!
   TBranch        *b_gpathLength;   //!
   TBranch        *b_gbeta;   //!
   TBranch        *b_gtdiff;   //!
   TBranch        *b_gdca;   //!
   TBranch        *b_gtrackindex;   //!
   TBranch        *b_gTrkMatchNum;   //!
   TBranch        *b_gprojMtdBackLeg;   //!
   TBranch        *b_gprojMtdModule;   //!
   TBranch        *b_gprojMtdCell;   //!
   TBranch        *b_gprojMtdPhi;   //!
   TBranch        *b_gprojMtdZ;   //!
   TBranch        *b_gprojMtdLength;   //!
   TBranch        *b_gtof2Mtd;   //!
   TBranch        *b_gnMatchMtdHits;   //!
   TBranch        *b_gmMtdHitIndex;   //!
   TBranch        *b_gmBackLeg;   //!
   TBranch        *b_gmModule;   //!
   TBranch        *b_gmCell;   //!
   TBranch        *b_gmLeTimeWest;   //!
   TBranch        *b_gmTotWest;   //!
   TBranch        *b_gmLeTimeEast;   //!
   TBranch        *b_gmTotEast;   //!
   TBranch        *b_gmLocalZ;   //!
   TBranch        *b_gmLocalY;   //!

   mtdEvent(TTree *tree=0);
   virtual ~mtdEvent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef mtdEvent_cxx
mtdEvent::mtdEvent(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../mtd_test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../mtd_test.root");
      }
      f->GetObject("mtdEvent",tree);

   }
   Init(tree);
}

mtdEvent::~mtdEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mtdEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mtdEvent::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void mtdEvent::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("nTrgIds", &nTrgIds, &b_nTrgIds);
   fChain->SetBranchAddress("trgId", trgId, &b_trgId);
   fChain->SetBranchAddress("bField", &bField, &b_bField);
   fChain->SetBranchAddress("vertexX", &vertexX, &b_vertexX);
   fChain->SetBranchAddress("vertexY", &vertexY, &b_vertexY);
   fChain->SetBranchAddress("vertexZ", &vertexZ, &b_vertexZ);
   fChain->SetBranchAddress("triggerTime", triggerTime, &b_triggerTime);
   fChain->SetBranchAddress("vpdVz", &vpdVz, &b_vpdVz);
   fChain->SetBranchAddress("tStart", &tStart, &b_tStart);
   fChain->SetBranchAddress("refMult", &refMult, &b_refMult);
   fChain->SetBranchAddress("nMtdRawHits", &nMtdRawHits, &b_nMtdRawHits);
   fChain->SetBranchAddress("flag", flag, &b_flag);
   fChain->SetBranchAddress("backlegRaw", backlegRaw, &b_backlegRaw);
   fChain->SetBranchAddress("chn", chn, &b_chn);
   fChain->SetBranchAddress("tdc", tdc, &b_tdc);
   fChain->SetBranchAddress("nMtdHits", &nMtdHits, &b_nMtdHits);
   fChain->SetBranchAddress("backleg", backleg, &b_backleg);
   fChain->SetBranchAddress("module", module, &b_module);
   fChain->SetBranchAddress("cell", cell, &b_cell);
   fChain->SetBranchAddress("leTimeWest", leTimeWest, &b_leTimeWest);
   fChain->SetBranchAddress("leTimeEast", leTimeEast, &b_leTimeEast);
   fChain->SetBranchAddress("totWest", totWest, &b_totWest);
   fChain->SetBranchAddress("totEast", totEast, &b_totEast);
   fChain->SetBranchAddress("TdiffWest", TdiffWest, &b_TdiffWest);
   fChain->SetBranchAddress("TdiffEast", TdiffEast, &b_TdiffEast);
   fChain->SetBranchAddress("ngTrks", &ngTrks, &b_ngTrks);
   fChain->SetBranchAddress("gpt", gpt, &b_gpt);
   fChain->SetBranchAddress("geta", geta, &b_geta);
   fChain->SetBranchAddress("gphi", gphi, &b_gphi);
   fChain->SetBranchAddress("ppt", ppt, &b_ppt);
   fChain->SetBranchAddress("peta", peta, &b_peta);
   fChain->SetBranchAddress("pphi", pphi, &b_pphi);
   fChain->SetBranchAddress("ghelixpx", ghelixpx, &b_ghelixpx);
   fChain->SetBranchAddress("ghelixpy", ghelixpy, &b_ghelixpy);
   fChain->SetBranchAddress("ghelixpz", ghelixpz, &b_ghelixpz);
   fChain->SetBranchAddress("ghelixox", ghelixox, &b_ghelixox);
   fChain->SetBranchAddress("ghelixoy", ghelixoy, &b_ghelixoy);
   fChain->SetBranchAddress("ghelixoz", ghelixoz, &b_ghelixoz);
   fChain->SetBranchAddress("gdedx", gdedx, &b_gdedx);
   fChain->SetBranchAddress("gnSigmaPi", gnSigmaPi, &b_gnSigmaPi);
   fChain->SetBranchAddress("gnSigmaK", gnSigmaK, &b_gnSigmaK);
   fChain->SetBranchAddress("gnSigmaP", gnSigmaP, &b_gnSigmaP);
   fChain->SetBranchAddress("gnSigmaE", gnSigmaE, &b_gnSigmaE);
   fChain->SetBranchAddress("gq", gq, &b_gq);
   fChain->SetBranchAddress("gtrackId", gtrackId, &b_gtrackId);
   fChain->SetBranchAddress("gIndex2Primary", gIndex2Primary, &b_gIndex2Primary);
   fChain->SetBranchAddress("gnFtPts", gnFtPts, &b_gnFtPts);
   fChain->SetBranchAddress("gnDedxPts", gnDedxPts, &b_gnDedxPts);
   fChain->SetBranchAddress("gchannel", gchannel, &b_gchannel);
   fChain->SetBranchAddress("gyLocal", gyLocal, &b_gyLocal);
   fChain->SetBranchAddress("gzLocal", gzLocal, &b_gzLocal);
   fChain->SetBranchAddress("gtdc", gtdc, &b_gtdc);
   fChain->SetBranchAddress("gtot", gtot, &b_gtot);
   fChain->SetBranchAddress("gtof", gtof, &b_gtof);
   fChain->SetBranchAddress("gpathLength", gpathLength, &b_gpathLength);
   fChain->SetBranchAddress("gbeta", gbeta, &b_gbeta);
   fChain->SetBranchAddress("gtdiff", gtdiff, &b_gtdiff);
   fChain->SetBranchAddress("gdca", gdca, &b_gdca);
   fChain->SetBranchAddress("gtrackindex", gtrackindex, &b_gtrackindex);
   fChain->SetBranchAddress("gTrkMatchNum", gTrkMatchNum, &b_gTrkMatchNum);
   fChain->SetBranchAddress("gprojMtdBackLeg", gprojMtdBackLeg, &b_gprojMtdBackLeg);
   fChain->SetBranchAddress("gprojMtdModule", gprojMtdModule, &b_gprojMtdModule);
   fChain->SetBranchAddress("gprojMtdCell", gprojMtdCell, &b_gprojMtdCell);
   fChain->SetBranchAddress("gprojMtdPhi", gprojMtdPhi, &b_gprojMtdPhi);
   fChain->SetBranchAddress("gprojMtdZ", gprojMtdZ, &b_gprojMtdZ);
   fChain->SetBranchAddress("gprojMtdLength", gprojMtdLength, &b_gprojMtdLength);
   fChain->SetBranchAddress("gtof2Mtd", gtof2Mtd, &b_gtof2Mtd);
   fChain->SetBranchAddress("gnMatchMtdHits", gnMatchMtdHits, &b_gnMatchMtdHits);
   fChain->SetBranchAddress("gmMtdHitIndex", gmMtdHitIndex, &b_gmMtdHitIndex);
   fChain->SetBranchAddress("gmBackLeg", gmBackLeg, &b_gmBackLeg);
   fChain->SetBranchAddress("gmModule", gmModule, &b_gmModule);
   fChain->SetBranchAddress("gmCell", gmCell, &b_gmCell);
   fChain->SetBranchAddress("gmLeTimeWest", gmLeTimeWest, &b_gmLeTimeWest);
   fChain->SetBranchAddress("gmTotWest", gmTotWest, &b_gmTotWest);
   fChain->SetBranchAddress("gmLeTimeEast", gmLeTimeEast, &b_gmLeTimeEast);
   fChain->SetBranchAddress("gmTotEast", gmTotEast, &b_gmTotEast);
   fChain->SetBranchAddress("gmLocalZ", gmLocalZ, &b_gmLocalZ);
   fChain->SetBranchAddress("gmLocalY", gmLocalY, &b_gmLocalY);
   Notify();
}

Bool_t mtdEvent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mtdEvent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mtdEvent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mtdEvent_cxx
