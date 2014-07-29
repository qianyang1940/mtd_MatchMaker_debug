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

#ifndef __CINT__
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TMath.h"
#include "TF1.h"
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
#include <TVirtualFitter.h>
#endif
using namespace std;
using std::cout;
using std::endl;
double Gaus(Double_t *x,Double_t *par)
{
	return par[0]*exp(-pow(x[0]-par[1],2)/(2*par[2]*par[2]));
}
void TimeWindoer()
{
	gStyle->SetOptTitle(0);
	gStyle->SetPalette(1,0);
	gStyle->SetOptLogz(1);
	gStyle->SetOptStat(000);
	gStyle->SetOptFit(1111);

	// ---------------------canvas set -----------
	TCanvas *c1 = new TCanvas("c1", "c1",600,600);
	c1->SetFillColor(10);
	c1->SetBorderMode(0);
	c1->SetBorderSize(2);
	c1->SetGrid(0,0);
	c1->SetFrameFillColor(0);
	c1->SetFrameBorderMode(0);
	c1->SetLeftMargin(0.10);
	c1->SetRightMargin(0.7);
	c1->SetBottomMargin(0.15);
	c1->SetTopMargin(0.025);
	c1->SetRightMargin(0.025);
cout<<"1"<<endl;	
	TFile *inputfile = new TFile("Timediff.histo.root"); 
	TH1D *hVzDiff = (TH1D*)inputfile->Get("hVzDiff");
	TH2D *hTdiffWestvsStrip= (TH2D*)inputfile->Get("hTdiffWestvsStrip");
	TH2D *hTdiffEastvsStrip= (TH2D*)inputfile->Get("hTdiffEastvsStrip");
	TH1D *TdiffWestdist[150];
	TH1D *TdiffEastdist[150];
cout<<"2"<<endl;
	char buff[512];
	TF1 *mygaus = new TF1("mygayyus",Gaus,2750,3000,3);
	TF1 *mygaus1 = new TF1("mygaus1",Gaus,2750,3000,3);
//	
//	mygaus->SetParameters(1000,2850,10);
/*
	ofstream outfile;
	outfile.open("timediff.dat",ios_base::trunc);
	cout<<"hah"<<endl;
	outfile<<setiosflags(ios_base::left)<<setw(12)<<"StripId";
	outfile<<setiosflags(ios_base::left)<<setw(12)<<"backleg";
	outfile<<setiosflags(ios_base::left)<<setw(12)<<"module";
	outfile<<setiosflags(ios_base::left)<<setw(12)<<"mean";
	outfile<<setiosflags(ios_base::left)<<setw(12)<<"sigma";
	outfile<<setiosflags(ios_base::left)<<setw(12)<<"chi2";
	outfile<<setiosflags(ios_base::left)<<setw(12)<<"NDF"<<endl;
	Int_t loop=0;
	int maxbin=-999;
	float mean = -999;
	
	for(Int_t bklegloop=1;bklegloop<31;bklegloop++)
	{
		for(Int_t moduleloop=1;moduleloop<6;moduleloop++)
		{
				sprintf(buff,"%d_ku",loop);
				Int_t projStartBin = (bklegloop-1)*60+(moduleloop-1)*12+2;
				Int_t projendBin = (bklegloop-1)*60+(moduleloop)*12+1;
				//cout<<"bin="<<projStartBin<<":"<<projendBin<<endl;
				TdiffWestdist[loop] = (TH1D*)hTdiffWestvsStrip->ProjectionY(buff,projStartBin,projendBin);
				maxbin = TdiffWestdist[loop]->GetMaximumBin();
				mean = TdiffWestdist[loop]->GetBinCenter(maxbin);
				mygaus->SetParameters(1000,mean,20);
				mygaus->SetParLimits(1,mean-10,mean+10);
				Float_t chi2=0;
				Int_t NDF =0;
				Float_t par[3]={0,0,0};
				//if(bklegloop==26&&(moduleloop==1||moduleloop==5))continue;
				TdiffWestdist[loop]->GetXaxis()->SetRangeUser(2600,3100);
				TdiffWestdist[loop]->Draw();
				if(mean<1000)continue;
				TdiffWestdist[loop]->Fit("mygaus","RV+");
				//mygaus->GetParameters(par);
				par[1]=	mygaus->GetParameter(1);
				par[2]=	mygaus->GetParameter(2);
				NDF=mygaus->GetNDF();
				chi2=mygaus->GetChisquare();
				outfile<<setiosflags(ios_base::left)<<setw(12)<<projStartBin;
				outfile<<setiosflags(ios_base::left)<<setw(12)<<bklegloop;
				outfile<<setiosflags(ios_base::left)<<setw(12)<<moduleloop;
				outfile<<setiosflags(ios_base::left)<<setw(12)<<par[1];
				outfile<<setiosflags(ios_base::left)<<setw(12)<<par[2];
				outfile<<setiosflags(ios_base::left)<<setw(12)<<chi2;
				outfile<<setiosflags(ios_base::left)<<setw(12)<<NDF<<endl;



				sprintf(buff,"%d_haha",loop);
				TdiffEastdist[loop] = (TH1D*)hTdiffEastvsStrip->ProjectionY(buff,projStartBin,projendBin);
				Float_t chi2=0;
				 NDF =0;
				TdiffEastdist[loop]->GetXaxis()->SetRangeUser(2600,3100);
				TdiffEastdist[loop]->Draw();
			//	TdiffEastdist[loop]->Fit("mygaus","RV+");
				maxbin = TdiffEastdist[loop]->GetMaximumBin();
				mean = TdiffEastdist[loop]->GetBinCenter(maxbin);
				mygaus1->SetParameters(1000,mean,par[2]);
				mygaus1->SetParLimits(1,mean-10,mean+10);
				mygaus1->SetParLimits(2,5,20);
				TdiffEastdist[loop]->Fit("mygaus1","RV+");
				par[1]=	mygaus1->GetParameter(1);
				par[2]=	mygaus1->GetParameter(2);
				NDF=mygaus1->GetNDF();
				chi2=mygaus1->GetChisquare();
				outfile<<setiosflags(ios_base::left)<<setw(12)<<" ";
				outfile<<setiosflags(ios_base::left)<<setw(12)<<" ";
				outfile<<setiosflags(ios_base::left)<<setw(12)<<" ";
				outfile<<setiosflags(ios_base::left)<<setw(12)<<par[1];
				outfile<<setiosflags(ios_base::left)<<setw(12)<<par[2];
				outfile<<setiosflags(ios_base::left)<<setw(12)<<chi2;
				outfile<<setiosflags(ios_base::left)<<setw(12)<<NDF<<endl;

				sprintf(buff,"./fig/%d_ku.gif",loop);
				c1->SaveAs(buff);

				loop++;
				cout<<"loop="<<loop<<endl;
				cout<<"+++++++++++++++++++++++++++++++++++++"<<endl;
		}
	}
*/
				TF1 *gAus = new TF1("gAus",Gaus,-5,5,3);
				gAus->SetParameters(1000,0,5);
				hVzDiff->GetXaxis()->SetRangeUser(-10,10);
				hVzDiff->Fit("gAus","RV+");
				sprintf(buff,"./fig/vzdiff.gif");
				c1->SaveAs(buff);
}
