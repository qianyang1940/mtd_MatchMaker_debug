//this .C is to plot histgram in mtd_*.root files
//No Tree analysis is done in this .C
//Qian Yang
//2014_3_27
//1.0

#include <iostream>//.h>
#include <fstream>//.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <math.h>
#include <cmath>

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TChain.h"
#include "TSystem.h"
#include "TTree.h"
#include "TBranch.h"

TLatex* drawLatex(Double_t x, Double_t y, char* text, Int_t textFont, Double_t textSize, Int_t colorIndex){
TLatex *latex = new TLatex(x,y,text);
latex->SetNDC();
latex->SetTextFont(textFont);
latex->SetTextSize(textSize);
latex->SetTextColor(colorIndex);
latex->Draw("same");
return latex;} 


plot()
{
	
	gStyle->SetOptStat(00000);
	gStyle->SetOptFit(00000);
	gStyle->SetTitleX(1.5);
	gStyle->SetStatX(0.92);
	gStyle->SetMarkerColor(2);
	gStyle->SetPalette(1);
	char inputfilename[256]="test.histo";
	char inputbuf[512];
	sprintf(inputbuf,"%s.root",inputfilename);
	TFile *inputfile = new TFile(inputbuf);
	
	if(!inputfile){
		cout<<"[e] cannot find file!!!!"<<endl;
		cout<<"[e] check filename"<<endl;
	}
	TCanvas *c1 = new TCanvas("c1","c1",800,600);
	c1->SetFillColor(10);
	c1->SetBorderMode(0);
	c1->SetBorderSize(2);
	c1->SetFrameFillColor(0);
	c1->SetFrameBorderMode(0);
	c1->SetGridx(0);
	c1->SetGridy(0);
	c1->SetLeftMargin(0.075);
	c1->SetBottomMargin(0.075);
	c1->SetTopMargin(0.025);
	c1->SetRightMargin(0.025);
	char buff[512];

	TH2F *hVertexXY = (TH2F*)inputfile->Get("hVertexXY");
	c1->SetLogz(1);
	hVerterxXY->Draw("colz")
	c1->SaveAs("vxvy.gif");

	TH1F *hVzDiff = (TH2F*)inputfile->Get("hVzDiff");
	c1->SetLogy(1);
	hVzDiff->Fit("gaus","R","",,-15,15);
	c1->SaveAs("vzdiff.gif");

	TH2D *hTdiffWestvsStrip = (TH2D*)inputfile->Get("hTdiffWestvsStrip");
	c1->SetLogy(0);
	hTdiffWestvsStrip->Draw("colz");
	c1->SaveAs("hTdiffWestvsStrip.gif");

	TH2D *hvzVpdvsvz = (TH2D*)inputfile->Get("hvzVpdvsvz");
	hvzVpdvsvz->Draw("colz");
	c1->SaveAs("vzVpdvsvz.gif");

	TH2D *HitsEtavsPhi = (TH2D*)inputfile->Get("HitsEtavsPhi");
	//c1->SetLogz(0);
	HitsEtavsPhi->GetYaxis()->SetRangeUser(23,30);
   	HitsEtavsPhi->Draw("colz");
	c1->SaveAs("hitsdensity.gif");

	TH2D *TrackEtavsPhi = (TH2D*)inputfile->Get("TrackEtavsPhi");
	TrackEtavsPhi->Draw("colz");
	c1->SaveAs("trackdensity.gif");

	TH2D *TrackEtavsPhi_1 = (TH1D*)inputfile->Get("TrackEtavsPhi_1");
	TH1D *TrackEta = (TH1D*)TrackEtavsPhi_1->ProjectionX("x1");
	TrackEta->Draw();
	c1->SaveAs("trackEta.gif");
	TH1D *TrackPhi = (TH1D*)TrackEtavsPhi_1->ProjectionY("Y1");
	TrackPhi->Draw();
	c1->SaveAs("trackPhi.gif");

	TH2D *MatchEtavsPhi = (TH2D*)inputfile->Get("MatchEtavsPhi");
	MatchEtavsPhi->Draw("colz");
	c1->SaveAs("matchEtavsPhi.gif");

	TH2D *hMatchRatio = new TH2D("hMatchRatio","hits matching Ratio;module;bkleg",5,1,6,30,1,31);
	for(Int_t i=1;i<6;i++)
	{
		for(Int_t j=1;j<31;j++)
		{
			float tothit= HitsEtavsPhi->GetBinContent(i,j);
			float mchhit = MatchEtavsPhi->GetBinContent(i,j);
		    float ratio =0;
			if(tothit==0)
			{
				ratio =0;
			};
			esle
			{
				ratio = mchhit/tothit;
			};
			hMatchRatio->SetBinContent(i,j,ratio);
		}
	}
	hMatchRatio->SetRangeUser(23,30);
	hMatchRatio->Draw("colz");
	c1->SaveAs("MatchRatio.gif");
	
	TH1D *hTrackMatchHits = (TH1D*)inputfile->Get("TrackMatchHits");
	hTrackMatchHits->Draw();
	c1->SaveAs("TrackMatchHits.gif");

	TH1D *hRefmult = (TH1D*)inputfile->Get("Refmult");
	c1->SetLogy(1);
   	hRefmult->Draw();
	c1->SaveAs("refmult.gif");
	c1->SetLogy(0);

	TH2D* hRefmultvsTracks = (TH2D*)inputfile->Get("RefmultvsTracks");
	hRefmultvsTracks->Draw("colz");
	c1->SaveAs("refmultvsTracks.gif");

	TH1D *hHits = (TH1D*)inputfile->Get("Hits");
	TH1D *hHits_cut = (TH1D*)inputfile->Get("Hits_cut");
	hHits->Draw();
	hHits_cut->SetLineColor(kRed);
	hHits_cut->Draw("same");
	sprintf(buff,before time diff cut);
	drawLatex(0.2,0.85,buff,62,0.06,1);
	sprintf(buff,after time diff cut);
	drawLatex(0.2,0.65,buff,62,0.06,2);
	c1->SaveAs("hit_distr.gif");

	TH1D* hTracks = (TH1D*)inputfile->Get("Tracks");
	hTracks->Draw("colz");
	c1->SaveAs("numb_track_cut.gif");

	TH2D* hHitDis = (TH2D*)inputfile->Get("hHitDis");
	hHitDis->Draw("colz");
	c1->SaveAs("hits_distruibution.gif");
	
	TH1D *hTimediff = (TH1D*)inputfile->Get("hTimediff");
	hTimediff->Draw();
	c1->SaveAs("Timediff_single.gif");

	TH1D *hTimediff_2tracks = (TH1D*)inputfile->Get("hTimediff_2tracks");
	hTimediff_2tracks->Draw();
	c1->SaveAs("Timediff_2track.gif");

	TH1D *hdeltaR_single = (TH1D*)inputfile->Get("hdetlaR_single");
	hdetlaR_single->Draw();
	c1->SaveAs("deltaR_single.gif");

	TH1D *hdeltaR_both = (TH1D*)inputfile->Get("hdetlaR_both");
	hdetlaR_both->Draw();
	c1->SaveAs("deltaR_both.gif");
	
	TH1D *hdeltaR_Timecut_l= (TH1D*)inputfile->Get("hdeltaR_Timecut_l");
	hdeltaR_Timecut_l->Draw();
	c1->SaveAs("deltaR_timecut_l.gif");

	TH1D *hdeltaR_Timecut_s= (TH1D*)inputfile->Get("hdeltaR_Timecut_s");
	hdeltaR_Timecut_s->Draw();
	c1->SaveAs("deltaR_timecut_s.gif");

	TH1D *hdeltaR_Ptcut_l= (TH1D*)inputfile->Get("hdeltaR_Ptcut_l");
	hdeltaR_Ptcut_l->Draw();
	c1->SaveAs("deltaR_Ptcut_l.gif");

	TH1D *hdeltaR_Ptcut_s= (TH1D*)inputfile->Get("hdeltaR_Ptcut_s");
	hdeltaR_Ptcut_s->Draw();
	c1->SaveAs("deltaR_Ptcut_s.gif");

	TH1D *hdeltaR_cut_l= (TH1D*)inputfile->Get("hdeltaR_cut_l");
	hdeltaR_cut_l->Draw();
	c1->SaveAs("deltaR_cut_l.gif");

	TH1D *hdeltaR_cut_s= (TH1D*)inputfile->Get("hdeltaR_cut_s");
	hdeltaR_cut_s->Draw();
	c1->SaveAs("deltaR_cut_s.gif");
}
