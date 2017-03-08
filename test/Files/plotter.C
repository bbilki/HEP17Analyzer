#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TRandom3.h"
#include "TH1F.h"
#include "TH1.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TColor.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TClonesArray.h"

#include <sys/types.h>
#include <dirent.h>

#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TVirtualFitter.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <TVector3.h>
#include <TMultiGraph.h>
#include <TH2Poly.h>
#include <TLine.h>

using namespace std;using namespace ROOT::Math;

struct edata
{
	vector <Int_t> *ieta;
	vector <Int_t> *iphi;
	vector <Int_t> *depth;
	vector <vector <Int_t>> *pulse;
	vector <vector <Int_t>> *tdc;
	vector <vector <Int_t>> *capid;
	vector <vector <Int_t>> *soi;
	vector <vector <Int_t>> *histo;
};
edata ed;

struct pedlist
{
	int ieta;
	int iphi;
	int depth;
	float pedmean;
	float pedsigma;
	float pedmeanerr;
	float pedsigmaerr;

	TH1F* h;
	TH1F* qtot; //total charge
	TH1F* p;//pulse
	TH1F* pn;//pulsenorm
};

struct ledlist
{
	int ieta;
	int iphi;
	int depth;
	float Qmean;
	float Qsigma;
	float Qmeanerr;
	float Qsigmaerr;
	float nbadcapid;
	TH1F* p;//pulse
	TH1F* pn;//pulsenorm
	TH1F* t;//tdc
	TH1F* q;//charge
};

double adc2fC_QIE11[256]={1.58, 4.73, 7.88, 11.0, 14.2, 17.3, 20.5, 23.6,26.8, 29.9, 33.1, 36.2, 39.4, 42.5, 45.7, 48.8,53.6, 60.1, 66.6, 73.0, 79.5, 86.0, 92.5, 98.9,105, 112, 118, 125, 131, 138, 144, 151,157, 164, 170, 177, 186, 199, 212, 225, 238, 251, 264, 277, 289, 302, 315, 328,341, 354, 367, 380, 393, 406, 418, 431,444, 464, 490, 516, 542, 568, 594, 620,569, 594, 619, 645, 670, 695, 720, 745,771, 796, 821, 846, 871, 897, 922, 947,960, 1010, 1060, 1120, 1170, 1220, 1270, 1320,1370, 1430, 1480, 1530, 1580, 1630, 1690, 1740,1790, 1840, 1890, 1940,2020, 2120, 2230, 2330,2430, 2540, 2640, 2740, 2850, 2950, 3050, 3150, 3260, 3360, 3460, 3570, 3670, 3770, 3880, 3980,4080, 4240, 4450, 4650, 4860, 5070, 5280, 5490,5080, 5280, 5480, 5680, 5880, 6080, 6280, 6480,6680, 6890, 7090, 7290, 7490, 7690, 7890, 8090,8400, 8810, 9220, 9630, 10000, 10400, 10900, 11300,11700, 12100, 12500, 12900, 13300, 13700, 14100, 14500,15000, 15400, 15800, 16200, 16800, 17600, 18400, 19300,20100, 20900, 21700, 22500, 23400, 24200, 25000, 25800,26600, 27500, 28300, 29100, 29900, 30700, 31600, 32400,33200, 34400, 36100, 37700, 39400, 41000, 42700, 44300,41100, 42700, 44300, 45900, 47600, 49200, 50800, 52500,54100, 55700, 57400, 59000, 60600, 62200, 63900, 65500,68000, 71300, 74700, 78000, 81400, 84700, 88000, 91400,94700, 98100, 101000, 105000, 108000,111000, 115000, 118000,121000, 125000, 128000, 131000, 137000, 145000, 152000, 160000,168000, 176000, 183000, 191000, 199000, 206000, 214000, 222000,230000, 237000, 245000, 253000, 261000, 268000, 276000, 284000,291000, 302000, 316000, 329000, 343000, 356000, 370000, 384000};

int HFMBoxMap[37]={0,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19};

int RunNo=0;int RunType=0;

int print()
{
	char hname[500];
	sprintf(hname,"../NTuples/N_%d.root",RunNo);
	TFile* inroot=new TFile(hname);
	TTree *tree = (TTree*)inroot->Get("Events");
	if(RunType==0 || RunType==1 || RunType==2 || RunType==3 || RunType==4)
	{
		tree->SetBranchAddress("ieta",&ed.ieta);
		tree->SetBranchAddress("iphi",&ed.iphi);
		tree->SetBranchAddress("depth",&ed.depth);
		tree->SetBranchAddress("pulse",&ed.pulse);
		tree->SetBranchAddress("tdc",&ed.tdc);
		tree->SetBranchAddress("capid",&ed.capid);
		tree->SetBranchAddress("soi",&ed.soi);
	}
	for(int i=0;i<tree->GetEntries();i++)
	{
		tree->GetEntry(i);
		for(int i1=0;i1<ed.ieta->size();i1++)
		{
			cout<<ed.ieta->at(i1)<<" "<<ed.iphi->at(i1)<<" "<<ed.depth->at(i1)<<endl;
			cout<<"                   "<<ed.pulse->at(i1)[0]<<" "<<ed.pulse->at(i1)[1]<<" "<<ed.pulse->at(i1)[2]<<" "<<ed.pulse->at(i1)[3]<<" "<<ed.pulse->at(i1)[4]<<" "<<ed.pulse->at(i1)[5]<<" "<<ed.pulse->at(i1)[6]<<" "<<ed.pulse->at(i1)[7]<<" "<<ed.pulse->at(i1)[8]<<" "<<ed.pulse->at(i1)[9]<<endl;
			cout<<"                   "<<ed.capid->at(i1)[0]<<" "<<ed.capid->at(i1)[1]<<" "<<ed.capid->at(i1)[2]<<" "<<ed.capid->at(i1)[3]<<" "<<ed.capid->at(i1)[4]<<" "<<ed.capid->at(i1)[5]<<" "<<ed.capid->at(i1)[6]<<" "<<ed.capid->at(i1)[7]<<" "<<ed.capid->at(i1)[8]<<" "<<ed.capid->at(i1)[9]<<endl;
			cout<<"                   "<<ed.tdc->at(i1)[0]<<" "<<ed.tdc->at(i1)[1]<<" "<<ed.tdc->at(i1)[2]<<" "<<ed.tdc->at(i1)[3]<<" "<<ed.tdc->at(i1)[4]<<" "<<ed.tdc->at(i1)[5]<<" "<<ed.tdc->at(i1)[6]<<" "<<ed.tdc->at(i1)[7]<<" "<<ed.tdc->at(i1)[8]<<" "<<ed.tdc->at(i1)[9]<<endl;
		}
	}
	inroot->Close();
}

int plotpeds()
{
	char hname[500];
	sprintf(hname,"../NTuples/N_%d.root",RunNo);
	TFile* inroot=new TFile(hname);
	TTree *tree = (TTree*)inroot->Get("Events");
	tree->SetBranchAddress("ieta",&ed.ieta);
	tree->SetBranchAddress("iphi",&ed.iphi);
	tree->SetBranchAddress("depth",&ed.depth);
	tree->SetBranchAddress("pulse",&ed.pulse);
	tree->SetBranchAddress("tdc",&ed.tdc);
	tree->SetBranchAddress("capid",&ed.capid);
	tree->SetBranchAddress("soi",&ed.soi);
	
	pedlist pl;vector <pedlist> PL;int plind=-1;
	tree->GetEntry(0);
	for(int i1=0;i1<ed.ieta->size();i1++)
	{
		plind=-1;
		for(int i2=0;i2<PL.size();i2++)
		{
			if(PL[i2].ieta==ed.ieta->at(i1) && PL[i2].iphi==ed.iphi->at(i1) && PL[i2].depth==ed.depth->at(i1))
			{
				plind=i2;break;
			}
		}
		if(plind==-1)
		{
			pl.ieta=ed.ieta->at(i1);
			pl.iphi=ed.iphi->at(i1);
			pl.depth=ed.depth->at(i1);
			pl.pedmean=0.;
			pl.pedsigma=0.;
			pl.pedmeanerr=0.;
			pl.pedsigmaerr=0.;
			sprintf(hname,"ADC Pulse %d %d %d",pl.ieta,pl.iphi,pl.depth);
			pl.p=new TH1F(hname,hname,10,-0.5,9.5); //pl.p is a pointer.
			pl.p->SetLineWidth(2);pl.p->SetLineColor(4);pl.p->SetFillColor(4);
			pl.p->GetXaxis()->SetTitle("TS");pl.p->GetXaxis()->CenterTitle();
			pl.p->GetYaxis()->SetTitle("Mean ADC / TS");pl.p->GetYaxis()->CenterTitle();
			sprintf(hname,"ADC Pulse %d %d %d norm",pl.ieta,pl.iphi,pl.depth);
			pl.pn=new TH1F(hname,hname,10,-0.5,9.5);
			sprintf(hname,"ADC Ped %d %d %d",pl.ieta,pl.iphi,pl.depth);
			pl.h=new TH1F(hname,hname,10,-0.5,9.5);
			pl.h->SetLineWidth(2);pl.h->SetLineColor(1);
			pl.h->GetXaxis()->SetTitle("ADC");pl.h->GetXaxis()->CenterTitle();
			pl.h->GetYaxis()->SetTitle("TS / ADC");pl.h->GetYaxis()->CenterTitle();
			sprintf(hname,"Qtot %d %d %d",pl.ieta,pl.iphi,pl.depth);
			pl.qtot=new TH1F(hname,hname,200,0,1000);
			pl.qtot->SetLineWidth(2);pl.qtot->SetLineColor(1);
			pl.qtot->GetXaxis()->SetTitle("Total Charge (fC)");pl.qtot->GetXaxis()->CenterTitle();
			pl.qtot->GetYaxis()->SetTitle("Events / 5fC");pl.qtot->GetYaxis()->CenterTitle();
			PL.push_back(pl);
		}
	}
	for(int i=0;i<tree->GetEntries();i++)
// 	for(int i=0;i<10;i++)
	{
		tree->GetEntry(i);
		if(i%200==0) cout<<i<<" / "<<tree->GetEntries()<<endl;
		for(int i1=0;i1<ed.ieta->size();i1++)
		{
			plind=-1;
			for(int i2=0;i2<PL.size();i2++)
			{
				if(PL[i2].ieta==ed.ieta->at(i1) && PL[i2].iphi==ed.iphi->at(i1) && PL[i2].depth==ed.depth->at(i1))
				{
					plind=i2;break;
				}
			}
			float Qtot=0; 
			for(int i2=0;i2<10;i2++)
			{
				PL[plind].h->Fill(ed.pulse->at(i1)[i2]);
				PL[plind].p->Fill(i2,ed.pulse->at(i1)[i2]);
				PL[plind].pn->Fill(i2);
				if(i2>0)
				{
					Qtot+=adc2fC_QIE11[ed.pulse->at(i1)[i2]];
				}
			}
			PL[plind].qtot->Fill(Qtot);
		}
	}
	double HMAX=0.;
	double HMAXP=0.;
	double QTXMIN=1000.;
	double QTXMAX=0.;
	double QTYMAX=0.;
	int iside=0;int iquad=0;
	TF1* tf=new TF1("gaus","gaus",0.,10.);
	for(int i2=0;i2<PL.size();i2++)
	{
		PL[i2].p->Divide(PL[i2].pn);
		if(PL[i2].p->GetBinContent(PL[i2].p->GetMaximumBin())>HMAXP){HMAXP=PL[i2].p->GetBinContent(PL[i2].p->GetMaximumBin());}
		if(PL[i2].qtot->GetBinContent(PL[i2].qtot->GetMaximumBin())>QTYMAX){QTYMAX=PL[i2].qtot->GetBinContent(PL[i2].qtot->GetMaximumBin());}
		if(PL[i2].qtot->GetBinCenter(PL[i2].qtot->FindFirstBinAbove(0))<QTXMIN){QTXMIN=PL[i2].qtot->GetBinCenter(PL[i2].qtot->FindFirstBinAbove(0));}
		if(PL[i2].qtot->GetBinCenter(PL[i2].qtot->FindLastBinAbove(1))>QTXMAX){QTXMAX=PL[i2].qtot->GetBinCenter(PL[i2].qtot->FindLastBinAbove(1));}
// 		if(PL[i2].h->GetMean()>0)
// 		{
			PL[i2].h->Fit(tf,"q","q",0.,10.);
			PL[i2].pedmean=tf->GetParameter(1);
			PL[i2].pedsigma=tf->GetParameter(2);
			PL[i2].pedmeanerr=tf->GetParError(1);
			PL[i2].pedsigmaerr=tf->GetParError(2);
			tf->SetLineStyle(2);
			if(PL[i2].h->GetBinContent(PL[i2].h->GetMaximumBin())>HMAX){HMAX=PL[i2].h->GetBinContent(PL[i2].h->GetMaximumBin());}
// 		}
// 		else
// 		{
// 			PL[i2].pedmean=0.;
// 			PL[i2].pedsigma=0.;
// 		}
	}
	HMAX+=20.;HMAXP+=20.;
	TFile* outPeds=new TFile("Pedestals.root","recreate");
	{
		TCanvas* cc[6];
		cc[0]=new TCanvas("cc1","cc1",4900,4900);//iphi=63, ieta=16-22
		gStyle->SetOptStat(0);
		gStyle->SetTitleFontSize(0.1);
		cc[0]->Divide(7,7,0,0);
		cc[1]=new TCanvas("cc2","cc2",4900,4900);//iphi=63, ieta=23-29
		cc[1]->Divide(7,7,0,0);
		cc[2]=new TCanvas("cc3","cc3",4900,3500);//iphi=64, ieta=16-20
		cc[2]->Divide(7,5,0,0);
		cc[3]=new TCanvas("cc4","cc4",4900,4900);//iphi=65, ieta=16-22
		cc[3]->Divide(7,7,0,0);
		cc[4]=new TCanvas("cc5","cc5",4900,4900);//iphi=65, ieta=23-29
		cc[4]->Divide(7,7,0,0);
		cc[5]=new TCanvas("cc6","cc6",4900,3500);//iphi=66, ieta=16-20
		cc[5]->Divide(7,5,0,0);
		int cci=0;
		
		for(int i2=0;i2<PL.size();i2++)
		{
			int ccind=-1; 
			if(PL[i2].iphi==63 && PL[i2].ieta>=16 && PL[i2].ieta<=22){ccind=0;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==63 && PL[i2].ieta>=23 && PL[i2].ieta<=29) {ccind=1;cci=(PL[i2].ieta-23)*7+PL[i2].depth;}
			else if(PL[i2].iphi==64) {ccind=2;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==65 && PL[i2].ieta>=16 && PL[i2].ieta<=22) {ccind=3;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==65 && PL[i2].ieta>=23 && PL[i2].ieta<=29) {ccind=4;cci=(PL[i2].ieta-23)*7+PL[i2].depth;}
			else if(PL[i2].iphi==66) {ccind=5;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			cc[ccind]->cd(cci);
			PL[i2].p->Draw("hist");
// 			PL[i2].p->GetYaxis()->SetRangeUser(0.,HMAXP);
			PL[i2].p->GetYaxis()->SetRangeUser(0.,20.);
			PL[i2].p->Write();
		}
		
		cc[0]->SaveAs("PulseShapes.pdf(");
		cc[1]->SaveAs("PulseShapes.pdf");
		cc[2]->SaveAs("PulseShapes.pdf");
		cc[3]->SaveAs("PulseShapes.pdf");
		cc[4]->SaveAs("PulseShapes.pdf");
		cc[5]->SaveAs("PulseShapes.pdf)");
		for(int i1=0;i1<6;i1++){delete cc[i1];}
	}
	{
		TCanvas* cc[6];
		cc[0]=new TCanvas("cc1","cc1",4900,4900);//iphi=63, ieta=16-22
		gStyle->SetOptStat(0);
		gStyle->SetTitleFontSize(0.1);
		cc[0]->Divide(7,7,0,0);
		cc[1]=new TCanvas("cc2","cc2",4900,4900);//iphi=63, ieta=23-29
		cc[1]->Divide(7,7,0,0);
		cc[2]=new TCanvas("cc3","cc3",4900,3500);//iphi=64, ieta=16-20
		cc[2]->Divide(7,5,0,0);
		cc[3]=new TCanvas("cc4","cc4",4900,4900);//iphi=65, ieta=16-22
		cc[3]->Divide(7,7,0,0);
		cc[4]=new TCanvas("cc5","cc5",4900,4900);//iphi=65, ieta=23-29
		cc[4]->Divide(7,7,0,0);
		cc[5]=new TCanvas("cc6","cc6",4900,3500);//iphi=66, ieta=16-20
		cc[5]->Divide(7,5,0,0);
		int cci=0;
		
		for(int i2=0;i2<PL.size();i2++)
		{
			int ccind=-1; 
			if(PL[i2].iphi==63 && PL[i2].ieta>=16 && PL[i2].ieta<=22){ccind=0;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==63 && PL[i2].ieta>=23 && PL[i2].ieta<=29) {ccind=1;cci=(PL[i2].ieta-23)*7+PL[i2].depth;}
			else if(PL[i2].iphi==64) {ccind=2;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==65 && PL[i2].ieta>=16 && PL[i2].ieta<=22) {ccind=3;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==65 && PL[i2].ieta>=23 && PL[i2].ieta<=29) {ccind=4;cci=(PL[i2].ieta-23)*7+PL[i2].depth;}
			else if(PL[i2].iphi==66) {ccind=5;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			cc[ccind]->cd(cci);
			PL[i2].h->Draw();PL[i2].h->GetYaxis()->SetRangeUser(0.,HMAX);
			PL[i2].h->Write();
		}
		
		cc[0]->SaveAs("PedDists.pdf(");
		cc[1]->SaveAs("PedDists.pdf");
		cc[2]->SaveAs("PedDists.pdf");
		cc[3]->SaveAs("PedDists.pdf");
		cc[4]->SaveAs("PedDists.pdf");
		cc[5]->SaveAs("PedDists.pdf)");
		for(int i1=0;i1<6;i1++){delete cc[i1];}
	}
	{
		TCanvas* cc[6];
		cc[0]=new TCanvas("cc1","cc1",4900,4900);//iphi=63, ieta=16-22
		gStyle->SetOptStat(0);
		gStyle->SetTitleFontSize(0.1);
		cc[0]->Divide(7,7,0,0);
		cc[1]=new TCanvas("cc2","cc2",4900,4900);//iphi=63, ieta=23-29
		cc[1]->Divide(7,7,0,0);
		cc[2]=new TCanvas("cc3","cc3",4900,3500);//iphi=64, ieta=16-20
		cc[2]->Divide(7,5,0,0);
		cc[3]=new TCanvas("cc4","cc4",4900,4900);//iphi=65, ieta=16-22
		cc[3]->Divide(7,7,0,0);
		cc[4]=new TCanvas("cc5","cc5",4900,4900);//iphi=65, ieta=23-29
		cc[4]->Divide(7,7,0,0);
		cc[5]=new TCanvas("cc6","cc6",4900,3500);//iphi=66, ieta=16-20
		cc[5]->Divide(7,5,0,0);
		int cci=0;
		gPad->SetLogy();		

		for(int i2=0;i2<PL.size();i2++)
		{
			int ccind=-1; 
			if(PL[i2].iphi==63 && PL[i2].ieta>=16 && PL[i2].ieta<=22){ccind=0;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==63 && PL[i2].ieta>=23 && PL[i2].ieta<=29) {ccind=1;cci=(PL[i2].ieta-23)*7+PL[i2].depth;}
			else if(PL[i2].iphi==64) {ccind=2;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==65 && PL[i2].ieta>=16 && PL[i2].ieta<=22) {ccind=3;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==65 && PL[i2].ieta>=23 && PL[i2].ieta<=29) {ccind=4;cci=(PL[i2].ieta-23)*7+PL[i2].depth;}
			else if(PL[i2].iphi==66) {ccind=5;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			cc[ccind]->cd(cci);
			gPad->SetLogy();
			PL[i2].qtot->Draw();PL[i2].qtot->GetYaxis()->SetRangeUser(.1,QTYMAX);
			PL[i2].qtot->GetXaxis()->SetRangeUser(QTXMIN,QTXMAX);
			PL[i2].qtot->Write();
		}
		
		cc[0]->SaveAs("Qtot.pdf(");
		cc[1]->SaveAs("Qtot.pdf");
		cc[2]->SaveAs("Qtot.pdf");
		cc[3]->SaveAs("Qtot.pdf");
		cc[4]->SaveAs("Qtot.pdf");
		cc[5]->SaveAs("Qtot.pdf)");
		for(int i1=0;i1<6;i1++){delete cc[i1];}
	}
	outPeds->Close();
	
	ofstream outfile("Pedestals.txt");
	for(int i2=0;i2<PL.size();i2++)
	{
		outfile<<PL[i2].ieta<<" "<<PL[i2].iphi<<" "<<PL[i2].depth<<" "<<PL[i2].pedmean<<" "<<PL[i2].pedmeanerr<<" "<<PL[i2].pedsigma<<" "<<PL[i2].pedsigmaerr<<endl;
	}
	outfile.close();
	sprintf(hname,"mv Pedestals.txt ../Plots/%d",RunNo);system(hname);
	sprintf(hname,"mv *.pdf ../Plots/%d",RunNo);system(hname);
	sprintf(hname,"mv Pedestals.root ../Histos/%d",RunNo);system(hname);
	inroot->Close();
}
int plotleds()
{
	char hname[500];
	sprintf(hname,"../NTuples/N_%d.root",RunNo);
	TFile* inroot=new TFile(hname);
	TTree *tree = (TTree*)inroot->Get("Events");
	tree->SetBranchAddress("ieta",&ed.ieta);
	tree->SetBranchAddress("iphi",&ed.iphi);
	tree->SetBranchAddress("depth",&ed.depth);
	tree->SetBranchAddress("pulse",&ed.pulse);
	tree->SetBranchAddress("tdc",&ed.tdc);
	tree->SetBranchAddress("capid",&ed.capid);
	tree->SetBranchAddress("soi",&ed.soi);
	
	ledlist pl;vector <ledlist> PL;int plind=-1;
	tree->GetEntry(0);
	for(int i1=0;i1<ed.ieta->size();i1++)
	{
		plind=-1;
		for(int i2=0;i2<PL.size();i2++)
		{
			if(PL[i2].ieta==ed.ieta->at(i1) && PL[i2].iphi==ed.iphi->at(i1) && PL[i2].depth==ed.depth->at(i1))
			{
				plind=i2;break;
			}
		}
		if(plind==-1)
		{
			pl.ieta=ed.ieta->at(i1);
			pl.iphi=ed.iphi->at(i1);
			pl.depth=ed.depth->at(i1);
			pl.Qmean=0.;
			pl.Qsigma=0.;
			pl.Qmeanerr=0.;
			pl.Qsigmaerr=0.;
			pl.nbadcapid=0.;
			sprintf(hname,"Q Pulse %d %d %d",pl.ieta,pl.iphi,pl.depth);
			pl.p=new TH1F(hname,hname,10,-0.5,9.5);
			pl.p->SetLineWidth(2);pl.p->SetLineColor(4);pl.p->SetFillColor(4);
			pl.p->GetXaxis()->SetTitle("TS");pl.p->GetXaxis()->CenterTitle();
			pl.p->GetYaxis()->SetTitle("Mean ADC / TS");pl.p->GetYaxis()->CenterTitle();
			sprintf(hname,"Q Pulse %d %d %d norm",pl.ieta,pl.iphi,pl.depth);
			pl.pn=new TH1F(hname,hname,10,-0.5,9.5);
			sprintf(hname,"TDC %d %d %d",pl.ieta,pl.iphi,pl.depth);
			pl.t=new TH1F(hname,hname,10,-0.5,9.5);
			pl.t->SetLineWidth(2);pl.t->SetLineColor(2);
			pl.t->GetXaxis()->SetTitle("TS");pl.t->GetXaxis()->CenterTitle();
			pl.t->GetYaxis()->SetTitle("Mean TDC / TS");pl.t->GetYaxis()->CenterTitle();
			sprintf(hname,"Q %d %d %d",pl.ieta,pl.iphi,pl.depth);
			pl.q=new TH1F(hname,hname,5000,0.,1000000.);
			pl.q->SetLineWidth(2);pl.q->SetLineColor(1);
			pl.q->GetXaxis()->SetTitle("Charge (fC)");pl.q->GetXaxis()->CenterTitle();
			pl.q->GetYaxis()->SetTitle("Events / 200 fC");pl.q->GetYaxis()->CenterTitle();
			
			PL.push_back(pl);
		}
	}
	double ped=0.;double sig=0.;bool capidOK=true;
	for(int i=0;i<tree->GetEntries();i++)
// 	for(int i=0;i<10;i++)
	{
		tree->GetEntry(i);
		if(i%200==0) cout<<i<<" / "<<tree->GetEntries()<<endl;
		for(int i1=0;i1<ed.ieta->size();i1++)
		{
			plind=-1;
			for(int i2=0;i2<PL.size();i2++)
			{
				if(PL[i2].ieta==ed.ieta->at(i1) && PL[i2].iphi==ed.iphi->at(i1) && PL[i2].depth==ed.depth->at(i1))
				{
					plind=i2;break;
				}
			}
			ped=0.;sig=0.;capidOK=true;
			for(int i2=0;i2<10;i2++)
			{
				PL[plind].p->Fill(i2,adc2fC_QIE11[ed.pulse->at(i1)[i2]]);
				PL[plind].pn->Fill(i2);
				PL[plind].t->Fill(i2,ed.tdc->at(i1)[i2]);
				if(i2==1 || i2==9) ped+=adc2fC_QIE11[ed.pulse->at(i1)[i2]];
				if(i2>=2 && i2<=8) sig+=adc2fC_QIE11[ed.pulse->at(i1)[i2]];
				if(i2>0)
				{
					if(!((ed.capid->at(i1)[i2]-ed.capid->at(i1)[i2-1])==1 || (ed.capid->at(i1)[i2]-ed.capid->at(i1)[i2-1])==-3))
					{
						capidOK=false;
					}
				}
			}
			ped/=2.;
			sig-=(ped*7.);
			PL[plind].q->Fill(sig);
			if(!capidOK) PL[plind].nbadcapid+=1.;
		}
	}
	double HMAX=0.;double HMAXX=0.;double HMIN=1000000.;
	double HMAXP=0.;double HMAXT=0.;
	int iside=0;int iquad=0;
	TF1* tf=new TF1("gaus","gaus",0.,1000000.);
	float fnevt=((float)tree->GetEntries());
	for(int i2=0;i2<PL.size();i2++)
	{
		PL[i2].p->Divide(PL[i2].pn);
		PL[i2].p->GetXaxis()->SetRangeUser(0.5,9.5);
		if(PL[i2].p->GetBinContent(PL[i2].p->GetMaximumBin())>HMAXP){HMAXP=PL[i2].p->GetBinContent(PL[i2].p->GetMaximumBin());}
		PL[i2].p->GetXaxis()->SetRangeUser(-0.5,9.5);
		PL[i2].t->Divide(PL[i2].pn);
		if(PL[i2].t->GetBinContent(PL[i2].t->GetMaximumBin())>HMAXT){HMAXT=PL[i2].t->GetBinContent(PL[i2].t->GetMaximumBin());}
		PL[i2].nbadcapid/=fnevt;
// 		if(PL[i2].h->GetMean()>0)
// 		{
			PL[i2].q->Fit(tf,"q","q",0.,1000000.);
			PL[i2].q->Fit(tf,"q","q",tf->GetParameter(1)-1.5*tf->GetParameter(2),tf->GetParameter(1)+1.5*tf->GetParameter(2));
			PL[i2].Qmean=tf->GetParameter(1);
			PL[i2].Qsigma=tf->GetParameter(2);
			PL[i2].Qmeanerr=tf->GetParError(1);
			PL[i2].Qsigmaerr=tf->GetParError(2);
			tf->SetLineStyle(2);
			if(PL[i2].q->GetBinContent(PL[i2].q->GetMaximumBin())>HMAX){HMAX=PL[i2].q->GetBinContent(PL[i2].q->GetMaximumBin());}
			if(PL[i2].q->GetBinCenter(PL[i2].q->FindLastBinAbove(0.))>HMAXX){HMAXX=PL[i2].q->GetBinCenter(PL[i2].q->FindLastBinAbove(0.));}
// 			PL[i2].q->GetXaxis()->SetRangeUser(10000.,1000000.);
			if(PL[i2].q->GetBinCenter(PL[i2].q->FindFirstBinAbove(10.))<HMIN){HMIN=PL[i2].q->GetBinCenter(PL[i2].q->FindFirstBinAbove(10.));}
// 		}
// 		else
// 		{
// 			PL[i2].pedmean=0.;
// 			PL[i2].pedsigma=0.;
// 		}
	}
	HMAX+=20.;HMAXP+=20.;
	HMAXT+=15.;HMAXX+=100.;HMIN-=100.;
	sprintf(hname,"LEDLaser.root");
// 	if(RunType==2){sprintf(hname,"LED.root");}
// 	else if(RunType==3){sprintf(hname,"Laser.root");}
	TFile* outPeds=new TFile(hname,"recreate");
	{
		TCanvas* cc[6];
		cc[0]=new TCanvas("cc1","cc1",4900,4900);//iphi=63, ieta=16-22
		gStyle->SetOptStat(0);
		gStyle->SetTitleFontSize(0.1);
		cc[0]->Divide(7,7,0,0);
		cc[1]=new TCanvas("cc2","cc2",4900,4900);//iphi=63, ieta=23-29
		cc[1]->Divide(7,7,0,0);
		cc[2]=new TCanvas("cc3","cc3",4900,3500);//iphi=64, ieta=16-20
		cc[2]->Divide(7,5,0,0);
		cc[3]=new TCanvas("cc4","cc4",4900,4900);//iphi=65, ieta=16-22
		cc[3]->Divide(7,7,0,0);
		cc[4]=new TCanvas("cc5","cc5",4900,4900);//iphi=65, ieta=23-29
		cc[4]->Divide(7,7,0,0);
		cc[5]=new TCanvas("cc6","cc6",4900,3500);//iphi=66, ieta=16-20
		cc[5]->Divide(7,5,0,0);
		int cci=0;
		
		for(int i2=0;i2<PL.size();i2++)
		{
			int ccind=-1; 
			if(PL[i2].iphi==63 && PL[i2].ieta>=16 && PL[i2].ieta<=22){ccind=0;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==63 && PL[i2].ieta>=23 && PL[i2].ieta<=29) {ccind=1;cci=(PL[i2].ieta-23)*7+PL[i2].depth;}
			else if(PL[i2].iphi==64) {ccind=2;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==65 && PL[i2].ieta>=16 && PL[i2].ieta<=22) {ccind=3;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==65 && PL[i2].ieta>=23 && PL[i2].ieta<=29) {ccind=4;cci=(PL[i2].ieta-23)*7+PL[i2].depth;}
			else if(PL[i2].iphi==66) {ccind=5;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			cc[ccind]->cd(cci);
			PL[i2].p->Draw("hist");
			PL[i2].p->GetYaxis()->SetRangeUser(0.,HMAXP);
// 			PL[i2].p->GetYaxis()->SetRangeUser(0.,20.);
			PL[i2].p->Write();
		}
		
		cc[0]->SaveAs("PulseShapes.pdf(");
		cc[1]->SaveAs("PulseShapes.pdf");
		cc[2]->SaveAs("PulseShapes.pdf");
		cc[3]->SaveAs("PulseShapes.pdf");
		cc[4]->SaveAs("PulseShapes.pdf");
		cc[5]->SaveAs("PulseShapes.pdf)");
		for(int i1=0;i1<6;i1++){delete cc[i1];}
	}
	{
		TCanvas* cc[6];
		cc[0]=new TCanvas("cc1","cc1",4900,4900);//iphi=63, ieta=16-22
		gStyle->SetOptStat(0);
		gStyle->SetTitleFontSize(0.1);
		cc[0]->Divide(7,7,0,0);
		cc[1]=new TCanvas("cc2","cc2",4900,4900);//iphi=63, ieta=23-29
		cc[1]->Divide(7,7,0,0);
		cc[2]=new TCanvas("cc3","cc3",4900,3500);//iphi=64, ieta=16-20
		cc[2]->Divide(7,5,0,0);
		cc[3]=new TCanvas("cc4","cc4",4900,4900);//iphi=65, ieta=16-22
		cc[3]->Divide(7,7,0,0);
		cc[4]=new TCanvas("cc5","cc5",4900,4900);//iphi=65, ieta=23-29
		cc[4]->Divide(7,7,0,0);
		cc[5]=new TCanvas("cc6","cc6",4900,3500);//iphi=66, ieta=16-20
		cc[5]->Divide(7,5,0,0);
		int cci=0;
		
		for(int i2=0;i2<PL.size();i2++)
		{
			int ccind=-1; 
			if(PL[i2].iphi==63 && PL[i2].ieta>=16 && PL[i2].ieta<=22){ccind=0;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==63 && PL[i2].ieta>=23 && PL[i2].ieta<=29) {ccind=1;cci=(PL[i2].ieta-23)*7+PL[i2].depth;}
			else if(PL[i2].iphi==64) {ccind=2;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==65 && PL[i2].ieta>=16 && PL[i2].ieta<=22) {ccind=3;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==65 && PL[i2].ieta>=23 && PL[i2].ieta<=29) {ccind=4;cci=(PL[i2].ieta-23)*7+PL[i2].depth;}
			else if(PL[i2].iphi==66) {ccind=5;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			cc[ccind]->cd(cci);
			PL[i2].q->Draw();PL[i2].q->GetYaxis()->SetRangeUser(0.,HMAX);PL[i2].q->GetXaxis()->SetRangeUser(HMIN,HMAXX);
			PL[i2].q->Write();
		}
		if(RunType==2)
		{
			cc[0]->SaveAs("QDists.pdf(");
			cc[1]->SaveAs("QDists.pdf");
			cc[2]->SaveAs("QDists.pdf");
			cc[3]->SaveAs("QDists.pdf");
			cc[4]->SaveAs("QDists.pdf");
			cc[5]->SaveAs("QDists.pdf)");
		}
		else if(RunType==3)
		{
			cc[0]->SaveAs("LQDists.pdf(");
			cc[1]->SaveAs("LQDists.pdf");
			cc[2]->SaveAs("LQDists.pdf");
			cc[3]->SaveAs("LQDists.pdf");
			cc[4]->SaveAs("LQDists.pdf");
			cc[5]->SaveAs("LQDists.pdf)");
		}
		for(int i1=0;i1<6;i1++){delete cc[i1];}
	}
	{
		TCanvas* cc[6];
		cc[0]=new TCanvas("cc1","cc1",4900,4900);//iphi=63, ieta=16-22
		gStyle->SetOptStat(0);
		gStyle->SetTitleFontSize(0.1);
		cc[0]->Divide(7,7,0,0);
		cc[1]=new TCanvas("cc2","cc2",4900,4900);//iphi=63, ieta=23-29
		cc[1]->Divide(7,7,0,0);
		cc[2]=new TCanvas("cc3","cc3",4900,3500);//iphi=64, ieta=16-20
		cc[2]->Divide(7,5,0,0);
		cc[3]=new TCanvas("cc4","cc4",4900,4900);//iphi=65, ieta=16-22
		cc[3]->Divide(7,7,0,0);
		cc[4]=new TCanvas("cc5","cc5",4900,4900);//iphi=65, ieta=23-29
		cc[4]->Divide(7,7,0,0);
		cc[5]=new TCanvas("cc6","cc6",4900,3500);//iphi=66, ieta=16-20
		cc[5]->Divide(7,5,0,0);
		int cci=0;
		
		for(int i2=0;i2<PL.size();i2++)
		{
			int ccind=-1; 
			if(PL[i2].iphi==63 && PL[i2].ieta>=16 && PL[i2].ieta<=22){ccind=0;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==63 && PL[i2].ieta>=23 && PL[i2].ieta<=29) {ccind=1;cci=(PL[i2].ieta-23)*7+PL[i2].depth;}
			else if(PL[i2].iphi==64) {ccind=2;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==65 && PL[i2].ieta>=16 && PL[i2].ieta<=22) {ccind=3;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			else if(PL[i2].iphi==65 && PL[i2].ieta>=23 && PL[i2].ieta<=29) {ccind=4;cci=(PL[i2].ieta-23)*7+PL[i2].depth;}
			else if(PL[i2].iphi==66) {ccind=5;cci=(PL[i2].ieta-16)*7+PL[i2].depth;}
			cc[ccind]->cd(cci);
			PL[i2].t->Draw();
			PL[i2].t->GetYaxis()->SetRangeUser(0.,80.);
			PL[i2].t->Write();
		}
		
		cc[0]->SaveAs("TDCShapes.pdf(");
		cc[1]->SaveAs("TDCShapes.pdf");
		cc[2]->SaveAs("TDCShapes.pdf");
		cc[3]->SaveAs("TDCShapes.pdf");
		cc[4]->SaveAs("TDCShapes.pdf");
		cc[5]->SaveAs("TDCShapes.pdf)");
		for(int i1=0;i1<6;i1++){delete cc[i1];}
	}
	outPeds->Close();
	ofstream outfile("LEDLaser.txt");
	for(int i2=0;i2<PL.size();i2++)
	{
		outfile<<PL[i2].ieta<<" "<<PL[i2].iphi<<" "<<PL[i2].depth<<" "<<PL[i2].Qmean<<" "<<PL[i2].Qmeanerr<<" "<<PL[i2].Qsigma<<" "<<PL[i2].Qsigmaerr<<" "<<PL[i2].nbadcapid<<endl;
	}
	outfile.close();
	sprintf(hname,"mv LEDLaser.txt ../Plots/%d",RunNo);system(hname);
	sprintf(hname,"mv *.pdf ../Plots/%d",RunNo);system(hname);
	sprintf(hname,"mv LEDLaser.root ../Histos/%d",RunNo);system(hname);
	inroot->Close();
}





// 
// int plotleds()
// {
// 	char hname[500];
// 	sprintf(hname,"../NTuples/N_%d.root",RunNo);
// 	TFile* inroot=new TFile(hname);
// 	TTree *tree = (TTree*)inroot->Get("Events");
// 	tree->SetBranchAddress("ieta",&ed.ieta);
// 	tree->SetBranchAddress("iphi",&ed.iphi);
// 	tree->SetBranchAddress("depth",&ed.depth);
// 	tree->SetBranchAddress("pulse",&ed.pulse);
// 	tree->SetBranchAddress("tdc",&ed.tdc);
// 	tree->SetBranchAddress("capid",&ed.capid);
// 	
// 	ledlist pl;vector <ledlist> PL;int plind=-1;
// 	tree->GetEntry(0);
// 	for(int i1=0;i1<ed.ieta->size();i1++)
// 	{
// 		plind=-1;
// 		for(int i2=0;i2<PL.size();i2++)
// 		{
// 			if(SEM[PL[i2].mapind].ieta==ed.ieta->at(i1) && SEM[PL[i2].mapind].iphi==ed.iphi->at(i1) && SEM[PL[i2].mapind].depth==ed.depth->at(i1))
// 			{
// 				plind=i2;break;
// 			}
// 		}
// 		if(plind==-1)
// 		{
// 			pl.mapind=-1;
// 			for(int i2=0;i2<SEM.size();i2++)
// 			{
// 				if(SEM[i2].ieta==ed.ieta->at(i1) && SEM[i2].iphi==ed.iphi->at(i1) && SEM[i2].depth==ed.depth->at(i1))
// 				{
// 					pl.mapind=i2;break;
// 				}
// 			}
// 			pl.Qmean=0.;
// 			pl.Qsigma=0.;
// 			pl.npemean=0.;
// 			pl.npeerr=0.;
// 			pl.gain=0.;
// 			pl.gainerr=0.;
// 			pl.nbadcapid=0.;
// 			sprintf(hname,"Pulse %d %s %s",SEM[pl.mapind].box,SEM[pl.mapind].PMTname.c_str(),SEM[pl.mapind].Channelname.c_str());
// 			pl.p=new TH1F(hname,hname,10,-0.5,9.5);
// 			pl.p->GetXaxis()->SetTitle("TS (x25 ns)");pl.p->GetXaxis()->CenterTitle();
// 			pl.p->GetYaxis()->SetTitle("Mean Charge per TS (fC)");pl.p->GetYaxis()->CenterTitle();
// 			pl.p->SetFillColor(4);pl.p->SetLineColor(4);
// 			sprintf(hname,"PulseNorm %d %s %s",SEM[pl.mapind].box,SEM[pl.mapind].PMTname.c_str(),SEM[pl.mapind].Channelname.c_str());
// 			pl.pn=new TH1F(hname,hname,10,-0.5,9.5);
// 			
// 			sprintf(hname,"TDC %d %s %s",SEM[pl.mapind].box,SEM[pl.mapind].PMTname.c_str(),SEM[pl.mapind].Channelname.c_str());
// 			pl.t=new TH1F(hname,hname,10,-0.5,9.5);
// 			pl.t->GetXaxis()->SetTitle("TS (x25 ns)");pl.t->GetXaxis()->CenterTitle();
// 			pl.t->GetYaxis()->SetTitle("Mean TDC per TS");pl.t->GetYaxis()->CenterTitle();
// 			pl.t->SetLineColor(2);
// 			
// 			sprintf(hname,"Charge %d %s %s",SEM[pl.mapind].box,SEM[pl.mapind].PMTname.c_str(),SEM[pl.mapind].Channelname.c_str());
// 			pl.q=new TH1F(hname,hname,10000,0.,250000.);
// 			pl.q->GetXaxis()->SetTitle("Charge (fC)");pl.q->GetXaxis()->CenterTitle();
// 			pl.q->GetYaxis()->SetTitle("Entries / 25 fC");pl.q->GetYaxis()->CenterTitle();
// 			pl.q->SetLineColor(1);pl.q->SetLineWidth(2);
// 			PL.push_back(pl);
// 		}
// 	}
// 	double ped=0.;double sig=0.;bool capidOK=true;int cap0=0;
// 	for(int i=0;i<tree->GetEntries();i++)
// 	for(int i=0;i<10;i++)
// 	{
// 		tree->GetEntry(i);
// 		if(i%200==0) cout<<i<<" / "<<tree->GetEntries()<<endl;
// 		for(int i1=0;i1<ed.ieta->size();i1++)
// 		{
// 			plind=-1;
// 			for(int i2=0;i2<PL.size();i2++)
// 			{
// 				if(SEM[PL[i2].mapind].ieta==ed.ieta->at(i1) && SEM[PL[i2].mapind].iphi==ed.iphi->at(i1) && SEM[PL[i2].mapind].depth==ed.depth->at(i1))
// 				{
// 					plind=i2;break;
// 				}
// 			}
// 			ped=0.;sig=0.;capidOK=true;
// 			for(int i2=0;i2<10;i2++)
// 			{
// 				cout<<i2<<" "<<ed.pulse->at(i1)[i2]<<" "<<adc2fC_QIE10[ed.pulse->at(i1)[i2]]<<endl;
// 				PL[plind].p->Fill(i2,adc2fC_QIE10[ed.pulse->at(i1)[i2]]);
// 				PL[plind].t->Fill(i2,ed.tdc->at(i1)[i2]);
// 				PL[plind].pn->Fill(i2);
// 				if(i2<=2) ped+=adc2fC_QIE10[ed.pulse->at(i1)[i2]];
// 				if(i2>=3 && i2<=7) sig+=adc2fC_QIE10[ed.pulse->at(i1)[i2]];
// 				if(i2>0)
// 				{
// 					if(!((ed.capid->at(i1)[i2]-ed.capid->at(i1)[i2-1])==1 || (ed.capid->at(i1)[i2]-ed.capid->at(i1)[i2-1])==-3))
// 					{
// 						capidOK=false;
// 					}
// 				}
// 			}
// 			ped/=3.;
// 			sig-=(5.*ped);
// 			PL[plind].q->Fill(sig);
// 			if(!capidOK) PL[plind].nbadcapid+=1.;
// 		}
// 	}
// 	double Qxmax[2][4]={{0.}};double Qymax[2][4]={{0.}};
// 	double Pmax[2][4]={{0.}};
// 	int iside=0;int iquad=0;float fnevt=((float)tree->GetEntries());
// 	TF1* tf=new TF1("gaus","gaus",0.,250000.);
// 	for(int i2=0;i2<PL.size();i2++)
// 	{
// 		iside=(SEM[PL[i2].mapind].box>0?0:1);
// 		iquad=(abs(SEM[PL[i2].mapind].box)-1)/9;
// 		PL[i2].p->Divide(PL[i2].pn);
// 		PL[i2].t->Divide(PL[i2].pn);
// 		PL[i2].nbadcapid/=fnevt;
// 		if(PL[i2].p->GetBinContent(PL[i2].p->GetMaximumBin())>Pmax[iside][iquad]){Pmax[iside][iquad]=PL[i2].p->GetBinContent(PL[i2].p->GetMaximumBin());}
// 		if(PL[i2].q->GetMean()>5)
// 		{
// 			PL[i2].q->Fit(tf,"q","q",0.,250000.);
// 			PL[i2].q->Fit(tf,"q","q",tf->GetParameter(1)-1.5*tf->GetParameter(2),tf->GetParameter(1)+1.5*tf->GetParameter(2));
// 			PL[i2].Qmean=tf->GetParameter(1);
// 			PL[i2].Qsigma=tf->GetParameter(2);
// 			PL[i2].npemean=1.15*pow(tf->GetParameter(1)/tf->GetParameter(2),2.);
// 			PL[i2].npeerr=2.*PL[i2].npemean*sqrt(pow(tf->GetParError(1)/tf->GetParameter(1),2.)+pow(tf->GetParError(2)/tf->GetParameter(2),2.));
// 			PL[i2].gain=((pow(tf->GetParameter(2),2.)/(tf->GetParameter(1)*1.15)));
// 			PL[i2].gainerr=sqrt(pow(tf->GetParError(1),2.)+pow(sqrt(2.)*tf->GetParError(2),2.));
// 			tf->SetLineStyle(2);
// 			if(PL[i2].q->GetBinContent(PL[i2].q->GetMaximumBin())>Qymax[iside][iquad]){Qymax[iside][iquad]=PL[i2].q->GetBinContent(PL[i2].q->GetMaximumBin());}
// 			if(PL[i2].q->GetBinCenter(PL[i2].q->FindLastBinAbove(0))>Qxmax[iside][iquad]){Qxmax[iside][iquad]=PL[i2].q->GetBinCenter(PL[i2].q->FindLastBinAbove(0));}
// 			if(PL[i2].p->GetBinContent(PL[i2].p->GetMaximumBin())>1000) cout<<SEM[PL[i2].mapind].boxname<<" "<<SEM[PL[i2].mapind].PMTname<<" "<<PL[i2].p->GetBinContent(PL[i2].p->GetMaximumBin())<<endl;
// 		}
// 	}
// 	for(int i1=0;i1<2;i1++)
// 	{
// 		for(int i2=0;i2<4;i2++)
// 		{
// 			Pmax[i1][i2]+=100.;Qxmax[i1][i2]+=50.;Qymax[i1][i2]+=20.;
// 		}
// 	}
// 	{
// 		TFile* outLeds=new TFile("LEDDetails.root","recreate");
// 		char cname[400];
// 		char cname1[400];
// 		char dname[400];
// 		char dname1[400];
// 		char ename[400];
// 		char ename1[400];
// 		char bname[100];
// 		string pmtnames[24]={"A2","A4","A6","A8","A1","A3","A5","A7","B2","B4","B6","B8","B1","B3","B5","B7","C2","C4","C6","C8","C1","C3","C5","C7"};
// 		string chnames[2]={"1+2","3+4"};
// 		string side[2]={"P","M"};
// 		for(int i1=0;i1<2;i1++)
// 		{
// 			for(int i2=0;i2<4;i2++)
// 			{
// 				sprintf(cname,"HF%s_Q%d_PulseShapes.pdf",side[i1].c_str(),(i2+1));
// 				sprintf(dname,"HF%s_Q%d_Charges.pdf",side[i1].c_str(),(i2+1));
// 				sprintf(ename,"HF%s_Q%d_TDCShapes.pdf",side[i1].c_str(),(i2+1));
// 				for(int i3=0;i3<9;i3++)
// 				{
// 					TCanvas* cc1=new TCanvas("cc1","cc1",4500,6000);
// 					TCanvas* cc2=new TCanvas("cc2","cc2",4500,6000);
// 					TCanvas* cc3=new TCanvas("cc3","cc3",4500,6000);
// 					gStyle->SetOptStat(0);
// 					gStyle->SetTitleFontSize(0.1);
// 					cc1->Divide(8,6,0,0);
// 					cc2->Divide(8,6,0,0);
// 					cc3->Divide(8,6,0,0);
// 					int cc1i=1;int cc2i=1;int cc3i=1;
// 					if((i2*9+i3+1)<10){sprintf(bname,"HF%s_0%d_Q%d",side[i1].c_str(),(i2*9+i3+1),(i2+1));}
// 					else{sprintf(bname,"HF%s_%d_Q%d",side[i1].c_str(),(i2*9+i3+1),(i2+1));}
// 					string bname1(bname);
// 					for(int is1=0;is1<24;is1++)
// 					{
// 						for(int is2=0;is2<2;is2++)
// 						{
// 							bool found=false;plind=-1;
// 							for(int ik1=0;ik1<PL.size();ik1++)
// 							{
// 								if(SEM[PL[ik1].mapind].PMTname==pmtnames[is1] && SEM[PL[ik1].mapind].Channelname==chnames[is2] && SEM[PL[ik1].mapind].boxname==bname1)
// 								{
// 									plind=ik1;break;
// 								}
// 							}
// 							
// 							iside=(SEM[PL[plind].mapind].box>0?0:1);
// 							iquad=(abs(SEM[PL[plind].mapind].box)-1)/9;
// 							
// 							cc1->cd(cc1i);
// 							PL[plind].p->Divide(PL[plind].pn);
// 							PL[plind].p->Draw("hist");
// 							PL[plind].p->GetYaxis()->SetRangeUser(-10.,Pmax[iside][iquad]);
// 							cc1i++;
// 							PL[plind].p->Write();
// 							
// 							cc2->cd(cc2i);
// 							PL[plind].q->Draw();
// 							PL[plind].q->GetYaxis()->SetRangeUser(0.,Qymax[iside][iquad]);PL[plind].q->GetXaxis()->SetRangeUser(0.,Qxmax[iside][iquad]);
// 							cc2i++;
// 							PL[plind].q->Write();
// 							
// 							cc3->cd(cc3i);
// 							PL[plind].t->Draw();
// 							PL[plind].t->GetYaxis()->SetRangeUser(0.,75.);
// 							cc3i++;
// 							PL[plind].t->Write();
// 						}
// 					}
// 					if(i3==0) sprintf(cname1,"%s(",cname);
// 					else if(i3<8) sprintf(cname1,"%s",cname);
// 					else sprintf(cname1,"%s)",cname);
// 					cc1->SaveAs(cname1);
// 					delete cc1;
// 					
// 					if(i3==0) sprintf(dname1,"%s(",dname);
// 					else if(i3<8) sprintf(dname1,"%s",dname);
// 					else sprintf(dname1,"%s)",dname);
// 					cc2->SaveAs(dname1);
// 					delete cc2;
// 					
// 					if(i3==0) sprintf(ename1,"%s(",ename);
// 					else if(i3<8) sprintf(ename1,"%s",ename);
// 					else sprintf(ename1,"%s)",ename);
// 					cc3->SaveAs(ename1);
// 					delete cc3;
// 				}
// 				sprintf(hname,"mv %s ../Plots/%d",cname,RunNo);system(hname);
// 				sprintf(hname,"mv %s ../Plots/%d",dname,RunNo);system(hname);
// 				sprintf(hname,"mv %s ../Plots/%d",ename,RunNo);system(hname);
// 			}
// 		}
// 		outLeds->Close();
// 		for(int ik1=0;ik1<PL.size();ik1++)
// 		{
// 			PL[ik1].p=0;PL[ik1].pn=0;PL[ik1].q=0;PL[ik1].t=0;
// 		}
// 		sprintf(hname,"mv LEDDetails.root ../Histos/%d",RunNo);system(hname);
// 	}
// 	ofstream outfile("log.txt");
// 	outfile<<"Run : "<<RunNo<<endl;
// 	{
// 		TFile* outGraphs=new TFile("LEDSummaries.root","recreate");
// 		char cname[400];
// 		char cname1[400];
// 		char bname[100];
// 		TGraphErrors* tg1[2][4][9][2];
// 		TGraphErrors* tg2[2][4][9][2];
// 		TGraphErrors* tg3[2][4][9][2];
// 		string pmtnames[24]={"A2","A4","A6","A8","A1","A3","A5","A7","B2","B4","B6","B8","B1","B3","B5","B7","C2","C4","C6","C8","C1","C3","C5","C7"};
// 		string chnames[2]={"1+2","3+4"};
// 		string side[2]={"P","M"};
// 		int np=0;
// 		double gdiff=0.;double gdifferr=0.;
// 		for(int i1=0;i1<2;i1++)
// 		{
// 			for(int i2=0;i2<4;i2++)
// 			{
// 				sprintf(cname,"HF%s_Q%d_LEDs.pdf",side[i1].c_str(),(i2+1));
// 				TCanvas* cc1=new TCanvas("cc1","cc1",900,900);
// 				TCanvas* cc2=new TCanvas("cc2","cc2",900,900);
// 				TCanvas* cc3=new TCanvas("cc3","cc3",900,900);
// 				gStyle->SetOptStat(0);
// 				cc1->Divide(3,3,0,0);
// 				cc2->Divide(3,3,0,0);
// 				cc3->Divide(3,3,0,0);
// 				int cc1i=1;int cc2i=1;int cc3i=1;
// 				double ymin1=1000.;double ymax1=0.;
// 				double ymin2=1000.;double ymax2=0.;
// 				double ymin3=1000.;double ymax3=0.;
// 				for(int i3=0;i3<9;i3++)
// 				{
// 					for(int i4=0;i4<2;i4++)
// 					{
// 						sprintf(hname,"HF%s Q%d %d %s NPE",side[i1].c_str(),i2+1,(i2*9+i3+1),chnames[i4].c_str());
// 						tg1[i1][i2][i3][i4]=new TGraphErrors();
// 						tg1[i1][i2][i3][i4]->SetName(hname);tg1[i1][i2][i3][i4]->SetTitle(hname);
// 						tg1[i1][i2][i3][i4]->SetMarkerColor(i4+1);tg1[i1][i2][i3][i4]->SetMarkerStyle(i4+24);tg1[i1][i2][i3][i4]->SetLineColor(i4+1);
// 						tg1[i1][i2][i3][i4]->GetXaxis()->SetTitle("PMT ID");tg1[i1][i2][i3][i4]->GetXaxis()->CenterTitle();
// 						tg1[i1][i2][i3][i4]->GetYaxis()->SetTitle("NPE");tg1[i1][i2][i3][i4]->GetYaxis()->CenterTitle();
// 						
// 						sprintf(hname,"HF%s Q%d %d %s Gains",side[i1].c_str(),i2+1,(i2*9+i3+1),chnames[i4].c_str());
// 						tg2[i1][i2][i3][i4]=new TGraphErrors();
// 						tg2[i1][i2][i3][i4]->SetName(hname);tg2[i1][i2][i3][i4]->SetTitle(hname);
// 						tg2[i1][i2][i3][i4]->SetMarkerColor(i4+1);tg2[i1][i2][i3][i4]->SetMarkerStyle(i4+24);tg2[i1][i2][i3][i4]->SetLineColor(i4+1);
// 						tg2[i1][i2][i3][i4]->GetXaxis()->SetTitle("PMT ID");tg2[i1][i2][i3][i4]->GetXaxis()->CenterTitle();
// 						tg2[i1][i2][i3][i4]->GetYaxis()->SetTitle("Gain (fC)");tg2[i1][i2][i3][i4]->GetYaxis()->CenterTitle();
// 						
// 						sprintf(hname,"HF%s Q%d %d %s GainDiff",side[i1].c_str(),i2+1,(i2*9+i3+1),chnames[i4].c_str());
// 						tg3[i1][i2][i3][i4]=new TGraphErrors();
// 						tg3[i1][i2][i3][i4]->SetName(hname);tg3[i1][i2][i3][i4]->SetTitle(hname);
// 						tg3[i1][i2][i3][i4]->SetMarkerColor(i4+1);tg3[i1][i2][i3][i4]->SetMarkerStyle(i4+24);tg3[i1][i2][i3][i4]->SetLineColor(i4+1);
// 						tg3[i1][i2][i3][i4]->GetXaxis()->SetTitle("PMT ID");tg3[i1][i2][i3][i4]->GetXaxis()->CenterTitle();
// 						tg3[i1][i2][i3][i4]->GetYaxis()->SetTitle("(G_{P5}-G_{SX5})/G_{SX5}");tg3[i1][i2][i3][i4]->GetYaxis()->CenterTitle();
// 					
// 						if((i2*9+i3+1)<10){sprintf(bname,"HF%s_0%d_Q%d",side[i1].c_str(),(i2*9+i3+1),(i2+1));}
// 						else{sprintf(bname,"HF%s_%d_Q%d",side[i1].c_str(),(i2*9+i3+1),(i2+1));}
// 						string bname1(bname);
// 						
// 						for(int is1=0;is1<24;is1++)
// 						{
// 							bool found=false;plind=-1;
// 							for(int ik1=0;ik1<PL.size();ik1++)
// 							{
// 								if(SEM[PL[ik1].mapind].PMTname==pmtnames[is1] && SEM[PL[ik1].mapind].Channelname==chnames[i4] && SEM[PL[ik1].mapind].boxname==bname1)
// 								{
// 									plind=ik1;break;
// 								}
// 							}
// 							tg1[i1][i2][i3][i4]->SetPoint(is1+1,((double)is1+1),PL[plind].npemean);
// 							tg1[i1][i2][i3][i4]->SetPointError(is1+1,0.,PL[plind].npeerr);
// 							if(PL[plind].npemean>ymax1) ymax1=PL[plind].npemean;
// 							if(PL[plind].npemean<ymin1) ymin1=PL[plind].npemean;
// 							
// 							tg2[i1][i2][i3][i4]->SetPoint(is1+1,((double)is1+1),PL[plind].gain);
// 							tg2[i1][i2][i3][i4]->SetPointError(is1+1,0.,PL[plind].gainerr);
// 							if(PL[plind].gain>ymax2) ymax2=PL[plind].gain;
// 							if(PL[plind].gain<ymin2) ymin2=PL[plind].gain;
// 							
// 							gdiff=(PL[plind].gain-SEM[PL[plind].mapind].SX5gain)/SEM[PL[plind].mapind].SX5gain;
// 							gdifferr=gdiff*sqrt(pow(PL[plind].gainerr/PL[plind].gain,2.)+pow(SEM[PL[plind].mapind].SX5gainerr/SEM[PL[plind].mapind].SX5gain,2.));
// 							tg3[i1][i2][i3][i4]->SetPoint(is1+1,((double)is1+1),gdiff);
// 							tg3[i1][i2][i3][i4]->SetPointError(is1+1,0.,gdifferr);
// 							if(gdiff>ymax3) ymax3=gdiff;
// 							if(gdiff<ymin3) ymin3=gdiff;
// 							
// 							if(PL[plind].npemean>35.)
// 							{
// 								outfile<<"NPE too high: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
// 							}
// 							if(PL[plind].npemean<5.)
// 							{
// 								outfile<<"NPE too low: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
// 							}
// 							if(PL[plind].gain>35.)
// 							{
// 								outfile<<"Gain too high: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
// 							}
// 							if(PL[plind].gain<15.)
// 							{
// 								outfile<<"Gain too low: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
// 							}
// 							if(fabs(gdiff)>0.2)
// 							{
// 								outfile<<"Gain difference too high: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
// 							}
// 						}
// 					}
// 				}
// 				ymax1+=5.;ymin1-=5.;
// 				ymax2+=5.;ymin2-=5.;
// 				for(int i3=0;i3<9;i3++)
// 				{
// 					cc1->cd(cc1i);
// 					sprintf(hname,"HF%s Q%d %d NPE",side[i1].c_str(),i2+1,(i2*9+i3+1));
// 					tg1[i1][i2][i3][0]->SetTitle(hname);
// 					tg1[i1][i2][i3][0]->Draw("AP");
// 					tg1[i1][i2][i3][0]->GetYaxis()->SetTitle("NPE");tg1[i1][i2][i3][0]->GetYaxis()->CenterTitle();
// 					tg1[i1][i2][i3][0]->GetXaxis()->SetTitle("PMT ID");tg1[i1][i2][i3][0]->GetXaxis()->CenterTitle();
// 					gPad->SetGridy(1);gPad->SetGridx(1);
// 					tg1[i1][i2][i3][0]->GetYaxis()->SetRangeUser(ymin1,ymax1);tg1[i1][i2][i3][0]->GetXaxis()->SetRangeUser(0.5,24.5);
// 					tg1[i1][i2][i3][1]->Draw("same P");
// 					cc1i++;
// 					tg1[i1][i2][i3][0]->Write();tg1[i1][i2][i3][1]->Write();
// 					
// 					cc2->cd(cc2i);
// 					sprintf(hname,"HF%s Q%d %d Gain",side[i1].c_str(),i2+1,(i2*9+i3+1));
// 					tg2[i1][i2][i3][0]->SetTitle(hname);
// 					tg2[i1][i2][i3][0]->Draw("AP");
// 					tg2[i1][i2][i3][0]->GetYaxis()->SetTitle("Gain (fC)");tg2[i1][i2][i3][0]->GetYaxis()->CenterTitle();
// 					tg2[i1][i2][i3][0]->GetXaxis()->SetTitle("PMT ID");tg2[i1][i2][i3][0]->GetXaxis()->CenterTitle();
// 					gPad->SetGridy(1);gPad->SetGridx(1);
// 					tg2[i1][i2][i3][0]->GetYaxis()->SetRangeUser(ymin2,ymax2);tg2[i1][i2][i3][0]->GetXaxis()->SetRangeUser(0.5,24.5);
// 					tg2[i1][i2][i3][1]->Draw("same P");
// 					cc2i++;
// 					tg2[i1][i2][i3][0]->Write();tg2[i1][i2][i3][1]->Write();
// 					
// 					cc3->cd(cc3i);
// 					sprintf(hname,"HF%s Q%d %d GainDiff",side[i1].c_str(),i2+1,(i2*9+i3+1));
// 					tg3[i1][i2][i3][0]->SetTitle(hname);
// 					tg3[i1][i2][i3][0]->Draw("AP");
// 					tg3[i1][i2][i3][0]->GetYaxis()->SetTitle("(G_{P5}-G_{SX5})/G_{SX5}");tg3[i1][i2][i3][0]->GetYaxis()->CenterTitle();
// 					tg3[i1][i2][i3][0]->GetXaxis()->SetTitle("PMT ID");tg3[i1][i2][i3][0]->GetXaxis()->CenterTitle();
// 					gPad->SetGridy(1);gPad->SetGridx(1);
// 					tg3[i1][i2][i3][0]->GetYaxis()->SetRangeUser(-1.,1.);tg3[i1][i2][i3][0]->GetXaxis()->SetRangeUser(0.5,24.5);
// 					tg3[i1][i2][i3][1]->Draw("same P");
// 					cc3i++;
// 					tg3[i1][i2][i3][0]->Write();tg3[i1][i2][i3][1]->Write();
// 				}
// 				sprintf(cname1,"%s(",cname);
// 				cc1->SaveAs(cname1);
// 				sprintf(cname1,"%s",cname);
// 				cc2->SaveAs(cname1);
// 				sprintf(cname1,"%s)",cname);
// 				cc3->SaveAs(cname1);
// 				delete cc1;delete cc2;delete cc3;
// 				sprintf(hname,"mv %s ../Plots/%d",cname,RunNo);system(hname);
// 			}
// 		}
// 		outGraphs->Close();
// 		sprintf(hname,"mv LEDSummaries.root ../Histos/%d",RunNo);system(hname);
// 	}
// 	{
// 		TFile* outGraphs=new TFile("CapIDSummaries.root","recreate");
// 		char cname[400];
// 		char cname1[400];
// 		char bname[100];
// 		TGraphErrors* tg1[2][4][9][2];
// 		string pmtnames[24]={"A2","A4","A6","A8","A1","A3","A5","A7","B2","B4","B6","B8","B1","B3","B5","B7","C2","C4","C6","C8","C1","C3","C5","C7"};
// 		string chnames[2]={"1+2","3+4"};
// 		string side[2]={"P","M"};
// 		int np=0;
// 		double gdiff=0.;double gdifferr=0.;
// 		for(int i1=0;i1<2;i1++)
// 		{
// 			for(int i2=0;i2<4;i2++)
// 			{
// 				sprintf(cname,"HF%s_Q%d_CapIDs.pdf",side[i1].c_str(),(i2+1));
// 				TCanvas* cc1=new TCanvas("cc1","cc1",900,900);
// 				gStyle->SetOptStat(0);
// 				cc1->Divide(3,3,0,0);
// 				int cc1i=1;
// 				for(int i3=0;i3<9;i3++)
// 				{
// 					for(int i4=0;i4<2;i4++)
// 					{
// 						sprintf(hname,"HF%s Q%d %d %s Bad CapID Fraction",side[i1].c_str(),i2+1,(i2*9+i3+1),chnames[i4].c_str());
// 						tg1[i1][i2][i3][i4]=new TGraphErrors();
// 						tg1[i1][i2][i3][i4]->SetName(hname);tg1[i1][i2][i3][i4]->SetTitle(hname);
// 						tg1[i1][i2][i3][i4]->SetMarkerColor(i4+1);tg1[i1][i2][i3][i4]->SetMarkerStyle(i4+24);tg1[i1][i2][i3][i4]->SetLineColor(i4+1);
// 						tg1[i1][i2][i3][i4]->GetXaxis()->SetTitle("PMT ID");tg1[i1][i2][i3][i4]->GetXaxis()->CenterTitle();
// 						tg1[i1][i2][i3][i4]->GetYaxis()->SetTitle("Bad CapID Fraction");tg1[i1][i2][i3][i4]->GetYaxis()->CenterTitle();
// 					
// 						if((i2*9+i3+1)<10){sprintf(bname,"HF%s_0%d_Q%d",side[i1].c_str(),(i2*9+i3+1),(i2+1));}
// 						else{sprintf(bname,"HF%s_%d_Q%d",side[i1].c_str(),(i2*9+i3+1),(i2+1));}
// 						string bname1(bname);
// 						
// 						for(int is1=0;is1<24;is1++)
// 						{
// 							bool found=false;plind=-1;
// 							for(int ik1=0;ik1<PL.size();ik1++)
// 							{
// 								if(SEM[PL[ik1].mapind].PMTname==pmtnames[is1] && SEM[PL[ik1].mapind].Channelname==chnames[i4] && SEM[PL[ik1].mapind].boxname==bname1)
// 								{
// 									plind=ik1;break;
// 								}
// 							}
// 							tg1[i1][i2][i3][i4]->SetPoint(is1+1,((double)is1+1),PL[plind].nbadcapid);
// 							tg1[i1][i2][i3][i4]->SetPointError(is1+1,0.,0.);
// 							
// 							if(PL[plind].nbadcapid>0.1)
// 							{
// 								outfile<<"CapID rotation failure too high: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
// 							}
// 						}
// 					}
// 				}
// 				for(int i3=0;i3<9;i3++)
// 				{
// 					cc1->cd(cc1i);
// 					sprintf(hname,"HF%s Q%d %d Bad CI Frac",side[i1].c_str(),i2+1,(i2*9+i3+1));
// 					tg1[i1][i2][i3][0]->SetTitle(hname);
// 					tg1[i1][i2][i3][0]->Draw("AP");
// 					tg1[i1][i2][i3][0]->GetYaxis()->SetTitle("Bad CapID Fraction");tg1[i1][i2][i3][0]->GetYaxis()->CenterTitle();
// 					tg1[i1][i2][i3][0]->GetXaxis()->SetTitle("PMT ID");tg1[i1][i2][i3][0]->GetXaxis()->CenterTitle();
// 					gPad->SetGridy(1);gPad->SetGridx(1);
// 					tg1[i1][i2][i3][0]->GetYaxis()->SetRangeUser(-0.1,1.1);tg1[i1][i2][i3][0]->GetXaxis()->SetRangeUser(0.5,24.5);
// 					tg1[i1][i2][i3][1]->Draw("same P");
// 					cc1i++;
// 					tg1[i1][i2][i3][0]->Write();tg1[i1][i2][i3][1]->Write();
// 				}
// 				sprintf(cname1,"%s",cname);
// 				cc1->SaveAs(cname1);
// 				delete cc1;
// 				sprintf(hname,"mv %s ../Plots/%d",cname,RunNo);system(hname);
// 			}
// 		}
// 		outGraphs->Close();
// 		sprintf(hname,"mv CapIDSummaries.root ../Histos/%d",RunNo);system(hname);
// 	}
// 	{
// 		string pnames[12][4]={{"A2_3+4","A4_3+4","A6_3+4","A8_3+4"},{"A2_1+2","A4_1+2","A6_1+2","A8_1+2"},{"A1_3+4","A3_3+4","A5_3+4","A7_3+4"},{"A1_1+2","A3_1+2","A5_1+2","A7_1+2"},{"B2_3+4","B4_3+4","B6_3+4","B8_3+4"},{"B2_1+2","B4_1+2","B6_1+2","B8_1+2"},{"B1_3+4","B3_3+4","B5_3+4","B7_3+4"},{"B1_1+2","B3_1+2","B5_1+2","B7_1+2"},{"C2_3+4","C4_3+4","C6_3+4","C8_3+4"},{"C2_1+2","C4_1+2","C6_1+2","C8_1+2"},{"C1_3+4","C3_3+4","C5_3+4","C7_3+4"},{"C1_1+2","C3_1+2","C5_1+2","C7_1+2"}};
// 		TFile *inrootP=new TFile("HFMG.root");
// 		TMultiGraph *mg;
// 		TFile* outPoly=new TFile("Poly.root","recreate");
// 		TCanvas* cc1=new TCanvas("cc1","cc1",600,600);
// 		int bin=0;
// 		gStyle->SetOptStat(0);
// 		gStyle->SetTitleFontSize(0.06);
// 		gStyle->SetPalette(kRainBow);
// 		TH2Poly *HFPNPE= new TH2Poly("HFP NPE","HFP NPE",-500,500,-500,500);
// 		TH2Poly *HFPGain= new TH2Poly("HFP Gains","HFP Gains",-500,500,-500,500);
// 		TH2Poly *HFPGainDiff= new TH2Poly("HFP (G_{P5}-G_{SX5})/G_{SX5}","HFP (G_{P5}-G_{SX5})/G_{SX5}",-500,500,-500,500);
// 		TH2Poly *HFMNPE= new TH2Poly("HFM NPE","HFM NPE",-500,500,-500,500);
// 		TH2Poly *HFMGain= new TH2Poly("HFM Gains","HFM Gains",-500,500,-500,500);
// 		TH2Poly *HFMGainDiff= new TH2Poly("HFM (G_{P5}-G_{SX5})/G_{SX5}","HFM (G_{P5}-G_{SX5})/G_{SX5}",-500,500,-500,500);
// 		inrootP->GetObject("HF",mg);
// 		TLine *line1 = new TLine(-500.,0.,500.,0.);line1->SetLineColor(1);line1->SetLineWidth(2);
// 		TLine *line2 = new TLine(0.,-500.,0.,500.);line2->SetLineColor(1);line2->SetLineWidth(2);
// 		for(int i1=0;i1<36;i1++)
// 		{
// 			for(int i2=0;i2<12;i2++)
// 			{
// 				for(int i3=0;i3<4;i3++)
// 				{
// 					sprintf(hname,"B%d_%s",(i1+1),pnames[i2][i3].c_str());
// 					bin=HFPNPE->AddBin(mg->GetListOfGraphs()->FindObject(hname));
// 					bin=HFPGain->AddBin(mg->GetListOfGraphs()->FindObject(hname));
// 					bin=HFPGainDiff->AddBin(mg->GetListOfGraphs()->FindObject(hname));
// 					bin=HFMNPE->AddBin(mg->GetListOfGraphs()->FindObject(hname));
// 					bin=HFMGain->AddBin(mg->GetListOfGraphs()->FindObject(hname));
// 					bin=HFMGainDiff->AddBin(mg->GetListOfGraphs()->FindObject(hname));
// 				}
// 			}
// 		}
// 		for(int i2=0;i2<PL.size();i2++)
// 		{
// 			if(SEM[PL[i2].mapind].box>0)
// 			{
// 				sprintf(hname,"B%d_%s_%s",SEM[PL[i2].mapind].box,SEM[PL[i2].mapind].PMTname.c_str(),SEM[PL[i2].mapind].Channelname.c_str());
// 				HFPNPE->Fill(hname,PL[i2].npemean);
// 				HFPGain->Fill(hname,PL[i2].gain);
// 				HFPGainDiff->Fill(hname,(PL[i2].gain-SEM[PL[i2].mapind].SX5gain)/SEM[PL[i2].mapind].SX5gain);
// 			}
// 			else
// 			{
// 				sprintf(hname,"B%d_%s_%s",HFMBoxMap[abs(SEM[PL[i2].mapind].box)],SEM[PL[i2].mapind].PMTname.c_str(),SEM[PL[i2].mapind].Channelname.c_str());
// 				HFMNPE->Fill(hname,PL[i2].npemean);
// 				HFMGain->Fill(hname,PL[i2].gain);
// 				HFMGainDiff->Fill(hname,(PL[i2].gain-SEM[PL[i2].mapind].SX5gain)/SEM[PL[i2].mapind].SX5gain);
// 			}
// 		}
// 		HFPNPE->SetMinimum(0.);HFPNPE->SetMaximum(30.);HFPNPE->Draw("colz");HFPNPE->GetXaxis()->SetLabelSize(0);HFPNPE->GetYaxis()->SetLabelSize(0);
// 		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFPNPE.png");
// 		HFPGain->SetMinimum(0.);HFPGain->SetMaximum(50.);HFPGain->Draw("colz");HFPGain->GetXaxis()->SetLabelSize(0);HFPGain->GetYaxis()->SetLabelSize(0);
// 		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFPGain.png");
// 		HFPGainDiff->SetMinimum(-1.);HFPGainDiff->SetMaximum(1.);HFPGainDiff->Draw("colz");HFPGainDiff->GetXaxis()->SetLabelSize(0);HFPGainDiff->GetYaxis()->SetLabelSize(0);
// 		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFPGainDiff.png");
// 		HFMNPE->SetMinimum(0.);HFMNPE->SetMaximum(30.);HFMNPE->Draw("colz");HFMNPE->GetXaxis()->SetLabelSize(0);HFMNPE->GetYaxis()->SetLabelSize(0);
// 		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFMNPE.png");
// 		HFMGain->SetMinimum(0.);HFMGain->SetMaximum(50.);HFMGain->Draw("colz");HFMGain->GetXaxis()->SetLabelSize(0);HFMGain->GetYaxis()->SetLabelSize(0);
// 		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFMGain.png");
// 		HFMGainDiff->SetMinimum(-1.);HFMGainDiff->SetMaximum(1.);HFMGainDiff->Draw("colz");HFMGainDiff->GetXaxis()->SetLabelSize(0);HFMGainDiff->GetYaxis()->SetLabelSize(0);
// 		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFMGainDiff.png");
// 		
// 		HFPNPE->Write();
// 		HFPGain->Write();
// 		HFPGainDiff->Write();
// 		HFMNPE->Write();
// 		HFMGain->Write();
// 		HFMGainDiff->Write();
// 		outPoly->Close();
// 		inrootP->Close();
// 		
// 		sprintf(hname,"mv Poly.root ../Histos/%d",RunNo);system(hname);
// 		sprintf(hname,"mv *.png ../Plots/%d",RunNo);system(hname);
// 	}
// 	outfile.close();
// 	sprintf(hname,"mv log.txt ../Plots/%d",RunNo);system(hname);
// 	inroot->Close();
// }

int main(int argc, char *argv[])
{
	int opt=atoi(argv[1]);
	RunNo=atoi(argv[2]);
	RunType=opt;
// 	RunType=atoi(argv[3]);
	
	if(opt==0)
	{
		print();
	}
	else if(opt==1)
	{
		plotpeds();
	}
	else if(opt==2 || opt==3)
	{
		plotleds();
	}
}

















