// -*- C++ -*-
//
// Package:    HCALCommissioning2017/HEP17Analyzer
// Class:      HEP17Analyzer
// 
/**\class HEP17Analyzer HEP17Analyzer.cc HCALCommissioning2017/HEP17Analyzer/plugins/HEP17Analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Burak Bilki
//         Created:  Mon, 27 Feb 2017 20:25:23 GMT
//
//


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "EventFilter/HcalRawToDigi/interface/HcalHTRData.h"
#include "EventFilter/HcalRawToDigi/interface/HcalDCCHeader.h"
#include "EventFilter/HcalRawToDigi/interface/HcalUnpacker.h"
#include "DataFormats/HcalDetId/interface/HcalOtherDetId.h"
#include "DataFormats/HcalDigi/interface/HcalQIESample.h"
#include "DataFormats/HcalDigi/interface/QIE11DataFrame.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalCalibDetId.h"
#include "EventFilter/HcalRawToDigi/interface/AMC13Header.h"
#include "EventFilter/HcalRawToDigi/interface/HcalUHTRData.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDHeader.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TProfile.h"
#include "TFile.h"
#include "TSystem.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TStyle.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

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

class HEP17Analyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
	public:
		explicit HEP17Analyzer(const edm::ParameterSet&);
		~HEP17Analyzer();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		TFile *_file;
		TTree *_tree;
		int histoFED;
		int EID;
		int numChannels;
		string outFileName;
		int runType;
		string digiCollection;
		vector <int> fedlist;
		
		edm::EDGetTokenT<HcalDataFrameContainer<QIE11DataFrame> > tok_QIE11DigiCollection_;
		edm::EDGetTokenT<HBHEDigiCollection> he_token;
		edm::EDGetTokenT<FEDRawDataCollection> raw_token;  
		edm::Handle<QIE11DigiCollection> qie11DigiCollection;
		edm::Handle<FEDRawDataCollection> raw_collection;  
};

HEP17Analyzer::HEP17Analyzer(const edm::ParameterSet& iConfig)
{
	tok_QIE11DigiCollection_ = consumes<HcalDataFrameContainer<QIE11DataFrame> >(edm::InputTag("hcalDigis"));
	he_token = consumes<HBHEDigiCollection>(edm::InputTag("hcalDigis"));
	raw_token = consumes<FEDRawDataCollection>(edm::InputTag("source"));
	runType = iConfig.getParameter<int>("RunType");//1:pedestal 2:LED 3:histogram
	outFileName=iConfig.getUntrackedParameter<string>("OutFileName");
	histoFED = iConfig.getParameter<int>("histoFED");
	_file = new TFile(outFileName.c_str(), "recreate");
	_tree = new TTree("Events", "Events");
	_tree->Branch("ieta", &ed.ieta);
	_tree->Branch("iphi", &ed.iphi);
	_tree->Branch("depth", &ed.depth);
	if(runType==1 || runType==2 || runType==3 || runType==4)//pedestal, LED, laser, global
	{
		_tree->Branch("pulse", &ed.pulse);
		_tree->Branch("tdc", &ed.tdc);
		_tree->Branch("capid", &ed.capid);
		_tree->Branch("soi", &ed.soi);
	}
	else if(runType==5)//histogram
	{
		_tree->Branch("histo", &ed.histo);
	}
	numChannels=0;
	EID=0;
}


HEP17Analyzer::~HEP17Analyzer()
{
	_file->cd();
	_file->Write();
	_file->Close();
}

void HEP17Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace hcal;
	
	iEvent.getByToken(tok_QIE11DigiCollection_,qie11DigiCollection);
	const QIE11DigiCollection& qie11dc=*(qie11DigiCollection);
	iEvent.getByToken(raw_token,raw_collection);  
	
	edm::ESHandle<HcalElectronicsMap> item;
	edm::ESHandle<HcalDbService> pSetup;
	iSetup.get<HcalDbRecord>().get(pSetup);
	iSetup.get<HcalElectronicsMapRcd>().get(item);
// 	if(EID==0)
// 	{
// 		for(int i1=10;i1<=1500;i1++)
// 		{
// 			if((raw_collection->FEDData(i1)).size()>0){fedlist.push_back(i1);}
// 		}
// 	}
// 	for(uint II=0;II<fedlist.size();II++)
// 	{
// 		const FEDRawData& raw = raw_collection->FEDData(fedlist[II]);
		
		vector <Int_t> vpulse;
		vector <Int_t> vtdc;
		vector <Int_t> vcapid;
		vector <Int_t> vsoi;
		if(runType==1 || runType==2 || runType==3 || runType==4)//pedestal, LED, laser, global
		{
			for (unsigned int j=0; j < qie11dc.size(); j++)
			{
				QIE11DataFrame qie11df = static_cast<QIE11DataFrame>(qie11dc[j]);
				DetId detid = qie11df.detid();
				HcalDetId hcaldetid = HcalDetId(detid);
				ed.ieta->push_back(hcaldetid.ieta());
				ed.iphi->push_back(hcaldetid.iphi());
				ed.depth->push_back(hcaldetid.depth());
				int nTS = qie11df.samples();
				vpulse.clear();vtdc.clear();vcapid.clear();vsoi.clear();
				for(int i=0; i<nTS; ++i)
				{
					vpulse.push_back(qie11df[i].adc());
					vtdc.push_back(qie11df[i].tdc());
					vcapid.push_back(qie11df[i].capid());
					vsoi.push_back(qie11df[i].soi());
				}
				ed.pulse->push_back(vpulse);
				ed.tdc->push_back(vtdc);
				ed.capid->push_back(vcapid);
				ed.soi->push_back(vsoi);
			}
		}
		else if(runType==5)//histogram -- need to edit for P5
		{
	// 		//the histos
	// 		const struct eventHeader* eh =(const struct eventHeader*)(raw.data());
	// 		const uint32_t* pData = (const uint32_t*) raw.data(); 
	// 		uint32_t numHistos  = ((eh->h3)>>16)&0xFFFF;
	// 		uint32_t numBins    = ((eh->h3)>>1)&0x0000FFFE; //includes overflow and header word
	// 		uint32_t fiber   = 0;
	// 		uint32_t channel = 0;
	// 		pData+=8;
	// 		int iz=0;
	// 		int nH=-1;
	// 		
	// 		for (unsigned int iHist = 0; iHist<numHistos; iHist++)
	// 		{
	// 			if(iHist==96) iz++;
	// 			fiber   = (*pData>>7)&0x1F;
	// 			channel = (*pData>>2)&0x1F;
	// 			
	// 			bool found=false;
	// 			for(int is1=0;is1<3;is1++)
	// 			{
	// 				for(int is2=0;is2<24;is2++)
	// 				{
	// 					for(int is3=0;is3<2;is3++)
	// 					{
	// 						if(MAP2Ch[is1][is2][is3][0]==((int)fiber) && MAP2Ch[is1][is2][is3][1]==((int)channel) && MAP2Ch[is1][is2][is3][2]==(iz+1)) {found=true;break;}
	// 					}
	// 					if(found){break;}
	// 				}
	// 				if(found){break;}
	// 			}
	// 			if(found)
	// 			{
	// 				nH++;
	// 				eda.ieta.push_back((int)fiber);
	// 				eda.iphi.push_back((int)channel);
	// 				eda.depth.push_back(iz+1);
	// 			}
	// 			for(unsigned int iBin = 0; iBin<numBins+1; iBin++)
	// 			{
	// 				if(iBin<60 && found)
	// 				{
	// 					eda.histo[iBin].push_back(pData[iBin+1]);
	// 				}
	// 			}
	// 			if(iHist<(numHistos-1))
	// 			{
	// 				pData+=(numBins+2);
	// 			}
	// 		}
		}
// 	}
	_tree->Fill();
	if(runType==1 || runType==2 || runType==3 || runType==4)//pedestal, LED, laser, global
	{
		ed.ieta->clear();
		ed.iphi->clear();
		ed.depth->clear();
		ed.pulse->clear();
		ed.tdc->clear();
		ed.capid->clear();
		ed.soi->clear();
	}
	else if(runType==5)//histogram -- need to edit for P5
	{
		
	}
	EID++;
}

void HEP17Analyzer::beginJob(){}

void HEP17Analyzer::endJob(){}

void HEP17Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HEP17Analyzer);
