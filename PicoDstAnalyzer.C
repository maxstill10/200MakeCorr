/**
 * \brief Example of how to read a file (list of files) using StFemtoEvent classes
 *
 * RunFemtoDstAnalyzer.C is an example of reading STAR FemtoDst format.
 * One can use either FemtoDst file or a list of femtoDst files (inFile.lis or
 * inFile.list) as an input, and preform physics analysis
 *
 * \author Grigory Nigmatkulov
 * \date May 29, 2018
 */

// This is needed for calling standalone classes (not needed on RACF)
#define _VANILLA_ROOT_

// C++ headers
#include <iostream>
#include <ctime>
#include <cmath>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom.h"


// FemtoDst headers
#include "StRoot/StPicoEvent/StPicoDstReader.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
/*
//We use our own Epd EP reconstruction
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StEpdUtil/StEpdEpInfo.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"
*/
//#include "StRoot/StFemtoEvent/BBCTile.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

// Load libraries (for ROOT_VERSTION_CODE >= 393215)
/*
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
R__LOAD_LIBRARY(libStPicoDst.so)
#endif
*/
// Forward declarations
bool EventCut(StPicoEvent *event);
Bool_t isGoodTrack(const TVector3& vect, const Int_t& nHits,
                   const Int_t& nHitsPoss);

const Float_t electron_mass = 0.0005485799;
const Float_t pion_mass = 0.13957061;
const Float_t kaon_mass = 0.493677;
const Float_t proton_mass = 0.9382720813;

const Float_t electron_mass_sqr = 0.000000301;
const Float_t pion_mass_sqr = 0.019479955;
const Float_t kaon_mass_sqr = 0.24371698;
const Float_t proton_mass_sqr = 0.880354499;

const Double_t lambda_mass = 1.115683;
const Double_t kaon0S_mass = 0.497611;

const Double_t alpha_Lam = 0.642;
const Double_t alpha_LamBar = -0.642;

//Func
double GetPsi(int iOrd, double Qx, double Qy);
bool ScipZeroWeight(double arr[], int size);
Float_t ZDCSMD( StPicoEvent *pEv, int eastwest, int verthori, int strip );
Float_t ZDCSMD_GetPosition( int eastwest, int verthori, int strip );
const int nSub = 3;

double const mZDCSMDCenterex = 0;
double const mZDCSMDCenterey = 0;
double const mZDCSMDCenterwx = 0;
double const mZDCSMDCenterwy = 0;

/*const Double_t mZDCSMDCenterex = 4.72466;
const Double_t mZDCSMDCenterey = 5.53629;
const Double_t mZDCSMDCenterwx = 4.39604;
const Double_t mZDCSMDCenterwy = 5.19968;*/



// inFile - is a name of name.FemtoDst.root file or a name
//          of a name.lis(t) files that contains a list of
//          name1.FemtoDst.root files
//_________________
void PicoDstAnalyzer(const Char_t *inFile, const Char_t *outputFile,
                      const string mode) {
//    gSystem->Load("StEpdUtil");
//    gSystem->Load("StRefMultCorr");
//    gSystem->Load("libStPicoDst");
    //mods
    /*
    1)6sem //To make profiles for recentring
    2)Centred //To make profiles for flatenning
    3)Flatt //To calculate Pz

    //Does not matter at Pz calculation
    4)Weighted
    5)WeightedCentred
    6)WeightedFlatt
    */

    //Histograms with corrections
    TH1F *Coef_A_n_TH_Psi1[10][nSub];
    TH1F *Coef_B_n_TH_Psi1[10][nSub];
    TH1F *Coef_A_n_TH_Psi2[10][nSub];
    TH1F *Coef_B_n_TH_Psi2[10][nSub];
    TH1F *Coef_A_n_TH_Psi3[10][nSub];
    TH1F *Coef_B_n_TH_Psi3[10][nSub];
    TH1F *Qvec1Prof_TH[nSub*2];
    TH1F *Qvec2Prof_TH[nSub*2];
    TH1F *Qvec3Prof_TH[nSub*2]; 
   

    if(mode == "Centred" || mode == "Flatt"){
        //Recentering
        TFile *input = new TFile("/star/u/mmorozov/14p5MakeCorrections/correctionsEP_output/output6sem_Psi1_2_3Corr.root", "read");
        for(int iSub=0; iSub!=2*nSub; iSub++){
            Qvec1Prof_TH[iSub] = (TH1F*)input->Get(Form("Qvec1Prof_%i", iSub));
            Qvec2Prof_TH[iSub] = (TH1F*)input->Get(Form("Qvec2Prof_%i", iSub));
            Qvec3Prof_TH[iSub] = (TH1F*)input->Get(Form("Qvec3Prof_%i", iSub));
        }

        if(mode == "Flatt"){

            //Flattening
            TFile *inputCentred = new TFile("/star/u/mmorozov/14p5MakeCorrections/correctionsEP_output/outputCentred_Psi1_2_3Corr.root", "read");
            
            for(int iProf=0; iProf!=10; iProf++){
                for(int iSub=0; iSub!=nSub; iSub++){
                    Coef_A_n_TH_Psi1[iProf][iSub] = (TH1F*)inputCentred->Get(Form("Coef_A_n_Psi1_%i_%i", iProf, iSub));
                    Coef_B_n_TH_Psi1[iProf][iSub] = (TH1F*)inputCentred->Get(Form("Coef_B_n_Psi1_%i_%i", iProf, iSub));
                    Coef_A_n_TH_Psi2[iProf][iSub] = (TH1F*)inputCentred->Get(Form("Coef_A_n_Psi2_%i_%i", iProf, iSub));
                    Coef_B_n_TH_Psi2[iProf][iSub] = (TH1F*)inputCentred->Get(Form("Coef_B_n_Psi2_%i_%i", iProf, iSub));
                    Coef_A_n_TH_Psi3[iProf][iSub] = (TH1F*)inputCentred->Get(Form("Coef_A_n_Psi3_%i_%i", iProf, iSub));
                    Coef_B_n_TH_Psi3[iProf][iSub] = (TH1F*)inputCentred->Get(Form("Coef_B_n_Psi3_%i_%i", iProf, iSub));
                }
            }      
        }
    }
    
    TFile *output = new TFile(outputFile, "RECREATE");
  
    //Reaction Plane  
    //0-West TPC
    //1-East TPC 
    TH1F *Qvec1Hist[9][2*nSub];
    TH1F *Psi1Hist[9][nSub];
    TH1F *Qvec2Hist[9][2*nSub];
    TH1F *Psi2Hist[9][nSub];
    TH1F *Qvec3Hist[9][2*nSub];
    TH1F *Psi3Hist[9][nSub];


    for(int iHist=0; iHist!=9; iHist++){
        for(int iSub=0; iSub!=nSub; iSub++){
            Qvec1Hist[iHist][2*iSub] = new TH1F(Form("Qvec1Hist_%i_%i", iHist, 2*iSub), "", 500, -1, 1);
            Qvec1Hist[iHist][2*iSub+1] = new TH1F(Form("Qvec1Hist_%i_%i", iHist, 2*iSub+1), "", 500, -1, 1);
            Qvec2Hist[iHist][2*iSub] = new TH1F(Form("Qvec2Hist_%i_%i", iHist, 2*iSub), "", 500, -1, 1);
            Qvec2Hist[iHist][2*iSub+1] = new TH1F(Form("Qvec2Hist_%i_%i", iHist, 2*iSub+1), "", 500, -1, 1);
            Qvec3Hist[iHist][2*iSub] = new TH1F(Form("Qvec3Hist_%i_%i", iHist, 2*iSub), "", 500, -1, 1);
            Qvec3Hist[iHist][2*iSub+1] = new TH1F(Form("Qvec3Hist_%i_%i", iHist, 2*iSub+1), "", 500, -1, 1);

            Psi1Hist[iHist][iSub] = new TH1F(Form("Psi1Hist_%i_%i", iHist, iSub), "", 140, 0, 7);
            Psi2Hist[iHist][iSub] = new TH1F(Form("Psi2Hist_%i_%i", iHist, iSub), "", 70, 0, 3.5);
            Psi3Hist[iHist][iSub] = new TH1F(Form("Psi3Hist_%i_%i", iHist, iSub), "", 42, 0, 2.1);
        }
    }

    //Profiles for centering
    TProfile *Qvec1Prof[2*nSub];
    TProfile *Qvec2Prof[2*nSub];
    TProfile *Qvec3Prof[2*nSub];
    for(int iSub=0; iSub!=2*nSub; iSub++){
        Qvec1Prof[iSub] = new TProfile(Form("Qvec1Prof_%i", iSub), "", 9, 0, 9);
        Qvec2Prof[iSub] = new TProfile(Form("Qvec2Prof_%i", iSub), "", 9, 0, 9);
        Qvec3Prof[iSub] = new TProfile(Form("Qvec3Prof_%i", iSub), "", 9, 0, 9);
    }
   
    
    //Profiles for flattening
    //A_n = (-2/n)*<sin(nPsi)>
    //B_n = (2/n)*<cos(nPsi)>
    TProfile *Coef_A_n_Psi1[10][nSub];
    TProfile *Coef_B_n_Psi1[10][nSub];
    TProfile *Coef_A_n_Psi2[10][nSub];
    TProfile *Coef_B_n_Psi2[10][nSub];
    TProfile *Coef_A_n_Psi3[10][nSub];
    TProfile *Coef_B_n_Psi3[10][nSub];


    for(int iProf=0; iProf!=10; iProf++){
        for(int iSub=0; iSub!=nSub; iSub++){
            Coef_A_n_Psi1[iProf][iSub] = new TProfile(Form("Coef_A_n_Psi1_%i_%i", iProf, iSub), "", 9, 0, 9);
            Coef_B_n_Psi1[iProf][iSub] = new TProfile(Form("Coef_B_n_Psi1_%i_%i", iProf, iSub), "", 9, 0, 9);
            Coef_A_n_Psi2[iProf][iSub] = new TProfile(Form("Coef_A_n_Psi2_%i_%i", iProf, iSub), "", 9, 0, 9);
            Coef_B_n_Psi2[iProf][iSub] = new TProfile(Form("Coef_B_n_Psi2_%i_%i", iProf, iSub), "", 9, 0, 9);
            Coef_A_n_Psi3[iProf][iSub] = new TProfile(Form("Coef_A_n_Psi3_%i_%i", iProf, iSub), "", 9, 0, 9);
            Coef_B_n_Psi3[iProf][iSub] = new TProfile(Form("Coef_B_n_Psi3_%i_%i", iProf, iSub), "", 9, 0, 9);
	    }
    }

    //Resolution
    TProfile *CosOfDiff_1 = new TProfile("CosOfDiff_1", "CosOfDiff for Psi_1", 9, 0, 9);
    TProfile *CosOfDiff_2 = new TProfile("CosOfDiff_2", "CosOfDiff for Psi_2", 9, 0, 9);
    TProfile *CosOfDiff_3 = new TProfile("CosOfDiff_3", "CosOfDiff for Psi_3", 9, 0, 9);
        

    //Variables
    //Reaction Plane
    double QWeight_1[2*nSub] = {}, QWeight_2[nSub] = {};
    double Qvec_1[2*nSub] = {}, Qvec_2[2*nSub] = {}, Qvec_3[2*nSub] = {};
    double Psi1[nSub] = {}, Psi2[nSub] = {}, Psi3[nSub] = {};
    double deltaPsi1[nSub] = {}, deltaPsi2[nSub] = {}, deltaPsi3[nSub] = {};
    double phi, w;
    int PP, TT, EW, iSide = 0, row;

   

    
    //-------------------------------Start analizing--------------------------------------------------
    std::cout << "Hi! Lets do some physics, Master!" << std::endl;

    Int_t runIdBins = 2000000;
    Int_t runIdRange[2] = { 11000000, 13000000 };
/*
    gSystem->Load("./libStFemtoDst.so");
    #if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
        gSystem->Load("./libStFemtoDst.so");
    #endif
*/
    StPicoDstReader* femtoReader = new StPicoDstReader(inFile);
    femtoReader->Init();

    //Long64_t events2read = femtoReader->chain()->GetEntries();

    // This is a way if you want to spead up IO
    std::cout << "Explicit read status for some branches" << std::endl;
    femtoReader->SetStatus("*",0);
    femtoReader->SetStatus("Event",1);
    femtoReader->SetStatus("Track",1);
    //femtoReader->SetStatus("KFP", 1);
    //femtoReader->SetStatus("EpdHit", 1);
    std::cout << "Status has been set" << std::endl;

    std::cout << "Now I know what to read, Master!" << std::endl;

    if( !femtoReader->chain() ) {
        std::cout << "No chain has been found." << std::endl;
    }
    Long64_t eventsInTree = femtoReader->tree()->GetEntries();
    std::cout << "eventsInTree: "  << eventsInTree << std::endl;
    Long64_t events2read = femtoReader->chain()->GetEntries();

    std::cout << "Number of events to read: " << events2read << std::endl;  

    
    // Loop over events
    for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

        std::cout << "Working on event #[" << (iEvent+1)
                << "/" << events2read << "]" << mode << std::endl;


        Bool_t readEvent = femtoReader->readPicoEvent(iEvent);
        if( !readEvent ) {
            std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
            break;
        }

        // Retrieve femtoDst
        StPicoDst *dst = femtoReader->picoDst();

        // Retrieve event information
        StPicoEvent *event = dst->event();
        if( !event ) {
            std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
            break;
        }


	    // Simple event cut
        if (EventCut(event) == false ) continue;
	
        
        //Centrality
        StRefMultCorr* refmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr_Run16_AuAu200_VpdMB5_P16ij() ;
        refmultCorrUtil -> init(event->runId());
        refmultCorrUtil -> initEvent(event->grefMult(), event->primaryVertex().z(), event->ZDCx());
        Bool_t isBadRun = refmultCorrUtil-> isBadRun(event->runId()); //reject bad runs
        //Bool_t isPileUpEvt = !refmultCorrUtil->passnTofMatchRefmultCut(1.*event->grefMult(), 1.*event->nBTOFMatch()); //reject pileup events

        int cent = refmultCorrUtil->getCentralityBin9() ;
        if (cent < 0 || isBadRun ) continue; 


        for(int iSub=0; iSub!=nSub; iSub++){
            Qvec_1[2*iSub] = 0.;
            Qvec_1[2*iSub+1] = 0.;
            QWeight_1[2*iSub] = 0.;
            Psi1[iSub] = 0.;
            Qvec_2[2*iSub] = 0.;
            Qvec_2[2*iSub+1] = 0.;
            QWeight_2[iSub] = 0.;
            Psi2[iSub] = 0.;
            Qvec_3[2*iSub] = 0.;
            Qvec_3[2*iSub+1] = 0.;
            Psi3[iSub] = 0.;
        }

        
        
        //-------------------------------Track analysis--------------------------------------------------
        Int_t nTracks = dst->numberOfTracks();

        // Track loop
        for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {
	   
        

            // Retrieve i-th femto track
            StPicoTrack *femtoTrack = dst->track(iTrk);

            if (!femtoTrack) continue;
            //std::cout << "Track #[" << (iTrk+1) << "/" << nTracks << "]"  << std::endl;
	        double eta = femtoTrack->pMom().Eta();
            double wEff = femtoTrack->pPt();

	        if (femtoTrack->nHits() < 15 || wEff < 0.15 || wEff > 2. || fabs(eta) > 1.5 || fabs(eta)<0.1) continue;
            if ((Float_t)femtoTrack->nHits()/femtoTrack->nHitsPoss() < 0.52 ) continue;

            
            double phi = femtoTrack->pMom().Phi(); //femtoTrack->isPrimary() ? femtoTrack->pMom().Phi() : 0.;
            if(eta>0) iSide = 0;//West
            else if(eta<0) iSide = 1;//East

            Qvec_2[2*iSide] += wEff*TMath::Cos(2*phi);
            Qvec_2[2*iSide+1] += wEff*TMath::Sin(2*phi);
            Qvec_3[2*iSide] += wEff*TMath::Cos(3*phi);
            Qvec_3[2*iSide+1] += wEff*TMath::Sin(3*phi);
            QWeight_2[iSide] += wEff;

            
        }//for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

		//.........................................................Qvec by ZDC calc..............................................................
        for( int iep=0; iep<2; iep++ ){ // east-west

                int nstrip;
                for( int ixy=0; ixy<2; ixy++ ){
                        if( ixy==0 ) nstrip = 8; // vertical strips (7 in x-direction)
                        else         nstrip = 9; // horizontal strips (8 in y-direction)
                        for( int is=1; is<nstrip; is++ ){
                            Float_t zdc_adc = ZDCSMD( event, iep, ixy, is );

                            if( zdc_adc<0 ) zdc_adc = 0;

                            Qvec_1[2*iep + ixy]   += ZDCSMD_GetPosition( iep, ixy, is ) * zdc_adc;
                            QWeight_1[2*iep + ixy] += zdc_adc;
                        }

                }
        }


            
        //-------------------------------reaction plane--------------------------------------------------
        if(!ScipZeroWeight(QWeight_1, 2*nSub-2)) continue;
        if(!ScipZeroWeight(QWeight_2, nSub-1)) continue;
            
        //Get Q vectors
        bool check = true;            
        for(int iSub=0; iSub!=nSub-1; iSub++){
            Qvec_1[2*iSub] = Qvec_1[2*iSub]/QWeight_1[2*iSub];
            Qvec_1[2*iSub+1] = Qvec_1[2*iSub+1]/QWeight_1[2*iSub+1];
            if(fabs(Qvec_1[2*iSub])>999 || fabs(Qvec_1[2*iSub+1])>999) check = false;

            Qvec_2[2*iSub] = Qvec_2[2*iSub]/QWeight_2[iSub];
            Qvec_2[2*iSub+1] = Qvec_2[2*iSub+1]/QWeight_2[iSub];
            if(fabs(Qvec_2[2*iSub])>999 || fabs(Qvec_2[2*iSub+1])>999) check = false;

            Qvec_3[2*iSub] = Qvec_3[2*iSub]/QWeight_2[iSub];
            Qvec_3[2*iSub+1] = Qvec_3[2*iSub+1]/QWeight_2[iSub];
            if(fabs(Qvec_3[2*iSub])>999 || fabs(Qvec_3[2*iSub+1])>999) check = false;
        }
        if(!check) continue;
        Qvec_1[4] = Qvec_1[0] - Qvec_1[2];
        Qvec_1[5] = Qvec_1[1] - Qvec_1[3];
        Qvec_2[4] = Qvec_2[0] + Qvec_2[2];
        Qvec_2[5] = Qvec_2[1] + Qvec_2[3];
        Qvec_3[4] = Qvec_3[0] + Qvec_3[2];
        Qvec_3[5] = Qvec_3[1] + Qvec_3[3];

        
        //Centring collibration
        if(mode == "Centred" || mode == "Flatt"){
            //Recentering       
            for(int iSub = 0; iSub!=2*nSub; iSub++){
                Qvec_1[iSub] -= Qvec1Prof_TH[iSub]->GetBinContent(cent+1);
                Qvec_2[iSub] -= Qvec2Prof_TH[iSub]->GetBinContent(cent+1);
                Qvec_3[iSub] -= Qvec3Prof_TH[iSub]->GetBinContent(cent+1);
            }
        }

        
        //Get Psi
        for(int iSub=0; iSub!=nSub; iSub++){
            Psi1[iSub] = GetPsi(1, Qvec_1[2*iSub], Qvec_1[2*iSub+1]);
            if(Psi1[iSub] > 7) check = false;

            Psi2[iSub] = GetPsi(2, Qvec_2[2*iSub], Qvec_2[2*iSub+1]);
            if(Psi2[iSub] > 3.5) check = false;

            Psi3[iSub] = GetPsi(3, Qvec_3[2*iSub], Qvec_3[2*iSub+1]);
            if(Psi3[iSub] > 2.1) check = false;
        }
        if(!check) continue;
                            

        //Flattening collibration
        if(mode == "Flatt"){
            //Flattaning
            for(int iProf=0; iProf!=10; iProf++){
                for(int iSub=0; iSub!=nSub; iSub++){
                    deltaPsi1[iSub] += Coef_A_n_TH_Psi1[iProf][iSub]->GetBinContent(cent+1)*TMath::Cos((iProf+1)*Psi1[iSub]) + \
                                        Coef_B_n_TH_Psi1[iProf][iSub]->GetBinContent(cent+1)*TMath::Sin((iProf+1)*Psi1[iSub]);
                    deltaPsi2[iSub] += Coef_A_n_TH_Psi2[iProf][iSub]->GetBinContent(cent+1)*TMath::Cos((iProf+1)*2*Psi2[iSub]) + \
                                        Coef_B_n_TH_Psi2[iProf][iSub]->GetBinContent(cent+1)*TMath::Sin((iProf+1)*2*Psi2[iSub]);
                    deltaPsi3[iSub] += Coef_A_n_TH_Psi3[iProf][iSub]->GetBinContent(cent+1)*TMath::Cos((iProf+1)*3*Psi3[iSub]) + \
                                        Coef_B_n_TH_Psi3[iProf][iSub]->GetBinContent(cent+1)*TMath::Sin((iProf+1)*3*Psi3[iSub]);
                }
            }
        
            for(int iSub=0; iSub!=nSub; iSub++){
                Psi1[iSub] += deltaPsi1[iSub];
                while(Psi1[iSub]>2*TMath::Pi()) Psi2[iSub] -= 2*TMath::Pi();
                while(Psi1[iSub]<0.0) Psi2[iSub] += 2*TMath::Pi();
                deltaPsi1[iSub] = 0.0;
                
                Psi2[iSub] += deltaPsi2[iSub]/2.;
                while(Psi2[iSub]>TMath::Pi()) Psi2[iSub] -= TMath::Pi();
                while(Psi2[iSub]<0.0) Psi2[iSub] += TMath::Pi();
                deltaPsi2[iSub] = 0.0;

                Psi3[iSub] += deltaPsi3[iSub]/3.;
                while(Psi3[iSub]>2*TMath::Pi()/3) Psi3[iSub] -= 2*TMath::Pi()/3;
                while(Psi3[iSub]<0.0) Psi3[iSub] += 2*TMath::Pi()/3;
                deltaPsi3[iSub] = 0.0;
            }
        }    
                        
                                    
        //Fill histograms
        for(int iSub = 0; iSub!=nSub; iSub++){
            Qvec1Hist[cent][2*iSub]->Fill(Qvec_1[2*iSub]);
            Qvec1Hist[cent][2*iSub+1]->Fill(Qvec_1[2*iSub+1]);
            Psi1Hist[cent][iSub]->Fill(Psi1[iSub]);

            Qvec2Hist[cent][2*iSub]->Fill(Qvec_2[2*iSub]);
            Qvec2Hist[cent][2*iSub+1]->Fill(Qvec_2[2*iSub+1]);
            Psi2Hist[cent][iSub]->Fill(Psi2[iSub]);

            Qvec3Hist[cent][2*iSub]->Fill(Qvec_3[2*iSub]);
            Qvec3Hist[cent][2*iSub+1]->Fill(Qvec_3[2*iSub+1]);
            Psi3Hist[cent][iSub]->Fill(Psi3[iSub]);
        }            


        if(mode == "6sem"){
            //Profiles for recentering
            for(int iSub=0; iSub!=nSub*2; iSub++){
                Qvec1Prof[iSub]->Fill(cent, Qvec_1[iSub]);
                Qvec2Prof[iSub]->Fill(cent, Qvec_2[iSub]);
                Qvec3Prof[iSub]->Fill(cent, Qvec_3[iSub]);
            }
        }    


        if(mode == "Centred"){
            //Profiles for flattening
            for(int iProf=0; iProf!=10; iProf++){
                for(int iSub=0; iSub!=nSub; iSub++){
                    Coef_A_n_Psi1[iProf][iSub]->Fill(cent, (-2.0/(iProf+1))*TMath::Sin((iProf+1)*Psi1[iSub]));
                    Coef_B_n_Psi1[iProf][iSub]->Fill(cent, (2.0/(iProf+1))*TMath::Cos((iProf+1)*Psi1[iSub]));

                    Coef_A_n_Psi2[iProf][iSub]->Fill(cent, (-2.0/(iProf+1))*TMath::Sin((iProf+1)*2*Psi2[iSub]));
                    Coef_B_n_Psi2[iProf][iSub]->Fill(cent, (2.0/(iProf+1))*TMath::Cos((iProf+1)*2*Psi2[iSub]));

                    Coef_A_n_Psi3[iProf][iSub]->Fill(cent, (-2.0/(iProf+1))*TMath::Sin((iProf+1)*3*Psi3[iSub]));
                    Coef_B_n_Psi3[iProf][iSub]->Fill(cent, (2.0/(iProf+1))*TMath::Cos((iProf+1)*3*Psi3[iSub]));
                }
            }        
        }

        CosOfDiff_1->Fill(cent, TMath::Cos(Psi1[1] - Psi1[0]));
        CosOfDiff_2->Fill(cent, TMath::Cos(2*(Psi2[1] - Psi2[0])));
        CosOfDiff_3->Fill(cent, TMath::Cos(3*(Psi3[1] - Psi3[0])));
        
        

    }//for(Long64_t iEvent=0; iEvent<events2read; iEvent++)


    //EP Distributions
    for(int iCent=0; iCent!=9; iCent++){
	    Qvec1Hist[iCent][0]->SetTitle(Form("Qx East EPD first harm (%i)", iCent));
        Qvec1Hist[iCent][1]->SetTitle(Form("Qy East EPD first harm (%i)", iCent));
        Qvec1Hist[iCent][2]->SetTitle(Form("Qx West EPD first harm (%i)", iCent));
        Qvec1Hist[iCent][3]->SetTitle(Form("Qy West EPD first harm (%i)", iCent));
        Qvec1Hist[iCent][4]->SetTitle(Form("Qx Comb EPD first harm (%i)", iCent));
        Qvec1Hist[iCent][5]->SetTitle(Form("Qy Comb EPD first harm (%i)", iCent));

        Qvec2Hist[iCent][0]->SetTitle(Form("Qx West TPC second harm (%i)", iCent));
        Qvec2Hist[iCent][1]->SetTitle(Form("Qy West TPC second harm (%i)", iCent));
        Qvec2Hist[iCent][2]->SetTitle(Form("Qx East TPC second harm (%i)", iCent));
        Qvec2Hist[iCent][3]->SetTitle(Form("Qy East TPC second harm (%i)", iCent));
        Qvec2Hist[iCent][4]->SetTitle(Form("Qx Comb TPC second harm (%i)", iCent));
        Qvec2Hist[iCent][5]->SetTitle(Form("Qy Comb TPC second harm (%i)", iCent));

        Qvec3Hist[iCent][0]->SetTitle(Form("Qx West TPC third harm (%i)", iCent));
        Qvec3Hist[iCent][1]->SetTitle(Form("Qy West TPC third harm (%i)", iCent));
        Qvec3Hist[iCent][2]->SetTitle(Form("Qx East TPC third harm (%i)", iCent));
        Qvec3Hist[iCent][3]->SetTitle(Form("Qy East TPC third harm (%i)", iCent));
        Qvec3Hist[iCent][4]->SetTitle(Form("Qx Comb TPC third harm (%i)", iCent));
        Qvec3Hist[iCent][5]->SetTitle(Form("Qy Comb TPC third harm (%i)", iCent));

        Psi1Hist[iCent][0]->SetTitle(Form("West Psi_1 EPD (%i)", iCent));
        Psi1Hist[iCent][1]->SetTitle(Form("East Psi_1 EPD (%i)", iCent));
        Psi1Hist[iCent][2]->SetTitle(Form("Comb Psi_1 EPD (%i)", iCent));
        
        Psi2Hist[iCent][0]->SetTitle(Form("West Psi_2 TPC (%i)", iCent));
        Psi2Hist[iCent][1]->SetTitle(Form("East Psi_2 TPC (%i)", iCent));
        Psi2Hist[iCent][2]->SetTitle(Form("Comb Psi_2 TPC (%i)", iCent));

        Psi3Hist[iCent][0]->SetTitle(Form("West Psi_3 TPC (%i)", iCent));
        Psi3Hist[iCent][1]->SetTitle(Form("East Psi_3 TPC (%i)", iCent));
        Psi3Hist[iCent][2]->SetTitle(Form("Comb Psi_3 TPC (%i)", iCent));

    }
    Qvec1Prof[0]->SetTitle("Profile Qx West EPD first harm");
    Qvec1Prof[1]->SetTitle("Profile Qy West EPD first harm");
    Qvec1Prof[2]->SetTitle("Profile Qx East EPD first harm");
    Qvec1Prof[3]->SetTitle("Profile Qy East EPD first harm");
    Qvec1Prof[4]->SetTitle("Profile Qx Comb EPD first harm");
    Qvec1Prof[5]->SetTitle("Profile Qy Comb EPD first harm");

    Qvec2Prof[0]->SetTitle("Profile Qx West TPC second harm");
    Qvec2Prof[1]->SetTitle("Profile Qy West TPC second harm");
    Qvec2Prof[2]->SetTitle("Profile Qx East TPC second harm");
    Qvec2Prof[3]->SetTitle("Profile Qy East TPC second harm");
    Qvec2Prof[4]->SetTitle("Profile Qx Comb TPC second harm");
    Qvec2Prof[5]->SetTitle("Profile Qy Comb TPC second harm");

    Qvec3Prof[0]->SetTitle("Profile Qx West TPC third harm");
    Qvec3Prof[1]->SetTitle("Profile Qy West TPC third harm");
    Qvec3Prof[2]->SetTitle("Profile Qx East TPC third harm");
    Qvec3Prof[3]->SetTitle("Profile Qy East TPC third harm");
    Qvec3Prof[4]->SetTitle("Profile Qx Comb TPC third harm");
    Qvec3Prof[5]->SetTitle("Profile Qy Comb TPC third harm");

    for(int i=0; i!=10; i++){
	    Coef_A_n_Psi1[i][0]->SetTitle(Form("A_n_%i West Psi_1 EPD", i));
	    Coef_B_n_Psi1[i][0]->SetTitle(Form("B_n_%i West Psi_1 EPD", i));
	    Coef_A_n_Psi1[i][1]->SetTitle(Form("A_n_%i East Psi_1 EPD", i));
        Coef_B_n_Psi1[i][1]->SetTitle(Form("B_n_%i East Psi_1 EPD", i));
	    Coef_A_n_Psi1[i][2]->SetTitle(Form("A_n_%i Comb Psi_1 EPD", i));
        Coef_B_n_Psi1[i][2]->SetTitle(Form("B_n_%i Comb Psi_1 EPD", i));
        
        Coef_A_n_Psi2[i][0]->SetTitle(Form("A_n_%i West Psi_2 TPC", i));
	    Coef_B_n_Psi2[i][0]->SetTitle(Form("B_n_%i West Psi_2 TPC", i));
	    Coef_A_n_Psi2[i][1]->SetTitle(Form("A_n_%i East Psi_2 TPC", i));
        Coef_B_n_Psi2[i][1]->SetTitle(Form("B_n_%i East Psi_2 TPC", i));
	    Coef_A_n_Psi2[i][2]->SetTitle(Form("A_n_%i Comb Psi_2 TPC", i));
        Coef_B_n_Psi2[i][2]->SetTitle(Form("B_n_%i Comb Psi_2 TPC", i));

        Coef_A_n_Psi3[i][0]->SetTitle(Form("A_n_%i West Psi_3 TPC", i));
	    Coef_B_n_Psi3[i][0]->SetTitle(Form("B_n_%i West Psi_3 TPC", i));
	    Coef_A_n_Psi3[i][1]->SetTitle(Form("A_n_%i East Psi_3 TPC", i));
        Coef_B_n_Psi3[i][1]->SetTitle(Form("B_n_%i East Psi_3 TPC", i));
	    Coef_A_n_Psi3[i][2]->SetTitle(Form("A_n_%i Comb Psi_3 TPC", i));
        Coef_B_n_Psi3[i][2]->SetTitle(Form("B_n_%i Comb Psi_3 TPC", i));
    }

    femtoReader->Finish();

    output->Write();

    std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
        << std::endl;

    
}

bool EventCut(StPicoEvent *event)
{
  bool cut = true;
  double vz = event->primaryVertex().Z(), vx = event->primaryVertex().X(), vy = event->primaryVertex().Y();
  double grefMult = event->grefMult(), tofMult = event->btofTrayMultiplicity();
  double vx_ave, vy_ave;

  if (fabs(vz) > 6. || fabs(vz-event->vzVpd()) > 3.) cut = false;
  if(fabs(vx)<1.e-5 && fabs(vy)<1.e-5 && fabs(vz)<1.e-5) cut = false;

  vx_ave = -0.205;
  vy_ave = -0.177;
  if (!event->isTrigger(520001) && !event->isTrigger(520011) && !event->isTrigger(520021) && !event->isTrigger(520031) &&
      !event->isTrigger(520041) && !event->isTrigger(520051)) cut = false;   
  if( tofMult<(-200+3.5*grefMult) ) cut = false;
  if( tofMult>( 180+5.8*grefMult) ) cut = false;

  double vxc = vx - vx_ave;
  double vyc = vy - vy_ave;

  if(( vxc*vxc + vyc*vyc) > 4) cut = false;
//  if (!GoodRun(event))  cut = false;
  return cut;


}

//________________
Bool_t isGoodTrack(const TVector3& mom, const Int_t& nHits,
                   const Int_t& nHitsPoss) {
  return ( mom.Mag() >= 0.15 &&
           mom.Mag() <= 3 &&
           nHits >= 14 &&
           (Float_t)nHits/nHitsPoss >= 0.51 );
}

//Get corner EPD
double GetPsi(int iOrd, double Qx, double Qy){
    double Psi;
    Psi = atan2(Qy, Qx)/iOrd;
    if(Psi<0.) Psi += 2*TMath::Pi() / iOrd;
    return Psi;

}


bool ScipZeroWeight(double arr[], int size) {
    for (size_t i = 0; i < size; ++i) {
        if (arr[i] == 0) {
            return false;
        }
    }
    return true;
}


Float_t ZDCSMD( StPicoEvent *pEv, int eastwest, int verthori, int strip ) {

        float val = 0;
        
        if     ( eastwest==0 && verthori==0 ) val = pEv->ZdcSmdEastVertical  (strip-1);
        else if( eastwest==0 && verthori==1 ) val = pEv->ZdcSmdEastHorizontal(strip-1);
        else if( eastwest==1 && verthori==0 ) val = pEv->ZdcSmdWestVertical  (strip-1);
        else if( eastwest==1 && verthori==1 ) val = pEv->ZdcSmdWestHorizontal(strip-1);

        
        return val;
}


Float_t ZDCSMD_GetPosition( int eastwest, int verthori, int strip ) {
// Get position of each slat;strip starts from 1

        Float_t zdcsmd_x[7] = {0.5,2,3.5,5,6.5,8,9.5};
        Float_t zdcsmd_y[8] = {1.25,3.25,5.25,7.25,9.25,11.25,13.25,15.25};

        if(eastwest==0 && verthori==0) return zdcsmd_x[strip-1]-mZDCSMDCenterex;
        if(eastwest==1 && verthori==0) return mZDCSMDCenterwx-zdcsmd_x[strip-1];
        if(eastwest==0 && verthori==1) return zdcsmd_y[strip-1]/sqrt(2.)-mZDCSMDCenterey;
        if(eastwest==1 && verthori==1) return zdcsmd_y[strip-1]/sqrt(2.)-mZDCSMDCenterwy;

  return 0;
}

