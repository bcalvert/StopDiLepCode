#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include <iostream>
#include <vector>
#include "TTree.h"
#include "TString.h"
#include "TProfile.h"
#include <TChain.h>
#include "TKey.h"
#include <cmath>
#include <sstream>
#include "../HeaderFiles/StopPlotHaddInfo.h"
void MergeRootfileHists( TDirectory *target, TList *sourcelist, vector<double> * weightVec, bool beVerbose) {
    TString path( (char*)strstr( target->GetPath(), ":" ) );
    path.Remove( 0, 2 );
    
    TFile *first_source = (TFile*)sourcelist->First();
    first_source->cd( path );
    TDirectory *current_sourcedir = gDirectory;
    //gain time, do not add the objects in the list in memory
    Bool_t status = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    
    // loop over all keys in this directory
    TChain *globChain = 0;
    TIter nextkey( current_sourcedir->GetListOfKeys() );
    TKey *key, *oldkey=0;
    int index;
    while ( (key = (TKey*)nextkey())) {
        //keep only the highest cycle number for each key
        if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;
        // read object from first source file
        first_source->cd( path );
        TObject *obj = key->ReadObj();
        
        if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
            // descendant of TH1 -> merge it
            index = 0;
            TH1 *h1 = (TH1*)obj;
            char * name = (char *) h1->GetName();
            if (beVerbose) cout << "h1 is " << name << endl;
            if (!strcmp(name, "h_qT_vs_SumEt_nVtx_3d") || !strcmp(name, "h_qT_vs_SumEt_nVtx_3dSmear")) continue;
            h1->Scale(weightVec->at(index));
            // loop over all source files and add the content of the
            // correspondant histogram to the one pointed to by "h1"
            TFile *nextsource = (TFile*)sourcelist->After( first_source );
            while ( nextsource ) {
                ++index;
                // make sure we are at the correct directory level by cd'ing to path
                nextsource->cd( path );
                TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
                if (key2) {
                    TH1 *h2 = (TH1*)key2->ReadObj();
                    if (!strcmp(name, "h_photonPt")) cout << "Integral of PhotonPt pre-scaling: " << h2->Integral() << endl;
                    h2->Scale(weightVec->at(index));
                    if (!strcmp(name, "h_photonPt")) cout << "Integral of PhotonPt post-scaling: " << h2->Integral() << endl;
                    h1->Add( h2 );
                    delete h2;
                }
                
                nextsource = (TFile*)sourcelist->After( nextsource );
            }
        }
        // now write the merged histogram (which is "in" obj) to the target file
        // note that this will just store obj in the current directory level,
        // which is not persistent until the complete directory itself is stored
        // by "target->Write()" below
        if ( obj ) {
            target->cd();
            
            //!!if the object is a tree, it is stored in globChain...
            if(obj->IsA()->InheritsFrom( TTree::Class() ))
                globChain->Merge(target->GetFile(),0,"keep");
            else
                obj->Write( key->GetName() );
        }
        
    } // while ( ( TKey *key = (TKey*)nextkey() ) )
    // save modifications to target file
    target->SaveSelf(kTRUE);
    TH1::AddDirectory(status);
}
int main( int argc, const char* argv[] ) {
    using namespace std;
    int whichNTuple   = 1;          //as with the plot making code, leave as 1 for now -- 0 is Oviedo, 1 is DESY    
    int whichTTBarSyst = 0;          // 0 is Madgraph, 1 is MC@NLO, 2 is Powheg
    bool doPURW       = 0;          // grab the nVtx reweighted MC files
    bool doSyst       = 1;          // look at systematics plots -- (6/25/13) don't turn off for now haven't validated code fully when not making systematics
    bool beVerbose    = 0;
    bool doSubLepCut  = 0;
    int  versNumber   = 1;
    int  subLepPtCut  = 10;
    bool doHardCodeNumParFiles = 0; // fix for an issue with Oviedo nTuples -- temporary (put in 7/18/13
    bool doExcSamp    = 0;          // For doing the exclusive DY and TTbar samples    
    bool doReReco     = 0;  
    for (int k = 0; k < argc; ++k) {
        cout << "argv[k] for k = " << k << " is: " << argv[k] << endl;
        if (strncmp (argv[k],"wNTuple", 7) == 0) {
            whichNTuple = strtol(argv[k+1], NULL, 10 );
        }
        else if (strncmp (argv[k],"wTTbarSyst", 10) == 0) {
            whichTTBarSyst = strtol(argv[k+1], NULL, 10 );
        }
        else if (strncmp (argv[k],"doPURW", 6) == 0) {
            doPURW = 1;
        }
        else if (strncmp (argv[k],"noSyst", 6) == 0) {
            doSyst = 0;
        }
        else if (strncmp (argv[k],"versNum", 7) == 0) {
            versNumber = strtol(argv[k+1], NULL, 10 );
        }
        else if (strncmp (argv[k],"beVerbose", 9) == 0) {
            beVerbose = 1;
        }
        else if (strncmp (argv[k],"doExcSamp", 11) == 0) {
            doExcSamp = 1;
        }
        else if (strncmp (argv[k],"doSubLepCut", 11) == 0) {
            doSubLepCut = 1;
            subLepPtCut = strtol(argv[k+1], NULL, 10 );
        }        
        else if (strncmp (argv[k],"doReReco", 8) == 0) {
            doReReco = 1;
        }
    }
    gROOT->ProcessLine("#include <vector>");
    
    ///Fix for Oviedo nTuples and number of par files
    if (whichNTuple == 0) doHardCodeNumParFiles = 1;    
    int hardCodeNumParFiles = 1;
    
    ///For doing exclusive samples, only applies for madgraph
    if (doExcSamp && whichTTBarSyst != 2) {
        cout << "Exclusive samples only for Madgraph!!" << endl;
        return 0;
    }
    
    // set up input/output strings
    vector<TString> * nameStringVec = new vector<TString>;
    TString nTupleString, PURWString, doSystString, outString = "";
    TString TTBarSystString = WhichTTBarString(whichTTBarSyst, whichNTuple);
    cout << "TTBarstring " << TTBarSystString << endl;    
    switch (whichNTuple) {
        case 0:
            if (versNumber == 2) {
                doSubLepCut = true;
                subLepPtCut = 20;
            }
            nTupleString += "_Oviedo";
            break;
        case 1:
            if (versNumber == 2) {
                nTupleString += "_DESY_SkimOutput";
            }
            nTupleString += "_DESY";
            break;
        default:
            cout << "whichNTuple not Ovi or DESY!!!" << endl;
            break;
    }
    PURWString = (doPURW) ? "_PURW": "";
    doSystString = (doSyst) ? "_wSyst" : "";    
    if (doSubLepCut) {
        doSystString += "_subLepPtCut";
        doSystString += subLepPtCut;
    }
    outString += "_Output.root";
    nameStringVec->push_back(TTBarSystString);
    nameStringVec->push_back(nTupleString);
    nameStringVec->push_back(PURWString);
    nameStringVec->push_back(doSystString);
    nameStringVec->push_back(outString);
    // set up input/output strings
    
    // boolean vector for which files to hadd
    vector<bool> * boolSampVec = new vector<bool>;
    bool doTTBar            = 1;
    bool doSingTop          = 1;
    bool doWLNu             = 1;
    bool doVV               = 1;
    bool doZDY              = 1;
    bool doQCD              = 1;
    bool doVG               = 1;
    bool doHiggs            = 1;
    bool doRare             = 1;
    if (whichTTBarSyst != 2) {
        doSingTop = 0;
        doZDY = 0;
        doVV = 0;
        doWLNu = 0;
        doQCD = 0;
        doVG = 0;
        doHiggs = 0;
        doRare = 0;
    }
    boolSampVec->push_back(doTTBar);
    boolSampVec->push_back(doTTBar); // two ttbar lists, one for diLep "Signal", and one for other decay modes mimicking diLep
    boolSampVec->push_back(doSingTop); // one sing top list
    boolSampVec->push_back(doZDY); // one ZDY list
    boolSampVec->push_back(doVV);
    boolSampVec->push_back(doVV);
    boolSampVec->push_back(doVV);  // three VV lists
    boolSampVec->push_back(doWLNu); // one WLNu list
    if (whichNTuple == 1) {
        boolSampVec->push_back(doQCD);
        boolSampVec->push_back(doQCD);
        boolSampVec->push_back(doQCD); // three QCD list
    }
    else {
        boolSampVec->push_back(doVG);
        boolSampVec->push_back(doVG);
        boolSampVec->push_back(doHiggs); //Three Higgs lists
        boolSampVec->push_back(doHiggs);
        boolSampVec->push_back(doHiggs);
        boolSampVec->push_back(doRare);
        boolSampVec->push_back(doRare);
    }
    vector<TList*> * fileListVec = FileListVec(whichNTuple, nameStringVec, boolSampVec, doExcSamp);
    vector<TFile*> * outFileVec = OutFileVec(whichNTuple, nameStringVec, boolSampVec, doExcSamp);
    // boolean vector for which files to hadd
    
    // vectors of vectors of doubles to contain weight factors
    vector<vector<double> *> * weightBasesVec = new vector<vector<double> *>;
    vector<double> * currWeightBaseVec;
    vector<vector<double> *> * weightVec = new vector<vector<double> *>;
    vector<double> * currWeightVec;
    TString nEventHistName = "weightedEvents";
    TString nParFileHistName = "h_numParFiles";
    float L_data;
    float indLumiDESY[4] = {892, 4404, 7032, 7274};
    float indLumiOvi[4] = {892, 4404, 7032, 7274};
    
    float indLumiOviReReco[4] = {876, 4404, 7016, 7360};
//    float L_data = 19602.901;
    if (whichNTuple == 0) {
        if (!doReReco) {
            L_data = indLumiOvi[0] + indLumiOvi[1] + indLumiOvi[2] + indLumiOvi[3];
        }
        else {
            L_data = indLumiOviReReco[0] + indLumiOviReReco[1] + indLumiOviReReco[2] + indLumiOviReReco[3];
        }
    }
    else {
        L_data = indLumiDESY[0] + indLumiDESY[1] + indLumiDESY[2] + indLumiDESY[3];
    }
    cout << "L_data " << L_data << endl;
    cout << "bool samp size " << boolSampVec->size() << endl;
    for (unsigned int i = 0; i < boolSampVec->size(); ++i) {
        vector<int> * numParFilesVec = new vector<int>;
        if (boolSampVec->at(i)) {
            currWeightBaseVec = WeightBaseVec(whichNTuple, fileListVec, i, nEventHistName, nParFileHistName, numParFilesVec);
            weightBasesVec->push_back(currWeightBaseVec);
            if (doHardCodeNumParFiles) {
                for (unsigned iHardCode = 0; iHardCode < numParFilesVec->size(); ++iHardCode) {
                    numParFilesVec->at(iHardCode) = hardCodeNumParFiles;
                }
            }
            currWeightVec = WeightVec(whichNTuple, L_data, currWeightBaseVec, i, numParFilesVec, doExcSamp);
            weightVec->push_back(currWeightVec);
        }
    }
    // vectors of vectors of doubles to contain weight factors
    
    for (unsigned int j = 0; j < boolSampVec->size(); ++j) {
        if (boolSampVec->at(j)) {
            switch (j) {
                case 0:
                    cout << "going to hadd TTBar" << endl;
                    break;
                case 2:
                    cout << "going to hadd SingTop" << endl;
                    break; 
                case 3:
                    cout << "going to hadd ZDY" << endl;
                    break; 
                case 4:
                    cout << "going to hadd VV" << endl;
                    break; 
                case 7:
                    cout << "going to hadd WLNu" << endl;
                    break;                 
                case 8:
                    if (whichNTuple == 1) {
                        cout << "going to hadd QCD" << endl;
                    }
                    else {
                        cout << "going to hadd VG" << endl;
                    }
                    break; 
                case 10:
                    if (whichNTuple == 0) {
                        cout << "going to hadd Higgs" << endl;
                    }
                case 13:
                    if (whichNTuple == 0) {
                        cout << "going to hadd Rare backgrounds" << endl;
                    }
                default:
                    break;
            }
            for (unsigned int iWeight = 0; iWeight < weightVec->at(j)->size(); ++iWeight) {
                cout << "weight for j = " << j << " and iWeight = " << iWeight << " is " << weightVec->at(j)->at(iWeight) << endl;
            }
            MergeRootfileHists(outFileVec->at(j), fileListVec->at(j), weightVec->at(j), beVerbose); //I think this needs some work..       
            outFileVec->at(j)->Close();
            switch (j) {
                case 1:
                    cout << "done hadd TTBar" << endl;
                    break;
                case 2:
                    cout << "done hadd SingTop" << endl;
                    break; 
                case 3:
                    cout << "done hadd ZDY" << endl;
                    break; 
                case 6:
                    cout << "done hadd VV" << endl;
                    break; 
                case 7:
                    cout << "done hadd WLNu" << endl;
                    break; 
                case 8:
                    if (whichNTuple == 1) {
                        cout << "done hadd QCD" << endl;
                    }
                    else {
                        cout << "done hadd VG" << endl;
                    }                
                    break; 
                case 10:
                    if (whichNTuple == 0) {
                        cout << "done hadd Higgs" << endl;
                    }
                case 13:
                    if (whichNTuple == 0) {
                        cout << "done hadd Rare backgrounds" << endl;
                    }
                default:
                    break;
            }
        }
    }
}