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
#include "TKey.h"
#include <cmath>
#include <sstream>
#include "../HeaderFiles/StopPlotHaddInfo.h"


/*
#include "TChain.h"
#include "TLorentzVector.h"

#include "TROOT.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TF1.h"
#include "TTree.h"
#include "TString.h"
#include "TRegexp.h"
#include "TProfile.h"
#include <fstream>
#include <string>
#include <iostream>
 */

//TList *FileList;
//TFile *Target;

//void MergeRootfile( TDirectory *target, TList *sourcelist );

/*
 void hadd() {
 // in an interactive ROOT session, edit the file names
 // Target and FileList, then
 // root > .L hadd.C
 // root > hadd()
 
 Target = TFile::Open( "result.root", "RECREATE" );
 
 FileList = new TList();
 FileList->Add( TFile::Open("hsimple1.root") );
 FileList->Add( TFile::Open("hsimple2.root") );
 
 MergeRootfile( Target, FileList );
 
 }
 */
void MergeRootfile( TDirectory *target, TList *sourcelist, vector<double> * weightVec, bool beVerbose) {
    //  cout << "Target path: " << target->GetPath() << endl;
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
//    double weight;
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
            //      cout << "Merging histogram " << obj->GetName() << endl;
            TH1 *h1 = (TH1*)obj;
            //            cout << "weight being used for h1: " << weights[weightsIndex] << endl;
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
                    //                    cout << "weight being used for h2: " << weights[weightsIndex] << endl;
                    h2->Scale(weightVec->at(index));
                    if (!strcmp(name, "h_photonPt")) cout << "Integral of PhotonPt post-scaling: " << h2->Integral() << endl;
                    h1->Add( h2 );
                    delete h2;
                }
                
                nextsource = (TFile*)sourcelist->After( nextsource );
            }
        }
        /*
        else if ( obj->IsA()->InheritsFrom( TTree::Class() ) ) {
            
            // loop over all source files create a chain of Trees "globChain"
            const char* obj_name= obj->GetName();
            
            globChain = new TChain(obj_name);
            globChain->Add(first_source->GetName());
            TFile *nextsource = (TFile*)sourcelist->After( first_source );
            //      const char* file_name = nextsource->GetName();
            // cout << "file name  " << file_name << endl;
            while ( nextsource ) {
                
                globChain->Add(nextsource->GetName());
                nextsource = (TFile*)sourcelist->After( nextsource );
            }
            
        } else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
            // it's a subdirectory
            
            cout << "Found subdirectory " << obj->GetName() << endl;
            
            // create a new subdir of same name and title in the target file
            target->cd();
            TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );
            
            // newdir is now the starting point of another round of merging
            // newdir still knows its depth within the target file via
            // GetPath(), so we can still figure out where we are in the recursion
            MergeRootfile( newdir, sourcelist );
            
        } else {
            
            // object is of no type that we know or can handle
            cout << "Unknown object type, name: "
            << obj->GetName() << " title: " << obj->GetTitle() << endl;
        }
        */
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
    //plotHaddersStop(bool doPURW = 1, int whichTTBarSyst = 0, bool doSyst)
    int whichNTuple   = 1;          //as with the plot making code, leave as 1 for now -- 0 is Oviedo, 1 is DESY    
    int whichTTBarSyst = 0;          // 0 is Madgraph, 1 is MC@NLO, 2 is Powheg
    bool doPURW       = 0;          // grab the nVtx reweighted MC files
    bool doSyst       = 1;          // look at systematics plots -- (6/25/13) don't turn off for now haven't validated code fully when not making systematics
    bool beVerbose    = 0;
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
        else if (strncmp (argv[k],"beVerbose", 9) == 0) {
            beVerbose = 1;
        }
    }
    gROOT->ProcessLine("#include <vector>");
//    TRint theApp("App", &argc, argv);
//    Bool_t retVal = kTRUE;
    
    
    // set up input/output strings
    vector<TString> * nameStringVec = new vector<TString>;
    TString nTupleString, PURWString, doSystString;
    TString TTBarSystString = WhichTTBarString(whichTTBarSyst);
    switch (whichNTuple) {
        case 0:
            nTupleString = "_Oviedo";
            break;
        case 1:
            nTupleString = "_DESY";
            break;
        default:
            cout << "whichNTuple not Ovi or DESY!!!" << endl;
            break;
    }
    PURWString = (doPURW) ? "_PURW": "";
    doSystString = (doSyst) ? "_wSyst" : "";
    nameStringVec->push_back(TTBarSystString);
    nameStringVec->push_back(nTupleString);
    nameStringVec->push_back(PURWString);
    nameStringVec->push_back(doSystString); 
    nameStringVec->push_back(TString("_Output.root"));
    // set up input/output strings
    
    // boolean vector for which files to hadd
    vector<bool> * boolSampVec = new vector<bool>;
    bool doTTBar =          1;
    bool doSingTop =        1;
    bool doWLNu =           1;
    bool doVV =             1;
    bool doZDY =            1;
    bool doQCD =            1;
    if (whichTTBarSyst != 2) {
        doSingTop = 0;
        doWLNu = 0;
        doVV = 0;
        doZDY = 0;
        doQCD = 0;
    }
    boolSampVec->push_back(doTTBar);
    boolSampVec->push_back(doTTBar); // two ttbar lists
    boolSampVec->push_back(doSingTop); // one sing top list
    boolSampVec->push_back(doZDY); // one ZDY list
    boolSampVec->push_back(doVV);
    boolSampVec->push_back(doVV);
    boolSampVec->push_back(doVV);  // three VV lists
    boolSampVec->push_back(doWLNu); // one WLNu list
    boolSampVec->push_back(doQCD);
    boolSampVec->push_back(doQCD);
    boolSampVec->push_back(doQCD); // one QCD list
    vector<TList*> * fileListVec = FileListVec(whichNTuple, nameStringVec, boolSampVec);
    vector<TFile*> * outFileVec = OutFileVec(nameStringVec, boolSampVec);
    // boolean vector for which files to hadd
    // vectors of vectors of doubles to contain weight factors
    vector<vector<double> *> * weightBasesVec = new vector<vector<double> *>;
    vector<double> * currWeightBaseVec;
    vector<vector<double> *> * weightVec = new vector<vector<double> *>;
    vector<double> * currWeightVec;
    TString nEventHistName = "weightedEvents";
    float L_data = 19602.901;
    for (unsigned int i = 0; i < boolSampVec->size(); ++i) {
        if (boolSampVec->at(i)) {
            currWeightBaseVec = WeightBaseVec(fileListVec, i, nEventHistName);
            weightBasesVec->push_back(currWeightBaseVec);
            currWeightVec = WeightVec(L_data, currWeightBaseVec, i);
            weightVec->push_back(currWeightVec);
        }
    }
    cout << "outFile Vec Size " << outFileVec->size() << endl;
    cout << "fileList Vec Size " << fileListVec->size() << endl;
    cout << "weight Vec Size " << weightVec->size() << endl;
    for (unsigned int j = 0; j < boolSampVec->size(); ++j) {
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
                cout << "going to hadd QCD" << endl;
                break; 
            default:
                break;
        }
        if (boolSampVec->at(j)) {
            cout << "j " << j <<  endl;
            cout << "weight Vec at j Size " << weightVec->at(j)->size() << endl;
            MergeRootfile(outFileVec->at(j), fileListVec->at(j), weightVec->at(j), beVerbose); //I think this needs some work..
        }
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
                cout << "done hadd QCD" << endl;
                break; 
            default:
                break;
        }
    }
//    theApp.Run(retVal);
}
