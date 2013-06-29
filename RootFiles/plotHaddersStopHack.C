#include "TROOT.h"
#include "TFile.h"
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
#include <cmath>
#include <sstream>

//TList *FileList;
//TFile *Target;

//void MergeRootfile( TDirectory *target, TList *sourcelist );


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
//double*   WeightArrayFiller(TList *sourcelist, double weights[], int arrayL) {
void  WeightArrayFiller(TList *sourcelist, double * weights) {
    TFile * first_source = (TFile*) sourcelist->First();
    cout << "first_source name = " << first_source->GetName() << endl;
    int index = 0;
    //    double * weights[ArrayLength]
    TH1F * h_eventCount = (TH1F*) first_source->Get("weightedEvents");
    float nEvents = h_eventCount->Integral();
    //    cout << "hi" << endl;
    cout << "nEvents " << nEvents << endl;
    cout << weights[index] << endl;
    weights[index] = nEvents;
    //    cout << "hi2" << endl;
    TFile *nextsource = (TFile*)sourcelist->After( first_source );
    //    cout << "hi3" << endl;
    while (nextsource) {
        //      cout << "hi4" << endl;
        ++index;
        cout << "nextsource " << nextsource->GetName() << endl;
        h_eventCount = (TH1F*) nextsource->Get("weightedEvents");
        //      cout << "hi5" << endl;
        nEvents = h_eventCount->Integral();
        cout << "nEvents " << nEvents << endl;
        //      cout << "hi6" << endl;
        weights[index] = nEvents;
        //      cout << "hi7" << endl;
        TFile *nextsource = (TFile*)sourcelist->After( nextsource );
        //      cout << "hi8" << endl;
    }
    //    return weights;
}
void MergeRootfile( TDirectory *target, TList *sourcelist, double weights[] ) {
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
    while ( (key = (TKey*)nextkey())) {
        //keep only the highest cycle number for each key
        if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;
        
        // read object from first source file
        first_source->cd( path );
        TObject *obj = key->ReadObj();
        
        if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
            // descendant of TH1 -> merge it
            weightsIndex = 0;
            //      cout << "Merging histogram " << obj->GetName() << endl;
            TH1 *h1 = (TH1*)obj;
            //            cout << "weight being used for h1: " << weights[weightsIndex] << endl;
            char * name = h1->GetName();
            cout << "h1 is " << name << endl;
            if (!strcmp(name, "h_qT_vs_SumEt_nVtx_3d") || !strcmp(name, "h_qT_vs_SumEt_nVtx_3dSmear")) continue;
            h1->Scale(weights[weightsIndex]);
            // loop over all source files and add the content of the
            // correspondant histogram to the one pointed to by "h1"
            TFile *nextsource = (TFile*)sourcelist->After( first_source );
            while ( nextsource ) {
                weightsIndex+=1;
                // make sure we are at the correct directory level by cd'ing to path
                nextsource->cd( path );
                TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
                if (key2) {
                    TH1 *h2 = (TH1*)key2->ReadObj();
                    if (!strcmp(name, "h_photonPt")) cout << "Integral of PhotonPt pre-scaling: " << h2->Integral() << endl;
                    //                    cout << "weight being used for h2: " << weights[weightsIndex] << endl;
                    h2->Scale(weights[weightsIndex]);
                    if (!strcmp(name, "h_photonPt")) cout << "Integral of PhotonPt post-scaling: " << h2->Integral() << endl;
                    h1->Add( h2 );
                    delete h2;
                }
                
                nextsource = (TFile*)sourcelist->After( nextsource );
            }
        }
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
/*
 TH1 * hadd(TFile** inputFileArray, float** weights, int ArrayLength, TString HistName) {
 
 for (int i = 0; i < ArrayLength; ++i) {
 
 }            
 }
 */
//int main( int argc, const char* argv[] ) {
void plotHaddersStopHack(bool doPURW = 1, int whichTTBarSyst = 0, bool doSyst) {
    bool FiveTwoX = 1;
    TList *FileListTTBarSig, *FileListTTBarBkg, *FileListSingTop, *FileListZDY, *FileListWLNu, *FileListWW, *FileListWZ, *FileListZZ, *FileListQCDMu, *FileListQCDEM, * FileListQCDBCEM;
    TFile *TargetTTBarSig, *TargetTTBarBkg, *TargetSingTop, *TargetZDY, *TargetWLNu, *TargetWW, *TargetWZ, *TargetZZ, *TargetQCDMu, *TargetQCDEM, * TargetQCDBCEM;
    gROOT->ProcessLine("#include <vector>");
    TString TTBarSystString, DESYString, PURWString;
    switch (whichTTBarSyst) {
        case 0:
            TTBarSystString = "";
            break;
        case 1:
            TTBarSystString = "_scaleup";
            break;
        case 2:
            TTBarSystString = "_scaledown";
            break;
        case 3:
            TTBarSystString = "_massup";
            break;
        case 4:
            TTBarSystString = "_massdown";
            break;
        case 5:
            TTBarSystString = "_matchingup";
            break;
        case 6:
            TTBarSystString = "_matchingdown";
            break;
        case 7:
            TTBarSystString = "_mcatnlo";
            break;
        case 8:
            TTBarSystString = "_powheg";
            break;
    }
    DESYString = "_DESY";
    doSystString = "";
    if (doPURW) PURWString = "_PURW";
    if (doSyst) doSystString = "_wSyst";
    bool doTTBar =          1;
    bool doSingTop =        1;
    bool doWLNu =           1;
    bool doVV =             1;
    bool doZDY =            1;
    bool doQCD =            1;
    if (whichTTBarSyst != 0) {
      doSingTop = 0;
      doWLNu = 0;
      doVV = 0;
      doZDY = 0;
      doQCD = 0;
    }

    if (doTTBar) {
        TargetTTBarSig = TFile::Open(TString("TTBarSig") + TTBarSystString + DESYString + PURWString + doSystString + TString("Haddplots.root"), "RECREATE");
        TargetTTBarBkg = TFile::Open(TString("TTBarBkg") + TTBarSystString + DESYString + PURWString + doSystString + TString("Haddplots.root"), "RECREATE");
    }
    if (doQCD) {
        TargetQCDMu = TFile::Open(TString("QCDMu") + DESYString + PURWString + doSystString + TString("Haddplots.root"), "RECREATE");   
        TargetQCDEM = TFile::Open(TString("QCDEM") + DESYString + PURWString + doSystString + TString("Haddplots.root"), "RECREATE");   
        TargetQCDBCEM = TFile::Open(TString("QCDBCEM") + DESYString + PURWString + doSystString + TString("Haddplots.root"), "RECREATE");   
    }
    if (doZDY) TargetZDY = TFile::Open(TString("ZDY") + DESYString + PURWString + doSystString + TString("Haddplots.root"), "RECREATE");
    if (doWLNu) TargetWLNu = TFile::Open(TString("WToLNu") + DESYString + PURWString + doSystString + TString("Haddplots.root"), "RECREATE");
    if (doSingTop) TargetSingTop = TFile::Open(TString("SingleTop") + DESYString + PURWString + doSystString + TString("Haddplots.root"), "RECREATE");
    if (doVV) {
        TargetWW = TFile::Open(TString("WW") + DESYString + PURWString + doSystString + TString("Haddplots.root"), "RECREATE");   
        TargetWZ = TFile::Open(TString("WZ") + DESYString + PURWString + doSystString + TString("Haddplots.root"), "RECREATE");   
        TargetZZ = TFile::Open(TString("ZZ") + DESYString + PURWString + doSystString + TString("Haddplots.root"), "RECREATE");   
    }    
    double WeightArrayTTBarSig[1];
    double WeightArrayTTBarSigBase[1];
    double WeightArrayTTBarBkg[1];
    double WeightArrayTTBarBkgBase[1];
    double WeightArraySingTop[2];
    double WeightArraySingTopBase[2];
    double WeightArrayWLNu[1];
    double WeightArrayWLNuBase[1];
    double WeightArrayWW[1];
    double WeightArrayWWBase[1];
    double WeightArrayWZ[1];
    double WeightArrayWZBase[1];
    double WeightArrayZZ[1];
    double WeightArrayZZBase[1];
    double WeightArrayZDY[2];
    double WeightArrayZDYBase[2];
    double WeightArrayQCDMu[6];
    double WeightArrayQCDMuBase[6];
    double WeightArrayQCDEM[3];
    double WeightArrayQCDEMBase[3];
    double WeightArrayQCDBCEM[3];
    double WeightArrayQCDBCEMBase[3];
    if (doTTBar) {
        FileListTTBarSig = new TList();
        FileListTTBarSig->Add(TFile::Open(TString("ttbarsignalplustau") + TTBarSystString + DESYString + PURWString + doSystString + TString("_Output.root")));
        WeightArrayFiller(FileListTTBarSig, WeightArrayTTBarSigBase);
        
        FileListTTBarBkg = new TList();
        FileListTTBarBkg->Add(TFile::Open(TString("ttbarbg") + TTBarSystString + DESYString + PURWString + doSystString + TString("_Output.root")));
        WeightArrayFiller(FileListTTBarBkg, WeightArrayTTBarBkgBase);
    }
    if (doSingTop) {
        FileListSingTop = new TList();
        FileListSingTop->Add(TFile::Open(TString("singletop_tw") + DESYString + PURWString + doSystString + TString("_Output.root")));
        FileListSingTop->Add(TFile::Open(TString("singleantitop_tw") + DESYString + PURWString + doSystString + TString("_Output.root")));
        WeightArrayFiller(FileListSingTop, WeightArraySingTopBase);
    }
    if (doZDY) {
        FileListZDY = new TList();
        FileListZDY->Add(TFile::Open(TString("dy1050") + DESYString + PURWString + doSystString + TString("_Output.root")));
        FileListZDY->Add(TFile::Open(TString("dy50inf") + DESYString + PURWString + doSystString + TString("_Output.root")));
        WeightArrayFiller(FileListZDY, WeightArrayZDYBase);
    }
    if (doVV) {
        FileListWW = new TList();
        FileListWW->Add(TFile::Open(TString("wwtoall") + DESYString + PURWString + doSystString + TString("_Output.root")));
        WeightArrayFiller(FileListWW, WeightArrayWWBase);
        FileListWZ = new TList();
        FileListWZ->Add(TFile::Open(TString("wztoall") + DESYString + PURWString + doSystString + TString("_Output.root")));
        WeightArrayFiller(FileListWZ, WeightArrayWZBase);
        FileListZZ = new TList();
        FileListZZ->Add(TFile::Open(TString("zztoall") + DESYString + PURWString + doSystString + TString("_Output.root")));
        WeightArrayFiller(FileListZZ, WeightArrayZZBase);
    }
    if (doWLNu) {
        FileListWLNu = new TList();
        FileListWLNu->Add(TFile::Open(TString("wtolnu") + DESYString + PURWString + doSystString + TString("_Output.root")));
        WeightArrayFiller(FileListWLNu, WeightArrayWLNuBase);
    }
    if (doQCD) {
        FileListQCDMu = new TList();
        FileListQCDMu->Add(TFile::Open(TString("qcdmu15") + DESYString + PURWString + doSystString + TString("_Output.root")));
        //        FileListQCDMu->Add(TFile::Open(TString("qcdmu15") + DESYString + PURWString + doSystString + TString("_Output.root")));
        WeightArrayFiller(FileListQCDMu, WeightArrayQCDMuBase);
        
        FileListQCDEM = new TList();
        FileListQCDEM->Add(TFile::Open(TString("qcdem2030") + DESYString + PURWString + doSystString + TString("_Output.root")));
        FileListQCDEM->Add(TFile::Open(TString("qcdem3080") + DESYString + PURWString + doSystString + TString("_Output.root")));
        FileListQCDEM->Add(TFile::Open(TString("qcdem80170") + DESYString + PURWString + doSystString + TString("_Output.root")));
        WeightArrayFiller(FileListQCDEM, WeightArrayQCDEMBase);
        
        FileListQCDBCEM = new TList();
        FileListQCDBCEM->Add(TFile::Open(TString("qcdbcem2030") + DESYString + PURWString + doSystString + TString("_Output.root")));
        FileListQCDBCEM->Add(TFile::Open(TString("qcdbcem3080") + DESYString + PURWString + doSystString + TString("_Output.root")));
        FileListQCDBCEM->Add(TFile::Open(TString("qcdbcem80170") + DESYString + PURWString + doSystString + TString("_Output.root")));
        WeightArrayFiller(FileListQCDBCEM, WeightArrayQCDBCEMBase);
    }
    
    //Cross Sections mostly taken from: https://twiki.cern.ch/twiki/bin/view/CMS/StandardModelCrossSectionsat8TeV
    //    float L_data = 12200;
    float L_data = 19602.901;
    double xsecTTbar = 244.849;
    double xsecSingTop = 11.1;
    double xsecWW = 54.838;
    double xsecWZ = 33.21;
    double xsecZZ = 17.654;
    double xsecZDY10to50 = 860.5;
    //double xsecZDY50toInf = 3532.8;
    double xsecZDY50toInf = 3503.71;
    double xsecWLNu = 36257.2;
    double xsecQCDMu15 = 3.640E8*3.7E-4;
    double xsecQCDMu20to30 = 2.870E8*6.500E-3;
    double xsecQCDMu30to50 = 6.609E7*12.20E-3;
    double xsecQCDMu50to80 = 8.802E6*21.80E-3;
    double xsecQCDMu80to120 = 1.024E6*39.50E-3;
    double xsecQCDMu120to170 = 1.578E5*47.30E-3;
    double xsecQCDEM20to30 = 2.886E8*10.10E-3;
    double xsecQCDEM30to80 = 7.433E7*62.10E-3;
    double xsecQCDEM80to170 = 1.191E6*153.9E-3;
    double xsecQCDBCEM20to30 = 2.886E8*5.800E-4;
    double xsecQCDBCEM30to80 = 7.424E7*2.250E-3;
    double xsecQCDBCEM80to170 = 1.191E6*10.90E-3;
    WeightArrayTTBarSig[0] = L_data * xsecTTbar / WeightArrayTTBarSigBase[0];
    WeightArrayTTBarBkg[0] = L_data * xsecTTbar / WeightArrayTTBarBkgBase[0];
    WeightArraySingTop[0] = L_data * xsecSingTop / WeightArraySingTopBase[0];
    WeightArrayWW[0]    = L_data * xsecWW / WeightArrayWWBase[0];
    WeightArrayWZ[0]    = L_data * xsecWZ / WeightArrayWZBase[0];
    WeightArrayZZ[0]    = L_data * xsecZZ / WeightArrayZZBase[0];
    WeightArrayWLNu[0]    = L_data * xsecWLNu / WeightArrayWLNuBase[0];
    WeightArrayZDY[0]    = L_data * xsecZDY10to50 / WeightArrayZDYBase[0];
    WeightArrayZDY[1]    = L_data * xsecZDY50toInf / WeightArrayZDYBase[1];
    WeightArrayQCDMu[0] = L_data * xsecQCDMu15 / WeightArrayQCDMuBase[0];
    WeightArrayQCDMu[1] = L_data * xsecQCDMu20to30 / WeightArrayQCDMuBase[1];
    WeightArrayQCDMu[2] = L_data * xsecQCDMu30to50 / WeightArrayQCDMuBase[2];
    WeightArrayQCDMu[3] = L_data * xsecQCDMu50to80 / WeightArrayQCDMuBase[3];
    WeightArrayQCDMu[4] = L_data * xsecQCDMu80to120 / WeightArrayQCDMuBase[4];
    WeightArrayQCDMu[5] = L_data * xsecQCDMu120to170 / WeightArrayQCDMuBase[5];
    WeightArrayQCDEM[0] = L_data * xsecQCDEM20to30 / WeightArrayQCDEMBase[0];
    WeightArrayQCDEM[1] = L_data * xsecQCDEM30to80 / WeightArrayQCDEMBase[1];
    WeightArrayQCDEM[2] = L_data * xsecQCDEM80to170 / WeightArrayQCDEMBase[2];
    WeightArrayQCDBCEM[0] = L_data * xsecQCDBCEM20to30 / WeightArrayQCDBCEMBase[0];
    WeightArrayQCDBCEM[1] = L_data * xsecQCDBCEM30to80 / WeightArrayQCDBCEMBase[1];
    WeightArrayQCDBCEM[2] = L_data * xsecQCDBCEM80to170 / WeightArrayQCDBCEMBase[2];
    if (doTTBar) {
        cout << "going to hadd TTBar" << endl;
        MergeRootfile( TargetTTBarSig, FileListTTBarSig, WeightArrayTTBarSig);
        MergeRootfile( TargetTTBarBkg, FileListTTBarBkg, WeightArrayTTBarBkg);
        cout << "done hadd TTBar" << endl;
    }
    if (doSingTop) {
        cout << "going to hadd SingTop" << endl;
        MergeRootfile( TargetSingTop, FileListSingTop, WeightArraySingTop);
        cout << "done hadd SingTop" << endl;
    }
    if (doVV) {
        cout << "going to hadd VV" << endl;
        MergeRootfile( TargetWW, FileListWW, WeightArrayWW);
        MergeRootfile( TargetWZ, FileListWZ, WeightArrayWZ);
        MergeRootfile( TargetZZ, FileListZZ, WeightArrayWW);
        cout << "done hadd VV" << endl;
    }
    if (doWLNu) {
        cout << "going to hadd WLNu" << endl;
        MergeRootfile( TargetWLNu, FileListWLNu, WeightArrayWLNu);
        cout << "done hadd WLNu" << endl;
    }
    if (doZDY) {
        cout << "going to hadd ZDY" << endl;
        MergeRootfile( TargetZDY, FileListZDY, WeightArrayZDY);
        cout << "done hadd ZDY" << endl;
    }
    if (doQCD) {
        cout << "going to hadd QCD" << endl;
        MergeRootfile( TargetQCDMu, FileListQCDMu, WeightArrayQCDMu);
        MergeRootfile( TargetQCDEM, FileListQCDEM, WeightArrayQCDEM);
        MergeRootfile( TargetQCDBCEM, FileListQCDBCEM, WeightArrayQCDBCEM);
        cout << "done hadd QCD" << endl;
    }
}
