R__LOAD_LIBRARY(libfastjet);
#include "fastjet/ClusterSequenceArea.hh"
#include"style.h"
#include <fastjet/PseudoJet.hh>
#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include <cmath>
using std::string;
using std::stringstream;
using namespace std;
using namespace fastjet;
TString flavor(int pdgid);
void sas()
{
  loadStyle();
  ofstream myfile;
  myfile.open ("jet.txt"); 
  myfile << "data begin\n"; 
  erhic::EventPythia *evt(NULL); 
  erhic::ParticleMC *particle(NULL); 
  TFile *f = new TFile("./ep_10_100_norad.root");
  TTree *tree = (TTree*)f->Get("EICTree");
  TDatabasePDG *pdgn = new TDatabasePDG();
  Int_t nEntries = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nEntries<<endl<<endl;
  tree->SetBranchAddress("event",&evt);
  tree->GetEntry(8523);
  myfile<<"#########"<<endl;
    for (int ipart = 0; ipart < evt->GetNTracks(); ++ipart)
{   int track_id = evt->GetTrack(ipart)->Id();
    int track_index=ipart+1;
    particle = evt->GetTrack(ipart);
   myfile<<" "<<track_index<<" "<<track_id<<" "<<"origin"<<particle->GetParentIndex()<<"  ks="<<particle->GetStatus()<<endl;
}
  tree->GetEntry(53196);
  myfile<<"#########"<<endl;
    for (int ipart = 0; ipart < evt->GetNTracks(); ++ipart)
{   int track_id = evt->GetTrack(ipart)->Id();
    int track_index=ipart+1;
    particle = evt->GetTrack(ipart);
   myfile<<" "<<track_index<<" "<<track_id<<" "<<"origin"<<particle->GetParentIndex()<<"  ks="<<particle->GetStatus()<<endl;
}
myfile.close();
}