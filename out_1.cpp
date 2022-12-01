R__LOAD_LIBRARY(libfastjet);
#include "fastjet/ClusterSequenceArea.hh"
#include"style.h"
#include <fastjet/PseudoJet.hh>
#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include <cmath>
#include "TLorentzVector.h"
#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>+;
#endif
using std::string;
using std::stringstream;
using namespace std;
using namespace fastjet;
using namespace std;
TString flavor(int pdgid);
double cal_angle(erhic::ParticleMC* part1,erhic::ParticleMC* part2);
void find_matching_index(erhic::EventPythia* evt, vector<PseudoJet>& constit, vector<int>& match_list);
int find_quark_origin(erhic::EventPythia* evt, erhic::ParticleMC* part);
static int jet_number=1;
int find_origin(erhic::EventPythia* evt, erhic::ParticleMC* part);
template <typename T, typename A>
int arg_max(std::vector<T, A> const& vec) {
  return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
}
void out_1()
{
  loadStyle();
  ofstream myfile;
  myfile.open ("jet_lepton_1.txt"); 
  //myfile << "data begin\n"; 
  erhic::EventPythia *evt(NULL); 
  erhic::ParticleMC *particle(NULL); 
  erhic::ParticleMC *quark(NULL); 
  TFile *f = new TFile("./ep_18_275_norad_2.root");
  TTree *tree = (TTree*)f->Get("EICTree");
  TDatabasePDG *pdgn = new TDatabasePDG();
  Int_t nEntries = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nEntries<<endl<<endl;

  tree->SetBranchAddress("event",&evt);
  TH2D *fraction[3][3];
  char *histname = new char[10];
  double etahigh[3]={-1,1,3.5};
  double etalow[3] = {-3.5,-1,1};
  double pthigh[3]={5.,10.,20.};
  double ptlow[3]={0.,5.,10.};
  for(int ihigh=0;ihigh<3;ihigh++)
    {
      for(int jhigh=0;jhigh<3;jhigh++)
       {
         sprintf(histname, "h_x_%d_%d",ihigh,jhigh); 
        fraction[ihigh][jhigh] = new TH2D(histname,"",50,-3,1.,100,1.,3.);
       }
    }  

  for(Int_t j=0;j<nEntries;j++){ 
   cout<<"***********************reading event"<<j<<"*********************"<<endl;
   tree->GetEntry(j);
   int process=evt->GetProcess();
   if(process!=99) continue;
  vector<PseudoJet> jet_constits;
  for (int ipart = 0; ipart < evt->GetNTracks(); ++ipart)
{  particle = evt->GetTrack(ipart);
   if (particle->GetStatus()==1 && fabs(particle->GetEta())<3.5 && particle->Id()!=11) 
   jet_constits.push_back(PseudoJet(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE()));
}

JetDefinition R1jetdef(antikt_algorithm, 1); // anti-kt and R = 1
//ClusterSequence cs(jet_constits, R1jetdef);
//vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
double maxrap = 7.0;
unsigned int n_repeat = 1; // default is 1
double ghost_area = 0.01; // this is the default
//ClusterSequenceArea clust_seq(jet_constits, R1jetdef, area_def);
const ClusterSequenceArea clust_seq(jet_constits, R1jetdef, VoronoiAreaSpec());
vector<PseudoJet> jets1 = sorted_by_pt(clust_seq.inclusive_jets());

//find scattered electron

particle = evt->GetTrack(2);
double phi_electron = particle->GetPhi();
int counter_quarks=0;
double phi_quark;
quark= evt->GetTrack(9);
double struck_pt = quark->GetPt();
double quark_eta = quark->GetEta();
double quark_phi = quark->GetPhi();


for (unsigned ijet = 0; ijet < jets1.size(); ijet++)
{  
    myfile<<"######"<<endl;
    int verbosity=0;
    double den=jets1[ijet].e();
    double jets_eta = jets1[ijet].eta();
    double jets_phi = jets1[ijet].phi();
    double deltaR =sqrt((jets_eta -quark_eta)*(jets_eta -quark_eta)+(jets_phi-quark_phi)*(jets_phi-quark_phi));
    double jets_pt = jets1[ijet].pt();
  //  if (verbosity>1) cout << "jet " << ijet << " (pt mass eta phi): "<< jets[ijet].pt() << " " << jets[ijet].m() << " "  << jets[ijet].rap() << " " << jets[ijet].phi() << endl;
    vector<PseudoJet> constituents = jets1[ijet].constituents(); // get jet constituents
   // match the constitutes to TTree indices
    vector<int> matching_index;
    find_matching_index(evt, constituents, matching_index);
    TParticlePDG *pdg;
    double other_energy=0;
    TLorentzVector neutral_hadron_vector = TLorentzVector(0,0,0,0);
    int counter=0; // at least two particles originating from struck quark
    vector<double> part_pt  = vector<double>();
    vector<int> index_list  = vector<int>();
    //
    for (unsigned iconstit = 0; iconstit < constituents.size(); iconstit++)
    {   
        int track_index = matching_index[iconstit];
       
       if (track_index>0)
        {    TString part;
            Int_t track_id = evt->GetTrack(track_index-1)->Id();     
            //TParticlePDG *
            index_list.push_back(track_index-1);
            part_pt.push_back(evt->GetTrack(track_index-1)->GetPt());
        }
      }
      int argmax = arg_max(part_pt);
      int marker = index_list[argmax];
      cout<<"argamx: "<<argmax<<endl;
      erhic::ParticleMC* max_part = evt->GetTrack(marker);
      if (max_part==NULL) continue;
      double max_eta = max_part->GetEta();
      double max_phi = max_part->GetPhi();
       Int_t max_id = evt->GetTrack(marker)->Id();  
      double deltaR_part = sqrt((jets_eta -max_eta)*(jets_eta -max_eta)+(jets_phi-max_phi)*(jets_phi-max_phi));

      myfile<<quark->Id()<<" "<<deltaR<<" "<<max_id<<" "<<deltaR_part<<endl;

   }
}
\
myfile.close();
}
/*TCanvas *C4 = new TCanvas();
C4->SetLogz();
thetajet->SetStats(0);
thetajet->GetXaxis()->SetTitle("eta");
thetajet->GetXaxis()->CenterTitle();
thetajet->GetYaxis()->SetTitle("Number of particles originating from struck quark");
thetajet->GetYaxis()->CenterTitle();
thetajet->Draw("colz");}*/

void find_matching_index(erhic::EventPythia* evt, vector<PseudoJet>& constit, vector<int>& match_list)
{
  vector<int> unmatch_list;
  for (int ipart = 0; ipart < evt->GetNTracks(); ++ipart)
  {
    erhic::ParticleMC* particle = evt->GetTrack(ipart);

    unmatch_list.push_back(particle->GetIndex());
  }

  for (int iconstit = 0; iconstit < constit.size(); ++iconstit)
  {cout<<"constit"<<constit[iconstit].px()<<"  "<<constit[iconstit].py()<<" "<<constit[iconstit].pz()<<endl;
    int match_counter = 0;
    for (int iunmatch = 0; iunmatch < unmatch_list.size(); ++iunmatch)
    {
      erhic::ParticleMC* particle = evt->GetTrack(unmatch_list[iunmatch]-1);
      if (fabs((particle->GetPx()-constit[iconstit].px())/particle->GetPx()) < 0.005&& fabs((particle->GetPy()-constit[iconstit].py())/particle->GetPy()) < 0.005 && fabs((particle->GetPz()-constit[iconstit].pz())/particle->GetPz()) < 0.005)
      {
        match_counter++;
        match_list.push_back(unmatch_list[iunmatch]);
        unmatch_list.erase(unmatch_list.begin()+iunmatch);

        int verbosity = 1; // to control how much info you want to print on the screen
        if (verbosity>1)
        {
          cout<<"match found!"<<endl;
         cout<<"constituent mom "<<constit[iconstit].px()<<", "<<constit[iconstit].py()<<", "<<constit[iconstit].pz()<<endl;
          cout<<"matched mom "<<particle->GetPx()<<", "<<particle->GetPy()<<", "<<particle->GetPz()<<endl;
        }

        break;
      }
     else continue;
    }
    if (match_counter==0) 
    {
      cout << "no match found" << endl;
      match_list.push_back(-1);
    }
  }
}

int find_quark_origin(erhic::EventPythia* evt, erhic::ParticleMC* part)
{
  // is current particle quark or anti-quark?
  if (part->Id()>=1 && part->Id()<=8) return part->GetIndex(); // current particle is quark
  if (part->Id()<=-1 && part->Id()>=-8) return part->GetIndex(); // current particle is anti-quark

  // if not, does this current particle have a parent?
  int parent_indx = part->GetParentIndex();
  if (parent_indx==0) return 0;
  int verbosity=2;
  erhic::ParticleMC* parent = evt->GetTrack(parent_indx-1);
  if (verbosity>1)
  {
    //if (part->GetParticle(parent->Id())) cout<< "parent id is " << part->GetParticle(parent->Id())->GetName() << endl;
    //else cout<< "parent id is " << parent->Id() << endl;
  }
  
  return find_quark_origin(evt, parent);
}

int find_origin(erhic::EventPythia* evt, erhic::ParticleMC* part)
{
  //
  int event_index= part->GetIndex();
  if (event_index>10)
    {
       cout<<event_index;
       int parent_indx = part->GetParentIndex();
       if (parent_indx==0) return 0;
       erhic::ParticleMC* parent = evt->GetTrack(parent_indx-1);
      return find_origin(evt, parent);
    }
  else if(event_index<=10)return event_index;
}

double cal_angle(erhic::ParticleMC* part1,erhic::ParticleMC* part2)
{
return sqrt((part1->GetPhi()-part2->GetPhi())*(part1->GetPhi()-part2->GetPhi())+(part2->GetEta()-part2->GetEta())*(part1->GetEta()-part2->GetEta()));
}
