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
void find_matching_index(erhic::EventPythia* evt, vector<PseudoJet>& constit, vector<int>& match_list);
int find_quark_origin(erhic::EventPythia* evt, erhic::ParticleMC* part);
static int jet_number=1;
int find_origin(erhic::EventPythia* evt, erhic::ParticleMC* part);
vector<int> list_info =vector<int>() ;
int find_quark(erhic::EventPythia* evt,vector<int> daughter);
void jets()
{
  loadStyle();
  ofstream myfile;
  myfile.open ("jets.txt"); 
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
  TH1D *fraction[3][3];
  char *histname = new char[10];
  double etahigh[3]={-1,1,3.5};
  double etalow[3] = {-3.5,-1,1};
  double pthigh[3]={5.,10.,20.};
  double ptlow[3]={1.,5.,10.};
  for(int ihigh=0;ihigh<3;ihigh++)
    {
      for(int jhigh=0;jhigh<3;jhigh++)
       {
         sprintf(histname, "h_x_%d_%d",ihigh,jhigh); 
        fraction[ihigh][jhigh] = new TH1D(histname,"",50,0.,1.1);
       }
    }  
  TH1F *eta = new TH1F("eta","",100,-5,5);
  TH1F *phi = new TH1F("phi","",50,0,2*3.1415); 
  TH1F *area = new TH1F("area","",50,0,10);
  TH1F *conts = new TH1F("conts","",25,0,25); 
  TH2F *h1=new TH2F("","",40,0,20,25,0,25);
     TH1F *ep = new TH1F("ep","",50,2,4); 
TH1F *ch = new TH1F("ch","",10,0,10); 
TH1F *charged_hadron = new TH1F("ch","",50,0,1.1); 
TH1F *neutral = new TH1F("neutral","",50,0,1.1); 
TH1F *photon = new TH1F("photon","",50,0,1.1); 
  TH1F *theta = new TH1F("theta","",50,-0.1,0.1); 
   TH2F *thetajet = new TH2F("thetajet","",50,-5,5,20,0,20); 
  TH1F *strucku= new TH1F("u","",50,0.,30.);
 TH1F *struckd=new TH1F("d","",50,0.,30.);
 TH2F *E = new TH2F("E","",50,0,1,50,0,200); 
 //int counter_K_plus = 0, counter_K_minus = 0, counter_K_L = 0,counter_pi_plus =0,counter_proton=0,counter_pi=0,counter_photon=0,counter_rho=0,counter_rho_plus=0,counter_neutron=0,counter_pi_minus=0;
  for(Int_t j=0;j<nEntries;j++)
  { 
   cout<<"***********************reading event"<<j<<"*********************"<<endl;
     tree->GetEntry(j);
    int process=evt->GetProcess();
   if(process!=99) continue;
  vector<PseudoJet> jet_constits;
    vector<int> first_daughter_quark = vector<int>();
  for (int ipart = 0; ipart < evt->GetNTracks(); ++ipart)
{  particle = evt->GetTrack(ipart);
   if (particle->GetStatus()==1 && fabs(particle->GetEta())<3.5 && particle->Id()!=11) 
     jet_constits.push_back( PseudoJet(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE()) );
     if(particle->GetParentIndex()==10)
     first_daughter_quark.push_back(particle->GetIndex());
}
  
JetDefinition R1jetdef(antikt_algorithm, 0.6); // anti-kt and R = 1
//ClusterSequence cs(jet_constits, R1jetdef);
//vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
double maxrap = 8.0;
unsigned int n_repeat = 1; // default is 1
double ghost_area = 0.01; // this is the default
//ClusterSequenceArea clust_seq(jet_constits, R1jetdef, area_def);
const ClusterSequenceArea clust_seq(jet_constits, R1jetdef, VoronoiAreaSpec());
vector<PseudoJet> jets1 = sorted_by_pt(clust_seq.inclusive_jets());

//find scattered electron

particle = evt->GetTrack(2);
double phi_electron = particle->GetPhi();
int counter_quarks=0;
TParticlePDG *pdg1;
double phi_quark;
particle= evt->GetTrack(9);
phi_quark = particle->GetPhi();
if (evt->GetTrack(9)->Id()==2){strucku->Fill(evt->GetTrueX()*evt->GetTrueQ2());E->Fill(evt->GetTrueX(),evt->GetTrueQ2());} 
else if (evt->GetTrack(9)->Id()==1)struckd->Fill(evt->GetTrueX()*evt->GetTrueQ2());

for (int ipart = 0;ipart< evt->GetNTracks(); ++ipart)
{ particle = evt->GetTrack(ipart);
 Int_t track_number = evt->GetTrack(ipart)->Id();  
 pdg1 =  pdgn->GetParticle(abs(track_number));
 if(pdg1==NULL) {counter_quarks++; continue;}
  
  TString n = pdg1->ParticleClass();
 if(n == "Quark") ++counter_quarks;
 if(counter_quarks==3)
  { phi_quark = particle->GetPhi();
    break;
  }
}
    erhic::ParticleMC* sq = evt->GetTrack(find_quark(evt,first_daughter_quark)-1);
    int sqdaughter = sq->GetChild1Index(); //string
    int daughter_begin = evt->GetTrack(sqdaughter-1)->GetChild1Index();
    int daughter_end = evt->GetTrack(sqdaughter-1)->GetChildNIndex();
    vector<int> string_daughter = vector<int>(); 
    for (int idaughter=daughter_begin;idaughter<=daughter_end;idaughter++)
    {
     erhic::ParticleMC* particle_candidates = evt->GetTrack(idaughter-1);
     if(particle_candidates->GetPz()<0)
        string_daughter.push_back(idaughter);
    if(particle_candidates->GetPz()>0&&particle_candidates->GetPz()<0.3&&particle_candidates->GetPt()>0.6)
      string_daughter.push_back(idaughter);
    }
for (unsigned ijet = 0; ijet < jets1.size(); ijet++)
{  
    myfile<<"######"<<endl;
    int verbosity=0;
  //  if (verbosity>1) cout << "jet " << ijet << " (pt mass eta phi): "<< jets[ijet].pt() << " " << jets[ijet].m() << " "  << jets[ijet].rap() << " " << jets[ijet].phi() << endl;
    vector<PseudoJet> constituents = jets1[ijet].constituents(); // get jet constituents
   // match the constitutes to TTree indices
    vector<int> matching_index;
    find_matching_index(evt, constituents, matching_index);
    double charged_hadrons=0, photons=0,neutral_hadrons=0;
    TParticlePDG *pdg;
    double strcuk_quark_energy=0;
     TLorentzVector remnant = TLorentzVector(0,0,0,0);
    double other_energy=0;
    TLorentzVector struck = TLorentzVector(0,0,0,0);
    int counter=0; // at least two particles originating from struck quark
    cout<<"size of con"<<constituents.size()<<endl;


    cout<<"get child"<<sqdaughter<<endl;
    for (unsigned iconstit = 0; iconstit < constituents.size(); iconstit++)
    {   
        int track_index = matching_index[iconstit];
        
        if (track_index>0)
        {    TString part;
            Int_t track_id = evt->GetTrack(track_index-1)->Id();     
            erhic::ParticleMC* particle = evt->GetTrack(track_index-1);
            TLorentzVector vec = particle->PxPyPzE();//get 4 momentum
            double energy = particle->GetE();    
            int origin = find_origin(evt,particle);
            cout<<"in the info list"<<list_info.size()<<endl;
            //if(std::find(list_info.begin(), list_info.end(), sqdaughter) != list_info.end())
            //origin =10;
          if(std::find(string_daughter.begin(), string_daughter.end(), track_index) != string_daughter.end())
            {
              origin = 10;
            } 
            for (int i_info=0;i_info<list_info.size();i_info++)
            {
              if(std::find(string_daughter.begin(), string_daughter.end(), list_info[i_info]) != string_daughter.end()) 
               {
               origin = 10;
               break;
               }
            } 
            cout<< origin;
            //myfile<<particle->GetPt()<<" "<<particle->GetRapidity()<<" "<<particle->GetPhi()<<" "<<track_id<<" "<<" "<<track_index<<" "<<origin<<"   energy"<<vec.E()<<"  "<<"jet pt"<<" "<<jets1[ijet].pt()<<"  "<<jets1[ijet].eta()<<" event "<<j<<endl;
            if(origin ==10) 
            { 
              struck+=vec;
              counter++;
            }
            else remnant+=vec;
            pdg =  pdgn->GetParticle(track_id);
           
         list_info.clear();
        }
        
        
    } 
     
   
      double den=jets1[ijet].e();
      double charged_hadron_ratio=struck.E()/den;
      double photon_ratio=remnant.E()/den;
      cout<<"struck energy = "<<struck.E()<<"total = "<<den<<endl;
      double neutral_hadron_ratio=other_energy/den; 
      if((jets1[ijet].pt()>=10&&jets1[ijet].pt()<20)&&(jets1[ijet].eta()>=-3.5&&jets1[ijet].eta()<-1))
        eta->Fill(fabs(jets1[ijet].phi()-phi_electron)-3.1415926);
     cout<<charged_hadron_ratio<<endl;
      for (int idex=0;idex<3;idex++)
        {for(int jdex=0;jdex<3;jdex++)
           {
             if((jets1[ijet].pt()>=ptlow[jdex]&&jets1[ijet].pt()<pthigh[jdex])&&(jets1[ijet].eta()>=etalow[idex]&&jets1[ijet].eta()<etahigh[idex]))
                fraction[idex][jdex]->Fill(photon_ratio);
           }
        }
   }
 }
myfile.close();
int k = 0;
TCanvas *c1 = new TCanvas("c1");
c1->Divide(3,3);
for (int ieta=0;ieta<3;ieta++) 
  {
    for (int ipt=0;ipt<3;ipt++) 
    {
      k++;
      c1->cd(k);
      c1->cd(k)->SetLogy();
      string ptcut ="<p_{t}<";
      string etacut = "<eta<";
      stringstream ss;
      ss<<ptlow[ipt]<<"GeV"<<ptcut<<pthigh[ipt]<<"GeV"<<"  "<<etalow[ieta]<<etacut<<etahigh[ieta];
      fraction[ieta][ipt]->SetTitle(ss.str().c_str());
      stringstream mean;
      
      fraction[ieta][ipt]->GetXaxis()->SetTitle("Energy");
      fraction[ieta][ipt]->GetYaxis()->SetTitle("Number of jets");
      fraction[ieta][ipt]->Draw();
      TLatex latex;
      latex.SetTextAlign(13);  //align at top
      mean<<"mean = "<<fraction[ieta][ipt]->GetMean();
      TPaveText *pave = new TPaveText(.8,.8,.95,.95,"brNDC");
      pave->SetTextAlign(12);
      char temp[20];
      sprintf(temp,"mean= %3.2f",fraction[ieta][ipt]->GetMean());
      TText *text = pave->AddText(temp);
      text->SetTextFont(72);
      text->SetTextAlign(22);
      
      pave->Draw();
      //latex.DrawLatex(.2,.9,mean.str().c_str());
      //stringstream mean;
      //mean<<"mean value = "<<fraction[ieta][ipt]->GetMean();
      //TLatex *l = new TLatex(10,40,mean.str().c_str());
      //l->Draw("same");
    }
  }
TCanvas *C4 = new TCanvas();
C4->SetLogz();
thetajet->SetStats(0);
thetajet->GetXaxis()->SetTitle("eta");
thetajet->GetXaxis()->CenterTitle();
thetajet->GetYaxis()->SetTitle("Number of particles originating from struck quark");
thetajet->GetYaxis()->CenterTitle();
thetajet->Draw("colz");
}

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
  //find the struck quark daughter
  //erhic::ParticleMC* sq = evt->GetTrack(14);
  //int sqdaughter = sq->GetChild1Index();
  int event_index= part->GetIndex();
  if (event_index>10)
    {
      //cout<<event_index;
       int parent_indx = part->GetParentIndex();
       list_info.push_back(parent_indx);
       cout<<"parent_indx"<<parent_indx<<endl;
       if (parent_indx==0) return 0;
       erhic::ParticleMC* parent = evt->GetTrack(parent_indx-1);
      return find_origin(evt, parent);
    }
  else if(event_index<=10)return event_index;
}
int find_quark(erhic::EventPythia* evt,vector<int> daughter)
{
  for (int i = 0;i<daughter.size();i++)
    {
      int number = evt->GetTrack(daughter[i]-1)->Id();
      if(number==evt->GetTrack(9)->Id())
        return daughter[i];
    }
}
