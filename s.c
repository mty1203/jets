R__LOAD_LIBRARY(libfastjet);
#include "fastjet/ClusterSequenceArea.hh"
#include"style.h"
#include <fastjet/PseudoJet.hh>
#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include <cmath>
#include "TLorentzVector.h"
#include <fstream>
#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>+;
#endif
void s(double threshold);
using std::string;
using std::stringstream;
using namespace std;
using namespace fastjet;
TString flavor(int pdgid);
double cal_angle(erhic::ParticleMC* part1,erhic::ParticleMC* part2);
void find_matching_index(erhic::EventPythia* evt, vector<PseudoJet>& constit, vector<int>& match_list);
int find_quark_origin(erhic::EventPythia* evt, erhic::ParticleMC* part);
static int jet_number=1;
int find_origin(erhic::EventPythia* evt, erhic::ParticleMC* part);
double result[3][3];

void plot()
{  loadStyle();
   vector<double> values = vector<double>();
   vector<double> replot[3][3];
   for (int q=0;q<3;q++)
     for (int c = 0;c<3;c++)
      replot[q][c]=vector<double>();
   for (double i=0.08;i<=0.5;i+=0.02)
   {
     values.push_back(i);
   }
   for (int j=0,j<values.size(),j++)
   {
     s(values.at(j));
     for (int q=0;q<3;q++)
     for (int c = 0;c<3;c++)
       replot[q][c].push_back(result[q][c]);
   }
   TGraph *re[3][3];
   for (int ihigh=0;ihigh<3;ihigh++)
     (int jhigh=0;jhigh<3;jhigh++)
       re[ihigh][jhigh] = new TGraph(values.size(),values,replot[ihigh][jhigh]);
int k=0;
TCanvas *c1 = new TCanvas("c1");
c1->Divide(3,3);
for (int ieta=0;ieta<3;ieta++) 
  {
    for (int ipt=0;ipt<3;ipt++) 
    {
      k++;
      c1->cd(k);
      c1->cd(k)->SetLog(1);
      string ptcut ="<p_{t}<";
      string etacut = "<eta<";
      stringstream ss;
      ss<<ptlow[ipt]<<"GeV"<<ptcut<<pthigh[ipt]<<"GeV"<<"  "<<etalow[ieta]<<etacut<<etahigh[ieta];
      re[ieta][ipt]->SetTitle(ss.str().c_str());
      stringstream mean;
      
      re[ieta][ipt]->GetXaxis()->SetTitle("Threshold");
      re[ieta][ipt]->GetYaxis()->SetTitle("Energy Fraction");
      re[ieta][ipt]->Draw();
    }
  }

}
void s(double threshold)
{
  
  
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
        fraction[ihigh][jhigh] = new TH1D(histname,"",55,0.,1.1);

       }
    }  
  
  //initialize every object
  string line;
  ifstream myfile;
  double pt=0,eta=0;
  TLorentzVector struck = TLorentzVector(0,0,0,0);
  TLorentzVector all_particle = TLorentzVector(0,0,0,0);
  myfile.open("jets_info1.txt"); 
  while(getline(myfile, line)) 
  {
    if(line=="######")
    {
      //calculate energy fraction and fill histogram accordingly

      double den=all_particle.E();
      double energy_arose_quark=struck.E()/den;
      //double photon_ratio=remnant.E()/den;
       for (int idex=0;idex<3;idex++)
        {
          for(int jdex=0;jdex<3;jdex++)
           {
             if((pt>=ptlow[jdex]&&pt<pthigh[jdex])&&(eta>=etalow[idex]&&eta<etahigh[idex]))
                fraction[idex][jdex]->Fill(energy_arose_quark);
           }
        }
      
      // clear all
       pt=0;
       eta=0;
       struck = TLorentzVector(0,0,0,0);
       all_particle = TLorentzVector(0,0,0,0);
    }
    else
    {
     double w, x, y, z,q;
     int m;
     sscanf(line.c_str(),"%d %lf %lf %lf %lf %lf %lf %lf",&m,&w,&x,&y,&z,&q,&pt,&eta);
     TLorentzVector vec = TLorentzVector(x,y,z,q);
     if (w<=threshold) struck+=vec;  //threshold is 0.3
     all_particle+=vec;
     cout<<m<<" "<<w<<endl;
    }   
  }
int k = 0;
for (int ieta=0;ieta<3;ieta++) 
  {
    for (int ipt=0;ipt<3;ipt++) 
    {
     double mean = fraction[ieta][ipt]->GetMean();
     result[ieta][ipt]=mean;
    }
  }


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