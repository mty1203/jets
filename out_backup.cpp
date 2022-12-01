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
TString flavor(int pdgid);
double cal_angle(erhic::ParticleMC* part1,erhic::ParticleMC* part2);
void find_matching_index(erhic::EventPythia* evt, vector<PseudoJet>& constit, vector<int>& match_list);
int find_quark_origin(erhic::EventPythia* evt, erhic::ParticleMC* part);
static int jet_number=1;
int find_origin(erhic::EventPythia* evt, erhic::ParticleMC* part);

void out_backup()
{
  loadStyle();
  TFile * theFile = new TFile("data.root","RECREATE");
  // TTree *T = new TTree("Tree","jets angle diff");
  // vector<double> angle = vector<double>();
  // vector<double> momentumx=vector<double>();
  // vector<double> momentumy= vector<double>();
  // vector<double> momentumz= vector<double>();
  // vector<double> Tenergy=  vector<double>();
  // vector<Int_t> pdgnumber=vector<Int_t>();
  // vector<Int_t> entrynumber=vector<Int_t>();
  // T->Branch("momentumx","vector<double>",&momentumx);
  // T->Branch("momentumy","vector<double>",&momentumy);
  // T->Branch("momentumz","vector<double>",&momentumz);
  // T->Branch("Tenergy","vector<double>",&Tenergy);
  // T->Branch("angle_diff","vector<double>",&angle);
  // T->Branch("pdg","vector<Int_t>",&pdgnumber);
  // T->Branch("entry","vector<Int_t>",&entrynumber);

  ofstream myfile;
  TRandom3 *rndm = new TRandom3(0);
  myfile.open ("./data/jet_nodet_10_100_1_heavy.txt"); 
  //myfile << "data begin\n"; 
  erhic::EventPythia *evt(NULL); 
  erhic::ParticleMC *particle(NULL); 
  erhic::ParticleMC *quark(NULL); 
  TFile *f = new TFile("./ep_10_100_norad_1million.root");
  TTree *tree = (TTree*)f->Get("EICTree");
  TDatabasePDG *pdgn = new TDatabasePDG();
  Int_t nEntries = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nEntries<<endl<<endl;

  tree->SetBranchAddress("event",&evt);
  // TH2D *fraction[3][3];
  // char *histname = new char[10];
  // double etahigh[3]={-1,1,3.5};
  // double etalow[3] = {-3.5,-1,1};
  // double pthigh[3]={5.,10.,20.};
  // double ptlow[3]={0.,5.,10.};
  // for(int ihigh=0;ihigh<3;ihigh++)
  //   {
  //     for(int jhigh=0;jhigh<3;jhigh++)
  //      {
  //        sprintf(histname, "h_x_%d_%d",ihigh,jhigh); 
  //       fraction[ihigh][jhigh] = new TH2D(histname,"",50,-3,1.,100,1.,3.);
  //      }
  //   }  
  // TH1F *eta = new TH1F("eta","",100,-5,5);
  // TH1F *phi = new TH1F("phi","",50,0,2*3.1415); 
  // TH1F *area = new TH1F("area","",50,0,10);
  // TH1F *conts = new TH1F("conts","",25,0,25); 
  // TH2F *h1=new TH2F("","",40,0,20,25,0,25);
  // TH1F *ep = new TH1F("ep","",50,2,4); 
  // TH1F *ch = new TH1F("ch","",10,0,10); 
  // TH1F *charged_hadron = new TH1F("ch","",50,0,1.1); 
  // TH1F *neutral = new TH1F("neutral","",50,0,1.1); 
  // TH1F *photon = new TH1F("photon","",50,0,1.1); 
  // TH1F *theta = new TH1F("theta","",50,-0.1,0.1); 
  // TH2F *thetajet = new TH2F("thetajet","",50,-5,5,20,0,20); 
  // TH1F *strucku= new TH1F("u","",50,0.,30.);
  // TH1F *struckd=new TH1F("d","",50,0.,30.);
  // TH2F *E = new TH2F("E","",50,0,1,50,0,200); 
 int counter_K_plus = 0, counter_K_minus = 0, counter_K_L = 0,counter_pi_plus =0,counter_proton=0,counter_pi=0,counter_photon=0,counter_rho=0,counter_rho_plus=0,counter_neutron=0,counter_pi_minus=0;
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
quark= evt->GetTrack(9);
double struck_pt = quark->GetPt();
double quark_eta = quark->GetEta();
double quark_phi = quark->GetPhi();
// double q2 = evt->GetTrueQ2();
// double x  =evt->GetTrueX();
//  for (int idex=0;idex<3;idex++)
//         {for(int jdex=0;jdex<3;jdex++)
//            {
//              if((struck_pt>=ptlow[jdex]&&struck_pt<pthigh[jdex])&&(quark_eta>=etalow[idex]&&quark_eta<etahigh[idex]))
//                 fraction[idex][jdex]->Fill(log10(x),log10(q2));
//            }
//         }



for (unsigned ijet = 0; ijet < jets1.size(); ijet++)
{  
    
    int verbosity=0;
    double counter_neutral_hadron = 0;
    double counter_charged_hadron = 0;
    double counter_photon = 0;
      double den=jets1[ijet].e();
      double jets_eta = jets1[ijet].eta();
      double jets_phi = jets1[ijet].phi();
      double jets_pt = jets1[ijet].pt();
      double deltaR =sqrt((jets_eta -quark_eta)*(jets_eta -quark_eta)+(jets_phi-quark_phi)*(jets_phi-quark_phi));
      if(deltaR<0.3&&(quark->Id()>3||quark->Id()<-3)&&jets_pt>=5)myfile<<"####"<<endl;
  //  if (verbosity>1) cout << "jet " << ijet << " (pt mass eta phi): "<< jets[ijet].pt() << " " << jets[ijet].m() << " "  << jets[ijet].rap() << " " << jets[ijet].phi() << endl;
    vector<PseudoJet> constituents = jets1[ijet].constituents(); // get jet constituents
   // match the constitutes to TTree indices
    vector<int> matching_index;
    find_matching_index(evt, constituents, matching_index);
    double charged_hadrons=0, photons=0,neutral_hadrons=0;
    double proton=0,kion=0,pion=0;
    TParticlePDG *pdg;
    double strcuk_quark_energy=0;
    TLorentzVector photon_vector = TLorentzVector(0,0,0,0);
    TLorentzVector pi_plus_vector = TLorentzVector(0,0,0,0);
    TLorentzVector pi_minus_vector = TLorentzVector(0,0,0,0);
    TLorentzVector kion_plus_vector = TLorentzVector(0,0,0,0);
    TLorentzVector kion_minus_vector = TLorentzVector(0,0,0,0);
    TLorentzVector proton_vector = TLorentzVector(0,0,0,0);
    TLorentzVector antiproton_vector = TLorentzVector(0,0,0,0);
     TLorentzVector charged_hadron_vector = TLorentzVector(0,0,0,0);
    double other_energy=0;
    TLorentzVector neutral_hadron_vector = TLorentzVector(0,0,0,0);
    int counter=0; // at least two particles originating from struck quark
    //
    for (unsigned iconstit = 0; iconstit < constituents.size(); iconstit++)
    {   
        int track_index = matching_index[iconstit];
        
      //   if (track_index>0)
      //   {    TString part;
      //       Int_t track_id = evt->GetTrack(track_index-1)->Id();   
      // //      pdgnumber.push_back(track_id);  
      //       erhic::ParticleMC* particle = evt->GetTrack(track_index-1);
      //       TLorentzVector vec = particle->PxPyPzE();//get 4 momentum
      //     //  momentumx.push_back(particle->GetPx());
      //    //   momentumy.push_back(particle->GetPy());
      //    //   momentumz.push_back(particle->GetPz());
      //  //     Tenergy.push_back(particle->GetE());
      //   //    entrynumber.push_back(j);
      //       double energy = particle->GetE();
      //       //int quark_origin = find_quark_origin(evt,particle);
      //       //int origin = find_origin(evt,particle);
      //       double angle_diff = cal_angle(particle,quark);
      //       angle.push_back(angle_diff);
      //       //if(jets1[ijet].eta()>=-3.5&&jets1[ijet].eta()<-1&&jets1[ijet].pt()>=5) 
             
      //       if(angle_diff<=0.2) 
      //       { 
      //         struck+=vec;
      //         counter++;
      //       }
      //       else remnant+=vec;
      //   }
       if (track_index>0)
        {    TString part;
            Int_t track_id = evt->GetTrack(track_index-1)->Id();     
            //TParticlePDG *
            
            //double ran1=rndm->Gaus(1,0.04);
            if(deltaR<0.3&&(quark->Id()>3||quark->Id()<-3)&&jets_pt>=5)
            {myfile<<quark->Id()<<" "<<evt->GetTrack(track_index-1)->GetPt()<<" "<<evt->GetTrack(track_index-1)->GetEta()-jets_eta<<" "<<evt->GetTrack(track_index-1)->GetPhi()-jets_phi<<" "<<track_id<<" "<<deltaR<<" "<<jets_pt<<" "<<jets_eta<<" "<<" "<<jets_phi<<" "<<den<<" "<<evt->GetTrack(track_index-1)->GetE()<<endl;}
            // pdg =  pdgn->GetParticle(track_id);
            // if(pdg== NULL) {cout<<"NULL pointer"<<endl; continue;}
            // if(track_id ==211 ||track_id==-211) {pi_plus_vector+=evt->GetTrack(track_index-1)->Get4Vector();pion++;}
            // if(track_id==321||track_id==-321){kion_plus_vector+=evt->GetTrack(track_index-1)->Get4Vector();kion++;}
            // if(track_id==2212||track_id==-2212){proton_vector+=evt->GetTrack(track_index-1)->Get4Vector();proton++;}
            // part = pdg->ParticleClass();
            // double charge = pdg->Charge();
            // if(part=="Meson"||part=="Baryon")
            //    {if(charge!=0) {counter_charged_hadron++;charged_hadron_vector+= evt->GetTrack(track_index-1)->Get4Vector();}
            //   else if(charge==0){counter_neutral_hadron++;neutral_hadron_vector+= evt->GetTrack(track_index-1)->Get4Vector();}}
            // if(track_id==22) {counter_photon++;photon_vector+=evt->GetTrack(track_index-1)->Get4Vector();}
            // cout<<"track_id"<<track_id<<endl;
        }
      }
  //  T->Fill(); 
  //  pdgnumber.clear();
  //  angle.clear();
 //   momentumx.clear();
  //  momentumy.clear();
  //  momentumz.clear();
    //Tenergy.clear();

      //den = counter_charged_hadron+counter_neutral_hadron
      //double charged_hadron_ratio=struck.E()/den;
     // double photon_ratio=remnant.E()/den;
      //cout<<"struck energy = "<<struck.E()<<"total = "<<den<<endl;
      //double neutral_hadron_ratio=other_energy/den; 
      //cout<<charged_hadron_ratio<<endl;
      // for (int idex=0;idex<3;idex++)
      //   {for(int jdex=0;jdex<3;jdex++)
      //      {
      //        if((jets1[ijet].pt()>=ptlow[jdex]&&jets1[ijet].pt()<pthigh[jdex])&&(jets1[ijet].eta()>=etalow[idex]&&jets1[ijet].eta()<etahigh[idex]))
      //           fraction[idex][jdex]->Fill(charged_hadron_ratio);
      //      }
      //   }
     // counter_charged_hadron = charged_hadron_vector.Pt();
      //counter_neutral_hadron = neutral_hadron_vector.Pt();
      //counter_photon = photon_vector.Pt();
     //myfile<<quark->Id()<<" "<<counter_charged_hadron<<" "<<counter_neutral_hadron<<" "<<counter_photon<<" "<<deltaR<<" "<<struck_pt<<" "<<quark_eta<<" "<<pi_plus_vector.Pt()<<" "<<kion_plus_vector.Pt()<<" "<<proton_vector.Pt()<<" "<<pion<<" "<<kion<<" "<<proton<<endl;

  //  }
}
// TFile *file = new TFile("./ep_10_100_norad_2.root");
//   TTree *tree1 = (TTree*)file->Get("EICTree");
//    Int_t nEntries1 = tree1->GetEntries();
//    for(Int_t j=0;j<nEntries1;j++){ 
//    cout<<"***********************reading event"<<j<<"*********************"<<endl;
//      tree1->GetEntry(j);
//     int process=evt->GetProcess();
//    if(process!=99) continue;
//   vector<PseudoJet> jet_constits;
//   for (int ipart = 0; ipart < evt->GetNTracks(); ++ipart)
// {  particle = evt->GetTrack(ipart);
//    if (particle->GetStatus()==1 && fabs(particle->GetEta())<3.5 && particle->Id()!=11) 
//    jet_constits.push_back(PseudoJet(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE()));
// }

// JetDefinition R1jetdef(antikt_algorithm, 0.6); // anti-kt and R = 1
// //ClusterSequence cs(jet_constits, R1jetdef);
// //vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
// double maxrap = 8.0;
// unsigned int n_repeat = 1; // default is 1
// double ghost_area = 0.01; // this is the default
// //ClusterSequenceArea clust_seq(jet_constits, R1jetdef, area_def);
// const ClusterSequenceArea clust_seq(jet_constits, R1jetdef, VoronoiAreaSpec());
// vector<PseudoJet> jets1 = sorted_by_pt(clust_seq.inclusive_jets());

// //find scattered electron

// particle = evt->GetTrack(2);
// double phi_electron = particle->GetPhi();
// int counter_quarks=0;
// TParticlePDG *pdg1;
// double phi_quark;
// quark= evt->GetTrack(9);
// double struck_pt = quark->GetPt();
// double quark_eta = quark->GetEta();
// double quark_phi = quark->GetPhi();
// double q2 = evt->GetTrueQ2();
// double x  =evt->GetTrueX();
//  for (int idex=0;idex<3;idex++)
//         {for(int jdex=0;jdex<3;jdex++)
//            {
//              if((struck_pt>=ptlow[jdex]&&struck_pt<pthigh[jdex])&&(quark_eta>=etalow[idex]&&quark_eta<etahigh[idex]))
//                 fraction[idex][jdex]->Fill(log10(x),log10(q2));
//            }
//         }


// for (unsigned ijet = 0; ijet < jets1.size(); ijet++)
// {  
//     myfile<<"######"<<endl;
//     int verbosity=0;
//     int counter_neutral_hadron = 0;
//     int counter_charged_hadron = 0;
//     int counter_photon = 0;
//       double den=jets1[ijet].e();
//       double jets_eta = jets1[ijet].eta();
//       double jets_phi = jets1[ijet].phi();
//       double deltaR =sqrt((jets_eta -quark_eta)*(jets_eta -quark_eta)+(jets_phi-quark_phi)*(jets_phi-quark_phi));
//   //  if (verbosity>1) cout << "jet " << ijet << " (pt mass eta phi): "<< jets[ijet].pt() << " " << jets[ijet].m() << " "  << jets[ijet].rap() << " " << jets[ijet].phi() << endl;
//     vector<PseudoJet> constituents = jets1[ijet].constituents(); // get jet constituents
//    // match the constitutes to TTree indices
//     vector<int> matching_index;
//     find_matching_index(evt, constituents, matching_index);
//     double charged_hadrons=0, photons=0,neutral_hadrons=0;
//     TParticlePDG *pdg;
//     double strcuk_quark_energy=0;
//     TLorentzVector photon_vector = TLorentzVector(0,0,0,0);
//      TLorentzVector charged_hadron_vector = TLorentzVector(0,0,0,0);
//     double other_energy=0;
//     TLorentzVector neutral_hadron_vector = TLorentzVector(0,0,0,0);
//     int counter=0; // at least two particles originating from struck quark
//     //
//     for (unsigned iconstit = 0; iconstit < constituents.size(); iconstit++)
//     {   
//         int track_index = matching_index[iconstit];
        
//       //   if (track_index>0)
//       //   {    TString part;
//       //       Int_t track_id = evt->GetTrack(track_index-1)->Id();   
//       // //      pdgnumber.push_back(track_id);  
//       //       erhic::ParticleMC* particle = evt->GetTrack(track_index-1);
//       //       TLorentzVector vec = particle->PxPyPzE();//get 4 momentum
//       //     //  momentumx.push_back(particle->GetPx());
//       //    //   momentumy.push_back(particle->GetPy());
//       //    //   momentumz.push_back(particle->GetPz());
//       //  //     Tenergy.push_back(particle->GetE());
//       //   //    entrynumber.push_back(j);
//       //       double energy = particle->GetE();
//       //       //int quark_origin = find_quark_origin(evt,particle);
//       //       //int origin = find_origin(evt,particle);
//       //       double angle_diff = cal_angle(particle,quark);
//       //       angle.push_back(angle_diff);
//       //       //if(jets1[ijet].eta()>=-3.5&&jets1[ijet].eta()<-1&&jets1[ijet].pt()>=5) 
             
//       //       if(angle_diff<=0.2) 
//       //       { 
//       //         struck+=vec;
//       //         counter++;
//       //       }
//       //       else remnant+=vec;
//       //   }
//        if (track_index>0)
//         {    TString part;
//             Int_t track_id = evt->GetTrack(track_index-1)->Id();     
//             //TParticlePDG *
//             pdg =  pdgn->GetParticle(track_id);
//             if(pdg== NULL) {cout<<"NULL pointer"<<endl; continue;}
//             part = pdg->ParticleClass();
//             double charge = pdg->Charge();
//             if(part=="Meson"||part=="Baryon")
//                {if(charge!=0) {counter_charged_hadron++;charged_hadron_vector+= evt->GetTrack(track_index-1)->Get4Vector();}
//               else if(charge==0){counter_neutral_hadron++;neutral_hadron_vector+= evt->GetTrack(track_index-1)->Get4Vector();}}
//             if(track_id==22) {counter_photon++;photon_vector+=evt->GetTrack(track_index-1)->Get4Vector();}
//             cout<<"track_id"<<track_id<<endl;
//         }
//       }
//   //  T->Fill(); 
//   //  pdgnumber.clear();
//   //  angle.clear();
//  //   momentumx.clear();
//   //  momentumy.clear();
//   //  momentumz.clear();
//     //Tenergy.clear();

//       //den = counter_charged_hadron+counter_neutral_hadron
//       //double charged_hadron_ratio=struck.E()/den;
//      // double photon_ratio=remnant.E()/den;
//       //cout<<"struck energy = "<<struck.E()<<"total = "<<den<<endl;
//       //double neutral_hadron_ratio=other_energy/den; 
//       //cout<<charged_hadron_ratio<<endl;
//       // for (int idex=0;idex<3;idex++)
//       //   {for(int jdex=0;jdex<3;jdex++)
//       //      {
//       //        if((jets1[ijet].pt()>=ptlow[jdex]&&jets1[ijet].pt()<pthigh[jdex])&&(jets1[ijet].eta()>=etalow[idex]&&jets1[ijet].eta()<etahigh[idex]))
//       //           fraction[idex][jdex]->Fill(charged_hadron_ratio);
//       //      }
//       //   }
//       counter_charged_hadron = charged_hadron_vector.Pt();
//       counter_neutral_hadron = neutral_hadron_vector.Pt();
//       counter_photon = photon_vector.Pt();
//       myfile<<quark->Id()<<" "<<counter_charged_hadron<<" "<<counter_neutral_hadron<<" "<<counter_photon<<" "<<deltaR<<" "<<struck_pt<<" "<<quark_eta<<endl;

//    }
// }

// TCanvas *c1 = new TCanvas("c1");
// c1->Divide(3,3);
//   int k = 0;
// for (int ieta=0;ieta<3;ieta++) 
//   {
//     for (int ipt=0;ipt<3;ipt++) 
//     {
//       k++;
//       c1->cd(k);
//       c1->cd(k)->SetLogz();
//       //c1->cd(k)->SetLogy();
//       //c1->cd(k)->SetLogx();
//       string ptcut ="<p_{t}^{sq}<";
//       string etacut = "<#eta^{sq}<";
//       stringstream ss;
//       ss<<ptlow[ipt]<<"GeV"<<ptcut<<pthigh[ipt]<<"GeV"<<"  "<<etalow[ieta]<<etacut<<etahigh[ieta];
//       fraction[ieta][ipt]->SetTitle(ss.str().c_str());
//       stringstream mean;
      
//       fraction[ieta][ipt]->GetXaxis()->SetTitle("lg(x)");
//       fraction[ieta][ipt]->GetYaxis()->SetTitle("lg(Q^{2})");
//       fraction[ieta][ipt]->Draw("colz");
      
//       TLatex latex;
//       latex.SetTextAlign(13);  //align at top
//       mean<<"mean = "<<fraction[ieta][ipt]->GetEntries();
//       TPaveText *pave = new TPaveText(.8,.8,.95,.95,"brNDC");
//       pave->SetTextAlign(12);
//       char temp[20];
//       sprintf(temp,"mean= %5.0f",fraction[ieta][ipt]->GetEntries());
//       TText *text = pave->AddText(temp);
//       text->SetTextFont(72);
//       text->SetTextAlign(22);
      
      //pave->Draw();
      //latex.DrawLatex(.2,.9,mean.str().c_str());
      //stringstream mean;
      //mean<<"mean value = "<<fraction[ieta][ipt]->GetMean();
      //TLatex *l = new TLatex(10,40,mean.str().c_str());
      //l->Draw("same");
    }
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
