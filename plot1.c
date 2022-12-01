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
template <typename T, typename A>
int arg_max(std::vector<T, A> const& vec) {
  return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
}
using std::string;
using std::stringstream;
using namespace std;
using namespace fastjet;
TString flavor(int pdgid);
void s(double threshold);
double cal_angle(erhic::ParticleMC* part1,erhic::ParticleMC* part2);
void find_matching_index(erhic::EventPythia* evt, vector<PseudoJet>& constit, vector<int>& match_list);
int find_quark_origin(erhic::EventPythia* evt, erhic::ParticleMC* part);
static int jet_number=1;
int find_origin(erhic::EventPythia* evt, erhic::ParticleMC* part);
static int counter;
TDatabasePDG *pdgn = new TDatabasePDG();
void plot1()
{  loadStyle();
s(0.3);
//   double etahigh[3]={-1,1,3.5};
//   double etalow[3] = {-3.5,-1,1};
//   double pthigh[3]={5.,10.,20.};
//   double ptlow[3]={1.,5.,10.};
//    vector<double> values = vector<double>();
//    vector<double> replot[3][3];
//    for (int q=0;q<3;q++)
//      for (int c = 0;c<3;c++)
//       replot[q][c]=vector<double>();
//    for (double i=0.06;i<=0.5;i+=0.02)
//    {
//      values.push_back(i);
//    }
//    for (int j=0;j<values.size();j++)
//    {
//      s(values.at(j));
//      for (int q=0;q<3;q++)
//      for (int c = 0;c<3;c++)
//        replot[q][c].push_back(result[q][c]);
//    }
//    TGraph *re[3][3];
//    for (int ihigh=0;ihigh<3;ihigh++)
//      for (int jhigh=0;jhigh<3;jhigh++)
//        re[ihigh][jhigh] = new TGraph(values.size(),&values[0],&replot[ihigh][jhigh][0]);

// TCanvas *c1 = new TCanvas("c1");
// c1->Divide(3,3);
// int k=0;
// for (int ieta=0;ieta<3;ieta++) 
//   {
//     for (int ipt=0;ipt<3;ipt++) 
//     {
//       k++;
//       c1->cd(k);
//       string ptcut ="<p_{t}<";
//       string etacut = "<eta<";
//       stringstream ss;
//       ss<<ptlow[ipt]<<"GeV"<<ptcut<<pthigh[ipt]<<"GeV"<<"  "<<etalow[ieta]<<etacut<<etahigh[ieta];
//       re[ieta][ipt]->SetTitle(ss.str().c_str());
//       stringstream mean;
//       re[ieta][ipt]->GetHistogram()->SetMaximum(1.);   // along          
//       re[ieta][ipt]->GetHistogram()->SetMinimum(0.);  
//       re[ieta][ipt]->GetXaxis()->SetTitle("Threshold");
//       re[ieta][ipt]->GetYaxis()->SetTitle("Energy Fraction");
//       re[ieta][ipt]->Draw();
//     }
//   }

}
void s(double threshold)
{
 //int flavor[2][3]={1,2,3,-1,-2,-3};
  int flavor[2][3]={4,5,6,-4,-5,-6};
  string os_X[8] = {"p^{-}","K^{-}","#pi^{-}","#pi^{+}","K^{+}","p^{+}","photon","neutral_hadrons"};
 //string os_X[4] = {"charged hadrons","neutral hadrons","leptons","photons"};
  //string os_y[6] = {"d","u","s","#bar{d}","#bar{u}","#bar{s}"};

  string os_y[6] = {"c","b","t","#bar{c}","#bar{b}","#bar{t}"};
  TH1D *fraction[2][3];
  TH1D *fraction1[2][3];
  TH1D *fraction2[2][3];
  TH1D *fraction3[2][3];
  TH1D *fraction4[2][3];
  TH1D *fraction5[2][3];
  TH1D *fraction6[2][3];
  TH1D *fraction7[2][3];
  char *histname = new char[10];
  double etahigh[3]={-1,1,3.5};
  double etalow[3] = {-3.5,-1,1};
  double pthigh[3]={5.,10.,20.};
  double ptlow[3]={0.,5.,10.};
  double high[3]={5,10,20};
  double a1=0, a2=0,a3=0,a4=0,a5=0,a6=0,a7=0;
  TFile * theFile = new TFile("dR.root","RECREATE");
  theFile->cd();
  for(int ihigh=0;ihigh<2;ihigh++)
    {
      for(int jhigh=0;jhigh<3;jhigh++)
       {
         sprintf(histname, "h_x_%d_%d",ihigh,jhigh); 
         cout<<flavor[ihigh][jhigh]<<endl;
        fraction[ihigh][jhigh] = new TH1D(histname,"",22,0.,1.1);

       }
    }  
   for(int ihigh=0;ihigh<2;ihigh++)
    {
      for(int jhigh=0;jhigh<3;jhigh++)
       {
         sprintf(histname, "h_x_%d_%d",ihigh,jhigh); 
        fraction1[ihigh][jhigh] = new TH1D(histname,"",22,0.,1.1);

       }
    }  
    for(int ihigh=0;ihigh<2;ihigh++)
    {
      for(int jhigh=0;jhigh<3;jhigh++)
       {
         sprintf(histname, "h_x_%d_%d",ihigh,jhigh); 
        fraction2[ihigh][jhigh] = new TH1D(histname,"",22,0.,1.1);

       }
    }  
      for(int ihigh=0;ihigh<2;ihigh++)
    {
      for(int jhigh=0;jhigh<3;jhigh++)
       {
         sprintf(histname, "h_x_%d_%d",ihigh,jhigh); 
        fraction3[ihigh][jhigh] = new TH1D(histname,"",22,0.,1.1);

       }
    }  
      for(int ihigh=0;ihigh<2;ihigh++)
    {
      for(int jhigh=0;jhigh<3;jhigh++)
       {
         sprintf(histname, "h_x_%d_%d",ihigh,jhigh); 
        fraction4[ihigh][jhigh] = new TH1D(histname,"",22,0.,1.1);

       }
    }  
      for(int ihigh=0;ihigh<2;ihigh++)
    {
      for(int jhigh=0;jhigh<3;jhigh++)
       {
         sprintf(histname, "h_x_%d_%d",ihigh,jhigh); 
        fraction5[ihigh][jhigh] = new TH1D(histname,"",22,0.,1.1);

       }
    }  
      for(int ihigh=0;ihigh<2;ihigh++)
    {
      for(int jhigh=0;jhigh<3;jhigh++)
       {
         sprintf(histname, "h_x_%d_%d",ihigh,jhigh); 
        fraction6[ihigh][jhigh] = new TH1D(histname,"",22,0.,1.1);

       }
    }  
  for(int ihigh=0;ihigh<3;ihigh++)
    {
      for(int jhigh=0;jhigh<3;jhigh++)
       {
         sprintf(histname, "h_x_%d_%d",ihigh,jhigh); 
        fraction7[ihigh][jhigh] = new TH1D(histname,"",22,0.,1.1);

       }
    }  
  //initialize every object
  string line;
  ifstream myfile;
  double pt=0,eta=0;
  //TLorentzVector struck = TLorentzVector(0,0,0,0);
  //TLorentzVector all_particle = TLorentzVector(0,0,0,0);
  myfile.open("jet_lepton.txt"); 
  int counter_all =0;
  while(getline(myfile, line)) 
  {
    if(line=="######") continue;
    else
    {
      int quark;
      double ch,nh,p;     double deltar, pt, eta,max_pt,jets_pt,delta,max_phi,max_eta; double pi_pt,kion_pt,proton_pt,kion_minus_pt,pi_minus_pt,anti_pt,neutral_pt,photon_pt,lepton_pt; int pi ,kion,proton;
     int max_id;
     //sscanf(line.c_str(),"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&quark,&pt,&kion_minus_pt,&pi_minus_pt,&deltar,&anti_pt,&eta,&pi_pt,&kion_pt,&proton_pt,&neutral_pt,&photon_pt,&lepton_pt);
     sscanf(line.c_str(),"%d %lf %d %lf %lf %lf %lf %lf",&quark,&deltar,&max_id,&max_pt,&jets_pt,&delta,&max_phi,&max_eta);//,&anti_pt,&eta,&pi_pt,&kion_pt,&proton_pt,&neutral_pt,&photon_pt,&lepton_pt);
    if(deltar<=threshold){
     if(quark ==1) a1++;
     else if(quark==2)a2++;
     else if(quark==3)a3++;
     else if(quark==-1)a4++;
     else if(quark==-2)a5++;
     else if(quark==-3)a6++;}
    //  double cmp_pt = pt+neutral_pt+photon_pt+lepton_pt;
    //  std::vector<double> my_vec = { anti_pt/cmp_pt,kion_minus_pt/cmp_pt,pi_minus_pt/cmp_pt,pi_pt/cmp_pt,kion_pt/cmp_pt,proton_pt/cmp_pt,photon_pt/cmp_pt,neutral_pt/cmp_pt};
     //pt here is charged hadron pt
     //std::vector<double> my_vec = { pt/cmp_pt,neutral_pt/cmp_pt,photon_pt/cmp_pt,lepton_pt/cmp_pt};
    //  int argmax = arg_max(my_vec);
    // TParticlePDG *pdg;
    //  pdg = pdgn->GetParticle(max_id);
    //  TString part = pdg->ParticleClass();
    //  double charge = pdg->Charge();
    //  if(part=="Meson"||part=="Baryon")
    //           {if(charge!=0) {max_id=0;}
    //           else if(charge==0)max_id=1;}
    //           if(part=="Lepton")max_id=2;
    //         if(max_id==22)max_id=3;
     switch (max_id){
       case -2212:
     if(deltar<=threshold)
      {counter_all++;
       for (int idex=0;idex<2;idex++)
        {for(int jdex=0;jdex<3;jdex++)
           {
             if(quark ==flavor[idex][jdex])
                if(max_pt/jets_pt>=0.)
                  fraction[idex][jdex]->Fill(delta);
           }
        }
      }
      break;

     case -321:
     if(deltar<=threshold)
      {counter_all++;
       for (int idex=0;idex<2;idex++)
        {for(int jdex=0;jdex<3;jdex++)
           {
             if(quark ==flavor[idex][jdex])
             if(max_pt/jets_pt>=0.)
                fraction1[idex][jdex]->Fill(delta);

           }
        }
      }
      break;
      case -211:
      if(deltar<=threshold)
      { counter_all++;
       for (int idex=0;idex<2;idex++)
        {for(int jdex=0;jdex<3;jdex++)
           {
             if(quark ==flavor[idex][jdex])
             if(max_pt/jets_pt>=0.)
                fraction2[idex][jdex]->Fill(delta);

           }
        }
      }
      break;
      case 211:
      if(deltar<=threshold)
      { counter_all++;
       for (int idex=0;idex<2;idex++)
        {for(int jdex=0;jdex<3;jdex++)
           {
             if(quark ==flavor[idex][jdex])
            if(max_pt/jets_pt>=0.)
               fraction3[idex][jdex]->Fill(delta);

           }
        }
      }
      break;
      case 321:
      if(deltar<=threshold)
      {counter_all++;
       for (int idex=0;idex<2;idex++)
        {for(int jdex=0;jdex<3;jdex++)
           {
             if(quark ==flavor[idex][jdex])
             if(max_pt/jets_pt>=0.)
                fraction4[idex][jdex]->Fill(delta);
       
                
           }
        }
      }
      break;
      case 2212:
      if(deltar<=threshold)
      { counter_all++;
       for (int idex=0;idex<2;idex++)
        {for(int jdex=0;jdex<3;jdex++)
           {
             if(quark ==flavor[idex][jdex])
              if(max_pt/jets_pt>=0.)
                fraction5[idex][jdex]->Fill(delta);
    
               
           }
        }
      }
      break;
      case 22:
      if(deltar<=threshold)
      { counter_all++;
       for (int idex=0;idex<2;idex++)
        {for(int jdex=0;jdex<3;jdex++)
           {
             if(quark ==flavor[idex][jdex])
              if(max_pt/jets_pt>=0.)
       fraction6[idex][jdex]->Fill(delta);

           }
        }
      }
      break;
// case 7:
//       if(deltar<=threshold)
//       { counter_all++;
//        for (int idex=0;idex<2;idex++)
//         {for(int jdex=0;jdex<3;jdex++)
//            {
//              if(quark ==flavor[idex][jdex])
//               if(my_vec[argmax]>=0.5)
//                { fraction7[idex][jdex]->Fill(neutral_pt/pt);
  
//                }
//            }
//         }
//       }
//       break;


    } 
    
  }
  }
vector<double> anti_proton  = vector<double>();
vector<double> kion_minus  = vector<double>();
vector<double> pi_minus  = vector<double>();
vector<double> pi_plus  = vector<double>();
vector<double> kion_plus  = vector<double>();
vector<double> proton  = vector<double>();
vector<double> photon  = vector<double>();
vector<double> neutral  = vector<double>();
int k=0;
for (int ieta=0;ieta<2;ieta++) 
  {
    for (int ipt=0;ipt<3;ipt++) 
    {
      k++;
      gStyle->SetPalette(1);
      //gPad->SetLogz();
      string ptcut ="<p_{t}^{sq}<";
      string etacut = "<#eta^{sq}<";
      stringstream ss;
      ss<<ptlow[ipt]<<"GeV"<<ptcut<<pthigh[ipt]<<"GeV"<<"  "<<etalow[ieta]<<etacut<<etahigh[ieta];
      fraction1[ieta][ipt]->SetTitle(ss.str().c_str());
      stringstream mean;
      //re[ieta][ipt]->GetHistogram()->SetMaximum(1.);   // along          
      //re[ieta][ipt]->GetHistogram()->SetMinimum(0.);  
      //fraction1[ieta][ipt]->GetXaxis()->SetTitle("#pi^{+}/#pi^{-} p_{t} fraction");//
      //fraction1[ieta][ipt]->GetXaxis()->SetTitle("K^{+}/K^{-} p_{t} fraction");
      //fraction1[ieta][ipt]->GetXaxis()->SetTitle("photon p_{t} fraction");
      fraction1[ieta][ipt]->GetXaxis()->SetTitle("number of photons");
      fraction1[ieta][ipt]->GetYaxis()->SetTitle("Number of Jets");//
      fraction[ieta][ipt]->SetLineColor(1);
      anti_proton.push_back(fraction[ieta][ipt]->GetMean());
      kion_minus.push_back(fraction1[ieta][ipt]->GetMean());
      pi_minus.push_back(fraction2[ieta][ipt]->GetMean());
      pi_plus.push_back(fraction3[ieta][ipt]->GetMean());
      kion_plus.push_back(fraction4[ieta][ipt]->GetMean());
      proton.push_back(fraction5[ieta][ipt]->GetMean());
      photon.push_back(fraction6[ieta][ipt]->GetMean());
      // neutral.push_back(fraction7[ieta][ipt]->GetEntries());
// fraction1[ieta][ipt]->SetMarkerSize(0.8);
// fraction1[ieta][ipt]->SetMarkerStyle(20);

//       fraction1[ieta][ipt]->Sumw2();
      
//       fraction1[ieta][ipt]->Scale(fraction[ieta][ipt]->Integral()/fraction1[ieta][ipt]->Integral());
//       fraction1[ieta][ipt]->Draw("E1");
//       fraction[ieta][ipt]->Draw("same");
      
//       TLatex  latex;
//       //latex.SetTextAlign(13);  //align at top
//       //mean<<"s_q entries = "<<fraction[ieta][ipt]->GetEntries();
//       mean<<fraction[ieta][ipt]->GetMean();
//       TPaveText *pave = new TPaveText(.85,.87,0.98,0.9,"brNDC");
//       pave->SetTextAlign(12);
//       char temp[20],temp1[20];
//       //sprintf(temp,"sq_entries= %5.0f",fraction[ieta][ipt]->GetEntries());
//       sprintf(temp,"Mean:%3.2f %3.2f" ,fraction[ieta][ipt]->GetMean(),fraction1[ieta][ipt]->GetMean());
//       sprintf(temp1,"Std:%3.2f %3.2f",fraction[ieta][ipt]->GetStdDev(),fraction1[ieta][ipt]->GetStdDev());
      

//       latex.DrawLatexNDC(0.51,0.85,Form("heavy flavor %.2f #pm %.2f",fraction1[ieta][ipt]->GetMean(),fraction1[ieta][ipt]->GetMeanError()));
//       latex.DrawLatexNDC(0.51,0.78,Form("light flavor %.2f #pm %.2f",fraction[ieta][ipt]->GetMean(),fraction[ieta][ipt]->GetMeanError()));
//       latex.Draw();
//       pave->AddText(temp);
//       TText *text = pave->AddText(temp1);
//       TLegend* leg = new TLegend(0.74,0.66,0.84,0.78);
//       leg->AddEntry(fraction[ieta][ipt],"light ","l");
      
//       leg->AddEntry(fraction1[ieta][ipt],"heavy","lep");
//       //leg->AddEntry(Form("heavy flavor %3.2f f#pm %3.2f",fraction1[ieta][ipt]->GetMean(),fraction1[ieta][ipt]->GetRMS()));
//       leg->Draw("same");
//       text->SetTextFont(72);
//       text->SetTextAlign(22);
      
//       //pave->Draw();
    }
  }
  double cnt[6] = {a1,a2,a3,a4,a5,a6};
  double q[6]={49612,321087,14751,24135,87199,15647};
  //double all =a1+a2+a3+a4+a5+a6+a7;
  TCanvas *c1 = new TCanvas("c1");
  c1->Divide(3,2);
  k=0;
  for (int ieta=0;ieta<2;ieta++) 
    {
      for (int ipt=0;ipt<3;ipt++) 
      {
        k++;
        c1->cd(k);
        gStyle->SetPalette(1);
        gStyle->SetHistMinimumZero(kTRUE);
        c1->cd(k)->SetGrid(); c1->cd(k)->SetBottomMargin(0.1);
 	c1->cd(k)->SetFillColor(10);  c1->cd(k)->SetBorderMode(0);
   fraction[ieta][ipt]->Sumw2();
     fraction1[ieta][ipt]->Sumw2();
     fraction2[ieta][ipt]->Sumw2();
     fraction3[ieta][ipt]->Sumw2();
     fraction4[ieta][ipt]->Sumw2();
     fraction5[ieta][ipt]->Sumw2();
     fraction6[ieta][ipt]->Sumw2();
     fraction[ieta][ipt]->Scale(fraction3[ieta][ipt]->Integral()/fraction[ieta][ipt]->Integral());
     fraction[ieta][ipt]->SetLineColor(1);
     fraction1[ieta][ipt]->Scale(fraction3[ieta][ipt]->Integral()/fraction1[ieta][ipt]->Integral());
     fraction1[ieta][ipt]->SetLineColor(2);
     fraction2[ieta][ipt]->Scale(fraction3[ieta][ipt]->Integral()/fraction2[ieta][ipt]->Integral());
     fraction2[ieta][ipt]->SetLineColor(28);
     fraction3[ieta][ipt]->SetLineColor(4);
     fraction3[ieta][ipt]->Scale(fraction3[ieta][ipt]->Integral()/fraction3[ieta][ipt]->Integral());
     fraction4[ieta][ipt]->SetLineColor(9);
     fraction4[ieta][ipt]->Scale(fraction3[ieta][ipt]->Integral()/fraction4[ieta][ipt]->Integral());
     fraction5[ieta][ipt]->SetLineColor(46);
     fraction5[ieta][ipt]->Scale(fraction3[ieta][ipt]->Integral()/fraction5[ieta][ipt]->Integral());
     fraction6[ieta][ipt]->SetLineColor(32);
     fraction6[ieta][ipt]->Scale(fraction3[ieta][ipt]->Integral()/fraction5[ieta][ipt]->Integral());
      fraction[ieta][ipt]->SetTitle(os_y[k-1].c_str());
    fraction[ieta][ipt]->Draw();
    fraction1[ieta][ipt]->Draw("same&&hist");
    fraction2[ieta][ipt]->Draw("same&&hist");
    fraction3[ieta][ipt]->Draw("same&&hist");
    fraction4[ieta][ipt]->Draw("same&&hist");
    fraction5[ieta][ipt]->Draw("same&&hist");
fraction6[ieta][ipt]->Draw("same&&hist");

  //  TH1F *h = new TH1F("h","",7,0,7);
  //  	h->SetFillColor(4);
	// h->SetBarWidth(0.4); h->SetBarOffset(0.1); h->SetStats(0);
	// h->SetMinimum(0.);
// h->SetMaximum(1.); 
//  cout<<cnt[k-1]<<endl;
// 		h->Fill(os_X[0].c_str(), anti_proton[k-1]/1.);
// 		h->GetXaxis()->SetBinLabel(1,os_X[0].c_str());
//     h->Fill(os_X[1].c_str(), kion_minus[k-1]/1.);
// 		h->GetXaxis()->SetBinLabel(2,os_X[1].c_str());
//     h->Fill(os_X[2].c_str(), pi_minus[k-1]/1.);
// 		h->GetXaxis()->SetBinLabel(3,os_X[2].c_str());  
//     h->Fill(os_X[3].c_str(), pi_plus[k-1]/1.);
// 	 	h->GetXaxis()->SetBinLabel(4,os_X[3].c_str());
//     h->Fill(os_X[4].c_str(), kion_plus[k-1]/1.);
// 		h->GetXaxis()->SetBinLabel(5,os_X[4].c_str());
//     h->Fill(os_X[5].c_str(), proton[k-1]/1.);
// 		h->GetXaxis()->SetBinLabel(6,os_X[5].c_str());
//     h->Fill(os_X[6].c_str(), photon[k-1]/1.);
// 		h->GetXaxis()->SetBinLabel(7,os_X[6].c_str());
//    // double eff = (anti_proton[k-1]+kion_minus[k-1]+pi_minus[k-1]+pi_plus[k-1]+kion_plus[k-1]+proton[k-1]+photon[k-1])/cnt[k-1];
//     //cout<<"eff:"<<eff<<endl;
//      h->SetTitle(os_y[k-1].c_str());

      TLatex latex;
      latex.SetTextAlign(13);  //align at top
      // mean<<"mean = "<<fraction[ieta][ipt]->GetEntries();
      TPaveText *pave = new TPaveText(.6,.75,.9,.9,"brNDC");
      // pave->SetTextAlign(12);
      // char temp[50];
      // sprintf(temp,"efficiency= %3.2f ",eff);
      // TText *text = pave->AddText(temp);
      // text->SetTextFont(72);
      // text->SetTextAlign(22);
    //   h->GetYaxis()->SetTitle("#DeltaR");
    //  pave->Draw();
    TLegend *leg = new TLegend(0.6,0.6,0.72,0.9);
	leg->SetTextSize(0.03); leg->SetBorderSize(2); leg->SetMargin(0.15);
	leg->SetFillColor(10);
	 leg->AddEntry(fraction6[ieta][ipt],"photon","l");
   leg->AddEntry(fraction[ieta][ipt],"p^{-}","lep");
   leg->AddEntry(fraction1[ieta][ipt],"K^{-}","l");
   leg->AddEntry(fraction2[ieta][ipt],"#pi^{-}","l");
   leg->AddEntry(fraction3[ieta][ipt],"#pi^{+}","l");
    leg->AddEntry(fraction4[ieta][ipt],"K^{+}","l");
    leg->AddEntry(fraction5[ieta][ipt],"p^{+}","l");
     leg->Draw();
      }
    }
    cout<<counter_all<<endl;
   
}
// void s(double threshold)
// {
  
  
//   TH1D *fraction[3][3];

//   char *histname = new char[10];
//   double etahigh[3]={-1,1,3.5};
//   double etalow[3] = {-3.5,-1,1};
//   double pthigh[3]={5.,10.,20.};
//   double ptlow[3]={1.,5.,10.};
//   for(int ihigh=0;ihigh<3;ihigh++)
//     {
//       for(int jhigh=0;jhigh<3;jhigh++)
//        {
//          sprintf(histname, "h_x_%d_%d",ihigh,jhigh); 
//         fraction[ihigh][jhigh] = new TH1D(histname,"",55,0.,1.1);

//        }
//     }  
  
//   //initialize every object
//   string line;
//   ifstream myfile;
//   double pt=0,eta=0;
//   TLorentzVector struck = TLorentzVector(0,0,0,0);
//   TLorentzVector all_particle = TLorentzVector(0,0,0,0);
//   myfile.open("jet.txt"); 
//   while(getline(myfile, line)) 
//   {
//     if(line=="######")
//     {
//       //calculate energy fraction and fill histogram accordingly

//       double den=all_particle.E();
//       double energy_arose_quark=struck.E()/den;
//       //double photon_ratio=remnant.E()/den;
//        for (int idex=0;idex<3;idex++)
//         {
//           for(int jdex=0;jdex<3;jdex++)
//            {
//              if((pt>=ptlow[jdex]&&pt<pthigh[jdex])&&(eta>=etalow[idex]&&eta<etahigh[idex]))
//                 fraction[idex][jdex]->Fill(energy_arose_quark);
//            }
//         }
      
//       // clear all
//        pt=0;
//        eta=0;
//        struck = TLorentzVector(0,0,0,0);
//        all_particle = TLorentzVector(0,0,0,0);
//     }
//     else
//     {
//      double w, x, y, z,q;
//      int m;
//      sscanf(line.c_str(),"%d %lf %lf %lf %lf %lf %lf %lf",&m,&w,&x,&y,&z,&q,&pt,&eta);
//      TLorentzVector vec = TLorentzVector(x,y,z,q);
//      if (w<=threshold) struck+=vec;  //threshold is 0.3
//      all_particle+=vec;
//      cout<<m<<" "<<w<<endl;
//     }   
//   }
// int k = 0;
// for (int ieta=0;ieta<3;ieta++) 
//   {
//     for (int ipt=0;ipt<3;ipt++) 
//     {
//      double mean = fraction[ieta][ipt]->GetMean();
//      result[ieta][ipt]=mean;
//     }
//   }


// }








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
return sqrt((part1->GetPhi()-part2->GetPhi())*(part1->GetPhi()-part2->GetPhi())+(part2->GetTheta()-part2->GetTheta())*(part1->GetTheta()-part2->GetTheta()));
}
