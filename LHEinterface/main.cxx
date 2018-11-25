#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "BucketofTops/BucketofTops.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "math.h"
#include "TLorentzVector.h"

using namespace std;

class BucketofTops; //forward declare

int main()
{
  //ifstream inFile("../tt_had_test_one.lhe");
  //ifstream inFile("../tt_had_test.lhe");
  //ifstream inFile("../tt_hadronic.lhe");
  //ifstream inFile("/afs/cern.ch/work/s/sosen/ChongbinTop/lhe/tt_hadronic.lhe");
  ifstream inFile("/afs/cern.ch/work/s/sosen/ChongbinTop/lhe/bbjjj.lhe");
  //ifstream inFile("../bbjjj_short.lhe");

  string line;
  bool event_flag = false; //switches on when finds an event
  bool event_meta = false; //event block readability switched off to skip the first event block line
  std::vector <TLorentzVector> specbjets;
  std::vector <TLorentzVector> specnonbjets;
  //mass
  TH1F htwmass("htwmass", "Mass of tw Buckets",150,0.0001,300); 
  TH1F htminmass("htminmass", "Mass of t- Buckets",150,0.0001,300); 
  TH1F ht0mass("ht0mass", "Mass of t0 Buckets",150,0.0001,300); 
  TH1F hXmass("hXmass", "Mass of the extra jets",110,-1,10); 
  // pT
  TH1F htwPt("htwPt", "Pt of tw Buckets",250,0,1200); 
  TH1F htminPt("htminPt", "Pt of t- Buckets",250,0,1200); 
  TH1F ht0Pt("ht0Pt", "Pt of t0 Buckets",250,0,1200); 
  TH1F hXPt("hXPt", "Pt of the extra jets",100,0,500); 
  // eta
  TH1F htweta("htweta", "#eta of tw Buckets",100,-10,10); 
  TH1F htmineta("htmineta", "#eta of t- Buckets",100,-10,10); 
  TH1F ht0eta("ht0eta", "#eta of t0 Buckets",100,-10,10);
  TH1F hXeta("hXeta", "#eta of the extra jets",100,-10,10);
  //W candidate mass
  TH1F hmW("hmW", "Mass of the (possible) W candidate",150,0.0001,300); 
  TH1F hmBucketPrim("hmBucketPrimitive", "Mass of the Entire Buckets before Recalculation",150,0,300); 
  TH1F hmratio("hmratio", "Mass Ratio Difference",120,-0.1,1.1); 
  

 
  int eventcounter = 0;
  int twcounter = 0;
  int tmincounter = 0;
  int t0counter = 0;
  int tXcounter = 0;
  while (getline(inFile, line)) 
  {
    if (line.find("<event>") != string::npos) 
    {
      event_flag = true;
      cout << "event: " << eventcounter << endl;
      eventcounter++;
      //cout << line << "\t" << event_flag << endl;
    }
    else if (line.find("</event>") != string::npos)
    {
      event_flag = false; //switch off the event block
      event_meta = false; //switch off the event block readability
      //discard events with less than two b jets
      if (specbjets.size() == 2)
      {
        BucketofTops *m_buckets = new BucketofTops(specbjets, specnonbjets);
        std::vector<bucketAlgo::bucket>& bucklist = *m_buckets->returnbucketlistptr();
	for(auto v: m_buckets->mWcand) {hmW.Fill(v);} 
	for(auto v: m_buckets->mBucketPrim) {hmBucketPrim.Fill(v);} 
	for(auto v: m_buckets->mratio) {hmratio.Fill(v);} 
	for(auto v: m_buckets->twmass) {htwmass.Fill(v);twcounter++;} 
	for(auto v: m_buckets->twPt) {htwPt.Fill(v);} 
	for(auto v: m_buckets->tweta) {htweta.Fill(v);} 
	for(auto v: m_buckets->tminmass) {htminmass.Fill(v);tmincounter++;} 
	for(auto v: m_buckets->tminPt) {htminPt.Fill(v);} 
	for(auto v: m_buckets->tmineta) {htmineta.Fill(v);} 
	for(auto v: m_buckets->t0mass) {ht0mass.Fill(v);t0counter++;} 
	for(auto v: m_buckets->t0Pt) {ht0Pt.Fill(v);} 
	for(auto v: m_buckets->t0eta) {ht0eta.Fill(v);} 
	for(auto v: m_buckets->Xmass) {hXmass.Fill(v);tXcounter++;} 
	for(auto v: m_buckets->XPt) {hXPt.Fill(v);} 
	for(auto v: m_buckets->Xeta) {hXeta.Fill(v);} 
      }
      specbjets.clear();
      specnonbjets.clear();
      //insert event operations before clearing the vector 
    }
    else{
      if (event_flag) 
      {
        if (event_meta)
        {
          istringstream iss(line);
	  std::vector<string> pinfo{istream_iterator<string>{iss},
            istream_iterator<string>{}};
	  std::string pxstr(pinfo[6]);
          double px = stof(pxstr);
	  std::string pystr(pinfo[7]);
          double py = stof(pystr);
	  std::string pzstr(pinfo[8]);
          double pz = stof(pzstr);
	  std::string Estr(pinfo[9]);
          double E = stof(Estr);
	  std::string pidstr(pinfo[0]);
          int pid = stoi(pidstr);
	  std::string statusstr(pinfo[1]);
          int status = stoi(statusstr);
	  
          if (status == 1) 
          {
            if (abs(pid) == 5)
	    {
	      TLorentzVector temb;
	      temb.SetPxPyPzE(px, py, pz, E);
	      specbjets.push_back(temb);
	    }
	    else
            {
	      TLorentzVector temnonb;
	      temnonb.SetPxPyPzE(px, py, pz, E);
	      specnonbjets.push_back(temnonb);
	    }

          }
        }
        else
        {
         event_meta = true; //to make the rest of the event block readable
        }
      }
    }
  }
  htwmass.Draw();
  htminmass.Draw();
  ht0mass.Draw();
  cout << "tw: " << twcounter << "\tt-: " << tmincounter << "\tt0: " << t0counter << "\ttX: "<< tXcounter << endl;
  TFile f("test.root", "RECREATE");
  TCanvas c ("c", "c", 800, 600);
  //mass
  htwmass.GetXaxis()->SetTitle("Mass (GeV)");
  htwmass.Write();
  htwmass.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_mass_tw.eps");
  c.Clear();
  htminmass.GetXaxis()->SetTitle("Mass (GeV)");
  htminmass.Write();
  htminmass.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_mass_t_.eps");
  c.Clear();
  ht0mass.GetXaxis()->SetTitle("Mass (GeV)");
  ht0mass.Write();
  ht0mass.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_mass_t0.eps");
  c.Clear();
  hXmass.GetXaxis()->SetTitle("Mass (GeV)");
  hXmass.Write();
  hXmass.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_mass_x.eps");
  c.Clear();

  //Pt
  htwPt.GetXaxis()->SetTitle("Pt (GeV)");
  htwPt.Write();
  htwPt.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_pt_tw.eps");
  c.Clear();
  htminPt.GetXaxis()->SetTitle("Pt (GeV)");
  htminPt.Write();
  htminPt.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_pt_t_.eps");
  c.Clear();
  ht0Pt.GetXaxis()->SetTitle("Pt (GeV)");
  ht0Pt.Write();
  ht0Pt.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_pt_t0.eps");
  c.Clear();
  hXPt.GetXaxis()->SetTitle("Pt (GeV)");
  hXPt.Write();
  hXPt.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_pt_x.eps");
  c.Clear();

  //Eta
  htweta.GetXaxis()->SetTitle("#eta");
  htweta.Write();
  htweta.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_eta_tw.eps");
  c.Clear();
  htmineta.GetXaxis()->SetTitle("#eta");
  htmineta.Write();
  htmineta.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_eta_t_.eps");
  c.Clear();
  ht0eta.GetXaxis()->SetTitle("#eta");
  ht0eta.Write();
  ht0eta.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_eta_t0.eps");
  c.Clear();
  hXeta.GetXaxis()->SetTitle("#eta");
  hXeta.Write();
  hXeta.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_eta_x.eps");
  c.Clear();

  //W candidate
  hmW.GetXaxis()->SetTitle("Mass (GeV)");
  hmW.Write();
  hmW.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_m_jk.eps");
  c.Clear();
  hmBucketPrim.GetXaxis()->SetTitle("Mass (GeV)");
  hmBucketPrim.Write();
  hmBucketPrim.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_m_b.eps");
  c.Clear();
  hmratio.GetXaxis()->SetTitle("Mass Ratio");
  hmratio.Write();
  hmratio.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_ratio.eps");
  c.Clear();

  f.Close();

  inFile.close();
  return 0;
  
}
