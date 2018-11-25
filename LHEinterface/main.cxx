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
  ifstream inFile("/afs/cern.ch/work/s/sosen/ChongbinTop/lhe/tt_hadronic.lhe");
  //ifstream inFile("../bbjjj_short.lhe");

  string line;
  bool event_flag = false; //switches on when finds an event
  bool event_meta = false; //event block readability switched off to skip the first event block line
  std::vector <TLorentzVector> specbjets;
  std::vector <TLorentzVector> specnonbjets;
  //vector <TLorentzVector> evt;
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
      //cout << line << "\t" << event_flag << endl;
      //finalstate::event ev1(evt);
      //cout << ">> " << ev1.nonbjet.size() << endl;
/*      for (int i = 0; i < evt.size(); ++i) {
        cout << evt[i].getPID() << "\t" << evt[i].getpX() << endl;
      }*/
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
      /*{
        //bucket algo PASS1 (to find tw buckets)//
        vector <bucketAlgo::bucket> B;
        double bucketP1massMax = 200; //GeV
        double bucketP1massMin = 155; //GeV
        double firstP1Bucketwt = 100;
        B = bucketAlgo::doublebucket(ev1, bucketP1massMax, bucketP1massMin, "tw", firstP1Bucketwt);
        //cout << "B1: " << B[0].getBucketLabel() << "\tB2: " << B[1].getBucketLabel() << endl;
        //pass2 find t- buckets//
        vector <bucketAlgo::bucket> tmincand; //fill tmin candidates
        int telseindex; //tw or t0 bucket
        for (int i = 0; i < B.size(); ++i)
        {
          if (B[i].getBucketMass() > -1) {hmW.Fill(B[i].WcandMnum());}
          if (B[i].getBucketMass() > -1) {hmBucketPrim.Fill(B[i].getBucketMass());}
          if (B[i].getBucketMass() > -1) {hmratio.Fill(B[i].WcandRatio());}
          if (B[i].getBucketLabel() == "t-") {tmincand.push_back(B[i]);}
          else 
          {
            telseindex = i;
          } 
        }
        //redo both buckets for t-
        double bucketP2massMax = 155; //GeV
        double bucketP2massMin = 75; //GeV
        double firstP2Bucketwt = 1;
        
        if (tmincand.size() == 2)
        {
          B = bucketAlgo::doublebucket(ev1, bucketP2massMax, bucketP2massMin, "t-", firstP2Bucketwt);
          //cout << "event: " << eventcounter << endl;
        }
        else if (tmincand.size() == 1)
        {
          B[1-telseindex] = bucketAlgo::singlebucket(ev1, B[telseindex], bucketP2massMax, bucketP2massMin);
        }
        for (int i = 0; i < B.size(); ++i)
        {
          if (B[i].getBucketLabel() == "tw") 
          {
              if (B[i].getBucketMass() > -1) {
        	  htwmass.Fill(B[i].getBucketMass());
        	  htwPt.Fill(B[i].getBucketPt());
        	  htweta.Fill(B[i].getBucketEta()); }
        	  ++twcounter;
          }
          else if (B[i].getBucketLabel() == "t-")
          {       
              if (B[i].getBucketMass() > -1) {
      	          htminmass.Fill(B[i].getBucketMass());
        	  htminPt.Fill(B[i].getBucketPt());
        	  htmineta.Fill(B[i].getBucketEta());}
                  ++tmincounter;
          }
          else 
          {
		  //cout << "LL: " << B[i].getBucketLabel() << "\tmass: " << B[i].getBucketMass() << endl;
        	  if (B[i].getBucketMass() > -1) {
		  ht0mass.Fill(B[i].getBucketMass());
        	  ht0Pt.Fill(B[i].getBucketPt());
        	  ht0eta.Fill(B[i].getBucketEta());}
		  if (B[i].getBucketEta() == 0) 
		  {
		          //vector<int> plll = B[i].getPIDlist();
		          vector<int> plll = B[i].getOrderlist();
	          }
                  ++t0counter;
          }
        }

          vector <finalstate::particle> Xtra = bucketAlgo::extra(ev1.EVT, B); // extra bucket
          //cout << "extrasize: "  << Xtra.size() << "\tBsize: " << B.size() << endl;
	  

          for (int nn=0; nn < Xtra.size(); ++nn)
          {
            ++tXcounter;
	    hXmass.Fill(Xtra[nn].getM());
            hXPt.Fill(Xtra[nn].getPt());
            hXeta.Fill(Xtra[nn].getEta());
          }
          
        //
      }*/
      specbjets.clear();
      specnonbjets.clear();
      //evt.clear(); //insert event operations before clearing the vector 
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

    
            //cout << p.getPID() << "\t status: " << p.getStatus() << "\t" << p.getpX() << endl;
            //evt.push_back(p);
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

/*{
    TLorentzVector b1, b2, nb1, nb2, nb3, nb4, nb5, nb6, nb7;
    b1.SetPxPyPzE(239.001, -44.1249, -150.152, 288.056);
    b2.SetPxPyPzE(158.135, 29.5485, -260.97, 307.221);
    std::vector<TLorentzVector> specbjets = {b1, b2};

    nb1.SetPxPyPzE(-61.5849 , 52.878     , -309.991   , 320.669);
    nb2.SetPxPyPzE(-67.7585 , 21.3895    , -93.6167   , 117.823);
    nb3.SetPxPyPzE(-60.8536 , -34.9778   , -63.2891   , 95.2948);
    nb4.SetPxPyPzE(46.7464  , -51.02     , 38.2181    , 79.4786);
    nb5.SetPxPyPzE(-41.4773 , -47.7591   , -132.374   , 147.095);
    nb6.SetPxPyPzE(42.4618  , -25.8559   , -51.3354   , 71.7133);
    nb7.SetPxPyPzE(10.5252  , 33.8559    , -44.5788   , 57.5742);
    std::vector<TLorentzVector> specnonbjets = {nb1, nb2, nb3, nb4, nb5, nb6, nb7};

    BucketofTops *m_buckets = new BucketofTops(specbjets, specnonbjets);
    std::vector<bucketAlgo::bucket> bucklist = m_buckets->Blist;
    std::cout << "buclet list size: " << bucklist.size() << std::endl;
    std::cout << "init bucket mass: " << bucklist[0].getBucketMass() << "\t" << bucklist[1].getBucketMass() << std::endl;
    std::cout << "init bucket label: " << bucklist[0].getBucketLabel() << "\t" << bucklist[1].getBucketLabel() << std::endl;
    std::cout << "init bucket pT: " << bucklist[0].getBucketPt() << "\t" << bucklist[1].getBucketPt() << std::endl;
    std::cout << "init bucket eta: " << bucklist[0].getBucketEta() << "\t" << bucklist[1].getBucketEta() << std::endl;
    std::cout << "init bucket WcandMnum(): " << bucklist[0].WcandMnum() << "\t" << bucklist[1].WcandMnum() << std::endl;
  return 0;
}*/


