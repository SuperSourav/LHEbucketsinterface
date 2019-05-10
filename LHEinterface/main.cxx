#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
//#include <random>
#include "BucketofTops/BucketofTops.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "math.h"
#include "TLorentzVector.h"
#include "TLine.h"
#include "TRandom.h"
#include "TLatex.h"
#include <TPaveStats.h>
#include "TStyle.h"
#include <map>
#include "TLegend.h"

using namespace std;

class BucketofTops; //forward declare
class topquark; //forward declare
TRandom *r = new TRandom();

TLorentzVector smear(TLorentzVector v)
{
  float E, eta, theta, phi, M, pt;
  //default_random_engine generator;
  //normal_distribution<double> distribution(E,sqrt(E));
  //float Enew = distribution(generator);
  E = v.E();
  //float Enew = r->Gaus(E, sqrt(E));
  float rand = r->Gaus(0, 1);
  float Enew = E + sqrt(E)*(rand);
  E = (Enew > 0) ? Enew : 0;
  M = v.M();
  eta = v.Eta();
  theta = 2*atan(exp(-1*eta));
  pt = sin(theta)*sqrt((E*E) - (M*M));
  phi = v.Phi();
  TLorentzVector vsmr;
  vsmr.SetPtEtaPhiE(pt, eta, phi, E);
  //if (v.Eta() >0) {cout << "+ve qudrant: " << theta << " Gaus " << rand << endl;}
  //else {cout << "-ve quadrant: "<< theta << " Gaus " << rand << endl;}
  //cout << "ori kin>>>" << " eta: " << v.Eta()    << " E: " << v.E()    << " pt: " << v.Pt()    << " pz: " << v.Pz()    << " M: " << v.M()    << endl;
  //cout << "smear kin>" << " eta: " << vsmr.Eta() << " E: " << vsmr.E() << " pt: " << vsmr.Pt() << " pz: " << vsmr.Pz() << " M: " << vsmr.M() << endl;
  return vsmr;
}


TLorentzVector LHEmomentumParser(string line){
  istringstream iss(line);
  std::vector<string> pinfo{istream_iterator<string>{iss},
  istream_iterator<string>{}};
  //cout << "line: " << line << endl;
  std::string pxstr(pinfo[6]);
  double px = stof(pxstr);
  std::string pystr(pinfo[7]);
  double py = stof(pystr);
  std::string pzstr(pinfo[8]);
  double pz = stof(pzstr);
  std::string Estr(pinfo[9]);
  double E = stof(Estr);
  TLorentzVector v;
  v.SetPxPyPzE(px, py, pz, E);
  return v;
}


class topquark
{
  public:
    int top_index, W_index, bottom_index;
    TLorentzVector top, W, bottom;
    int Wdkq1_index, Wdkq2_index;
    TLorentzVector Wdkq1, Wdkq2;


  ~topquark() {}

  topquark()
  {
    top_index     =  -9999;
    W_index       =  -9999;
    bottom_index  =  -9999;
    Wdkq1_index   =  -9999;
    Wdkq2_index   =  -9999;
  }

  topquark(int tindex, 
	   TLorentzVector vTop, 
	   int Windex, 
	   TLorentzVector vW, 
	   int bindex, 
	   TLorentzVector vBottom, 
	   int q1index, 
	   TLorentzVector vq1, 
	   int q2index, 
	   TLorentzVector vq2)
  {
    top_index     =  tindex;
    W_index       =  Windex;
    bottom_index  =  bindex;
    Wdkq1_index   =  q1index;
    Wdkq2_index   =  q2index;
    top           =  vTop;
    W             =  vW;
    bottom        =  vBottom;
    Wdkq1         =  vq1;
    Wdkq2         =  vq2;
  }

  void settopindex(int tindex){
    top_index = tindex;
  }

  void settopP(TLorentzVector v){
    top = v;
  }

  void setWindex(int Windex){
    W_index = Windex;
  }

  void setWP(TLorentzVector v){
    W = v;
  }

  void setbottomindex(int bindex){
    bottom_index = bindex;
  }

  void setbottomP(TLorentzVector v){
    bottom = v;
  }

  void setWdkProducts(map<int,int>parentWdict, map<int, string> eventDict){
    vector<int> wdk;
    for(auto w: parentWdict){
      if (w.second == W_index) {wdk.push_back(w.first);}
    }
    Wdkq1_index = wdk[0];
    Wdkq1 = LHEmomentumParser(eventDict[Wdkq1_index]);
    Wdkq2_index = wdk[1];
    Wdkq2 = LHEmomentumParser(eventDict[Wdkq2_index]);
  }

  void printMembers(){
    cout << "top: \t" << top_index << endl;
    cout << "top: \t px: " << top.Px() << " \t py: " << top.Py() << " \t pz: " << top.Pz() << " \t E :" << top.E() << endl;
    cout << "W: \t" << W_index << "--> " << Wdkq1_index << "+" << Wdkq2_index << endl;
    cout << "W: \t px: " << W.Px() << " \t py: " << W.Py() << " \t pz: " << W.Pz() << " \t E :" << W.E() << endl;
    cout << "Wdkq1: \t px: " << Wdkq1.Px() << " \t py: " << Wdkq1.Py() << " \t pz: " << Wdkq1.Pz() << " \t E :" << Wdkq1.E() << endl;
    cout << "Wdkq2: \t px: " << Wdkq2.Px() << " \t py: " << Wdkq2.Py() << " \t pz: " << Wdkq2.Pz() << " \t E :" << Wdkq2.E() << endl;
    cout << "bottom: \t" << bottom_index << endl;
    cout << "bottom: \t px: " << bottom.Px() << " \t py: " << bottom.Py() << " \t pz: " << bottom.Pz() << " \t E :" << bottom.E() << endl;
  }

  double validateMass(){
    return (top.M() - (bottom + Wdkq1 + Wdkq2).M());
  }
};



map<int, topquark> topParser(vector<string> lines){
  map<int, string> eventDict;
  map<int, topquark> topDict;
  //cout << "inside the topParser" << endl;
  int partindex = 1;

  for (auto line: lines){
    //cout << line << endl;
    eventDict[partindex] = line;
    istringstream iss(line);
    std::vector<string> pinfo{istream_iterator<string>{iss},
      istream_iterator<string>{}};
    std::string pidstr(pinfo[0]);
    int pid = stoi(pidstr);
    std::string statusstr(pinfo[1]);
    if (abs(pid) == 6) {
      topquark t;
      t.settopindex(partindex);
      t.settopP(LHEmomentumParser(eventDict[partindex]));
      topDict[partindex] = t;
    }
    partindex++;
  }

  vector<int> W_index;
  vector<int> b_index;
  map <int, int> parentWdict;
  
  for (auto p: eventDict){
    string l = p.second;
    istringstream iss(l);
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
    std::string parent1str(pinfo[2]);
    int parent1index = stoi(parent1str);
    std::string parent2str(pinfo[3]);
    int parent2index = stoi(parent2str);
    if (status == 1){
      if (abs(pid) == 5) {
        b_index.push_back(parent1index);
        topDict[parent1index].setbottomindex(p.first);
        topDict[parent1index].setbottomP(LHEmomentumParser(eventDict[p.first]));
        //cout << p.first << "-b- " << eventDict[parent1index] << endl;
      }
      else {
         parentWdict[p.first] = parent1index;
      }
    }
    else {
      if (abs(pid) == 24) {
        W_index.push_back(parent1index);
        topDict[parent1index].setWindex(p.first);
        topDict[parent1index].setWP(LHEmomentumParser(eventDict[p.first]));
        //cout << p.first << "-W- " << eventDict[parent1index] << endl;
      }
    }
  }

  for (auto& t: topDict){
    t.second.setWdkProducts(parentWdict, eventDict);
  }
  
  //validation
  /*for (auto t: topDict){
    cout << "top:\t" << eventDict[t.second.top_index] << endl;
    t.second.printMembers();
    cout << "Validation:   " << t.second.validateMass() << "\t" << t.second.top.M() << endl;
  }*/
  return topDict;
}


int particleIndexfinder(vector<string> lines, TLorentzVector p){
  int partindex = 1;
  float dRmin = pow(10,10); //arbit large number
  int matchindex = -9999;
  for (auto line: lines){
    if (p.DeltaR(LHEmomentumParser(line)) < dRmin){matchindex = partindex; dRmin = p.DeltaR(LHEmomentumParser(line));}
    partindex++;
  }
  //cout << "test matchindex: " << matchindex << "\tpx: " << p.Px() << "\tpy: " << p.Py() << "\tpz: " << p.Pz() << "\tE: " << p.E() << endl;
  //cout << "test line: dR " << dRmin << lines.at(matchindex-1) << endl;
  return matchindex;
}

int bucketmatcher(vector<string> lines, TLorentzVector b, vector <TLorentzVector> nonbs, map<int, topquark>& topDict){
  //1 -> contamination; //2 -> subset; // 0 -> total match 
  int match = 1; //
  if (nonbs.size() == 2) {
    int bindex = particleIndexfinder(lines, b);
    int nb1index = particleIndexfinder(lines, nonbs[0]);
    int nb2index = particleIndexfinder(lines, nonbs[1]);
    //cout << ">>>>bindex: " << bindex << "\tnb1index: " << nb1index << "\tnb2index: " << nb2index << endl;
    //for (auto t: topDict){
    //cout << "B>>>bindex: " << t.second.bottom_index << "\tnb1index: " << t.second.Wdkq1_index << "\tnb2index: " << t.second.Wdkq2_index << endl;
    //}

    for (auto t: topDict){
      if (bindex == t.second.bottom_index){
        if ((nb1index == t.second.Wdkq1_index) && (nb2index == t.second.Wdkq2_index)){ match = 0; }
        else if ((nb1index == t.second.Wdkq2_index) && (nb2index == t.second.Wdkq1_index)){ match = 0; }
      }
    }
  }
  else {
    int bindex = particleIndexfinder(lines, b);
    int nb1index = particleIndexfinder(lines, nonbs[0]);
    for (auto t: topDict){
      if (bindex == t.second.bottom_index){
        if ((nb1index == t.second.Wdkq1_index) || (nb1index == t.second.Wdkq2_index)){ match = 2; }
      }
    }
  }
  return match;
}



int main()
{
  gStyle->SetOptStat("irmen");
		  //001001111);
  //ifstream inFile("../tt_had_test_one.lhe");
  //ifstream inFile("../tt_had_test.lhe");
  //ifstream inFile("../tt_hadronic.lhe");
  ifstream inFile("/afs/cern.ch/work/s/sosen/ChongbinTop/lhe/tt_hadronic.lhe");
  //ifstream inFile("/afs/cern.ch/work/s/sosen/ChongbinTop/lhe/bbjjj.lhe");
  //ifstream inFile("../bbjjj_short.lhe");
  
  string line;
  std::vector <string> lines;
  int Wcount = 0;
  bool event_flag = false; //switches on when finds an event
  bool event_meta = false; //event block readability switched off to skip the first event block line
  std::vector <TLorentzVector> specbjets;
  std::vector <TLorentzVector> specnonbjets;
  //mass
  TH1F htwt0mass("htwt0mass", "Mass of tw and t0 Buckets superposed",150,0.0001,300); 
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
  //TH1F hmBucketPrim("hmBucketPrimitive", "Mass of the Entire Buckets before Recalculation",150,0,300); 
  //TH1F hmBucketPrim0("hmBucketPrimitiveB1", "Mass of the B1 Bucket before Recalculation",150,0,300); 
  //TH1F hmBucketPrim1("hmBucketPrimitiveB2", "Mass of the B2 Bucket before Recalculation",150,0,300); 
  TH1F hmBucketPrim("hmBucketPrimitive", "Mass of the Entire Buckets before Recalculation",75,100,250); 
  TH1F hmBucketPrim0("hmBucketPrimitiveB1", "Mass of the B1 Bucket before Recalculation",75,100,250); 
  TH1F hmBucketPrim1("hmBucketPrimitiveB2", "Mass of the B2 Bucket before Recalculation",75,100,250); 
  TH1F hmratio("hmratio", "Mass Ratio Difference",120,-0.1,1.1); 
  
  //njets
  TH1F hnonbjetinit("hnonbjetinit", "non b jets before buckets", 16, -0.5, 15.5);
  TH1F hnonbjetB1("hnonbjetB1", "non b jets in B1", 16, -0.5, 15.5);
  TH1F hnonbjetB2("hnonbjetB2", "non b jets in B2", 16, -0.5, 15.5);
  TH1F hnonbjetBISR("hnonbjetBISR", "non b jets in BISR", 16, -0.5, 15.5);

  //truthtops
  TH1F hMVal("hMVal", "mass difference between top and its dk product", 200, -10, 10);

  //truthmatched buckets
  //
  TH1F hmBucketPrim0C("hmBucketPrimitiveB1C", "Mass of the B1 Bucket before Recalculation",75,100,250); 
  TH1F hmBucketPrim1C("hmBucketPrimitiveB2C", "Mass of the B2 Bucket before Recalculation",75,100,250); 

  //njets
  TH1F hnonbjetB1C("hnonbjetB1C", "non b jets in B1", 16, -0.5, 15.5);
  TH1F hnonbjetB2C("hnonbjetB2C", "non b jets in B2", 16, -0.5, 15.5);

  //wrong buckets
  //
  TH1F hmBucketPrim0W1("hmBucketPrimitiveB1W1", "Mass of the B1 Bucket before Recalculation",75,100,250); 
  TH1F hmBucketPrim0W2("hmBucketPrimitiveB1W2", "Mass of the B1 Bucket before Recalculation",75,100,250); 
  TH1F hmBucketPrim1W1("hmBucketPrimitiveB2W1", "Mass of the B2 Bucket before Recalculation",75,100,250); 
  TH1F hmBucketPrim1W2("hmBucketPrimitiveB2W2", "Mass of the B2 Bucket before Recalculation",75,100,250); 

  //njets
  TH1F hnonbjetB1W1("hnonbjetB1W1", "non b jets in B1", 16, -0.5, 15.5);
  TH1F hnonbjetB1W2("hnonbjetB1W2", "non b jets in B1", 16, -0.5, 15.5);
  TH1F hnonbjetB2W1("hnonbjetB2W1", "non b jets in B2", 16, -0.5, 15.5);
  TH1F hnonbjetB2W2("hnonbjetB2W2", "non b jets in B2", 16, -0.5, 15.5);


//  
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
      //if (eventcounter == 1000) {break;}
      cout << "event: " << eventcounter << endl;
      eventcounter++; 
      //cout << line << "\t" << event_flag << endl;
    }
    else if (line.find("</event>") != string::npos)
    {
      event_flag = false; //switch off the event block
      event_meta = false; //switch off the event block readability
      //discard events with less than two b jets
      if ((specbjets.size() == 2) && (Wcount == 2))
      {
        map<int, topquark> topDict = topParser(lines);
        for (auto t: topDict){
          hMVal.Fill(t.second.validateMass());
        }
        BucketofTops *m_buckets = new BucketofTops(specbjets, specnonbjets);
        std::vector<bucketAlgo::bucket>& bucklist = *m_buckets->returnbucketlistptr();
	//prebuckets are buckets before labels are added to them
	TLorentzVector prebucket1bjet = m_buckets->B1bjet;
	vector<TLorentzVector> prebucket1nonbjets = m_buckets->B1nonbjets;
        int B1truthmatchFlag = bucketmatcher(lines, prebucket1bjet, prebucket1nonbjets, topDict);
	TLorentzVector prebucket2bjet = m_buckets->B2bjet;
	vector<TLorentzVector> prebucket2nonbjets = m_buckets->B2nonbjets; 
        int B2truthmatchFlag = bucketmatcher(lines, prebucket2bjet, prebucket2nonbjets, topDict);
        //for (auto t: topDict){
    //cout << "B>bindex: " << t.second.bottom_index << "\tnb1index: " << t.second.Wdkq1_index << "\tnb2index: " << t.second.Wdkq2_index << endl;}
	//cout << "B1flag: " << B1truthmatchFlag << "\tB2flag: " << B2truthmatchFlag << endl;
	//
        for(auto v: m_buckets->mWcand) {hmW.Fill(v);} 
        for(auto v: m_buckets->mBucketPrim) {hmBucketPrim.Fill(v);}
        std::vector<float> m_mBucketPrim = m_buckets->mBucketPrim;
        hmBucketPrim0.Fill(m_mBucketPrim.at(0));
        if (B1truthmatchFlag == 0) {hmBucketPrim0C.Fill(m_mBucketPrim.at(0));} //truthmatched
	else if (B1truthmatchFlag == 1) {hmBucketPrim0W1.Fill(m_mBucketPrim.at(0));} //subset
	else {hmBucketPrim0W2.Fill(m_mBucketPrim.at(0));} //wrong combination
	/*cout << "B1 mass: " << m_mBucketPrim.at(0) <<endl;
        for(auto v: bucklist[0].nonBJETS) {
	  cout << "B1 nonb-jet:" << "\tpx: " << v.Px() << "\tpy: " << v.Py() << "\tpz: " << v.Pz() << "\tE: " << v.E() << endl;
	}
        TLorentzVector bjetB1 = bucklist[0].BJET;
	cout << "B1 b-jet:" << "\tpx: " << bjetB1.Px() << "\tpy: " << bjetB1.Py() << "\tpz: " << bjetB1.Pz() << "\tE: " << bjetB1.E() << endl;
        */

        hmBucketPrim1.Fill(m_mBucketPrim.at(1)); 
        if (B2truthmatchFlag == 0) {hmBucketPrim1C.Fill(m_mBucketPrim.at(1));} //truthmatched
	else if (B2truthmatchFlag == 1) {hmBucketPrim1W1.Fill(m_mBucketPrim.at(1));} //subset
	else {hmBucketPrim1W2.Fill(m_mBucketPrim.at(1));} //wrong combination
	/*cout << "B2 mass: " << m_mBucketPrim.at(1) <<endl;
        for(auto v: bucklist[1].nonBJETS) {
	  cout << "B1 nonb-jet:" << "\tpx: " << v.Px() << "\tpy: " << v.Py() << "\tpz: " << v.Pz() << "\tE: " << v.E() << endl;
	}
        TLorentzVector bjetB2 = bucklist[1].BJET;
	cout << "B2 b-jet:" << "\tpx: " << bjetB2.Px() << "\tpy: " << bjetB2.Py() << "\tpz: " << bjetB2.Pz() << "\tE: " << bjetB2.E() << endl;
	*/

	for(auto v: m_buckets->mratio) {hmratio.Fill(v);} 
	for(auto v: m_buckets->twmass) {htwmass.Fill(v);twcounter++;htwt0mass.Fill(v);} 
	for(auto v: m_buckets->twPt) {htwPt.Fill(v);} 
	for(auto v: m_buckets->tweta) {htweta.Fill(v);} 
	for(auto v: m_buckets->tminmass) {htminmass.Fill(v);tmincounter++;} 
	for(auto v: m_buckets->tminPt) {htminPt.Fill(v);} 
	for(auto v: m_buckets->tmineta) {htmineta.Fill(v);} 
	for(auto v: m_buckets->t0mass) {ht0mass.Fill(v);t0counter++;htwt0mass.Fill(v);} 
	for(auto v: m_buckets->t0Pt) {ht0Pt.Fill(v);} 
	for(auto v: m_buckets->t0eta) {ht0eta.Fill(v);} 
	for(auto v: m_buckets->Xmass) {hXmass.Fill(v);tXcounter++;} 
	for(auto v: m_buckets->XPt) {hXPt.Fill(v);} 
	for(auto v: m_buckets->Xeta) {hXeta.Fill(v);} 
	hnonbjetinit.Fill(m_buckets->nonbinitcount); 
	hnonbjetB1.Fill(m_buckets->nonb1count); 
        if (B1truthmatchFlag == 0) {hnonbjetB1C.Fill(m_buckets->nonb1count);} //truthmatched
	else if (B1truthmatchFlag == 1) {hnonbjetB1W1.Fill(m_buckets->nonb1count);} //subset
	else {hnonbjetB1W2.Fill(m_buckets->nonb1count);} //contamination
	hnonbjetB2.Fill(m_buckets->nonb2count); 
        if (B2truthmatchFlag == 0) {hnonbjetB2C.Fill(m_buckets->nonb2count);} //truthmatched
	else if (B2truthmatchFlag == 1) {hnonbjetB2W1.Fill(m_buckets->nonb2count);} //subset
	else {hnonbjetB2W2.Fill(m_buckets->nonb2count);} //contamination
	hnonbjetBISR.Fill(m_buckets->nonbISRcount); 
      }
      specbjets.clear();
      specnonbjets.clear();
      lines.clear();
      Wcount = 0;
      //insert event operations before clearing the vector 
    }
    else{
      if (event_flag) 
      {
        if (event_meta)
        {
	  lines.push_back(line);
          istringstream iss(line);
	  std::vector<string> pinfo{istream_iterator<string>{iss},
            istream_iterator<string>{}};
	  //cout << "line: " << line << endl;
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
	      //specbjets.push_back(smear(temb));
	    }
	    else
            {
	      TLorentzVector temnonb;
	      temnonb.SetPxPyPzE(px, py, pz, E);
	      specnonbjets.push_back(temnonb);
	      //specnonbjets.push_back(smear(temnonb));
	    }
          }
	  else if (abs(pid) == 24){Wcount++;}
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
  htwt0mass.Draw();
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
  hMVal.GetXaxis()->SetTitle("Mass (GeV)");
  hMVal.Write();
  hMVal.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_massdiff_truth_.eps");
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
  htwt0mass.GetXaxis()->SetTitle("Mass (GeV)");
  htwt0mass.Write();
  htwt0mass.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_mass_twt0.eps");
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
  /*// Retrieve the stat box
  TPaveStats *psmBucketPrim = (TPaveStats*)c.GetPrimitive("stats");
  psmBucketPrim->SetName("mystats");
  TList *listOfLinesmBucketPrim = psmBucketPrim->GetListOfLines();
  // Add a new line in the stat box.
  // Note that "=" is a control character
  int bminmBucketPrim = hmBucketPrim.GetXaxis()->FindBin(100);
  int bmaxmBucketPrim = hmBucketPrim.GetXaxis()->FindBin(250);
  TLatex *mytmBucketPrim = new TLatex(0,0,Form("Integral: (%.f)",hmBucketPrim.Integral(bminmBucketPrim, bmaxmBucketPrim)));
  mytmBucketPrim->SetTextFont(42);
  mytmBucketPrim->SetTextSize(0.04);
  listOfLinesmBucketPrim->Add(mytmBucketPrim);
  // the following line is needed to avoid that the automatic redrawing of stats
  hmBucketPrim.SetStats(0);
  c.Modified();*/
  c.Print("cpp_reconstructed_top_m_b.eps");
  c.Clear();
  hmBucketPrim0.GetXaxis()->SetTitle("Mass (GeV)");
  hmBucketPrim0.Write();
  hmBucketPrim0.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_m_b1.eps");
  c.Clear();
  auto legM0 = TLegend( 0.65, 0.75, 0.88, 0.88);
  legM0.SetFillColor(0);
  legM0.SetLineColor(0);
  hmBucketPrim0W1.GetXaxis()->SetTitle("Mass (GeV)");
  hmBucketPrim0W1.SetLineColor(2);
  hmBucketPrim0W1.Write();
  hmBucketPrim0W1.SetMaximum(max(max(hmBucketPrim0W1.GetMaximum(), hmBucketPrim0W2.GetMaximum()), hmBucketPrim0C.GetMaximum())*1.1);
  hmBucketPrim0W1.Draw();
  legM0.AddEntry(&hmBucketPrim0W1, Form("subset (%.f)", hmBucketPrim0W1.Integral()), "l");
  hmBucketPrim0W2.GetXaxis()->SetTitle("Mass (GeV)");
  hmBucketPrim0W2.SetLineColor(3);
  hmBucketPrim0W2.Write();
  hmBucketPrim0W2.Draw("same");
  legM0.AddEntry(&hmBucketPrim0W2, Form("contamination (%.f)", hmBucketPrim0W2.Integral()), "l");
  hmBucketPrim0C.GetXaxis()->SetTitle("Mass (GeV)");
  hmBucketPrim0C.SetLineColor(1);
  hmBucketPrim0C.Write();
  hmBucketPrim0C.Draw("same");
  legM0.AddEntry(&hmBucketPrim0C, Form("correct (%.f)", hmBucketPrim0C.Integral()), "l");
  legM0.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_m_b1_truthmatched.eps");
  c.Clear();
  hmBucketPrim1.GetXaxis()->SetTitle("Mass (GeV)");
  hmBucketPrim1.Write();
  hmBucketPrim1.Draw();
  TLine *marker = new TLine(173,0,173,50000);
  marker->SetLineColor(kRed);
  marker->Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_m_b2.eps");
  c.Clear();
  auto legM1 = TLegend( 0.65, 0.75, 0.88, 0.88);
  legM1.SetFillColor(0);
  legM1.SetLineColor(0);
  hmBucketPrim1W1.GetXaxis()->SetTitle("Mass (GeV)");
  hmBucketPrim1W1.SetLineColor(2);
  hmBucketPrim1W1.Write();
  hmBucketPrim1W1.SetMaximum(max(max(hmBucketPrim1W1.GetMaximum(), hmBucketPrim1W2.GetMaximum()), hmBucketPrim1C.GetMaximum())*1.1);
  hmBucketPrim1W1.Draw();
  legM1.AddEntry(&hmBucketPrim1W1, Form("subset (%.f)", hmBucketPrim1W1.Integral()), "l");
  hmBucketPrim1W2.GetXaxis()->SetTitle("Mass (GeV)");
  hmBucketPrim1W2.SetLineColor(3);
  hmBucketPrim1W2.Write();
  hmBucketPrim1W2.Draw("same");
  legM1.AddEntry(&hmBucketPrim1W2, Form("contamination (%.f)", hmBucketPrim1W2.Integral()), "l");
  hmBucketPrim1C.GetXaxis()->SetTitle("Mass (GeV)");
  hmBucketPrim1C.SetLineColor(1);
  hmBucketPrim1C.Write();
  hmBucketPrim1C.Draw("same");
  legM1.AddEntry(&hmBucketPrim1C, Form("correct (%.f)", hmBucketPrim1C.Integral()), "l");
  legM1.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_m_b2_truthmatched.eps");
  c.Clear();
  hmratio.GetXaxis()->SetTitle("Mass Ratio");
  hmratio.Write();
  hmratio.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_ratio.eps");
  c.Clear();
  hnonbjetinit.GetXaxis()->SetTitle("nonb jets per event before buckets");
  hnonbjetinit.Write();
  hnonbjetinit.Draw();
  c.Update();
  c.Print("cpp_reconstructed_nonbjetcountinit.eps");
  c.Clear();
  hnonbjetB1.GetXaxis()->SetTitle("nonb jets per event in B1");
  hnonbjetB1.Write();
  hnonbjetB1.Draw();
  c.Update();
  c.Print("cpp_reconstructed_nonbjetcountB1.eps");
  c.Clear();
  auto legnb1 = TLegend( 0.65, 0.75, 0.88, 0.88);
  legnb1.SetFillColor(0);
  legnb1.SetLineColor(0);
  hnonbjetB1W1.GetXaxis()->SetTitle("nonb jets per event in B1");
  hnonbjetB1W1.SetLineColor(2);
  hnonbjetB1W1.Write();
  hnonbjetB1W1.SetMaximum(max(max(hnonbjetB1W1.GetMaximum(), hnonbjetB1W2.GetMaximum()), hnonbjetB1C.GetMaximum())*1.1);
  hnonbjetB1W1.Draw();
  legnb1.AddEntry(&hnonbjetB1W1, Form("subset (%.f)", hnonbjetB1W1.Integral()), "l");
  hnonbjetB1W2.GetXaxis()->SetTitle("nonb jets per event in B1");
  hnonbjetB1W2.SetLineColor(3);
  hnonbjetB1W2.Write();
  hnonbjetB1W2.Draw("same");
  legnb1.AddEntry(&hnonbjetB1W2, Form("contamination (%.f)", hnonbjetB1W2.Integral()), "l");
  hnonbjetB1C.GetXaxis()->SetTitle("nonb jets per event in B1");
  hnonbjetB1C.SetLineColor(1);
  hnonbjetB1C.Write();
  hnonbjetB1C.Draw("same");
  legnb1.AddEntry(&hnonbjetB1C, Form("correct (%.f)", hnonbjetB1C.Integral()), "l");
  legnb1.Draw();
  c.Update();
  c.Print("cpp_reconstructed_nonbjetcountB1truthmatched.eps");
  c.Clear();
  hnonbjetB2.GetXaxis()->SetTitle("nonb jets per event in B2");
  hnonbjetB2.Write();
  hnonbjetB2.Draw();
  c.Update();
  c.Print("cpp_reconstructed_nonbjetcountB2.eps");
  c.Clear();
  auto legnb2 = TLegend( 0.65, 0.75, 0.88, 0.88);
  legnb2.SetFillColor(0);
  legnb2.SetLineColor(0);
  hnonbjetB2W1.GetXaxis()->SetTitle("nonb jets per event in B2");
  hnonbjetB2W1.SetLineColor(2);
  hnonbjetB2W1.Write();
  hnonbjetB2W1.SetMaximum(max(max(hnonbjetB2W1.GetMaximum(), hnonbjetB2W2.GetMaximum()), hnonbjetB2C.GetMaximum())*1.1);
  hnonbjetB2W1.Draw();
  legnb2.AddEntry(&hnonbjetB2W1, Form("subset (%.f)", hnonbjetB2W1.Integral()), "l");
  hnonbjetB2W2.GetXaxis()->SetTitle("nonb jets per event in B2");
  hnonbjetB2W2.SetLineColor(3);
  hnonbjetB2W2.Write();
  hnonbjetB2W2.Draw("same");
  legnb2.AddEntry(&hnonbjetB2W2, Form("contamination (%.f)", hnonbjetB2W2.Integral()), "l");
  hnonbjetB2C.GetXaxis()->SetTitle("nonb jets per event in B2");
  hnonbjetB2C.SetLineColor(1);
  hnonbjetB2C.Write();
  hnonbjetB2C.Draw("same");
  legnb2.AddEntry(&hnonbjetB2C, Form("correct (%.f)", hnonbjetB2C.Integral()), "l");
  legnb2.Draw();
  c.Update();
  c.Print("cpp_reconstructed_nonbjetcountB2truthmatched.eps");
  c.Clear();
  hnonbjetBISR.GetXaxis()->SetTitle("nonb jets per event in BISR");
  hnonbjetBISR.Write();
  hnonbjetBISR.Draw();
  c.Update();
  c.Print("cpp_reconstructed_nonbjetcountBISR.eps");
  c.Clear();

  f.Close();

  inFile.close();
  return 0;
  
}
