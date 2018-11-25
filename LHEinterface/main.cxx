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
}


