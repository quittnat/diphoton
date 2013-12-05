#include <vector>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include "TTree.h"
#include "TString.h"
#include "TFile.h"

const int nbins = 26;
const bool exclude_EEEE = 0;
const bool exclude_EBEE = 0;

using namespace std;

int min_integral(vector<TTree*> &t, TString var, float a, float b){
  vector<int> bla;
  for (int i=0; i<3; i++) {
    if (exclude_EEEE && i==2) continue;
    if (exclude_EBEE && i==1) continue;
    bla.push_back(t.at(i)->GetEntries(Form("%s>=%f && %s<%f",var.Data(),a,var.Data(),b)));
  }
  sort(bla.begin(),bla.end());
  return bla[0];
};

void determine_binning(const TString filename, const TString var, const float min, const float max, const int splits, const int min_events){

  TFile *f = new TFile(filename.Data(),"read");
  TDirectoryFile *dir; 
  f->GetObject("roofit",dir);

  TString title[3];
  title[0] = Form("obs_roodset_EBEB_%s_b%d",var.Data(),nbins);
  title[1] = Form("obs_roodset_EBEE_%s_b%d",var.Data(),nbins);
  title[2] = Form("obs_roodset_EEEE_%s_b%d",var.Data(),nbins);

  vector<TTree*> t;
  for (int i=0; i<3; i++) { t.push_back(NULL); dir->GetObject(title[i].Data(),t[i]); assert (t[i]);}

  vector<float> b;
  b.push_back(min);

  float split = (max-min)/splits;
  float low = min;
  float high = min;
  do {
    do {
      high+=split;
    }
    while (high<max && min_integral(t,var,low,high)<min_events);
    if (high<max) b.push_back(high); else b.at(b.size()-1)=max;    
    low=high;
  } 
  while (high<max);

  cout << var.Data() << ": " << b.size()-1 << " bins" << endl;
  cout << "{"; for (unsigned int i=0; i<b.size(); i++) cout << b.at(i) << (i<b.size()-1 ? "," : ""); cout << "}" << endl;
  for (unsigned int i=0; i<b.size()-1; i++) cout << min_integral(t,var,b.at(i),b.at(i+1)) << (i<b.size()-2 ? "," : ""); cout << endl;


}
