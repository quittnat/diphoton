#ifndef __BINSDEF__
#define __BINSDEF__

#include <map>
#include <vector>
#include "TString.h"
#include "TMath.h"
#include <iostream>
#include <assert.h>
const float Pi = TMath::Pi();
const float MaxDrExperiment = sqrt(pow(5.,2)+pow(Pi,2)); // = 5.90

TString __variables__[] = {
TString("invmass"),\
TString("diphotonpt"),\
TString("costhetastar"),\
TString("dphi"),\
TString("dR"),\
TString("njets"),\
TString("1jet_jpt"),\
TString("1jet_dR_lead_j"),\
TString("1jet_dR_trail_j"),\
TString("1jet_dR_close_j"),\
TString("1jet_dR_far_j"),\
TString("2jet_j1pt"),\
TString("2jet_j2pt"),\
TString("2jet_deta_jj"),\
TString("2jet_dphi_jj"),\
TString("2jet_dR_jj"),\
TString("2jet_mjj"),\
TString("2jet_zeppen"),\
TString("2jet_dphi_gg_jj")\
};
TString __names__[] = {
TString("m_{#gamma #gamma}, incl."),\
TString("p_{T, #gamma #gamma}, incl."),\
TString("|cos #theta^{*}|, incl."),\
TString("#Delta #phi_{#gamma #gamma}, incl."),\
TString("#Delta R_{#gamma #gamma}, incl."),\
TString("N_{jets}"),\
TString("p_{T,j}, #geq 1j"),\
TString("p_{T,#gamma ,j}^{lead}, #geq 1j"),\
TString("p_{T,#gamma ,j}^{trail}, #geq 1j"),\
TString("p_{T,#gamma ,j}^{close}, #geq 1j"),\
TString("p_{T,#gamma ,j}^{far}, #geq 1j"),\
TString("p_{T,j}^{lead}, #geq 2j"),\
TString("p_{T,j}^{trail}, #geq 2j"),\
TString("#Delta #eta_{jj}, #geq 2j"),\
TString("#Delta #phi_{jj}, #geq 2j"),\
TString("#Delta R_{jj}, #geq 2j"),\
TString("m_{jj}, #geq 2j"),\
TString("Zeppenfeld var., #geq 2j"),\
TString("#Delta #phi_{#gamma #gamma,jj}, #geq 2j")\
};
TString __units__[] = {
TString("GeV"),\
TString("GeV"),\
TString(""),\
TString(""),\
TString(""),\
TString(""),\
TString("GeV"),\
TString(""),\
TString(""),\
TString(""),\
TString(""),\
TString("GeV"),\
TString("GeV"),\
TString(""),\
TString(""),\
TString(""),\
TString("GeV"),\
TString(""),\
TString("")\
};


static const int n_bins=30;
int __nbins__[] = { // this should always be the effective length of mybinsdef_ array minus 1 (last number there is for overflow)
  18,25,13,26,22,5,10,13,14,14,13,5,4,5,5,5,5,5,5
};

float mybinsdef_invmass[n_bins+1]={0,40,60,68,72,76,80,84,88,92,96,100,108,116,132,168,228,800,800.01};
float mybinsdef_diphotonpt[n_bins+1]={0,5,7,9,11,13,15,16,17,19,21,23,25,27,30,33,36,40,45,51,59,67,75,85,200,200.01};
float mybinsdef_costhetastar[n_bins+1]={0,0.15,0.185,0.215,0.24,0.265,0.29,0.32,0.36,0.415,0.784999,0.944999,1,1.01};
float mybinsdef_dphi[n_bins+1]={0,0.518363,0.84823,1.27235,1.68075,1.9635,2.1677,2.34049,2.46615,2.5604,2.63894,2.70177,2.7646,2.81173,2.85885,2.89027,2.92168,2.9531,2.98451,3.01593,3.04735,3.06305,3.07876,3.09447,3.11018,3.14159,3.15159};
float mybinsdef_dR[n_bins+1]={0,0.885757,1.38769,1.86009,2.18487,2.39155,2.53917,2.65727,2.74585,2.8049,2.86395,2.923,2.95253,2.98205,3.01158,3.0411,3.07063,3.10015,3.12968,3.1592,3.69066,5.90505,5.91505};
float mybinsdef_njets[n_bins+1]={0,1,2,3,5,5.01};
float mybinsdef_1jet_jpt[n_bins+1]={0,25,30,35,40,50,60,70,85,600,600.01};
float mybinsdef_1jet_dR_lead_j[n_bins+1]={0,1.91914,2.3325,2.5687,2.74585,2.923,3.07063,3.21825,3.39541,3.60208,3.86781,4.25164,5.90505,5.91505};
float mybinsdef_1jet_dR_trail_j[n_bins+1]={0,1.32864,1.65341,1.91914,2.18487,2.42107,2.65727,2.86395,3.07063,3.30683,3.57256,3.92686,4.39926,5.90505,5.91505};
float mybinsdef_1jet_dR_close_j[n_bins+1]={0,1.26959,1.53531,1.74199,1.94867,2.15534,2.36202,2.5687,2.74585,2.923,3.12968,3.45446,3.98591,5.90505,5.91505};
float mybinsdef_1jet_dR_far_j[n_bins+1]={0,2.36202,2.62775,2.83443,2.98205,3.12968,3.2773,3.45446,3.63161,3.83828,4.10401,4.45832,5.90505,5.91505};
float mybinsdef_2jet_j1pt[n_bins+1]={0,40,55,75,600,600.01};
float mybinsdef_2jet_j2pt[n_bins+1]={0,25,35,600,600.01};
float mybinsdef_2jet_deta_jj[n_bins+1]={0,0.611,1.222,2.068,9.4,9.41};
float mybinsdef_2jet_dphi_jj[n_bins+1]={0,0.691151,1.33518,2.07345,3.14159,3.15159};
float mybinsdef_2jet_dR_jj[n_bins+1]={0,1.50579,2.36202,3.01158,5.90505,5.91505};
float mybinsdef_2jet_mjj[n_bins+1]={0,60,100,160,1000,1000.01};
float mybinsdef_2jet_zeppen[n_bins+1]={0,0.98,1.82,2.625,7,7.01};
float mybinsdef_2jet_dphi_gg_jj[n_bins+1]={0,2.46615,2.82743,2.98451,3.14159,3.15159};


float*  __binsdef__[] = {
mybinsdef_invmass,\
mybinsdef_diphotonpt,\
mybinsdef_costhetastar,\
mybinsdef_dphi,\
mybinsdef_dR,\
mybinsdef_njets,\
mybinsdef_1jet_jpt,\
mybinsdef_1jet_dR_lead_j,\
mybinsdef_1jet_dR_trail_j,\
mybinsdef_1jet_dR_close_j,\
mybinsdef_1jet_dR_far_j,\
mybinsdef_2jet_j1pt,\
mybinsdef_2jet_j2pt,\
mybinsdef_2jet_deta_jj,\
mybinsdef_2jet_dphi_jj,\
mybinsdef_2jet_dR_jj,\
mybinsdef_2jet_mjj,\
mybinsdef_2jet_zeppen,\
mybinsdef_2jet_dphi_gg_jj\
};

std::vector<TString> diffvariables_list (__variables__, __variables__ + sizeof(__variables__) / sizeof(TString) );
std::vector<TString> diffvariables_names_list_help_ (__names__, __names__ + sizeof(__names__) / sizeof(TString) );
std::vector<TString> diffvariables_units_list_help_ (__units__, __units__ + sizeof(__units__) / sizeof(TString) );
std::vector<int> diffvariables_nbins_list_help_ (__nbins__, __nbins__ + sizeof(__nbins__) / sizeof(int) );
std::vector<float*> diffvariables_binsdef_list_help_ (__binsdef__, __binsdef__ + sizeof(__binsdef__) / sizeof(float*) );

std::map<TString,TString> diffvariables_names_list_;
TString diffvariables_names_list(TString diffvariable){
  if (diffvariables_names_list_.size()==0) {
    assert(diffvariables_list.size()==diffvariables_names_list_help_.size());
    for (size_t i=0; i<diffvariables_names_list_help_.size(); i++){
      diffvariables_names_list_[diffvariables_list.at(i)]=diffvariables_names_list_help_.at(i);
    }
  }
  return diffvariables_names_list_[diffvariable];
};
std::map<TString,TString> diffvariables_units_list_;
TString diffvariables_units_list(TString diffvariable){
  if (diffvariables_units_list_.size()==0) {
    assert(diffvariables_list.size()==diffvariables_units_list_help_.size());
    for (size_t i=0; i<diffvariables_units_list_help_.size(); i++){
      diffvariables_units_list_[diffvariables_list.at(i)]=diffvariables_units_list_help_.at(i);
    }
  }
  return diffvariables_units_list_[diffvariable];
};
std::map<TString,int> diffvariables_nbins_list_;
int diffvariables_nbins_list(TString diffvariable){
  if (diffvariables_nbins_list_.size()==0) {
    assert(diffvariables_list.size()==diffvariables_nbins_list_help_.size());
    for (size_t i=0; i<diffvariables_nbins_list_help_.size(); i++){
      diffvariables_nbins_list_[diffvariables_list.at(i)]=diffvariables_nbins_list_help_.at(i);
    }
  }
  return diffvariables_nbins_list_[diffvariable];
};
std::map<TString,float*> diffvariables_binsdef_list_;
float* diffvariables_binsdef_list(TString diffvariable){
  if (diffvariables_binsdef_list_.size()==0) {
    assert(diffvariables_list.size()==diffvariables_binsdef_list_help_.size());
    for (size_t i=0; i<diffvariables_binsdef_list_help_.size(); i++){
      diffvariables_binsdef_list_[diffvariables_list.at(i)]=diffvariables_binsdef_list_help_.at(i);
    }
  }
  return diffvariables_binsdef_list_[diffvariable];
};

const int nclosest = 5;
const int nclosestmore = 40;

const Int_t n_histobins = 96;
const Float_t leftrange = -3;
const Float_t rightrange = 9;

const float default_threshold_adaptive_binning = -999;

const int dy_dataset_id = 30;

const float beam_energy = 7000;
const float pass_veto_closejets_dRcut = 1.0;

const int n_templatebins_max = 1000; 
int n_templatebins = 0;
Double_t templatebinsboundaries[n_templatebins_max+1];

int n_templates_pt=4;
float binsdef_single_gamma_pt[n_bins+1]={20,35,50,80,150};
int n_templates_EB_eta=7;
int n_templates_EE_eta=5;
float binsdef_single_gamma_EB_eta[n_bins+1]={0,0.2,0.4,0.6,0.8,1,1.2,1.4442};
float binsdef_single_gamma_EE_eta[n_bins+1]={1.566,1.653,1.8,2,2.2,2.5};

float AbsDeltaPhi(double phi1, double phi2){
  // From cmssw reco::deltaPhi()
  double result = phi1 - phi2;
  while( result >   TMath::Pi() ) result -= TMath::TwoPi();
  while( result <= -TMath::Pi() ) result += TMath::TwoPi();
  return TMath::Abs(result);
};

const int n_eta_cats = n_templates_EB_eta;
int n_eta1eta2_cats = n_eta_cats*n_eta_cats;
float *etabins = binsdef_single_gamma_EB_eta+0;

// FOR PHOTON COMPONENT
// 030903p1 2011 dataset
float eff_areas_EB_data[n_bins] = {2.703034e-01,2.678859e-01,2.722684e-01,2.720999e-01,2.643882e-01,2.480913e-01,1.706292e-01};
float eff_areas_EE_data[n_bins] = {5.329403e-02,7.733851e-02,1.091783e-01,1.339074e-01,1.068975e-01};
float eff_areas_EB_mc[n_bins] = {2.840231e-01,2.859291e-01,2.825974e-01,2.930248e-01,2.766801e-01,2.567621e-01,1.797386e-01};
float eff_areas_EE_mc[n_bins] = {5.101126e-02,7.539486e-02,1.088317e-01,1.385660e-01,1.077761e-01};

// // 030903p1 2012 dataset
// float eff_areas_EB_data[n_bins] = {3.070820e-01,3.078518e-01,3.089162e-01,3.103081e-01,3.048204e-01,2.865540e-01,1.999668e-01};
// float eff_areas_EE_data[n_bins] = {5.330422e-02,7.396740e-02,9.679689e-02,1.125413e-01,8.631587e-02};



const int n_ptbins_forreweighting = 4;
Float_t ptbins_forreweighting[n_ptbins_forreweighting+1]={20,35,50,80,999};
//const int n_ptbins_forreweighting = 1;
//Float_t ptbins_forreweighting[n_ptbins_forreweighting+1]={0,300};

typedef struct {
  TString name;
  TString title;
  bool is_on_raw;
  bool is_on_effunf;
  bool is_uncorrelated;
  bool is_1catcorrelated;
  bool is_allcatcorrelated;
} source_systematic_struct;

source_systematic_struct ConstructSystematic(TString name_, TString title_, bool is_on_raw_, bool is_on_effunf_, bool is_uncorrelated_, bool is_1catcorrelated_, bool is_allcatcorrelated_) {
  source_systematic_struct a;
  a.name=name_;
  a.title=title_;
  a.is_on_raw=is_on_raw_;
  a.is_on_effunf=is_on_effunf_;
  a.is_uncorrelated=is_uncorrelated_;
  a.is_1catcorrelated=is_1catcorrelated_;
  a.is_allcatcorrelated=is_allcatcorrelated_;
  assert((int)(a.is_uncorrelated)+(int)(a.is_1catcorrelated)+(int)(a.is_allcatcorrelated)==1);
  assert((int)(a.is_on_raw)+(int)(a.is_on_effunf)==1);
  return a;
};

source_systematic_struct __systematics__[]={
  ConstructSystematic("purefitbias","Fit bias",1,0,1,0,0);
  ConstructSystematic("zee","Zee subtraction",0,1,0,1,0);
  ConstructSystematic("templatestatistics","Template statistics",1,0,0,1,0);
  ConstructSystematic("efficiency","Efficiency uncertainty",0,1,0,1,0); // to be splitted in different scale factors?
  ConstructSystematic("unfolding","Unfolding uncertainty",0,1,0,1,0);
  ConstructSystematic("templateshapeMCpromptdrivenEB","Prompt template shape EB",1,0,0,0,1);
  ConstructSystematic("templateshapeMCfakedrivenEB","Fake template shape EB",1,0,0,0,1);
  ConstructSystematic("templateshapeMCpromptdrivenEE","Prompt template shape EE",1,0,0,0,1);
  ConstructSystematic("templateshapeMCfakedrivenEE","Fake template shape EE",1,0,0,0,1);
  ConstructSystematic("templateshape2frag","Fragmentation effect on template",1,0,0,0,1);
  ConstructSystematic("noise_mixing","Event mixing effect on template",1,0,0,0,1);
};
std::vector<source_systematic_struct> systematics_list (__systematics__, __systematics__ + sizeof(__systematics__) / sizeof(source_systematic_struct) );



#endif
