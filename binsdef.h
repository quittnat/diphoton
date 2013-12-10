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


static const int n_bins=26;
int __nbins__[] = { // this should always be the effective length of mybinsdef_ array minus 1 (last number there is for overflow)
16,\
21,\
8,\
14,\
19,\
5,\
6,\
8,\
8,\
8,\
7,\
4,\
4,\
5,\
5,\
5,\
5,\
5,\
5
};

float mybinsdef_invmass[n_bins+1]={0,40,60,70,75,80,85,90,95,100,110,120,150,250,400,800,801};
float mybinsdef_diphotonpt[n_bins+1]={0,6,10,12,14,16,18,20,22,24,28,34,40,50,60,70,80,90,100,120,200,201};
float mybinsdef_costhetastar[n_bins+1]={0,0.20,0.28,0.36,0.44,0.60,0.90,1.00,1.01};
float mybinsdef_dphi[n_bins+1]={0,0.2*Pi,0.4*Pi,0.6*Pi,0.7*Pi,0.8*Pi,0.84*Pi,0.88*Pi,0.90*Pi,0.92*Pi,0.94*Pi,0.96*Pi,0.98*Pi,1.0*Pi,1.01*Pi};
float mybinsdef_dR[n_bins+1]={0,0.9,1.4,1.85,2.15,2.35,2.5,2.6,2.7,2.8,2.85,2.9,2.95,3,3.05,3.1,3.15,3.2,MaxDrExperiment,MaxDrExperiment+0.01};
float mybinsdef_njets[n_bins+1]={0,1,2,3,5,5.01};
float mybinsdef_1jet_jpt[n_bins+1]={0,40,50,60,80,1000,1001};
float mybinsdef_1jet_dR_lead_j[n_bins+1]={0,2.35,2.7,2.95,3.15,3.4,3.8,MaxDrExperiment,MaxDrExperiment+0.01};
float mybinsdef_1jet_dR_trail_j[n_bins+1]={0,1.25,1.85,2.35,2.75,3.1,3.55,MaxDrExperiment,MaxDrExperiment+0.01};
float mybinsdef_1jet_dR_close_j[n_bins+1]={0,1.15,1.7,2.15,2.55,2.85,3.2,MaxDrExperiment,MaxDrExperiment+0.01};
float mybinsdef_1jet_dR_far_j[n_bins+1]={0,2.65,2.95,3.15,3.4,3.7,MaxDrExperiment,MaxDrExperiment+0.01};
float mybinsdef_2jet_j1pt[n_bins+1]={0,50,70,1000,1001};
float mybinsdef_2jet_j2pt[n_bins+1]={0,40,50,1000,1001};
float mybinsdef_2jet_deta_jj[n_bins+1]={0,0.55,1.15,1.95,9.4,9.41};
float mybinsdef_2jet_dphi_jj[n_bins+1]={0,0.95,1.9,2.6,Pi,1.01*Pi};
float mybinsdef_2jet_dR_jj[n_bins+1]={0,1.7,2.6,3.1,MaxDrExperiment,MaxDrExperiment+0.01};
float mybinsdef_2jet_mjj[n_bins+1]={0,80,120,180,2000,2001};
float mybinsdef_2jet_zeppen[n_bins+1]={0,1.25,1.95,2.2,7,7.01};
float mybinsdef_2jet_dphi_gg_jj[n_bins+1]={0,2.55,2.9,3.05,Pi,1.01*Pi};

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
float binsdef_single_gamma_EE_eta[n_bins+1]={1.56,1.653,1.8,2,2.2,2.5};

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


#endif
