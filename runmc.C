{

  
  gROOT->ProcessLine(".L template_production.C+O");
  TString dir;

  //    dir = TString("~/nobackup/ntuples/gg_minitree_mc_030903p1_28dec/");
  // dir = TString("~/nobackup/ntuples/gg_minitree_mc_030903p2_21jan_8TeV_test2/");
  dir = TString("~/nobackup/ntuples/gg_minitree_mc_030903p2_29gen_testRCone3020/");
  //    dir = TString("~/nobackup/ntuples/gg_minitree_mc_030903p2_30gen_test53v7a/");

 gen_templates(dir+TString("allmc.root"),"signal",0,"outphoton_allmc_sig.root","photoniso",5e6);
 gen_templates(dir+TString("allmc.root"),"randomcone",0,"outphoton_allmc_rcone.root","photoniso",5e6);

 return;

// gen_templates(dir+TString("allmc.root"),"sigsig",0,"outphoton_allmc_sigsig.root","photoniso",5e5);
// gen_templates(dir+TString("allmc.root"),"sigbkg",0,"outphoton_allmc_sigbkg.root","photoniso",-1);
// gen_templates(dir+TString("allmc.root"),"bkgbkg",0,"outphoton_allmc_bkgbkg.root","photoniso",-1);
// gen_templates(dir+TString("allmc.root"),"standard",0,"outphoton_allmc_standard.root","photoniso",-1);
//
// gen_templates(dir+TString("allmc.root"),"standard_2pgen",0,"outphoton_allmc_standard_2pgen.root","photoniso",-1);
// gen_templates(dir+TString("allmc.root"),"standard_2fgen",0,"outphoton_allmc_standard_2fgen.root","photoniso",-1);



 gen_templates(dir+TString("allmc.root"),"signal",0,"outphoton_allmc_sig.root","photoniso",5e6);
// gen_templates(dir+TString("allmc.root"),"background",0,"outphoton_allmc_bkg.root","photoniso",5e5);
// gen_templates(dir+TString("allmc.root"),"2pgen",0,"outphoton_allmc_2pgen.root","photoniso",-1);
// gen_templates(dir+TString("allmc.root"),"1p1fbothgen",0,"outphoton_allmc_1p1fbothgen.root","photoniso",-1);
// gen_templates(dir+TString("allmc.root"),"1prcone1fgen",0,"outphoton_allmc_1prcone1fgen.root","photoniso",-1);
// gen_templates(dir+TString("allmc.root"),"1pgen1fside",0,"outphoton_allmc_1pgen1fside.root","photoniso",-1);
// gen_templates(dir+TString("allmc.root"),"2fgen",0,"outphoton_allmc_2fgen.root","photoniso",-1);
 gen_templates(dir+TString("allmc.root"),"randomcone",0,"outphoton_allmc_rcone.root","photoniso",5e6);
// gen_templates(dir+TString("allmc.root"),"sieiesideband",0,"outphoton_allmc_sieiesideband.root","photoniso",-1);
//
//
//
// gen_templates(dir+TString("DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"standard",0,"outphoton_dy_standard.root","photoniso",-1);
// gen_templates(dir+TString("DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"standard_pixelrev",0,"outphoton_dy_pixelrev.root","photoniso",-1);



 // gen_templates(dir+TString("allmc.root"),"standard_2frag",0,"outphoton_allmc_standard_2frag.root","photoniso",-1);
// gen_templates(dir+TString("allmc.root"),"signal_2frag",0,"outphoton_allmc_sig_2frag.root","photoniso",5e5);
// gen_templates(dir+TString("allmc.root"),"fragmentation",0,"outphoton_allmc_frag.root","photoniso",5e5);
// gen_templates(dir+TString("allmc.root"),"nofragmentation",0,"outphoton_allmc_nofrag.root","photoniso",5e5);


}


