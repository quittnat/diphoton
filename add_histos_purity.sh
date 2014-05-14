for i in invmass diphotonpt costhetastar dphi dR njets 1jet_jpt 1jet_dR_lead_j 1jet_dR_trail_j 1jet_dR_close_j 1jet_dR_far_j 2jet_j1pt 2jet_j2pt 2jet_deta_jj 2jet_dphi_jj 2jet_dR_jj 2jet_mjj 2jet_zeppen 2jet_dphi_gg_jj; do for j in EBEB EBEE EEEE; do hadd histo_purity_${i}_${j}_allbins.root histo_purity_${i}_${j}_b*; done; done


