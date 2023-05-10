void implementantiisosf() {
	TFile *filein=new TFile("/scratchnvme/bruschin/Newtest/WRemnants/wremnants/data/testMuonSF/allSmooth_GtoH3Dsfiso2.root");
	TFile *filein2=new TFile("/scratchnvme/bruschin/Newtest/smoothLeptonScaleFactorsredofullspline/GtoH/allSmooth_GtoH.root");
	TFile *fileout=new TFile("/scratchnvme/bruschin/Newtest/WRemnants/wremnants/data/testMuonSF/allSmooth_GtoH3Dsfiso3_.root","RECREATE");
	fileout->cd();
	for (auto&& keyAsObj : *(filein->GetListOfKeys())){
		auto key = (TKey*) keyAsObj;
		if (std::string(key->GetClassName())==std::string("TH2D")) {
			TH2D* histo=(TH2D*)filein->Get(key->GetName())->Clone();
			histo->SetName(key->GetName());
			histo->Write();
		}
		if (std::string(key->GetClassName())==std::string("TH3D")) {
			TH3D* histo=(TH3D*)filein->Get(key->GetName())->Clone();
			if (std::string(key->GetName())==std::string("SF_nomiAndAlt_onlyDataVar_GtoH_antiiso_both")) {
				TH3D* histo2=(TH3D*)filein->Get("SF_nomiAndAlt_onlyDataVar_GtoH_iso_both");
				TH3D* histo3=(TH3D*)filein2->Get("effMC_nomiAndAlt_GtoH_iso_both");
				for (unsigned int i=0; i!=histo->GetXaxis()->GetNbins(); i++) {
					for (unsigned int j=0; j!=histo->GetYaxis()->GetNbins(); j++) {
						for (unsigned int h=0; h!=histo->GetZaxis()->GetNbins(); h++) {
							double antiisosf=(1.-histo2->GetBinContent(i+1,j+1,h+1)*histo3->GetBinContent(i+1,j+1,1))/(1.-histo3->GetBinContent(i+1,j+1,1));
							histo->SetBinContent(i+1,j+1,h+1,antiisosf);
						}
					}
				}
			}
			if (std::string(key->GetName())==std::string("SF_nomiAndAlt_onlyMCVar_GtoH_antiiso_both")) {
				TH3D* histo2=(TH3D*)filein->Get("SF_nomiAndAlt_onlyDataVar_GtoH_iso_both");
				TH3D* histo3=(TH3D*)filein2->Get("effMC_nomiAndAlt_GtoH_iso_both");
				for (unsigned int i=0; i!=histo->GetXaxis()->GetNbins(); i++) {
					for (unsigned int j=0; j!=histo->GetYaxis()->GetNbins(); j++) {
						for (unsigned int h=0; h!=histo->GetZaxis()->GetNbins(); h++) {
							double antiisosf=(1.-histo2->GetBinContent(i+1,j+1,1)*histo3->GetBinContent(i+1,j+1,1))/(1.-histo3->GetBinContent(i+1,j+1,1));
							histo->SetBinContent(i+1,j+1,h+1,antiisosf);
						}
					}
				}
			}
			histo->SetName(key->GetName());
			histo->Write();
		}
	}
	fileout->Close();
}
