void implementisosf() {
	TFile *filein=new TFile("/gpfs/ddn/cms/user/bruschin/Newtest/WRemnants/wremnants/data/testMuonSF/allSmooth_GtoH3Dsfisoabcd.root");
	TFile *fileout=new TFile("/gpfs/ddn/cms/user/bruschin/Newtest/WRemnants/wremnants/data/testMuonSF/allSmooth_GtoH3Dsfisohighpt.root","RECREATE");
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
			if (std::string(key->GetName())==std::string("SF_nomiAndAlt_onlyDataVar_GtoH_iso_both")) {
				TH3D* histo2=(TH3D*)filein->Get("SF_nomiAndAlt_GtoH_iso_both");
				int npar=(histo->GetZaxis()->GetNbins()-2)/2, npar2=(histo2->GetZaxis()->GetNbins()-2)/2;
				for (unsigned int i=0; i!=histo2->GetXaxis()->GetNbins(); i++) {
					for (unsigned int j=0; j!=histo2->GetYaxis()->GetNbins(); j++) {
						histo->SetBinContent(i+1,j+1,1,histo2->GetBinContent(i+1,j+1,1));
						histo->SetBinContent(i+1,j+1,histo->GetZaxis()->GetNbins(),histo2->GetBinContent(i+1,j+1,histo2->GetZaxis()->GetNbins()));
						for (int h=0; h!=npar2; h++) {
							histo->SetBinContent(i+1,j+1,2+h,histo2->GetBinContent(i+1,j+1,2+h));
							histo->SetBinContent(i+1,j+1,2+h+npar,histo2->GetBinContent(i+1,j+1,2+h+npar2));
						}
						for (int h=npar2; h!=npar; h++) {
							histo->SetBinContent(i+1,j+1,2+h,histo2->GetBinContent(i+1,j+1,1));
							histo->SetBinContent(i+1,j+1,2+h+npar,histo2->GetBinContent(i+1,j+1,1));
						}
					}
				}
			}
			if (std::string(key->GetName())==std::string("SF_nomiAndAlt_onlyMCVar_GtoH_iso_both")) {
				TH3D* histo2=(TH3D*)filein->Get("SF_nomiAndAlt_GtoH_iso_both");
				for (unsigned int i=0; i!=histo->GetXaxis()->GetNbins(); i++) {
					for (unsigned int j=0; j!=histo->GetYaxis()->GetNbins(); j++) {
						for (unsigned int h=0; h!=histo->GetZaxis()->GetNbins(); h++) {
							histo->SetBinContent(i+1,j+1,h+1,histo2->GetBinContent(i+1,j+1,1));
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
