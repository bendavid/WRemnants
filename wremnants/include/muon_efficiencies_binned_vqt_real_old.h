#ifndef WREMNANTS_MUON_EFFICIENCIES_BINNED_VQT_REAL_H
#define WREMNANTS_MUON_EFFICIENCIES_BINNED_VQT_REAL_H

#include <boost/histogram/axis.hpp>
#include "TFile.h"
#include "TH2F.h"
#include <vector>
#include <array>
#include <string>

std::vector<TH2F*> VectorZQT;
std::vector<std::vector<TH2F*> > VectorZQTTrigger;
std::vector<float> vqtbinning;

void initializeDifferentialSFsNormal(std::vector<std::string> names) {
  for (unsigned int i=0; i!=names.size(); i++) {
    TFile *tfile=new TFile(names[i].c_str());
    TH2F* Histo=((TH2F*)tfile->Get("SF2D_nominal")->Clone());
	std::string name("SF2D_nominal_");
	name+=std::to_string(i);
    Histo->SetName(name.c_str());
    VectorZQT.push_back(Histo);
  }
}
void initializeDifferentialSFsNormalTrigger(std::vector<std::string> namesplus, std::vector<std::string> namesminus) {
  std::vector<TH2F*> VectorZQTTriggerPlus;
  for (unsigned int i=0; i!=namesplus.size(); i++) {
    TFile *tfile=new TFile(namesplus[i].c_str());
    TH2F* Histo=((TH2F*)tfile->Get("SF2D_nominal")->Clone());
	std::string name("SF2D_nominal_");
	name+=std::to_string(i);
    Histo->SetName(name.c_str());
    VectorZQTTriggerPlus.push_back(Histo);
  }
  std::vector<TH2F*> VectorZQTTriggerMinus;
  for (unsigned int i=0; i!=namesminus.size(); i++) {
    TFile *tfile=new TFile(namesminus[i].c_str());
    TH2F* Histo=((TH2F*)tfile->Get("SF2D_nominal")->Clone());
	std::string name("SF2D_nominal_");
	name+=std::to_string(i);
    Histo->SetName(name.c_str());
    VectorZQTTriggerMinus.push_back(Histo);
  }
  VectorZQTTrigger.push_back(VectorZQTTriggerMinus);
  VectorZQTTrigger.push_back(VectorZQTTriggerPlus);
}

void initializeVQTBinning(std::vector<float> binning) {
  vqtbinning=binning;
}

namespace wrem_vqt_real {

    template<typename HIST_IDIPTRIGISO, typename HIST_TRACKING, typename HIST_RECO>
    class muon_efficiency_binned_helper_base {
    public:

        muon_efficiency_binned_helper_base(HIST_IDIPTRIGISO &&sf_idip_trig_iso, HIST_TRACKING &&sf_tracking, HIST_RECO &&sf_reco, bool includeTrigger) :
            sf_idip_trig_iso_(std::make_shared<const HIST_IDIPTRIGISO>(std::move(sf_idip_trig_iso))),
            sf_tracking_(std::make_shared<const HIST_TRACKING>(std::move(sf_tracking))),
            sf_reco_(std::make_shared<const HIST_RECO>(std::move(sf_reco))) {includeTrigger_=includeTrigger;}

        std::array<double,5> scale_factor_array(int pt_idx, int pt_idx_reco, int sapt_idx,
                                                int eta_idx, int saeta_idx,
                                                int charge_idx,
                                                bool pass_iso, bool with_trigger,
                                                int idx_nom_alt) const {

            auto const eff_type_idx_reco = idx_reco_;
            auto const eff_type_idx_tracking = idx_tracking_;
            auto const eff_type_idx_idip = idx_idip_;
            auto const eff_type_idx_trig = idx_trig_;
            auto const eff_type_idx_iso_pass = with_trigger ? idx_iso_triggering_ : idx_iso_nontriggering_;
            auto const eff_type_idx_iso_fail = with_trigger ? idx_antiiso_triggering_ : idx_antiiso_nontriggering_;
            auto const eff_type_idx_iso = pass_iso ? eff_type_idx_iso_pass : eff_type_idx_iso_fail;

            const double reco      = sf_reco_->at(eta_idx, pt_idx_reco, charge_idx, eff_type_idx_reco,      idx_nom_alt).value();
            const double tracking  = sf_tracking_->at(saeta_idx, sapt_idx,    charge_idx, eff_type_idx_tracking,  idx_nom_alt).value();
            const double idip      = sf_idip_trig_iso_->at(  eta_idx,   pt_idx,    charge_idx, eff_type_idx_idip,      idx_nom_alt).value();
            double trig = 1.0;
            if (with_trigger) trig = sf_idip_trig_iso_->at(  eta_idx,   pt_idx,    charge_idx, eff_type_idx_trig,      idx_nom_alt).value();
            const double iso       = sf_idip_trig_iso_->at(  eta_idx,   pt_idx,    charge_idx, eff_type_idx_iso,       idx_nom_alt).value();

            std::array<double,5> ret = {reco, tracking, idip, trig, iso};
            // for(int i = 0; i < ret.size(); i++)
            //     std::cout << "Scale factor i = " << i << " --> " << ret[i] << std::endl;

            return ret;
            
        }

        double scale_factor_product(float pt, float eta, float sapt, float saeta, int charge, bool pass_iso, bool with_trigger, float vqt, int idx_nom_alt) const {

            auto const eta_idx = sf_idip_trig_iso_->template axis<0>().index(eta);
            auto const pt_idx  = sf_idip_trig_iso_->template axis<1>().index(pt);
            auto const charge_idx = sf_idip_trig_iso_->template axis<2>().index(charge);
            auto const pt_idx_reco = sf_reco_->template axis<1>().index(pt);
            auto const saeta_idx   = sf_tracking_->template axis<0>().index(saeta);
            auto const sapt_idx    = sf_tracking_->template axis<1>().index(sapt);

            std::array<double,5> allSF = scale_factor_array(pt_idx, pt_idx_reco, sapt_idx, eta_idx, saeta_idx, charge_idx, pass_iso, with_trigger, idx_nom_alt);
            double sf = 1.0;
            if (includeTrigger_) {
              for(int i = 0; i < (allSF.size()-2); i++) {
                  // std::cout << "Scale factor i = " << i << " --> " << allSF[i] << std::endl;
                  sf *= allSF[i];
              }
            }
            else {
              for(int i = 0; i < (allSF.size()-1); i++) {
                  // std::cout << "Scale factor i = " << i << " --> " << allSF[i] << std::endl;
                  sf *= allSF[i];
              }
            }
            bool cond=false;
            int i1=-1;
            for (unsigned int i=0; i!=(vqtbinning.size()-1); i++) {
              if ((vqt>=vqtbinning[i])&&(vqt<vqtbinning[i+1])) {
                i1=i;
		        cond=true;
              }
            }
            if (cond) {
              if ((vqt>=vqtbinning[0])&&(vqt<=vqtbinning[vqtbinning.size()-1])) sf *= VectorZQT[i1]->GetBinContent(VectorZQT[i1]->FindBin(eta,pt));
              else {
                if (vqt<vqtbinning[0]) sf *= VectorZQT[0]->GetBinContent(VectorZQT[0]->FindBin(eta,pt));
                else sf *= VectorZQT[vqtbinning.size()-2]->GetBinContent(VectorZQT[vqtbinning.size()-2]->FindBin(eta,pt));
              }
              if (includeTrigger_) {
                unsigned int chargeidx = (1+charge)/2;
                if ((vqt>=vqtbinning[0])&&(vqt<=vqtbinning[vqtbinning.size()-1])) sf *= VectorZQTTrigger[chargeidx][i1]->GetBinContent(VectorZQTTrigger[chargeidx][i1]->FindBin(eta,pt));
                else {
                  if (vqt<vqtbinning[0]) sf *= VectorZQTTrigger[chargeidx][0]->GetBinContent(VectorZQTTrigger[chargeidx][0]->FindBin(eta,pt));
                  else sf *= VectorZQTTrigger[chargeidx][vqtbinning.size()-2]->GetBinContent(VectorZQTTrigger[chargeidx][vqtbinning.size()-2]->FindBin(eta,pt));
                }
              }
	        }
            //std::cout << "Scale factor product " << sf << std::endl;
            return sf;

        }

        using syst_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<5>>; // 5 bins for reco, tracking, idip, trigger, iso(notrig) in this order

        syst_tensor_t sf_syst_var(float pt, float eta, float sapt, float saeta, int charge, bool pass_iso, bool with_trigger) const {

            syst_tensor_t res;

            auto const eta_idx = sf_idip_trig_iso_->template axis<0>().index(eta);
            auto const pt_idx  = sf_idip_trig_iso_->template axis<1>().index(pt);
            auto const charge_idx = sf_idip_trig_iso_->template axis<2>().index(charge);
            auto const pt_idx_reco = sf_reco_->template axis<1>().index(pt);
            auto const saeta_idx   = sf_tracking_->template axis<0>().index(saeta);
            auto const sapt_idx    = sf_tracking_->template axis<1>().index(sapt);

            std::array<double,5> allSF_nomi = scale_factor_array(pt_idx, pt_idx_reco, sapt_idx, eta_idx, saeta_idx, charge_idx, pass_iso, with_trigger, idx_nom_);
            std::array<double,5> allSF_alt  = scale_factor_array(pt_idx, pt_idx_reco, sapt_idx, eta_idx, saeta_idx, charge_idx, pass_iso, with_trigger, idx_alt_);

            // anticorrelation between iso and antiiso already embedded in the numbers stored in the histograms
            // also the alternate comes from data efficiency variation only, so the anticorrelation in the efficiencies is preserved in the scale factors
            
            // order is reco-tracking-idip-trigger-iso
            for(int i = 0; i < allSF_nomi.size(); i++) {
                res(i) = allSF_alt[i] / allSF_nomi[i]; 
            }
                
            return res;
     
        }

    protected:

        std::shared_ptr<const HIST_IDIPTRIGISO> sf_idip_trig_iso_;
        std::shared_ptr<const HIST_TRACKING> sf_tracking_;
        std::shared_ptr<const HIST_RECO> sf_reco_;

        // hardcoded things to keep original pt binning for tnp when using smoothed efficiencies
        //boost::histogram::axis::variable<double> originalTnpPtBins{24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 47., 50., 55., 60., 65.};
        
        // cache the bin indices since the string category lookup is slow
        int idx_reco_ = sf_reco_->template axis<3>().index("reco");
        int idx_tracking_ = sf_tracking_->template axis<3>().index("tracking");
        int idx_idip_ = sf_idip_trig_iso_->template axis<3>().index("idip");
        int idx_trig_ = sf_idip_trig_iso_->template axis<3>().index("trigger");
        int idx_iso_triggering_        = sf_idip_trig_iso_->template axis<3>().index("iso");
        int idx_antiiso_triggering_    = sf_idip_trig_iso_->template axis<3>().index("antiiso");
        int idx_iso_nontriggering_     = sf_idip_trig_iso_->template axis<3>().index("isonotrig");
        int idx_antiiso_nontriggering_ = sf_idip_trig_iso_->template axis<3>().index("antiisonotrig");

        int idx_nom_ = sf_idip_trig_iso_->template axis<4>().index(0);
        int idx_alt_ = sf_idip_trig_iso_->template axis<4>().index(1);

        bool includeTrigger_;
    };

    // base template for one-lepton case
    template<bool do_other, typename HIST_IDIPTRIGISO, typename HIST_TRACKING, typename HIST_RECO>
    class muon_efficiency_binned_helper:
        public muon_efficiency_binned_helper_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> {

    public:

        using base_t = muon_efficiency_binned_helper_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO>;
        // inherit constructor
        using base_t::base_t;

        muon_efficiency_binned_helper(const base_t &other) : base_t(other) {}
        
        double operator() (float pt, float eta, float sapt, float saeta, int charge, bool pass_iso, float vqt) {
            constexpr bool with_trigger = true;
            return base_t::scale_factor_product(pt, eta, sapt, saeta, charge, pass_iso, with_trigger, vqt, base_t::idx_nom_);
        }

    };

    // specialization for two-lepton case
    template<typename HIST_IDIPTRIGISO, typename HIST_TRACKING, typename HIST_RECO>
    class muon_efficiency_binned_helper<true, HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> :
        public muon_efficiency_binned_helper_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> {

    public:

        using base_t = muon_efficiency_binned_helper_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO>;
        // inherit constructor
        using base_t::base_t;
        
        muon_efficiency_binned_helper(const base_t &other) : base_t(other) {}

        double operator() (float trig_pt, float trig_eta, float trig_sapt, float trig_saeta, int trig_charge,
                           float nontrig_pt, float nontrig_eta, float nontrig_sapt, float nontrig_saeta, int nontrig_charge) {
            constexpr bool with_trigger = true;
            constexpr bool without_trigger = false;
            constexpr bool pass_iso = true;
            const double sftrig = base_t::scale_factor_product(trig_pt, trig_eta, trig_sapt, trig_saeta, trig_charge, pass_iso, with_trigger, base_t::idx_nom_);
            const double sfnontrig = base_t::scale_factor_product(nontrig_pt, nontrig_eta, nontrig_sapt, nontrig_saeta, nontrig_charge, pass_iso, without_trigger, base_t::idx_nom_);
            return sftrig*sfnontrig;
        }

    };
        
    // base template for one lepton case
    template<bool do_other, typename HIST_IDIPTRIGISO, typename HIST_TRACKING, typename HIST_RECO>
    class muon_efficiency_binned_helper_syst :
        public muon_efficiency_binned_helper_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> {

    public:

        using base_t = muon_efficiency_binned_helper_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO>;

        using tensor_t = typename base_t::syst_tensor_t;
        
        // inherit constructor
        using base_t::base_t;
        
        muon_efficiency_binned_helper_syst(const base_t &other) : base_t(other) {}
        
        tensor_t operator() (float pt, float eta, float sapt, float saeta, int charge, bool pass_iso, double nominal_weight = 1.0) {
            constexpr bool with_trigger = true;
            return nominal_weight*base_t::sf_syst_var(pt, eta, sapt, saeta, charge, pass_iso, with_trigger);
        }

    };

    // specialization for two-lepton case
    template<typename HIST_IDIPTRIGISO, typename HIST_TRACKING, typename HIST_RECO>
    class muon_efficiency_binned_helper_syst<true, HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> :
        public muon_efficiency_binned_helper_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> {

    public:

        using base_t = muon_efficiency_binned_helper_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO>;
        using tensor_t = typename base_t::syst_tensor_t;

        // inherit constructor
        using base_t::base_t;
        
        muon_efficiency_binned_helper_syst(const base_t &other) : base_t(other) {}

        tensor_t operator() (float trig_pt, float trig_eta, float trig_sapt, float trig_saeta, int trig_charge,
                             float nontrig_pt, float nontrig_eta, float nontrig_sapt, float nontrig_saeta, int nontrig_charge, double nominal_weight = 1.0) {
            constexpr bool with_trigger = true;
            constexpr bool without_trigger = false;
            constexpr bool pass_iso = true;
            const tensor_t variation_trig = base_t::sf_syst_var(trig_pt, trig_eta, trig_sapt, trig_saeta, trig_charge, pass_iso, with_trigger);
            const tensor_t variation_nontrig = base_t::sf_syst_var(nontrig_pt, nontrig_eta, nontrig_sapt, nontrig_saeta, nontrig_charge, pass_iso, without_trigger);
            return nominal_weight*variation_trig*variation_nontrig;
        }

    };

    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////

    
    ////
    // STAT UNCERTAINTY
    // this is now an independent class with respect to the previous one which only deals with nominal and statistical variations
    ////
    //// BASE CLASS FOR HELPER_STAT
    template<int NEtaBins, int NPtBins, typename HIST_SF>
    class muon_efficiency_binned_helper_stat_base {
    public:

        muon_efficiency_binned_helper_stat_base(HIST_SF &&sf_type) :
            sf_type_(std::make_shared<const HIST_SF>(std::move(sf_type))) {}

        // number of eta bins, number of eigen variations for pt axis, then 2 charges
        using stat_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<NEtaBins, NPtBins, 2>>;
        
        int checkEffTypeInAxis(boost::histogram::axis::category<std::string> axis, const std::string& match = "match") {
            int ret = -1;
            for (Int_t i = 0; i < axis.size(); i++) {
                if (match == axis.value(i)) {
                    ret = i;
                    break;
                }
            }
            return ret;
        }
        
        // general case with no isolation
        stat_tensor_t sf_stat_var(float pt, float eta, int charge) const {
            stat_tensor_t res;
            res.setConstant(1.0);

            auto const eta_idx = sf_type_->template axis<0>().index(eta);
            auto const pt_idx = sf_type_->template axis<1>().index(pt);
            auto const charge_idx = sf_type_->template axis<2>().index(charge);
            auto const eff_type_idx = 0; // TODO FIXME, use first (and only existent) bin

            // overflow/underflow are attributed to adjacent bin
            auto const tensor_eta_idx = std::clamp(eta_idx, 0, NEtaBins - 1);
            auto const tensor_pt_idx =  std::clamp(pt_idx, 0, NPtBins - 1);

            auto const &cell_nomi = sf_type_->at(eta_idx,
                                                 pt_idx,
                                                 charge_idx,
                                                 eff_type_idx);

            const double sf_value   = cell_nomi.value();
            const double sf_statunc = std::sqrt(cell_nomi.variance());
            const double sf_stat_variation = (sf_value + sf_statunc) / sf_value;
                    
            res(tensor_eta_idx, tensor_pt_idx, charge_idx) *= sf_stat_variation;
                
            return res;
        }

        // special case for isolation
        stat_tensor_t sf_stat_var_iso(float pt, float eta, int charge, bool pass_iso, bool with_trigger) const {
            stat_tensor_t res;
            res.setConstant(1.0);

            auto const eta_idx = sf_type_->template axis<0>().index(eta);
            auto const pt_idx = sf_type_->template axis<1>().index(pt);
            auto const charge_idx = sf_type_->template axis<2>().index(charge);

            auto const eff_type_idx_iso_pass = with_trigger ? idx_iso_triggering_ : idx_iso_nontriggering_;
            auto const eff_type_idx_iso_fail = with_trigger ? idx_antiiso_triggering_ : idx_antiiso_nontriggering_;
            auto const eff_type_idx_iso = pass_iso ? eff_type_idx_iso_pass : eff_type_idx_iso_fail;

            // overflow/underflow are attributed to adjacent bin
            auto const tensor_eta_idx = std::clamp(eta_idx, 0, NEtaBins - 1);
            auto const tensor_pt_idx =  std::clamp(pt_idx, 0, NPtBins - 1);

            auto const &cell_nomi = sf_type_->at(eta_idx,
                                                 pt_idx,
                                                 charge_idx,
                                                 eff_type_idx_iso);
                
            const double sf_value   = cell_nomi.value();
            const double sf_statunc = std::sqrt(cell_nomi.variance());
            // anti-correlation between iso and anti-iso SF's is not exact, but an excellent approximation
            double sf_shifted =  sf_value + sf_statunc;
            if (not pass_iso) {
                sf_shifted = sf_value - sf_statunc;
                sf_shifted = std::clamp(sf_shifted, 0.0, sf_shifted);
            }
            const double sf_stat_variation = sf_shifted / sf_value;
                    
            res(tensor_eta_idx, tensor_pt_idx, charge_idx) *= sf_stat_variation;
                
            return res;
        }

    protected:

        std::shared_ptr<const HIST_SF> sf_type_;
        // cache the bin indices since the string category lookup is slow
        int isTriggerStep_ = sf_type_->template axis<3>().value(0) == "trigger"; // special treatment for the stat variation in 2 lepton case
        // check if axis name exists in histogram, return -1 (invalid index) if not found
        int idx_iso_triggering_        = checkEffTypeInAxis(sf_type_->template axis<3>(), "iso");
        int idx_antiiso_triggering_    = checkEffTypeInAxis(sf_type_->template axis<3>(), "antiiso");
        int idx_iso_nontriggering_     = checkEffTypeInAxis(sf_type_->template axis<3>(), "isonotrig");
        int idx_antiiso_nontriggering_ = checkEffTypeInAxis(sf_type_->template axis<3>(), "antiisonotrig");

    };

    ////
    //// General case (isolation is treated separately)
    ////

    // base template for one lepton case
    template<bool do_other, int NEtaBins, int NPtBins, typename HIST_SF>
    class muon_efficiency_binned_helper_stat :
        public muon_efficiency_binned_helper_stat_base<NEtaBins, NPtBins, HIST_SF> {
        
    public:
        
        using stat_base_t = muon_efficiency_binned_helper_stat_base<NEtaBins, NPtBins, HIST_SF>;
        using tensor_t = typename stat_base_t::stat_tensor_t;
  
        using stat_base_t::stat_base_t;
        
        tensor_t operator() (float pt, float eta, int charge, double nominal_weight = 1.0) {
            return nominal_weight*stat_base_t::sf_stat_var(pt, eta, charge);
        }

    };

    // specialization for two-lepton case
    template<int NEtaBins, int NPtBins, typename HIST_SF>
    class muon_efficiency_binned_helper_stat<true, NEtaBins, NPtBins, HIST_SF> :
        public muon_efficiency_binned_helper_stat_base<NEtaBins, NPtBins, HIST_SF> {

    public:

        using stat_base_t = muon_efficiency_binned_helper_stat_base<NEtaBins, NPtBins, HIST_SF>;
        using tensor_t = typename stat_base_t::stat_tensor_t;

        using stat_base_t::stat_base_t;

        tensor_t operator() (float trig_pt, float trig_eta, int trig_charge,
                             float nontrig_pt, float nontrig_eta, int nontrig_charge, double nominal_weight = 1.0) {

            const tensor_t variation_trig = stat_base_t::sf_stat_var(trig_pt, trig_eta, trig_charge);
            if (stat_base_t::isTriggerStep_) {
                return nominal_weight * variation_trig;
            } else {
                const tensor_t variation_nontrig = stat_base_t::sf_stat_var(nontrig_pt, nontrig_eta, nontrig_charge);
                return nominal_weight * variation_trig * variation_nontrig;
            }
        }

    };

    ////
    //// Isolation
    ////

    // base template for one lepton case
    template<bool do_other, int NEtaBins, int NPtBins, typename HIST_SF>
    class muon_efficiency_binned_helper_stat_iso :
        public muon_efficiency_binned_helper_stat_base<NEtaBins, NPtBins, HIST_SF> {
        
    public:
        
        using stat_base_t = muon_efficiency_binned_helper_stat_base<NEtaBins, NPtBins, HIST_SF>;
        using tensor_t = typename stat_base_t::stat_tensor_t;
  
        using stat_base_t::stat_base_t;
        
        tensor_t operator() (float pt, float eta, int charge, bool pass_iso, double nominal_weight = 1.0) {
            constexpr bool with_trigger = true;
            return nominal_weight*stat_base_t::sf_stat_var_iso(pt, eta, charge, pass_iso, with_trigger);
        }

    };

    // specialization for two-lepton case
    template<int NEtaBins, int NPtBins, typename HIST_SF>
    class muon_efficiency_binned_helper_stat_iso<true, NEtaBins, NPtBins, HIST_SF> :
        public muon_efficiency_binned_helper_stat_base<NEtaBins, NPtBins, HIST_SF> {

    public:

        using stat_base_t = muon_efficiency_binned_helper_stat_base<NEtaBins, NPtBins, HIST_SF>;
        using tensor_t = typename stat_base_t::stat_tensor_t;

        using stat_base_t::stat_base_t;

        tensor_t operator() (float trig_pt, float trig_eta, int trig_charge,
                             float nontrig_pt, float nontrig_eta, int nontrig_charge, double nominal_weight = 1.0) {
            constexpr bool with_trigger = true;
            constexpr bool without_trigger = false;
            constexpr bool pass_iso = true;

            const tensor_t variation_trig    = stat_base_t::sf_stat_var_iso(   trig_pt,    trig_eta,    trig_charge, pass_iso, with_trigger);
            const tensor_t variation_nontrig = stat_base_t::sf_stat_var_iso(nontrig_pt, nontrig_eta, nontrig_charge, pass_iso, without_trigger);
            return nominal_weight * variation_trig * variation_nontrig;
            
        }

    };


}

#endif
