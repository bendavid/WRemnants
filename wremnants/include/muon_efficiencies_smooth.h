#ifndef WREMNANTS_MUON_EFFICIENCIES_SMOOTH_H
#define WREMNANTS_MUON_EFFICIENCIES_SMOOTH_H

#include <boost/histogram/axis.hpp>
#include <array>

namespace wrem {

    // TODO use enums for integer/boolean/category axes so that the code is less error-prone?

    template<int NSysts, typename HIST_SF>
    class muon_efficiency_smooth_helper_base {
    public:

        muon_efficiency_smooth_helper_base(HIST_SF &&sf_all) :
            sf_all_(std::make_shared<const HIST_SF>(std::move(sf_all))) {
        }
    
        std::array<double,5> scale_factor_array(int pt_idx, int eta_idx, int sapt_idx, int saeta_idx, int charge_idx, bool pass_iso, bool with_trigger, int idx_nom_alt) const {

            auto const eff_type_idx_reco = idx_reco_;
            auto const eff_type_idx_tracking = idx_tracking_;
            auto const eff_type_idx_idip = idx_idip_;
            auto const eff_type_idx_trig = idx_trig_;
            auto const eff_type_idx_iso_pass = with_trigger ? idx_iso_triggering_ : idx_iso_nontriggering_;
            auto const eff_type_idx_iso_fail = with_trigger ? idx_antiiso_triggering_ : idx_antiiso_nontriggering_;
            auto const eff_type_idx_iso = pass_iso ? eff_type_idx_iso_pass : eff_type_idx_iso_fail;

            const double reco      = sf_all_->at(  eta_idx,   pt_idx, charge_idx, eff_type_idx_reco,      idx_nom_alt).value();
            const double tracking  = sf_all_->at(saeta_idx, sapt_idx, charge_idx, eff_type_idx_tracking,  idx_nom_alt).value();
            const double idip      = sf_all_->at(  eta_idx,   pt_idx, charge_idx, eff_type_idx_idip,      idx_nom_alt).value();
            double trig = 1.0;
            if (with_trigger) trig = sf_all_->at(  eta_idx,   pt_idx, charge_idx, eff_type_idx_trig,      idx_nom_alt).value();
            const double iso       = sf_all_->at(  eta_idx,   pt_idx, charge_idx, eff_type_idx_iso,       idx_nom_alt).value();

            std::array<double,5> ret = {reco, tracking, idip, trig, iso};
            // for(int i = 0; i < ret.size(); i++)
            //     std::cout << "Scale factor i = " << i << " --> " << ret[i] << std::endl;

            return ret;
            
        }

        double scale_factor_product(float pt, float eta, float sapt, float saeta, int charge, bool pass_iso, bool with_trigger, int idx_nom_alt) const {

            auto const eta_idx = sf_all_->template axis<0>().index(eta);
            auto const pt_idx = sf_all_->template axis<1>().index(pt);
            auto const charge_idx = sf_all_->template axis<2>().index(charge);
            auto const saeta_idx = sf_all_->template axis<0>().index(saeta);
            auto const sapt_idx = sf_all_->template axis<1>().index(sapt);

            std::array<double,5> allSF = scale_factor_array(pt_idx, eta_idx, sapt_idx, saeta_idx, charge_idx, pass_iso, with_trigger, idx_nom_alt);
            double sf = 1.0;
            for(int i = 0; i < allSF.size(); i++) {
                // std::cout << "Scale factor i = " << i << " --> " << allSF[i] << std::endl;
                sf *= allSF[i];
            }
            // std::cout << "Scale factor product " << sf << std::endl;
            return sf;

        }

        using syst_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<5, NSysts>>; // 5 bins for reco, tracking, idip, trigger, iso(notrig) in this order

        syst_tensor_t sf_syst_var(float pt, float eta, float sapt, float saeta, int charge, bool pass_iso, bool with_trigger) const {

            syst_tensor_t res;

            auto const eta_idx =     sf_all_->template axis<0>().index(eta);
            auto const pt_idx =      sf_all_->template axis<1>().index(pt);
            auto const saeta_idx =   sf_all_->template axis<0>().index(saeta);
            auto const sapt_idx =    sf_all_->template axis<1>().index(sapt);
            auto const charge_idx =  sf_all_->template axis<2>().index(charge);

            std::array<double,5> allSF_nomi = scale_factor_array(pt_idx, eta_idx, sapt_idx, saeta_idx, charge_idx, pass_iso, with_trigger, idx_nom_);

            for(int ns = 0; ns < NSysts; ns++) {
                
                std::array<double,5> allSF_alt  = scale_factor_array(pt_idx, eta_idx, sapt_idx, saeta_idx, charge_idx, pass_iso, with_trigger, sf_all_->template axis<4>().index(ns+1) ); // 0 is the nominal, systs starts from 1
                
                // anticorrelation between iso and antiiso already embedded in the numbers stored in the histograms
                // also the alternate comes from data efficiency variation only, so the anticorrelation in the efficiencies is preserved in the scale factors
            
                // order is reco-tracking-idip-trigger-iso
                for(int i = 0; i < allSF_nomi.size(); i++) {
                    // if (allSF_nomi[i] <= 0)
                    //     std::cout << "allSF_nomi/alt[" << i << "] = " << allSF_nomi[i] << "/" << allSF_alt[i] << " --> pt/eta/charge = " << pt << "/" << eta << "/" << charge << std::endl;
                    res(i, ns) = allSF_alt[i] / allSF_nomi[i]; 
                }

            }
            
            return res;
     
        }

    protected:

        std::shared_ptr<const HIST_SF> sf_all_;
        // cache the bin indices since the string category lookup is slow
        int idx_reco_ = sf_all_->template axis<3>().index("reco");
        int idx_tracking_ = sf_all_->template axis<3>().index("tracking");
        int idx_idip_ = sf_all_->template axis<3>().index("idip");
        int idx_trig_ = sf_all_->template axis<3>().index("trigger");
        int idx_iso_triggering_ = sf_all_->template axis<3>().index("iso");
        int idx_antiiso_triggering_ = sf_all_->template axis<3>().index("antiiso");
        int idx_iso_nontriggering_ = sf_all_->template axis<3>().index("isonotrig");
        int idx_antiiso_nontriggering_ = sf_all_->template axis<3>().index("antiisonotrig");

        int idx_nom_ = sf_all_->template axis<4>().index(0);
        // int idx_alt_ = sf_all_->template axis<4>().index(NSysts);
        
    };

    // base template for one-lepton case
    template<bool do_other, int NSysts, typename HIST_SF>
    class muon_efficiency_smooth_helper:
        public muon_efficiency_smooth_helper_base<NSysts, HIST_SF> {

    public:

        using base_t = muon_efficiency_smooth_helper_base<NSysts, HIST_SF>;
        // inherit constructor
        using base_t::base_t;

        muon_efficiency_smooth_helper(const base_t &other) : base_t(other) {}
        
        double operator() (float pt, float eta, float sapt, float saeta, int charge, bool pass_iso) {
            constexpr bool with_trigger = true;
            return base_t::scale_factor_product(pt, eta, sapt, saeta, charge, pass_iso, with_trigger, base_t::idx_nom_);
        }

    };

    // specialization for two-lepton case
    template<typename HIST_SF, int NSysts>
    class muon_efficiency_smooth_helper<true, NSysts, HIST_SF> :
        public muon_efficiency_smooth_helper_base<NSysts, HIST_SF> {

    public:

        using base_t = muon_efficiency_smooth_helper_base<NSysts, HIST_SF>;
        // inherit constructor
        using base_t::base_t;
        
        muon_efficiency_smooth_helper(const base_t &other) : base_t(other) {}

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
    template<bool do_other, int NSysts, typename HIST_SF>
    class muon_efficiency_smooth_helper_syst :
        public muon_efficiency_smooth_helper_base<NSysts, HIST_SF> {

    public:

        using base_t = muon_efficiency_smooth_helper_base<NSysts, HIST_SF>;
        using tensor_t = typename base_t::syst_tensor_t;
        
        // inherit constructor
        using base_t::base_t;
        
        muon_efficiency_smooth_helper_syst(const base_t &other) : base_t(other) {}
        
        tensor_t operator() (float pt, float eta, float sapt, float saeta, int charge, bool pass_iso, double nominal_weight = 1.0) {
            constexpr bool with_trigger = true;
            return nominal_weight*base_t::sf_syst_var(pt, eta, sapt, saeta, charge, pass_iso, with_trigger);
        }

    };

    // specialization for two-lepton case
    template<int NSysts, typename HIST_SF>
    class muon_efficiency_smooth_helper_syst<true, NSysts, HIST_SF> :
        public muon_efficiency_smooth_helper_base<NSysts, HIST_SF> {

    public:

        using base_t = muon_efficiency_smooth_helper_base<NSysts, HIST_SF>;
        using tensor_t = typename base_t::syst_tensor_t;

        // inherit constructor
        using base_t::base_t;
        
        muon_efficiency_smooth_helper_syst(const base_t &other) : base_t(other) {}

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
    template<int NEtaBins, int NPtEigenBins, int NCharges, typename HIST_SF>
    class muon_efficiency_smooth_helper_stat_base {
    public:

        muon_efficiency_smooth_helper_stat_base(HIST_SF &&sf_type) :
            sf_type_(std::make_shared<const HIST_SF>(std::move(sf_type))) {}

        // number of eta bins, number of eigen variations for pt axis, then 2 charges
        using stat_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<NEtaBins, NPtEigenBins, NCharges>>;

        //using base_t = muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins, HIST_SF>;        
        //muon_efficiency_smooth_helper_stat_base(const base_t &other) : base_t(other) {}

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
            auto const eigen_axis = sf_type_->template axis<4>();

            // overflow/underflow are attributed to adjacent bin
            auto const tensor_eta_idx = std::clamp(eta_idx, 0, NEtaBins - 1);
            //auto const tensor_pt_idx =  std::clamp(pt_idx, 0, NPtBins_ - 1);

            auto const &cell_nomi = sf_type_->at(eta_idx,
                                                 pt_idx,
                                                 charge_idx,
                                                 eff_type_idx,                                                      
                                                 idx_nom_);

            const double sf_nomi = cell_nomi.value();
            // loop on dimension with the parameters variations (the histogram already has the alternate SF) 
            // start from 1 because first bin contains the nominal SF
            for (int tensor_eigen_idx = 1; tensor_eigen_idx <= NPtEigenBins; tensor_eigen_idx++) {

                auto const eigen_axis_idx = eigen_axis.index(tensor_eigen_idx);
                
                auto const &cell_stat = sf_type_->at(eta_idx,
                                                     pt_idx,
                                                     charge_idx,
                                                     eff_type_idx,
                                                     eigen_axis_idx);

                const double sf_stat = cell_stat.value();
                const double sf_stat_variation = sf_stat / sf_nomi;
                    
                res(tensor_eta_idx, tensor_eigen_idx-1, charge_idx) *= sf_stat_variation;

            }
                
            return res;
        }

        // special case for isolation
        stat_tensor_t sf_stat_var_iso(float pt, float eta, int charge, bool pass_iso, bool with_trigger) const {
            stat_tensor_t res;
            res.setConstant(1.0);

            auto const eta_idx = sf_type_->template axis<0>().index(eta);
            auto const pt_idx = sf_type_->template axis<1>().index(pt);
            auto const charge_idx = sf_type_->template axis<2>().index(charge);
            auto const eigen_axis = sf_type_->template axis<4>();

            auto const eff_type_idx_iso_pass = with_trigger ? idx_iso_triggering_ : idx_iso_nontriggering_;
            auto const eff_type_idx_iso_fail = with_trigger ? idx_antiiso_triggering_ : idx_antiiso_nontriggering_;
            auto const eff_type_idx_iso = pass_iso ? eff_type_idx_iso_pass : eff_type_idx_iso_fail;

            // overflow/underflow are attributed to adjacent bin
            auto const tensor_eta_idx = std::clamp(eta_idx, 0, NEtaBins - 1);
            //auto const tensor_pt_idx =  std::clamp(pt_idx, 0, NPtBins_ - 1);

            auto const &cell_nomi = sf_type_->at(eta_idx,
                                                 pt_idx,
                                                 charge_idx,
                                                 eff_type_idx_iso,  
                                                 idx_nom_);
                
            const double sf_nomi = cell_nomi.value();
            // loop on dimension with the parameters variations (the histogram already has the alternate SF) 
            // start from 1 because first bin contains the nominal SF
            for (int tensor_eigen_idx = 1; tensor_eigen_idx <= NPtEigenBins; tensor_eigen_idx++) {

                auto const eigen_axis_idx = eigen_axis.index(tensor_eigen_idx);
                
                auto const &cell_stat = sf_type_->at(eta_idx,
                                                     pt_idx,
                                                     charge_idx,
                                                     eff_type_idx_iso,
                                                     eigen_axis_idx);

                const double sf_stat = cell_stat.value();
                const double sf_stat_variation = sf_stat / sf_nomi;
                res(tensor_eta_idx, tensor_eigen_idx-1, charge_idx) *= sf_stat_variation;

            }
                
            return res;
        }

    protected:

        std::shared_ptr<const HIST_SF> sf_type_;
        // int NPtBins_ = sf_type_->template axis<1>().size(); // there are no overflows, so extent == size 
        // cache the bin indices since the string category lookup is slow
        int idx_nom_ = sf_type_->template axis<4>().index(0); // input effStat axis is organized as DownVar - nomi - UpVar, with nomi centered at 0
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
    template<bool do_other, int NEtaBins, int NPtEigenBins, int NCharges, typename HIST_SF>
    class muon_efficiency_smooth_helper_stat :
        public muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins, NCharges, HIST_SF> {
        
    public:
        
        using stat_base_t = muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins, NCharges, HIST_SF>;
        using tensor_t = typename stat_base_t::stat_tensor_t;
  
        using stat_base_t::stat_base_t;
        
        tensor_t operator() (float pt, float eta, int charge, double nominal_weight = 1.0) {
            return nominal_weight*stat_base_t::sf_stat_var(pt, eta, charge);
        }

    };

    // specialization for two-lepton case
    template<int NEtaBins, int NPtEigenBins, int NCharges, typename HIST_SF>
    class muon_efficiency_smooth_helper_stat<true, NEtaBins, NPtEigenBins, NCharges, HIST_SF> :
        public muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins, NCharges, HIST_SF> {

    public:

        using stat_base_t = muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins, NCharges, HIST_SF>;
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
    template<bool do_other, int NEtaBins, int NPtEigenBins, int NCharges, typename HIST_SF>
    class muon_efficiency_smooth_helper_stat_iso :
        public muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins, NCharges, HIST_SF> {
        
    public:
        
        using stat_base_t = muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins, NCharges, HIST_SF>;
        using tensor_t = typename stat_base_t::stat_tensor_t;
  
        using stat_base_t::stat_base_t;
        
        tensor_t operator() (float pt, float eta, int charge, bool pass_iso, double nominal_weight = 1.0) {
            constexpr bool with_trigger = true;
            return nominal_weight*stat_base_t::sf_stat_var_iso(pt, eta, charge, pass_iso, with_trigger);
        }

    };

    // specialization for two-lepton case
    template<int NEtaBins, int NPtEigenBins, int NCharges, typename HIST_SF>
    class muon_efficiency_smooth_helper_stat_iso<true, NEtaBins, NPtEigenBins, NCharges, HIST_SF> :
        public muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins, NCharges, HIST_SF> {

    public:

        using stat_base_t = muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins, NCharges, HIST_SF>;
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
