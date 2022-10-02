#ifndef WREMNANTS_MUON_EFFICIENCIES_H
#define WREMNANTS_MUON_EFFICIENCIES_H

namespace wrem {

    // TODO use enums for integer/boolean/category axes so that the code is less error-prone?

    template<typename HIST_IDIPTRIGISO, typename HIST_TRACKING, typename HIST_RECO>
    class muon_efficiency_helper_base {
    public:

        muon_efficiency_helper_base(HIST_IDIPTRIGISO &&sf_idip_trig_iso, HIST_TRACKING &&sf_tracking, HIST_RECO &&sf_reco) :
            sf_idip_trig_iso_(std::make_shared<const HIST_IDIPTRIGISO>(std::move(sf_idip_trig_iso))),
            sf_tracking_(std::make_shared<const HIST_TRACKING>(std::move(sf_tracking))),
            sf_reco_(std::make_shared<const HIST_RECO>(std::move(sf_reco))) {}

        // FIXME: why not using a double directly as in scale_factor_product?
        auto const &sf_tracking(float saeta, float sapt, int charge, int idx_nom_alt) const {
            //TODO index for nom_alt are in principle known at compile time

            ////
            // CENTRAL VALUES
            ////

            return sf_tracking_->at(sf_tracking_->template axis<0>().index(saeta),
                                    sf_tracking_->template axis<1>().index(sapt),
                                    sf_tracking_->template axis<2>().index(charge),
                                    idx_nom_alt);
        }
    
        double scale_factor_product(float pt, float eta, float sapt, float saeta, int charge, bool pass_iso, bool with_trigger, int idx_nom_alt) const {

            auto const eta_idx = sf_idip_trig_iso_->template axis<0>().index(eta);
            auto const pt_idx = sf_idip_trig_iso_->template axis<1>().index(pt);
            auto const pt_idx_reco = sf_reco_->template axis<1>().index(pt);
            auto const charge_idx = sf_idip_trig_iso_->template axis<2>().index(charge);
            auto const eff_type_idx_idip_trig = with_trigger ? idx_idip_trig_ : idip_idx_;
            auto const eff_type_idx_iso_pass = with_trigger ? idx_iso_triggering_ : idx_iso_nontriggering_;
            auto const eff_type_idx_iso_fail = with_trigger ? idx_antiiso_triggering_ : idx_antiiso_nontriggering_;
            auto const eff_type_idx_iso = pass_iso ? eff_type_idx_iso_pass : eff_type_idx_iso_fail;

            const double idip = sf_idip_trig_iso_->at(eta_idx, pt_idx, charge_idx, eff_type_idx_idip_trig, idx_nom_alt).value();
            const double iso = sf_idip_trig_iso_->at(eta_idx, pt_idx, charge_idx, eff_type_idx_iso, idx_nom_alt).value();
            const double reco = sf_reco_->at(eta_idx, pt_idx_reco, charge_idx, idx_nom_alt).value();
            const double tracking = sf_tracking(saeta, sapt, charge, idx_nom_alt).value();

            return idip*iso*tracking*reco;
        }

    protected:

        std::shared_ptr<const HIST_IDIPTRIGISO> sf_idip_trig_iso_;
        std::shared_ptr<const HIST_TRACKING> sf_tracking_;
        std::shared_ptr<const HIST_RECO> sf_reco_;

        // cache the bin indices since the string category lookup is slow
        int idip_idx_ = sf_idip_trig_iso_->template axis<3>().index("idip");
        int idx_idip_trig_ = sf_idip_trig_iso_->template axis<3>().index("idip_trig");
        int idx_iso_triggering_ = sf_idip_trig_iso_->template axis<3>().index("iso_triggering");
        int idx_antiiso_triggering_ = sf_idip_trig_iso_->template axis<3>().index("antiiso_triggering");
        int idx_iso_nontriggering_ = sf_idip_trig_iso_->template axis<3>().index("iso_nontriggering");
        int idx_antiiso_nontriggering_ = sf_idip_trig_iso_->template axis<3>().index("antiiso_nontriggering");

        int idx_nom_ = sf_idip_trig_iso_->template axis<4>().index(0);
        int idx_alt_ = sf_idip_trig_iso_->template axis<4>().index(1);

    };

    // base template for one-lepton case
    template<bool do_other, typename HIST_IDIPTRIGISO, typename HIST_TRACKING, typename HIST_RECO>
    class muon_efficiency_helper:
        public muon_efficiency_helper_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> {

    public:

        using base_t = muon_efficiency_helper_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO>;

        // inherit constructor
        using base_t::base_t;

        double operator() (float pt, float eta, float sapt, float saeta, int charge, bool pass_iso) {
            constexpr bool with_trigger = true;
            return base_t::scale_factor_product(pt, eta, sapt, saeta, charge, pass_iso, with_trigger, base_t::idx_nom_);
        }

    };

    // specialization for two-lepton case
    template<typename HIST_IDIPTRIGISO, typename HIST_TRACKING, typename HIST_RECO>
    class muon_efficiency_helper<true, HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> :
        public muon_efficiency_helper_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> {

    public:

        using base_t = muon_efficiency_helper_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO>;

        // inherit constructor
        using base_t::base_t;

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

    ////
    // SYST FOR EVERYTHING
    ////
    //// BASE CLASS FOR HELPER_SYST
    template<typename HIST_IDIPTRIGISO, typename HIST_TRACKING, typename HIST_RECO>
    class muon_efficiency_helper_syst_base:
        public muon_efficiency_helper_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> {
    public:

        using base_t = muon_efficiency_helper_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO>;
        // inherit constructor
        using base_t::base_t;

        muon_efficiency_helper_syst_base(const base_t &other) : base_t(other) {}

        using syst_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<4>>; // 4 bins for idip(*trigger), iso(notrig), reco, tracking, in this order

        syst_tensor_t sf_syst_var(float pt, float eta, float sapt, float saeta, int charge, bool pass_iso, bool with_trigger) const {
            syst_tensor_t res;

            auto const eta_idx =     base_t::sf_idip_trig_iso_->template axis<0>().index(eta);
            auto const pt_idx =      base_t::sf_idip_trig_iso_->template axis<1>().index(pt);
            auto const pt_idx_reco = base_t::sf_reco_->template axis<1>().index(pt);
            auto const saeta_idx =   base_t::sf_tracking_->template axis<0>().index(saeta);
            auto const sapt_idx =    base_t::sf_tracking_->template axis<0>().index(sapt);
            auto const charge_idx =  base_t::sf_idip_trig_iso_->template axis<2>().index(charge);
            auto const eff_type_idx_idip_trig = with_trigger ? base_t::idx_idip_trig_          : base_t::idip_idx_;
            auto const eff_type_idx_iso_pass =  with_trigger ? base_t::idx_iso_triggering_     : base_t::idx_iso_nontriggering_;
            auto const eff_type_idx_iso_fail =  with_trigger ? base_t::idx_antiiso_triggering_ : base_t::idx_antiiso_nontriggering_;
            auto const eff_type_idx_iso = pass_iso ? eff_type_idx_iso_pass : eff_type_idx_iso_fail;

            const double sf_reco = base_t::sf_reco_->at(eta_idx,
                                                        pt_idx_reco,
                                                        charge_idx,
                                                        base_t::idx_nom_).value();
            
            const double sf_reco_alt = base_t::sf_reco_->at(eta_idx,
                                                            pt_idx_reco,
                                                            charge_idx,
                                                            base_t::idx_alt_).value();
            
            const double variation_factor_reco = sf_reco_alt/sf_reco;
            res(0) = variation_factor_reco;

            const double sf_tracking = base_t::sf_tracking_->at(saeta_idx,
                                                                sapt_idx,
                                                                charge_idx,
                                                                base_t::idx_nom_).value();
    
            const double sf_tracking_alt = base_t::sf_tracking_->at(saeta_idx,
                                                                    sapt_idx,
                                                                    charge_idx,
                                                                    base_t::idx_alt_).value();

            const double variation_factor_tracking = sf_tracking_alt/sf_tracking;
            res(1) = variation_factor_tracking;

            const double sf_idip_trig = base_t::sf_idip_trig_iso_->at(eta_idx,
                                                                      pt_idx,
                                                                      charge_idx,
                                                                      eff_type_idx_idip_trig,
                                                                      base_t::idx_nom_).value();

            const double sf_idip_trig_alt = base_t::sf_idip_trig_iso_->at(eta_idx,
                                                                          pt_idx,
                                                                          charge_idx,
                                                                          eff_type_idx_idip_trig,
                                                                          base_t::idx_alt_).value();

            const double variation_factor_idip_trig = sf_idip_trig_alt/sf_idip_trig;
            res(2) = variation_factor_idip_trig;

            const double sf_iso = base_t::sf_idip_trig_iso_->at(eta_idx,
                                                                pt_idx,
                                                                charge_idx,
                                                                eff_type_idx_iso,
                                                                base_t::idx_nom_).value();

            const double sf_iso_alt = base_t::sf_idip_trig_iso_->at(eta_idx,
                                                                    pt_idx,
                                                                    charge_idx,
                                                                    eff_type_idx_iso,
                                                                    base_t::idx_alt_).value();

            // anti-correlation between iso and anti-iso SF's is not exact, but an excellent
            // approximation
            const double variation_factor_iso = pass_iso ? sf_iso_alt/sf_iso : (sf_iso - (sf_iso_alt - sf_iso))/sf_iso;
            res(3) = variation_factor_iso;

            return res;

        }

    };
        
    // base template for one lepton case
    template<bool do_other, typename HIST_IDIPTRIGISO, typename HIST_TRACKING, typename HIST_RECO>
    class muon_efficiency_helper_syst :
        public muon_efficiency_helper_syst_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> {

    public:

        using base_t = muon_efficiency_helper_syst_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO>;
        using tensor_t = typename base_t::syst_tensor_t;

        muon_efficiency_helper_syst(const base_t &other) : base_t(other) {}

        tensor_t operator() (float pt, float eta, float sapt, float saeta, int charge, bool pass_iso, double nominal_weight = 1.0) {
            constexpr bool with_trigger = true;
            return nominal_weight*base_t::sf_syst_var(pt, eta, sapt, saeta, charge, pass_iso, with_trigger);
        }

    };

    // specialization for two-lepton case
    template<typename HIST_IDIPTRIGISO, typename HIST_TRACKING, typename HIST_RECO>
    class muon_efficiency_helper_syst<true, HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> :
        public muon_efficiency_helper_syst_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> {

    public:

        using base_t = muon_efficiency_helper_syst_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO>;
        using tensor_t = typename base_t::syst_tensor_t;

        muon_efficiency_helper_syst(const base_t &other) : base_t(other) {}

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

    ////
    // STAT UNCERTAINTY
    ////
    //// BASE CLASS FOR HELPER_STAT
    template<int NEtaBins, int NPtBins, typename HIST_IDIPTRIGISO, typename HIST_TRACKING, typename HIST_RECO>
    class muon_efficiency_helper_stat_base:
        public muon_efficiency_helper_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> {
    public:

        using base_t = muon_efficiency_helper_base<HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO>;
        // inherit constructor
        using base_t::base_t;

        muon_efficiency_helper_stat_base(const base_t &other) : base_t(other) {}
        
        using stat_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<NEtaBins, NPtBins, 2, 2>>;

        stat_tensor_t sf_idip_trig_iso_stat_var(float pt, float eta, int charge, bool pass_iso, bool with_trigger) const {
            stat_tensor_t res;
            res.setConstant(1.0);

            auto const eta_idx = base_t::sf_idip_trig_iso_->template axis<0>().index(eta);
            auto const pt_idx = base_t::sf_idip_trig_iso_->template axis<1>().index(pt);
            auto const charge_idx = base_t::sf_idip_trig_iso_->template axis<2>().index(charge);
            auto const eff_type_idx_idip_trig = with_trigger ? base_t::idx_idip_trig_ : base_t::idip_idx_;
            auto const eff_type_idx_iso_pass = with_trigger ? base_t::idx_iso_triggering_ : base_t::idx_iso_nontriggering_;
            auto const eff_type_idx_iso_fail = with_trigger ? base_t::idx_antiiso_triggering_ : base_t::idx_antiiso_nontriggering_;
            auto const eff_type_idx_iso = pass_iso ? eff_type_idx_iso_pass : eff_type_idx_iso_fail;

            // overflow/underflow are attributed to adjacent bin
            auto const tensor_eta_idx = std::clamp(eta_idx, 0, NEtaBins - 1);
            auto const tensor_pt_idx = std::clamp(pt_idx, 0, NPtBins - 1);

            auto const &cell_idip_trig = base_t::sf_idip_trig_iso_->at(eta_idx,
                                                                       pt_idx,
                                                                       charge_idx,
                                                                       eff_type_idx_idip_trig,
                                                                       base_t::idx_nom_);

            const double sf_idip_trig = cell_idip_trig.value();
            const double err_idip_trig = std::sqrt(cell_idip_trig.variance());

            const double variation_factor_idip_trig = (sf_idip_trig + err_idip_trig)/sf_idip_trig;

            res(tensor_eta_idx, tensor_pt_idx, charge_idx, 0) *= variation_factor_idip_trig;

            auto const &cell_iso = base_t::sf_idip_trig_iso_->at(eta_idx,
                                                                 pt_idx,
                                                                 charge_idx,
                                                                 eff_type_idx_iso,
                                                                 base_t::idx_nom_);

            const double sf_iso = cell_iso.value();
            const double err_iso = std::sqrt(cell_iso.variance());

            // anti-correlation between iso and anti-iso SF's is not exact, but an excellent approximation
            const double variation_factor_iso = pass_iso ? (sf_iso + err_iso)/sf_iso : (sf_iso - err_iso)/sf_iso;

            res(tensor_eta_idx, tensor_pt_idx, charge_idx, 1) *= variation_factor_iso;

            return res;
        }

        using stat_tensor_singleStep_t = Eigen::TensorFixedSize<double, Eigen::Sizes<NEtaBins, NPtBins, 2>>;

        //TODO avoid duplication of code above
        
        stat_tensor_singleStep_t sf_singleStep_stat_var(float pt, float eta, int charge, int step) const {
            stat_tensor_singleStep_t res;
            res.setConstant(1.0);

            auto const sf_singleStep_ = (step == 0) ? base_t::sf_reco_ : base_t::sf_tracking_; // might be more general to read other stuff too
            
            auto const eta_idx = sf_singleStep_->template axis<0>().index(eta);
            auto const pt_idx = sf_singleStep_->template axis<1>().index(pt);
            auto const charge_idx = sf_singleStep_->template axis<2>().index(charge);

            // overflow/underflow are attributed to adjacent bin
            auto const tensor_eta_idx = std::clamp(eta_idx, 0, NEtaBins - 1);
            auto const tensor_pt_idx = std::clamp(pt_idx, 0, NPtBins - 1);

            auto const &cell_singleStep = sf_singleStep_->at(eta_idx,
                                                             pt_idx,
                                                             charge_idx,
                                                             base_t::idx_nom_);

            const double sf_singleStep = cell_singleStep.value();
            const double err_singleStep = std::sqrt(cell_singleStep.variance());

            const double variation_factor_singleStep = (sf_singleStep + err_singleStep)/sf_singleStep;

            res(tensor_eta_idx, tensor_pt_idx, charge_idx) *= variation_factor_singleStep;
            return res;
        }

    };


    ////
    // STAT FOR EVERYTHING EXCEPT RECO AND TRACKING
    ////
    
    // base template for one lepton case
    template<bool do_other, int NEtaBins, int NPtBins, typename HIST_IDIPTRIGISO, typename HIST_TRACKING, typename HIST_RECO>
    class muon_efficiency_helper_stat :
        public muon_efficiency_helper_stat_base<NEtaBins, NPtBins, HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> {
        
    public:
        
        using base_t = muon_efficiency_helper_stat_base<NEtaBins, NPtBins, HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO>;
  
        muon_efficiency_helper_stat(const base_t &other) : base_t(other) {}

        using tensor_t = typename base_t::stat_tensor_t;
        
        tensor_t operator() (float pt, float eta, int charge, bool pass_iso, double nominal_weight = 1.0) {
            constexpr bool with_trigger = true;
            return nominal_weight*base_t::sf_idip_trig_iso_stat_var(pt, eta, charge, pass_iso, with_trigger);
        }

    };

    // specialization for two-lepton case
    template<int NEtaBins, int NPtBins, typename HIST_IDIPTRIGISO, typename HIST_TRACKING, typename HIST_RECO>
    class muon_efficiency_helper_stat<true, NEtaBins, NPtBins, HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> :
        public muon_efficiency_helper_stat_base<NEtaBins, NPtBins, HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> {

    public:

        using base_t = muon_efficiency_helper_stat_base<NEtaBins, NPtBins, HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO>;
        using tensor_t = typename base_t::stat_tensor_t;

        muon_efficiency_helper_stat(const base_t &other) : base_t(other) {}

        tensor_t operator() (float trig_pt, float trig_eta, int trig_charge,
                             float nontrig_pt, float nontrig_eta, int nontrig_charge, double nominal_weight = 1.0) {
            constexpr bool with_trigger = true;
            constexpr bool without_trigger = false;
            constexpr bool pass_iso = true;

            const tensor_t variation_trig = base_t::sf_idip_trig_iso_stat_var(trig_pt, trig_eta, trig_charge, pass_iso, with_trigger);

            const tensor_t variation_nontrig = base_t::sf_idip_trig_iso_stat_var(nontrig_pt, nontrig_eta, nontrig_charge, pass_iso, without_trigger);

            return nominal_weight*variation_trig*variation_nontrig;
        }

    };

    ////
    // STAT FOR RECO AND TRACKING
    ////
    
    // base template for one lepton case
    template<bool do_other, int NEtaBins, int NPtBins, typename HIST_IDIPTRIGISO, typename HIST_TRACKING, typename HIST_RECO>
    class muon_efficiency_helper_singleStep_stat :
        public muon_efficiency_helper_stat_base<NEtaBins, NPtBins, HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> {

    public:

        using base_t = muon_efficiency_helper_stat_base<NEtaBins, NPtBins, HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO>;
        using tensor_t = typename base_t::stat_tensor_singleStep_t;

        muon_efficiency_helper_singleStep_stat(const base_t &other) : base_t(other) {}

        tensor_t operator() (float pt, float eta, int charge, int step, double nominal_weight = 1.0) {
            return nominal_weight*base_t::sf_singleStep_stat_var(pt, eta, charge, step);
        }

    };

    // specialization for two-lepton case
    template<int NEtaBins, int NPtBins, typename HIST_IDIPTRIGISO, typename HIST_TRACKING, typename HIST_RECO>
    class muon_efficiency_helper_singleStep_stat<true, NEtaBins, NPtBins, HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> :
        public muon_efficiency_helper_stat_base<NEtaBins, NPtBins, HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO> {

    public:

        using base_t = muon_efficiency_helper_stat_base<NEtaBins, NPtBins, HIST_IDIPTRIGISO, HIST_TRACKING, HIST_RECO>;
        using tensor_t = typename base_t::stat_tensor_singleStep_t;

        muon_efficiency_helper_singleStep_stat(const base_t &other) : base_t(other) {}

        tensor_t operator() (float trig_pt, float trig_eta, int trig_charge,
                             float nontrig_pt, float nontrig_eta, int nontrig_charge,
                             int step, double nominal_weight = 1.0) {

            const tensor_t variation_trig = base_t::sf_singleStep_stat_var(trig_pt, trig_eta, trig_charge, step);

            const tensor_t variation_nontrig = base_t::sf_singleStep_stat_var(nontrig_pt, nontrig_eta, nontrig_charge, step);

            return nominal_weight*variation_trig*variation_nontrig;
        }

    };

}

#endif
