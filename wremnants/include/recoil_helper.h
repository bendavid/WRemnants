#include <ROOT/RVec.hxx>
#include <iostream>
#include <vector>
#include <thread>
#include <chrono>
#include <cmath>
#include "defines.h"
#include "tfliteutils.h"


namespace wrem {

using ROOT::VecOps::RVec;


double print_val(double in, double out) {

    cout << in << " " << out << endl;
        return 0;
}


auto testFunc() {
    Eigen::TensorFixedSize<double, Eigen::Sizes<500>> tmp;
    tmp.setConstant(100.);

    auto ret = std::tuple(20, tmp);
    return ret;
}


class TestHelper {
public:
    using out_t = std::tuple<Eigen::TensorFixedSize<double, Eigen::Sizes<500>>, double>;
    TestHelper() {}

    out_t operator() () {

        out_t ret;
        std::get<1>(ret) = 20;
        std::get<0>(ret).setConstant(100.);
        return ret;
    }
};


class RecoilCalibrationHelperOld {

public:

    RecoilCalibrationHelperOld(const std::string &filename, const std::string &func_name, double ut_min_, double ut_max_) :
        helper_(std::make_shared<narf::tflite_helper>(filename, func_name, ROOT::GetThreadPoolSize())) {
            ut_max = ut_max_;
            ut_min = ut_min_;
    }

    double operator() (const double &pt, const double &ut) {

        if(ut < ut_max and ut > ut_min) {
            Eigen::TensorFixedSize<double, Eigen::Sizes<>> pt_tensor;
            Eigen::TensorFixedSize<double, Eigen::Sizes<>> ut_tensor;
            Eigen::TensorFixedSize<double, Eigen::Sizes<>> ut_corr_tensor;
            pt_tensor(0) = pt;
            ut_tensor(0) = ut;

            auto const inputs = std::tie(pt_tensor, ut_tensor);
            auto outputs = std::tie(ut_corr_tensor);
            (*helper_)(inputs, outputs);
            return ut_corr_tensor(0);
        }
        else {
            return ut;
        }
    }


private:
    std::shared_ptr<narf::tflite_helper> helper_;
    double ut_max;
    double ut_min;
};





class ResponseCorrector {

public:

    ResponseCorrector(const std::string &filename, const std::string &func_name) :
        helper_(std::make_shared<narf::tflite_helper>(filename, func_name, ROOT::GetThreadPoolSize())) {}

    double operator() (const double &pt) {

        Eigen::TensorFixedSize<double, Eigen::Sizes<>> pt_tensor;
        Eigen::TensorFixedSize<double, Eigen::Sizes<>> pt_corr_tensor;
        pt_tensor(0) = pt;

        auto const inputs = std::tie(pt_tensor);
        auto outputs = std::tie(pt_corr_tensor);
        (*helper_)(inputs, outputs);
        return pt_corr_tensor(0);
    }


private:
    std::shared_ptr<narf::tflite_helper> helper_;
    double ut_max;
    double ut_min;
};












template <size_t N_UNC>
class RecoilCalibrationHelper {

public:

    struct scalar_tensor {
        using out_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<N_UNC>>;
        using out_scalar_t = Eigen::TensorFixedSize<double, Eigen::Sizes<>>;

        out_scalar_t ut_corr;
        out_tensor_t unc_weights;
    };

    using out_t = scalar_tensor;

    RecoilCalibrationHelper(const std::string &filename, const std::string &func_name, const double ut_min_, const double ut_max_) :
        helper_(std::make_shared<narf::tflite_helper>(filename, func_name, ROOT::GetThreadPoolSize())) {
            ut_max = ut_max_;
            ut_min = ut_min_;
    }

    out_t operator() (const double &pt, const double &ut) {

        Eigen::TensorFixedSize<double, Eigen::Sizes<>> pt_tensor;
        Eigen::TensorFixedSize<double, Eigen::Sizes<>> ut_tensor;
        pt_tensor(0) = pt;
        ut_tensor(0) = ut;

        out_t ret;
        if(ut < ut_max and ut > ut_min) {
            auto const inputs = std::tie(pt_tensor, ut_tensor);
            auto outputs = std::tie(ret.ut_corr, ret.unc_weights);
            (*helper_)(inputs, outputs);
        }
        else {
            ret.ut_corr = ut_tensor;
            ret.unc_weights.setConstant(1.);
        }
        return ret;
    }


private:
    std::shared_ptr<narf::tflite_helper> helper_;
    double ut_max;
    double ut_min;
};



template <size_t N_UNC>
class RecoilCalibrationHelperNoUnc {

public:

    struct scalar_tensor {
        using out_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<N_UNC>>;
        using out_scalar_t = Eigen::TensorFixedSize<double, Eigen::Sizes<>>;

        out_scalar_t ut_corr;
        out_tensor_t unc_weights;
    };

    using out_t = scalar_tensor;

    RecoilCalibrationHelperNoUnc(const std::string &filename, const std::string &func_name, const double ut_min_, const double ut_max_) :
        helper_(std::make_shared<narf::tflite_helper>(filename, func_name, ROOT::GetThreadPoolSize())) {
            ut_max = ut_max_;
            ut_min = ut_min_;
    }

    out_t operator() (const double &pt, const double &ut) {

        Eigen::TensorFixedSize<double, Eigen::Sizes<>> pt_tensor;
        Eigen::TensorFixedSize<double, Eigen::Sizes<>> ut_tensor;
        pt_tensor(0) = pt;
        ut_tensor(0) = ut;

        out_t ret;
        ret.unc_weights.setConstant(1.);
        if(ut < ut_max and ut > ut_min) {
            auto const inputs = std::tie(pt_tensor, ut_tensor);
            auto outputs = std::tie(ret.ut_corr);
            (*helper_)(inputs, outputs);
        }
        else {
            ret.ut_corr = ut_tensor;
        }
        return ret;
    }


private:
    std::shared_ptr<narf::tflite_helper> helper_;
    double ut_max;
    double ut_min;
};


}


