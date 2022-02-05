#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

double scaleWeight(double weight, double scale) {
    return std::exp(scale*std::log(std::abs(weight))*std::copysign(1., weight));
}

template <size_t ETABINS, size_t T>
Eigen::TensorFixedSize<double, Eigen::Sizes<2, ETABINS>> dummyScaleFromMassWeights(
    Eigen::TensorFixedSize<double, Eigen::Sizes<T>>& weights, double eta, double scale, bool isW=true) {
    const double refMass = isW ? 80351.81229 : 91153.50974;
    const size_t centralIdx = 10;
    const double scaleMeV = refMass*scale;
    const int step10MeV = static_cast<int>(scaleMeV/10)+1;
    if (centralIdx+step10MeV > T-1)
        throw std::out_of_range("Maximum allowed range for momentum scale uncertainty is 100 MeV!");
    // Find weight (10 MeV steps) closest, but smaller than desired scale, then scale down
    const double scaleFac = scaleMeV/(10*step10MeV);

    Eigen::TensorFixedSize<double, Eigen::Sizes<2, ETABINS>> outWeights;
    outWeights.setConstant(1.);

    const double etaStep = 2*2.4/ETABINS;
    size_t ieta = (std::clamp(eta, -2.4, 2.4)+2.4)/etaStep;

    // Down weight, then up weight
    outWeights(0, ieta) = scaleWeight(weights[centralIdx-step10MeV], scaleFac);
    outWeights(1, ieta) = scaleWeight(weights[centralIdx+step10MeV], scaleFac);
    return outWeights;
}
