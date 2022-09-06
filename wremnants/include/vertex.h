#ifndef WREMNANTS_VERTEX_H
#define WREMNANTS_VERTEX_H

#include "TH2D.h"
#include <memory>

namespace wrem {

class vertex_helper {
public:

  vertex_helper(const TH2D &weights) :
    vertexweights_(make_shared_TH1<const TH2D>(weights)),
        nBinsX_(make_shared<const int>(weights.GetNbinsX())),
        nBinsY_(make_shared<const int>(weights.GetNbinsY())) {
        // store these to avoid call to TH2 methods everytime
        // nBinsX_ = vertexweights_.GetNbinsX();
        // nBinsY_ = vertexweights_.GetNbinsY();
    }

    // returns the vertex weight
    double operator() (float genVtx_z, float nTrueInt) const {
        int xbin = std::clamp(vertexweights_->GetXaxis()->FindFixBin(genVtx_z), 1, *nBinsX_);
        int ybin = std::clamp(vertexweights_->GetYaxis()->FindFixBin(nTrueInt), 1, *nBinsY_);
        return vertexweights_->GetBinContent(xbin, ybin);
    }


private:
    std::shared_ptr<const TH2D> vertexweights_;
    std::shared_ptr<const int> nBinsX_, nBinsY_;
};

}

#endif
