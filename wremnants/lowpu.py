
import narf

# load lowPU specific libs
#ROOT.gInterpreter.AddIncludePath(f"{pathlib.Path(__file__).parent}/include/")
narf.clingutils.Declare('#include "lowpu_utils.h"')
narf.clingutils.Declare('#include "lowpu_efficiencies.h"')
narf.clingutils.Declare('#include "lowpu_prefire.h"')
narf.clingutils.Declare('#include "lowpu_rochester.h"')
narf.clingutils.Declare('#include "electron_selections.h"')


def lepSF_systs(df, results, sName, sVars, defineExpr, baseName, baseAxes, baseCols):

    if sName not in df.GetColumnNames():
        df = df.Define(sName, defineExpr)
        df = df.Define(f"{sName}_tensor", f"Eigen::TensorFixedSize<double, Eigen::Sizes<{sVars}>> res; auto w = nominal_weight*{sName}; std::copy(std::begin(w), std::end(w), res.data()); return res;")
    results.append(df.HistoBoost(f"{baseName}_{sName}", [*baseAxes], [*baseCols, f"{sName}_tensor"]))
    return df