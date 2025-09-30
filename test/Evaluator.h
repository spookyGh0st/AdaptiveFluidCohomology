#pragma once

#include <eider/AdaptiveFluidSolver.h>
#include "eider/cfd.h"
#include "eider/homotopy.h"
#include "eider/util.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/meshio.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <memory>
#include <string>
#include <vector>

using namespace geometrycentral::surface;


struct EvData {
    ManifoldSurfaceMesh& mesh;
    IntrinsicGeometryInterface& geom;
    velocity_wrapper vel;
    wc_wrapper rhs;
    wc_wrapper wc;
    Harmonic_basis h;
    DOPRI5_sample dp5s;
    std::vector<FaceData<double>> dc;
    double time_per_sim_sec;
    double poison_residual_error;
};

struct EvVector {
    std::string name;
    std::vector<float> data;
    std::function<double(EvData)> f;
    void evaluate(const EvData& ev){
        data.push_back(float(f(ev)));
    }
    EvVector(const std::string& name, const std::function<double(EvData)>& f) :name(name), f(f) {}
};


struct Evaluator {
    std::vector<EvVector> y;
    std::vector<float> x;
    void reg(const std::string& name, const std::function<float(EvData)>& f){
        y.emplace_back(name,f);
    }

    void onStep(const EvData& data, float x_value){
        x.push_back(x_value);
        for (EvVector& v: y) {
            v.evaluate(data);
        }
    }
    void saveCSV_T(const std::string &filename) const;

    std::map<std::string, std::vector<size_t>> group_columns() {
        std::map<std::string, std::vector<size_t>> groups;

        for (size_t i = 0; i < y.size(); ++i) {
            const std::string& name = y[i].name;

            // Find numeric suffix by scanning backwards
            size_t pos = name.size();
            while (pos > 0 && isdigit(name[pos - 1])) --pos;

            std::string prefix = (pos == 0 || pos == name.size()) ? name : name.substr(0, pos);

            groups[prefix].push_back(i);
        }

        return groups;
    }
};

void to_csv(const DoeflerConf& conf, const std::filesystem::path& filename);
void to_csv(const DOPRI5_conf& conf, const std::filesystem::path& filename);



