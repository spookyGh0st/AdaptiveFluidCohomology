#include "cfd.h"
#include "util.h"

#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/utilities/vector2.h"
#include "poisson.h"
#include "refine.h"

namespace geometrycentral::surface {

wc_wrapper operator+(const wc_wrapper &a, const wc_wrapper &b) {
    assert(a.c.size() == b.c.size());

    wc_wrapper result;
    result.w = a.w + b.w;
    result.c.resize(a.c.size());
    for (size_t i = 0; i < a.c.size(); ++i)
        result.c[i] = a.c[i] + b.c[i];
    return result;
}

wc_wrapper operator*(double s, const wc_wrapper &wc) {
    wc_wrapper result;
    result.w = wc.w * s;
    result.c.resize(wc.c.size());
    for (size_t i = 0; i < wc.c.size(); ++i)
        result.c[i] = wc.c[i] * s;
    return result;
}

velocity_wrapper velocity(
    ManifoldSurfaceMesh &mesh,
    IntrinsicGeometryInterface &geom,
    const wc_wrapper &wc,
    const std::vector<FaceData<Vector2>> &h,
    const StreamFunctionSolver &S) {
    geom.requireHalfedgeVectorsInFace();
    VertexData<double> f(mesh, 0);
    S.solve(mesh, geom, f, wc.w);
    FaceData<Vector2> u(mesh, Vector2::zero());
    for (Face face : mesh.faces()) {
        u[face] = -grad(geom, face, f).rotate90();
    }

    assert(wc.c.size() == h.size());
    for (std::size_t i = 0; i < h.size(); i++)
        for (Face face : mesh.faces())
            u[face] += wc.c[i] * h[i][face];

    geom.unrequireHalfedgeVectorsInFace();
    return {f, u };
}

wc_wrapper evalRHS(
    ManifoldSurfaceMesh &mesh,
    IntrinsicGeometryInterface &geom,
    const wc_wrapper &wc,
    const std::vector<FaceData<Vector2>> &h,
    const StreamFunctionSolver &S,
    std::vector<FaceData<double>>* face_dc
    )
{
    geom.requireFaceAreas();
    geom.requireHalfedgeVectorsInFace();
    auto u = velocity(mesh, geom, wc, h, S);

    auto l = FaceData<Vector2>(mesh, Vector2::zero());
    for (Face f : mesh.faces()) {
        l[f] = Lamb(f, wc.w, u.u);
    }

    VertexData<double> dw(mesh, 0);
    for (Vertex v : mesh.vertices()) {
        dw[v] = -derive(v, u.u, geom, wc.w);
    }
    std::vector<double> dc(wc.c.size(), 0);
    for (std::size_t i = 0; i < wc.c.size(); i++){
        for (Face f : mesh.faces()){
            auto dcf = dot(l[f], h[i][f]) * geom.faceAreas[f];
            if(face_dc) face_dc->at(i)[f] = dcf;
            dc[i] += dcf;
        }
    }

    geom.unrequireFaceAreas();
    geom.unrequireHalfedgeVectorsInFace();
    return {dw, dc};
}

wc_wrapper RK4Step(
    ManifoldSurfaceMesh &mesh,
    IntrinsicGeometryInterface &geom,
    const std::vector<FaceData<Vector2>> &h,
    const wc_wrapper &x,
    double dt,
    const StreamFunctionSolver &S) {
    auto F = [&mesh, &geom, &h, &S](const wc_wrapper &wc) -> wc_wrapper { return evalRHS(mesh, geom, wc, h, S); };
    wc_wrapper k1 = F(x), k2 = F(x + (dt / 2) * k1);
    wc_wrapper k3 = F(x + (dt / 2) * k2), k4 = F(x + dt * k3);
    return x + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
}

using F_type = std::function<wc_wrapper(const wc_wrapper &)>;

std::array<wc_wrapper, 7> RK6_k(const std::array<std::array<double, 7>, 7> &a, const wc_wrapper &y0, double h, const F_type &F) {
    constexpr std::size_t n = 7;
    std::array<wc_wrapper, n> k{};
    for (int i = 0; i < n; ++i) {
        wc_wrapper x = y0;
        for (int j = 0; j < i; ++j) {
            if (a[i][j] != 0)
                x = x + h * a[i][j] * k[j];
        }
        k[i] = F(x);
    }
    return k;
}
wc_wrapper RK6_2(const std::array<double, 7> &b, const wc_wrapper &y0, double h, const std::array<wc_wrapper, 7> &k) {
    constexpr std::size_t n = 7;
    wc_wrapper y1 = y0;
    for (int i = 0; i < n; ++i) {
        if (b[i] != 0)
            y1 = y1 + (h * b[i]) * k[i];
    }
    return y1;
}

double error(ManifoldSurfaceMesh &mesh, const wc_wrapper &y0, wc_wrapper &y1, wc_wrapper &y_hat, const double &Atol_i, const double &Rtol_i) {
    double sum = 0;
    for (Vertex v : mesh.vertices()) {
        double sc_i = Atol_i + std::max(std::abs(y0.w[v]), std::abs(y1.w[v])) * Rtol_i;
        assert(sc_i > 0);
        sum += std::pow((y1.w[v] - y_hat.w[v]) / sc_i, 2);
    }
    for (std::size_t i = 0; i < y0.c.size(); i++) {
        double sc_i = Atol_i + std::max(std::abs(y0.c[i]), std::abs(y1.c[i])) * Rtol_i;
        assert(sc_i > 0);
        sum += std::pow((y1.c[i] - y_hat.c[i]) / sc_i, 2);
    }
    double n = (mesh.nVertices() + y0.c.size());
    return std::sqrt(sum / n);
}

inline double compute_h_new(double h, double err, double q, double facmin, double facmax) {
    double h_o = std::pow(1. / err, 1. / (q + 1));
    double fac = std::pow(0.38, 1. / (q + 1));
    return h * std::min(facmax, std::max(facmin, fac * h_o));
}

inline bool accept_step(double err) { return err < 1; }

std::array<wc_wrapper, 2> DOPRI5(const wc_wrapper &y0, double h, const F_type &F) {
    constexpr std::size_t n = 7;
    std::array<std::array<double, n>, n> a = {{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                               {1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                               {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                               {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0, 0.0},
                                               {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0, 0.0},
                                               {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0, 0.0},
                                               {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0}}};
    std::array<double, 7> b = {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0};
    std::array<double, 7> b_hat = {5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0};
    std::array<wc_wrapper, 7> k = RK6_k(a, y0, h, F);

    return {RK6_2(b, y0, h, k), RK6_2(b_hat, y0, h, k)};
}

DOPRI5_sample DOPRI5_step(
    ManifoldSurfaceMesh &mesh,
    const wc_wrapper &y0,
    double h,
    const F_type &F,
    const DOPRI5_conf &conf) {
    double h_past;
    DOPRI5_sample sample;
    bool accepted = false;
    double facMax = conf.faxmax;
    std::array<wc_wrapper, 2> y;
    int attempts = -1;
    while (!accepted) {
        attempts++;
        constexpr double q = 5; // min of order y and y_hat
        y = DOPRI5(y0, h, F);
        h_past = h; // store time step
        double err = error(mesh, y0, y[0], y[1], conf.Atol_i, conf.Rtol_i);
        accepted = accept_step(err);
        h = compute_h_new(h, err, q, conf.facmin, facMax);
        facMax = 1; // addvised to override after step-rejection
    }
    return {y[0], h_past, h, attempts};
}

DOPRI5_sample adaptive_step(
    ManifoldSurfaceMesh &mesh,
    IntrinsicGeometryInterface &geom,
    const std::vector<FaceData<Vector2>> &h,
    const wc_wrapper &x,
    double dt,
    const StreamFunctionSolver &S,
    const DOPRI5_conf &conf) {
    F_type F = [&mesh, &geom, &h, &S](const wc_wrapper &wc) -> wc_wrapper { return evalRHS(mesh, geom, wc, h, S); };
    return DOPRI5_step(mesh, x, dt, F, conf);
}
DOPRI5_conf DOPRI5Preset(DOPRI5PresetConf preset) {
    DOPRI5_conf conf;
    switch (preset) {
    case DOPRI5PresetConf::LOW: // Low precision (fast)
        conf.Rtol_i = 1e-3;
        conf.Atol_i = 1e-6;
        conf.facmin = 0.2;
        conf.faxmax = 10.0;
        break;
    case DOPRI5PresetConf::MEDIUM:
        // Medium precision (balanced)
        conf.Rtol_i = 1e-6;
        conf.Atol_i = 1e-9;
        conf.facmin = 0.2;
        conf.faxmax = 5.0;
        break;

    case DOPRI5PresetConf::HIGH:
        conf.Rtol_i = 1e-8;
        conf.Atol_i = 1e-11;
        conf.facmin = 0.25;
        conf.faxmax = 2.0;
        break;
    case DOPRI5PresetConf::VERY_HIGH:
        conf.Rtol_i = 1e-10;
        conf.Atol_i = 1e-13;
        conf.facmin = 0.25;
        conf.faxmax = 1.5;
        break;
    }
    return conf;
}

} // namespace geometrycentral::surface
