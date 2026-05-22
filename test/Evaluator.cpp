#include "Evaluator.h"

template <typename T>
void exportCSV(std::ostream &out, const std::string &headerName, const std::vector<T> &baseColumn, const std::vector<EvVector> &evectors) {
    out << headerName;
    for (auto const &ev : evectors)
        out << ",\"" << ev.name << "\"";
    out << "\n";

    size_t rows = baseColumn.size();
    for (auto const &ev : evectors)
        rows = std::max(rows, ev.data.size());

    for (size_t i = 0; i < rows; ++i) {
        if (i < baseColumn.size())
            out << baseColumn[i];
        out << "";
        for (auto const &ev : evectors) {
            out << ",";
            if (i < ev.data.size())
                out << ev.data[i];
        }
        out << "\n";
    }
}
void Evaluator::saveCSV_T(const std::string &filename) const {
    std::ofstream file(filename);
    exportCSV(file, "time", x, y);
}

void to_csv(const DoeflerConf &conf, const std::filesystem::path &filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filename.string());
    }

    // Write header
    file << ",Doerfler" << std::endl;

    // Refine row
    file << "$\\theta_r$," << std::fixed << std::setprecision(3) << conf.theta_refine << std::endl;
    file << "$\\epsilon_r$," << std::scientific << std::setprecision(1) << conf.threshold_refine << std::endl; // scientific for threshold
    file << "$\\theta_c$," << std::fixed << std::setprecision(3) << conf.theta_coarse << std::endl;
    file << "$\\epsilon_c$," << std::scientific << std::setprecision(1) << conf.threshold_coarse << std::endl; // scientific for threshold

    file.close();
}

// Function to generate CSV as string
void to_csv(const DOPRI5_conf &conf, const std::filesystem::path &filename) {
    std::ofstream ss(filename);
    if (!ss.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filename.string());
    }

    // Header for pgfplotstable
    ss << ",Dopri56\n";

    // Absolute tolerance
    ss << "$\\mathrm{Atol}$,"
       << std::scientific << std::setprecision(1) << conf.Atol_i << "\n";

    // Relative tolerance
    ss << "$\\mathrm{Rtol}$,"
       << std::scientific << std::setprecision(1) << conf.Rtol_i << "\n";

    // Max step factor
    ss << "$\\mathrm{f_{max}}$, "
       << std::fixed << std::setprecision(3) << conf.faxmax << "\n";

    // Min step factor
    ss << "$\\mathrm{f_{min}}$, "
       << std::fixed << std::setprecision(3) << conf.facmin << "\n";
    // Header
    ss.close();
}

void registerProperties(Evaluator &ev, ExportProperty p, int h_size) {
    if (p & EXPORT_DT)
        ev.reg("dt (s)", [](EvData d) { return d.dp5s.t_past; });

    if (p & EXPORT_WTST)
        ev.reg("wall time per simulation time (s)", [](EvData d) { return d.time_per_sim_sec; });

    if (p & EXPORT_RESIDUAL)
        ev.reg(R"($\eta_R$)", [](EvData d) { return d.poison_residual_error; });

    if (p & EXPORT_int_w)
        ev.reg(R"($\int_M w$)", [](EvData d) { return integral(d.wc.w, d.geom); });

    if (p & EXPORT_int_dwdt)
        ev.reg(R"($\int_M \frac{d}{dt} w$)", [](EvData d) { return integral(d.rhs.w, d.geom); });

    if (p & EXPORT_C) {
        for (int i = 0; i < h_size; ++i)
            ev.reg("$c_IDX(" + std::to_string(i) + ")$",
                   [i](EvData d) { return d.wc.c[i]; });
    }

    if (p & EXPORT_dCdW) {
        for (int i = 0; i < h_size; ++i)
            ev.reg(R"($\frac{d}{dt} c_IDX()" + std::to_string(i) + ")$",
                   [i](EvData d) { return d.rhs.c[i]; });
    }

    if (p & EXPORT_ATTEMPTS)
        ev.reg("attempts", [](EvData d) { return d.dp5s.attempts; });

    if (p & EXPORT_int_psi)
        ev.reg(R"($\int_M \psi$)", [](EvData d) { return integral(d.vel.stream_function, d.geom); });

    if (p & EXPORT_velocity)
        ev.reg(R"($\|u\|_2$)", [](EvData d) { return L2Norm(d.vel.u, d.geom); });

    if (p & EXPORT_nF)
        ev.reg(R"($\#F$)", [](EvData d) { return d.mesh.nFaces(); });

    if (p & EXPORT_nV)
        ev.reg(R"($\#V$)", [](EvData d) { return d.mesh.nVertices(); });

    if (p & EXPORT_DHC_HI) {
        for (int i = 0; i < h_size; ++i)
            ev.reg(R"($\| h_IDX()" + std::to_string(i) + R"() - \widehat{h}_IDX($)" + std::to_string(i) + R"()  \|_2)",
                   [i](EvData d) { return L2Norm(d.h[i] - d.h_interpol[i], d.geom); });
    }

    if (p & EXPORT_Vort_Y) {
        ev.reg("vorticity center",
               [](EvData d) {
                   Vertex high_v;
                   double max = std::numeric_limits<double>::lowest();
                   for (Vertex v : d.mesh.vertices()) {
                       if (d.wc.w[v] > max) {
                           high_v = v;
                           max = d.wc.w[v];
                       }
                   }
                   return d.vertexPosition[high_v].y;
               });
    }
}

ExportProperty defaultTimeProperties() {
    ExportProperty flags = (EXPORT_DT |
                            // EXPORT_ATTEMPTS |
                            EXPORT_velocity |
                            EXPORT_int_psi |
                            EXPORT_int_w |
                            EXPORT_int_dwdt |
                            EXPORT_WTST |
                            EXPORT_RESIDUAL |
                            EXPORT_C |
                            EXPORT_dCdW
                            // EXPORT_nF |
                            // EXPORT_nV |
                            // EXPORT_DHC_HI
    );
    return flags;
}
