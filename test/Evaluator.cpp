#include "Evaluator.h"

template <typename T>
void exportCSV( std::ostream& out, const std::string& headerName, const std::vector<T>& baseColumn, const std::vector<EvVector>& evectors)
{
    out << headerName;
    for (auto const& ev : evectors) out << ",\"" << ev.name<< "\"";
    out << "\n";

    size_t rows = baseColumn.size();
    for (auto const& ev : evectors) rows = std::max(rows, ev.data.size());

    for (size_t i = 0; i < rows; ++i) {
        if (i < baseColumn.size()) out << baseColumn[i];
        out << "";
        for (auto const& ev : evectors) {
            out << ",";
            if (i < ev.data.size()) out << ev.data[i];
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
    file << "$\\epsilon_r$," << std::scientific << std::setprecision(1) << conf.threshold_refine << std::endl;  // scientific for threshold
    file << "$\\theta_c$," << std::fixed << std::setprecision(3) << conf.theta_coarse << std::endl;
    file << "$\\epsilon_c$," << std::scientific << std::setprecision(1) << conf.threshold_coarse << std::endl;  // scientific for threshold

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
