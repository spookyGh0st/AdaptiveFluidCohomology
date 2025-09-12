#include "Evaluator.h"

template <typename T>
void exportCSV( std::ostream& out, const std::string& headerName, const std::vector<T>& baseColumn, const std::vector<EvVector>& evectors)
{
    out << headerName;
    for (auto const& ev : evectors) out << "," << ev.name;
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

