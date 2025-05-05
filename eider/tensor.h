#ifndef SL_TENSOR_H
#define SL_TENSOR_H


#include <array>
#include <cstddef>
#include <stdexcept>
#include <functional>
#include <span>

//---------------------------------------------------------------------
// Helper Functions
//---------------------------------------------------------------------

// Compute integer power: returns base^exp.
constexpr std::size_t ipow(std::size_t base, std::size_t exp) {
    return exp == 0 ? 1 : base * ipow(base, exp - 1);
}

// Compute binomial coefficient.
constexpr std::size_t binom(std::size_t n, std::size_t k) {
    return k == 0 ? 1 : (n * binom(n - 1, k - 1)) / k;
}

// Given a multi–index (stored in an array) with 'count' entries (each in [0, n)),
// compute the corresponding flat index (assuming lexicographic order).
template<std::size_t count>
constexpr std::size_t multi_index_to_flat(const std::array<std::size_t, count>& indices, std::size_t n) {
    std::size_t flat = 0;
    for (std::size_t i = 0; i < count; ++i) {
        flat = flat * n + indices[i];
    }
    return flat;
}

// Helper to iterate over all multi–indices of a given length 'count' with entries in [0, n).
template<std::size_t count>
void iterate_indices(std::size_t n,
                     std::array<std::size_t, count>& indices,
                     std::size_t pos,
                     const std::function<void(const std::array<std::size_t, count>&)>& f)
{
    if (pos == count) {
        f(indices);
        return;
    }
    for (std::size_t i = 0; i < n; ++i) {
        indices[pos] = i;
        iterate_indices<count>(n, indices, pos + 1, f);
    }
}

//---------------------------------------------------------------------
// Vector Space Tags
//---------------------------------------------------------------------

template<typename T, std::size_t n>
struct vector_space {
    typedef T scalar_t;
    static constexpr std::size_t dimension = n;
};

using R3 = vector_space<double,3>;
using R2 = vector_space<double,3>;

// Tag used to indicate that a tensor is totally antisymmetric.
struct alternating_tag {};
struct symmetric_tag {};

//---------------------------------------------------------------------
// General Tensor Template
//---------------------------------------------------------------------
// A general tensor with k contravariant (upper) and l covariant (lower) indices
// defined over a vector space V. Its components are stored in a flat array of size
// (V::dimension)^(k+l).
template<typename V, std::size_t k, std::size_t l, typename Symmetry = void>
struct Tensor {
    typedef typename V::scalar_t T;
    static constexpr std::size_t n = V::dimension;
    static constexpr std::size_t dim = ipow(n, k + l);
    std::array<T, dim> data{};
    typedef Tensor<V,1,0,void> vector_t;
    typedef Tensor<T,0,1,void> covector_t;

    // Access tensor components via an operator() that takes r+s indices.
    template<typename... Indices>
    T& operator()(Indices... indices) {
        static_assert(sizeof...(Indices) == k + l, "Wrong number of indices.");
        std::array<std::size_t, k + l> idx = { static_cast<std::size_t>(indices)... };
        std::size_t flat = multi_index_to_flat(idx, n);
        return data[flat];
    }
    template<typename... Indices>
    const T& operator()(Indices... indices) const {
        static_assert(sizeof...(Indices) == k + l, "Wrong number of indices.");
        std::array<std::size_t, k + l> idx = { static_cast<std::size_t>(indices)... };
        std::size_t flat = multi_index_to_flat(idx, n);
        return data[flat];
    }
    T operator()(std::span<covector_t, k> covecs, std::span<vector_t, l> vecs) const {
        throw std::runtime_error("not not implemented");
    }

    vector_t operator()(std::span<covector_t, k-1> covecs, std::span<vector_t , l> vecs) const {
        throw std::runtime_error("not not implemented");
    }
};

template <typename T,std::size_t n>
using Vector = Tensor<vector_space<T,n>,1,0,void>;


        // When s == 0 and we use alternating, we store only the independent components,
// which are indexed by sorted k–subsets (with count = binom(n, k)).
template<typename V, std::size_t k>
struct Tensor<V, k, 0, alternating_tag> {
    using T = V::scalar_t;
    static constexpr std::size_t n = V::dimension;
    static constexpr std::size_t dim = binom(n, k);
    std::array<T, dim> data{};

    constexpr T& operator()(const std::size_t& i) {
        return data[i];
    }

    constexpr const T& operator()(const std::size_t& i) const {
        return data[i];
    }
};


// Alias: A KVector is a totally antisymmetric contravariant tensor of rank k.
template<typename V, std::size_t k>
using KVector = Tensor<V, k, 0, alternating_tag>;

//---------------------------------------------------------------------
// Exterior Operations: Wedge and Hodge Star for KVectors
//---------------------------------------------------------------------

template<typename V, std::size_t k, std::size_t l>
KVector<V, k+l> wedge(const KVector<V, k>& a, const KVector<V, l>& b) {
    using T = typename V::scalar_t;
    static_assert(k + l <= V::dimension, "Wedge product: k+l must be at most the dimension of V.");
    KVector<V, k+l> result{};
    constexpr std::size_t n = V::dimension;
    constexpr std::size_t dim_k = binom(n, k);
    constexpr std::size_t dim_l = binom(n, l);
    constexpr std::size_t dim_kl = binom(n, k+l);

    // Precompute sorted combinations for each degree.
    std::array<std::array<std::size_t, k>, dim_k> combos_k{};
    std::array<std::array<std::size_t, l>, dim_l> combos_l{};
    std::array<std::array<std::size_t, k+l>, dim_kl> combos_kl{};
    {
        std::array<std::size_t, k> curr{};
        std::size_t idx = 0;
        std::function<void(std::size_t, std::size_t)> gen = [&](std::size_t start, std::size_t pos) {
            if (pos == k) { combos_k[idx++] = curr; return; }
            for (std::size_t i = start; i < n; ++i) {
                curr[pos] = i;
                gen(i + 1, pos + 1);
            }
        };
        gen(0, 0);
    }
    {
        std::array<std::size_t, l> curr{};
        std::size_t idx = 0;
        std::function<void(std::size_t, std::size_t)> gen = [&](std::size_t start, std::size_t pos) {
            if (pos == l) { combos_l[idx++] = curr; return; }
            for (std::size_t i = start; i < n; ++i) {
                curr[pos] = i;
                gen(i + 1, pos + 1);
            }
        };
        gen(0, 0);
    }
    {
        std::array<std::size_t, k+l> curr{};
        std::size_t idx = 0;
        std::function<void(std::size_t, std::size_t)> gen = [&](std::size_t start, std::size_t pos) {
            if (pos == k+l) { combos_kl[idx++] = curr; return; }
            for (std::size_t i = start; i < n; ++i) {
                curr[pos] = i;
                gen(i + 1, pos + 1);
            }
        };
        gen(0, 0);
    }
    // Loop over basis elements of a and b.
    for (std::size_t i = 0; i < dim_k; ++i) {
        T a_coeff = a.data[i];
        if (a_coeff == T{}) continue;
        const auto& I = combos_k[i];
        for (std::size_t j = 0; j < dim_l; ++j) {
            T b_coeff = b.data[j];
            if (b_coeff == T{}) continue;
            const auto& J = combos_l[j];
            // Check that I and J are disjoint.
            bool disjoint = true;
            for (std::size_t p = 0; p < k; ++p)
                for (std::size_t q = 0; q < l; ++q)
                    if (I[p] == J[q]) { disjoint = false; break; }
            if (!disjoint) continue;
            // Merge I and J into a sorted (k+l)-tuple.
            std::array<std::size_t, k+l> merged{};
            std::size_t p = 0, q = 0, r = 0;
            while (p < k && q < l) {
                if (I[p] < J[q])
                    merged[r++] = I[p++];
                else
                    merged[r++] = J[q++];
            }
            while (p < k) merged[r++] = I[p++];
            while (q < l) merged[r++] = J[q++];
            // Compute sign from the shuffle.
            int sign = 1;
            for (std::size_t p1 = 0; p1 < k; ++p1)
                for (std::size_t q1 = 0; q1 < l; ++q1)
                    if (I[p1] > J[q1]) sign = -sign;
            // Find the index in combos_kl corresponding to merged.
            for (std::size_t m = 0; m < dim_kl; ++m) {
                bool match = true;
                for (std::size_t t = 0; t < k+l; ++t) {
                    if (combos_kl[m][t] != merged[t]) { match = false; break; }
                }
                if (match) {
                    result.data[m] += sign * (a_coeff * b_coeff);
                    break;
                }
            }
        }
    }
    return result;
}

template<typename V, std::size_t k>
KVector<V,V::dimension - k> hodge(const KVector<V, k>& a) {
    using T = V::scalar_t;
    constexpr std::size_t n = V::dimension;
    KVector<V, n - k> result{};
    constexpr std::size_t dim_k = binom(n, k);
    constexpr std::size_t dim_nk = binom(n, n - k);
    std::array<std::array<std::size_t, k>, dim_k> combos_k{};
    std::array<std::array<std::size_t, n - k>, dim_nk> combos_nk{};
    {
        std::array<std::size_t, k> curr{};
        std::size_t idx = 0;
        std::function<void(std::size_t, std::size_t)> gen = [&](std::size_t start, std::size_t pos) {
            if (pos == k) { combos_k[idx++] = curr; return; }
            for (std::size_t i = start; i < n; ++i) {
                curr[pos] = i;
                gen(i + 1, pos + 1);
            }
        };
        gen(0, 0);
    }
    {
        std::array<std::size_t, n - k> curr{};
        std::size_t idx = 0;
        std::function<void(std::size_t, std::size_t)> gen = [&](std::size_t start, std::size_t pos) {
            if (pos == n - k) { combos_nk[idx++] = curr; return; }
            for (std::size_t i = start; i < n; ++i) {
                curr[pos] = i;
                gen(i + 1, pos + 1);
            }
        };
        gen(0, 0);
    }
    for (std::size_t i = 0; i < dim_k; ++i) {
        T coeff = a.data[i];
        if (coeff == T{}) continue;
        const auto& I = combos_k[i];
        std::array<std::size_t, n - k> comp{};
        std::size_t pos = 0;
        for (std::size_t j = 0; j < n; ++j) {
            bool found = false;
            for (std::size_t t = 0; t < k; ++t) {
                if (I[t] == j) { found = true; break; }
            }
            if (!found) comp[pos++] = j;
        }
        int sign = 1;
        for (std::size_t p = 0; p < k; ++p)
            for (std::size_t q = 0; q < n - k; ++q)
                if (I[p] > comp[q]) sign = -sign;
        for (std::size_t m = 0; m < dim_nk; ++m) {
            bool match = true;
            for (std::size_t t = 0; t < n - k; ++t) {
                if (combos_nk[m][t] != comp[t]) { match = false; break; }
            }
            if (match) { result.data[m] = sign * coeff; break; }
        }
    }
    return result;
}

//---------------------------------------------------------------------
// T^1_1 Tensor (Matrix) and Matrix Multiplication
//---------------------------------------------------------------------

// A T^1_1 tensor is a general tensor with one upper and one lower index.
template<typename T, typename V>
using Tensor1_1 = Tensor<T, 1, 1, V>;

// TODO: refactor into general function
template<typename T, typename V>
Tensor1_1<T, V> operator*(const Tensor1_1<T, V>& A, const Tensor1_1<T, V>& B) {
    constexpr std::size_t n = V::dimension;
    Tensor1_1<T, V> C{};
    for (std::size_t i1 = 0; i1 < n; ++i1)
        for (std::size_t i2 = 0; i2 < n; ++i2) {
            T sum = T{};
            for (std::size_t i3 = 0; i3 < n; ++i3)
                sum += A(i1, i3) * B(i3, i2);
            C(i1, i2) = sum;
        }
    return C;
}

//---------------------------------------------------------------------
// Change-of-Basis for General (Unsymmetric) Tensors
//---------------------------------------------------------------------

template<typename T, std::size_t r, std::size_t s, typename V>
Tensor<T, r, s, V> change_basis(const Tensor<T, r, s, V>& X,
                                const Tensor1_1<T, V>& P,
                                const Tensor1_1<T, V>& Pinv)
{
    constexpr std::size_t n = V::dimension;
    Tensor<T, r, s, V> result{};
    std::array<std::size_t, r> newU{};
    std::array<std::size_t, s> newL{};
    iterate_indices<r>(n, newU, 0, [&](const std::array<std::size_t, r>& newUpper) {
        iterate_indices<s>(n, newL, 0, [&](const std::array<std::size_t, s>& newLower) {
            std::size_t flat_new = multi_index_to_flat(newUpper, n) * ipow(n, s)
                                   + multi_index_to_flat(newL, n);
            T sum = T{};
            std::array<std::size_t, r> oldU{};
            iterate_indices<r>(n, oldU, 0, [&](const std::array<std::size_t, r>& oldUpper) {
                std::array<std::size_t, s> oldL{};
                iterate_indices<s>(n, oldL, 0, [&](const std::array<std::size_t, s>& oldLower) {
                    std::size_t flat_old = multi_index_to_flat(oldUpper, n) * ipow(n, s)
                                           + multi_index_to_flat(oldL, n);
                    T prod = T{1};
                    for (std::size_t i = 0; i < r; ++i)
                        prod *= P(newUpper[i], oldUpper[i]);
                    for (std::size_t j = 0; j < s; ++j)
                        prod *= Pinv(oldL[j], newL[j]);
                    sum += prod * X.data[flat_old];
                });
            });
            result.data[flat_new] = sum;
        });
    });
    return result;
}







#endif//SL_TENSOR_H
