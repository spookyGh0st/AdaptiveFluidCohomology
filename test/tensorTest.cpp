#include <gtest/gtest.h>
#include <sl/tensor.h>
TEST(kvecTest,tmp) {
    using E3 = vector_space<double, 3>;
    std::cout << "=== Euclidean Space Example (R^3) ===\n";

    /*
    // Create a general (unsymmetric) tensor T^2_0 over E3.
    Tensor<E3, 2, 0> T20{};
    for (std::size_t i = 0; i < E3::dimension; ++i)
        for (std::size_t j = 0; j < E3::dimension; ++j)
            T20(i, j) = static_cast<double>(i * E3::dimension + j + 1);
    std::cout << "Original T^2_0 tensor (unsymmetric):\n";
    for (std::size_t i = 0; i < E3::dimension; ++i) {
        for (std::size_t j = 0; j < E3::dimension; ++j)
            std::cout << T20(i, j) << " ";
        std::cout << "\n";
    }
     */
    std::cout << "\n";
}

TEST(ExteriorAlgebraTest,wedge_product) {
    // Exterior algebra over E3:
    // Let’s build a 1–vector and compute its Hodge star.
    KVector<R3, 1> e0,e1,e2;
    e0(0) = 1; e1(1) = 1; e2(2) = 1;
    auto e01 = wedge(e0,e1), e12 = wedge(e1,e2), e20=wedge(e2,e0);
    for (auto c : e01.data)
        std::cout << c << " ";
    auto e123 = wedge(e01,e2);
    // TODO:
    // ASSERT_EQ(wedge(e01,e2),wedge(e20,e1));
}

TEST(ExteriorAlgebraTest,hodgeStar) {
    // Exterior algebra over E3:
    // Let’s build a 1–vector and compute its Hodge star.
    KVector<R3, 1> v1{};
    // The basis elements for a 1–vector over E3 correspond to {0}, {1}, {2}.
    v1.data[0] = 0.0; v1.data[1] = 0.0; v1.data[2] = 1.0;
    auto star_v1 = hodge(v1);
    std::cout << "Hodge star of v1 (a 1–vector in R^3 becomes a 2–vector):\n";
    for (auto c : star_v1.data)
        std::cout << c << " ";
    ASSERT_NEAR(star_v1(0),1,0.01);
}
