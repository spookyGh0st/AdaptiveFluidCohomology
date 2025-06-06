#pragma once

#include <gtest/gtest.h>
#include <tbb/tbb.h>

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
    // Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
}


using namespace oneapi::tbb;

class ApplyFoo {
    float *const my_a;
public:
    void operator()( const blocked_range<size_t>& r ) const {
        float *a = my_a;
        for( size_t i=r.begin(); i!=r.end(); ++i )
            a[i] *= 2;
    }
    ApplyFoo( float a[] ) :
        my_a(a)
    {}
};
// Demonstrate some basic assertions.
TEST(HelloTest, testTBB) {
    std::vector<float> a( 1000000000);
    parallel_for( blocked_range<size_t>(0,a.size()),
    [=,&a](const blocked_range<size_t>& r) {
        for(size_t i=r.begin(); i!=r.end(); ++i)
            a[i] *= 2;
    }
);
}
