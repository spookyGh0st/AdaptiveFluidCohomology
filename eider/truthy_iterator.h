#pragma once
#include <iterator>

template<typename T, size_t N>
class truthy_iterator {
    T* ptr;
    T* end_ptr;

    void skip_false() {
        while (ptr != end_ptr && !static_cast<bool>(*ptr))
            ++ptr;
    }

public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;

    truthy_iterator(T* start, T* end) : ptr(start), end_ptr(end) { skip_false(); }

    reference operator*() const { return *ptr; }
    pointer operator->() const { return ptr; }

    truthy_iterator& operator++() {
        ++ptr;
        skip_false();
        return *this;
    }

    bool operator==(const truthy_iterator& other) const { return ptr == other.ptr; }
    bool operator!=(const truthy_iterator& other) const { return ptr != other.ptr; }
};

template<typename T, size_t N>
class truthy_range {
    T* data;
public:
    truthy_range(T (&arr)[N]) : data(arr) {}
    auto begin() { return truthy_iterator<T, N>(data, data + N); }
    auto end() { return truthy_iterator<T, N>(data + N, data + N); }
};
