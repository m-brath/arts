#pragma once

#include <debug.h>

#include <algorithm>
#include <concepts>
#include <functional>

#include "matpack_mdspan.h"

namespace matpack {
template <class Compare>
class grid_t;

template <class Compare>
class grid_view_t {
  grid_t<Compare>& x;

 public:
  [[nodiscard]] Index size() const { return x.x.size(); }

  grid_view_t(grid_t<Compare>& x_) : x(x_) {}

  [[nodiscard]] const Vector& vec() const { return x.x; }

  grid_view_t& operator=(const grid_t<Compare>& y) {
    if (&y.x != &x.x) VectorView{x.x} = y.x;
    return *this;
  }

  grid_view_t& operator=(const ConstVectorView& y) {
    *this = grid_t<Compare>{y};
    return *this;
  }
};

template <class Compare>
class grid_t {
  Vector x;

  friend class grid_view_t<Compare>;

 public:
  [[nodiscard]] const Vector& vec() const { return x; }
  [[nodiscard]] Vector vec() && { return std::move(x); }

  using value_type = Numeric;

  static constexpr bool is_sorted(const exact_md<Numeric, 1> auto& x) {
    return stdr::is_sorted(x, Compare{});
  }

  static void assert_sorted(const Vector& x) {
    ARTS_USER_ERROR_IF(not is_sorted(x), "Wrong sorting:\n{:B,}", x);
  }

  //! Unsafe constructor
  grid_t(Index N) : x(N) {}

  grid_t(std::initializer_list<Numeric> il) : x(il) { assert_sorted(x); }

  template <typename... Ts>
  explicit grid_t(Ts&&... ts)
    requires std::constructible_from<Vector, Ts...>
      : x(Vector{std::forward<Ts>(ts)...}) {
    assert_sorted(x);
  }

  template <class Iterator, class Func>
  grid_t(const Iterator& begin, const Iterator& end, Func fun)
      : x(std::distance(begin, end)) {
    std::transform(begin, end, x.begin(), fun);
    assert_sorted(x);
  }

  grid_t(Vector&& in) : x(std::move(in)) { assert_sorted(x); }

  grid_t& operator=(Vector&& in) {
    assert_sorted(in);
    x = std::move(in);
    return *this;
  }

  template <typename T>
  grid_t& operator=(const T& in) {
    assert_sorted(in);
    x = in;
    return *this;
  }

  grid_t()                         = default;
  grid_t(grid_t&&)                 = default;
  grid_t(const grid_t&)            = default;
  grid_t& operator=(grid_t&&)      = default;
  grid_t& operator=(const grid_t&) = default;

  operator Vector() && { return Vector{std::move(x)}; }
  operator const Vector&() const { return vec(); }
  operator ConstVectorView() const { return vec(); }
  operator StridedConstVectorView() const { return vec(); }

  template <access_operator Op>
  [[nodiscard]] constexpr auto operator[](const Op& op) const {
    return x[op];
  }

  [[nodiscard]] constexpr auto size() const { return x.size(); }
  [[nodiscard]] constexpr auto begin() const { return x.begin(); }
  [[nodiscard]] constexpr auto end() const { return x.end(); }
  [[nodiscard]] constexpr auto empty() const { return x.empty(); }
  [[nodiscard]] constexpr Numeric front() const { return x.front(); }
  [[nodiscard]] constexpr Numeric back() const { return x.back(); }

  [[nodiscard]] constexpr std::array<Index, 1> shape() const {
    return x.shape();
  }
};

template <class C1, class C2>
[[nodiscard]] constexpr bool operator==(const grid_t<C1>& x,
                                        const grid_t<C2>& y) {
  return x.vec() == y.vec();
}

template <class C1, class C2>
[[nodiscard]] constexpr bool operator!=(const grid_t<C1>& x,
                                        const grid_t<C2>& y) {
  return x.vec() != y.vec();
}

template <class Compare, exact_md<Numeric, 1> md>
[[nodiscard]] constexpr bool operator==(const grid_t<Compare>& x, const md& y) {
  return x.vec() == y;
}

template <class Compare, exact_md<Numeric, 1> md>
[[nodiscard]] constexpr bool operator!=(const grid_t<Compare>& x, const md& y) {
  return x.vec() != y;
}

template <class Compare, exact_md<Numeric, 1> md>
[[nodiscard]] constexpr bool operator==(const md& y, const grid_t<Compare>& x) {
  return y == x.vec();
}

template <class Compare, exact_md<Numeric, 1> md>
[[nodiscard]] constexpr bool operator!=(const md& y, const grid_t<Compare>& x) {
  return y != x.vec();
}
}  // namespace matpack

using AscendingGrid        = matpack::grid_t<std::less_equal<>>;
using AscendingGridView    = matpack::grid_view_t<std::less_equal<>>;
using DescendingGrid       = matpack::grid_t<std::greater_equal<>>;
using ArrayOfAscendingGrid = std::vector<AscendingGrid>;

template <class Compare>
struct std::formatter<matpack::grid_t<Compare>> {
  std::formatter<matpack::strided_view_t<const Numeric, 1>> fmt;

  [[nodiscard]] constexpr auto& inner_fmt() { return fmt.inner_fmt(); }
  [[nodiscard]] constexpr auto& inner_fmt() const { return fmt.inner_fmt(); }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return fmt.parse(ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const matpack::grid_t<Compare>& v,
                              FmtContext& ctx) const {
    return fmt.format(v, ctx);
  }
};

template <class Compare>
struct std::formatter<matpack::grid_view_t<Compare>> {
  std::formatter<Vector> fmt;

  [[nodiscard]] constexpr auto& inner_fmt() { return fmt.inner_fmt(); }
  [[nodiscard]] constexpr auto& inner_fmt() const { return fmt.inner_fmt(); }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return fmt.parse(ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const matpack::grid_view_t<Compare>& v,
                              FmtContext& ctx) const {
    return fmt.format(v.vec(), ctx);
  }
};
