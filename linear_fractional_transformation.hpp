#pragma once

#include <cmath>

template<class T>
class linear_fractional_transformation
{
public:
	T a, b, c, d;

	linear_fractional_transformation() : a(1), b(0), c(0), d(1) {}

	template<class T1, class T2, class T3, class T4>
	linear_fractional_transformation(
		const T1 aa,
		const T2 bb,
		const T3 cc,
		const T4 dd
	) : a(aa), b(bb), c(cc), d(dd) {}
	
	template<class U>
	inline auto operator*(const linear_fractional_transformation<U> &other) const
	{
		return linear_fractional_transformation(
			a * other.a + b * other.c,
			a * other.b + b * other.d,
			c * other.a + d * other.c,
			c * other.b + d * other.d
		);
	}

	template<class U>
	inline auto operator()(const U z) const
	{
		return (a * z + b) / (c * z + d);
	}

	template<class U>
	inline auto map(const U &zz) const
	{
		U ww = zz;
		for (auto &e : ww) e = operator()(e);
		return ww;
	}

	inline auto inverse() const
	{
		return linear_fractional_transformation(d, -b, -c, a);
	}

	inline auto inverse_matrix() const
	{
		const auto k = det();
		return linear_fractional_transformation(d / k, -b / k, -c / k, a / k);
	}
	
	inline auto det() const
	{
		return a * d - b * c;
	}

	inline auto cancel_out() const
	{
		const auto k = std::sqrt(det());
		return linear_fractional_transformation(a / k, b / k, c / k, d / k);
	}
	
	inline auto Euclidean_norm() const
	{
		const auto
			aa = std::abs(a),
			bb = std::abs(b),
			cc = std::abs(c),
			dd = std::abs(d);
		return std::sqrt(aa * aa + bb * bb + cc * cc + dd * dd);
	}
};

template<class T, class U>
inline bool equal_coefficients(
	const linear_fractional_transformation<T> &f,
	const linear_fractional_transformation<U> &g
) {
	return f.a == g.a && f.b == g.b && f.c == g.c && f.d == g.d;
}

template<class T, class U>
inline bool operator==(
	const linear_fractional_transformation<T> &f,
	const linear_fractional_transformation<U> &g
) {
	//const auto ff = f.cancel_out(), gg = g.cancel_out();
	return equal_coefficients(f.cancel_out(), g.cancel_out());
}

template<class T>
inline auto inverse(const linear_fractional_transformation<T> &f)
{
	return f.inverse();
}

template<class T>
inline auto inverse_matrix(const linear_fractional_transformation<T> &f)
{
	return f.inverse_matrix();
}

template<class T>
inline auto det(const linear_fractional_transformation<T> &f)
{
	return f.det();
}

template<class T>
inline auto cancel_out(const linear_fractional_transformation<T> &f)
{
	return f.cancel_out();
}

template<class T>
inline auto Euclidean_norm(const linear_fractional_transformation<T> &f)
{
	return f.Euclidean_norm();
}

template<class T, class U>
inline auto Euclidean_distance(
	const linear_fractional_transformation<T> &f,
	const linear_fractional_transformation<U> &g
) {
	const auto
		aa = std::abs(f.a - g.a),
		bb = std::abs(f.b - g.b),
		cc = std::abs(f.c - g.c),
		dd = std::abs(f.d - g.d);
	return std::sqrt(aa * aa + bb * bb + cc * cc + dd * dd);
}