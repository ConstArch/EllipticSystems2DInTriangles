#pragma once

#include <vector>
#include "linear_fractional_transformation.hpp"
#include "numeric_tools.hpp"

template<class T>
class theta_series
{
private:
	struct member_t
	{
		linear_fractional_transformation<T> f;
		T c, d;
		member_t() : f(), c(), d() {}
		member_t(
			const linear_fractional_transformation<T> &ff,
			const T cc, const T dd
		) : f(ff), c(cc), d(dd) {}
	};

	std::vector<member_t> members;

	unsigned int m;

public:
	theta_series() : m(0), members(0) {}

	template<class U, class V>
	theta_series(
		const unsigned int mm,
		const linear_fractional_transformation<U> &h,
		const std::vector<linear_fractional_transformation<V>> &G
	) : m(mm), members(G.size())
	{
		for (const auto g : G)
			members.emplace_back(h * g, g.c, g.d);
	}

	template<class U, class V>
	void build(
		const unsigned int mm,
		const linear_fractional_transformation<U> &h,
		const std::vector<linear_fractional_transformation<V>> &G
	) {
		m = mm;
		members.reserve(G.size());
		for (const auto g : G)
			members.emplace_back(h * g, g.c, g.d);
	}

	template<class U>
	inline auto operator()(const U z) const
	{
		decltype(members.front().c * z) w = 0;
		for (const auto member : members)
			w += member.f(z) * int_pow(member.c * z + member.d, -2 * m);
		return w;
	}

	template<class U>
	inline auto map(const U &zz) const
	{
		U ww = zz;
		for (auto &e : ww) e = operator()(e);
		return ww;
	}
	
	template<class U>
	inline auto transform(U &zz) const
	{
		for (auto &z : zz) z = operator()(z);
		return zz;
	}
	
	/*template<class iterator_t>
	inline auto transform(const iterator_t begin, const iterator_t end) const
	{
		for (auto it = begin; it != end; ++it)
			*it = operator()(*it);
		return zz;
	}//*/
	
	inline std::size_t members_count() const
	{
		return members.size();
	}
};