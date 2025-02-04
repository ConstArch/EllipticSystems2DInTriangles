#pragma once

#include <complex>
#include <list>
#include <vector>
#include <cmath>
#include <thread>
#include "linear_fractional_transformation.hpp"
#include "theta_series.hpp"
#include "numeric_tools.hpp"
#include "io_tools.hpp"

using namespace std::complex_literals;

template<class T> const auto I = std::complex<T>(0, 1);

template<class real>
struct triangle
{
private:
	using complex = std::complex<real>;
	
public:
	complex A, B, C;
	
	triangle() {}
	
	triangle(const complex AA, const complex BB, const complex CC) : A(AA), B(BB), C(CC) {}
};

template<class T>
std::istream &operator>>(std::istream &in, triangle<T> &tr)
{
	return in >> tr.A >> tr.B >> tr.C;
}

template<class real>
class solution
{
private:
	enum flag_t { GEN_1, GEN_2, INV_1, INV_2 };

	using complex = std::complex<real>;
	using transform = linear_fractional_transformation<complex>;
	
	complex a, b;
	real tau;
	transform P;
	theta_series<complex> th1, th2;

	void __build__(
		const complex zeta,
		const unsigned int level,
		const unsigned int m,
		const transform &h1,
		const transform &h2
	) {
		const auto zeta_c = std::conj(zeta);

		const auto A1 = (zeta - (real)1) / (zeta_c - (real)1);
		const auto B1 = -(real)2 * I<real> * std::imag(zeta) / (zeta_c - (real)1);
		const auto A2 = zeta / zeta_c;
		const auto B2 = 0;

		const auto a1 = (A1 - tau) / (1 - tau);
		const auto b1 = B1 / (1 - tau);
		const auto a2 = (A2 - tau) / (1 - tau);
		const auto b2 = B2 / (1 - tau);

		const transform Id, T1(a1, b1, 0, 1), T2(a2, b2, 0, 1);

		const auto invP = inverse_matrix(P);

		const auto S1 = cancel_out(invP * T1 * P);
		const auto S2 = cancel_out(invP * T2 * P);
		const auto I1 = inverse(S1), I2 = inverse(S2);

		std::list<std::vector<transform>> data;

		std::size_t len = 4;
		data.push_back(std::vector<transform>({ Id, S1, S2, I1, I2 }));
		std::vector<flag_t> flags = { GEN_1, GEN_2, INV_1, INV_2 };
		std::vector<transform> dt;
		std::vector<flag_t> fg;
		for (unsigned int i = 0; i < level; ++i)
		{
			len *= 3;
			dt.reserve(len);
			fg.reserve(len);
			std::size_t j = 0;
			for (const auto e : data.back())
			{
				switch (flags[j++])
				{
				case GEN_1:
					dt.push_back(e * S1); fg.push_back(GEN_1);
					dt.push_back(e * S2); fg.push_back(GEN_2);
					dt.push_back(e * I2); fg.push_back(INV_2);
					break;
				
				case GEN_2:
					dt.push_back(e * S1); fg.push_back(GEN_1);
					dt.push_back(e * S2); fg.push_back(GEN_2);
					dt.push_back(e * I1); fg.push_back(INV_1);
					break;
				
				case INV_1:
					dt.push_back(e * S2); fg.push_back(GEN_2);
					dt.push_back(e * I1); fg.push_back(INV_1);
					dt.push_back(e * I2); fg.push_back(INV_2);
					break;
				
				case INV_2:
					dt.push_back(e * S1); fg.push_back(GEN_1);
					dt.push_back(e * I1); fg.push_back(INV_1);
					dt.push_back(e * I2); fg.push_back(INV_2);
					break;
				}
			}
			data.push_back(std::move(dt));
			flags = std::move(fg);
			dt.clear(); fg.clear();
		}
		dt.reserve(1 + 4 * pow_sum(3, 0, level));
		for (const auto &d : data)
			dt.insert(dt.end(), d.begin(), d.end());
		delete_dublicates(dt, equal_vectors_fp(Euclidean_distance<complex, complex>));
		
		/*typed_ofstream<transform> fout(
			"group.dat", std::ios::out | std::ios::binary | std::ios::trunc
		);
		fout.write_vector(dt);
		fout.close();//*/
		
		th1.build(m, h1, dt);
		th2.build(m, h2, dt);
	}
	
public:
	solution() {}
	
	template<class T>
	solution(
		const real tt,
		const triangle<T> &tr,
		const transform &PP,
		const unsigned int level,
		const unsigned int m,
		const transform &h1,
		const transform &h2
	) : a((T)1 / (tr.B - tr.A)), b(tr.A / (tr.A - tr.B)), tau(tt), P(PP)
	{
		__build__(a * tr.C + b, level, m, h1, h2);
	}

	template<class T>
	void build(
		const real tt,
		const triangle<T> &tr,
		const transform &PP,
		const unsigned int level,
		const unsigned int m,
		const transform &h1,
		const transform &h2
	) {
		a = (T)1 / (tr.B - tr.A);
		b = tr.A / (tr.A - tr.B);
		tau = tt; P = PP;
		__build__(a * tr.C + b, level, m, h1, h2);
	}
	
	inline std::size_t members_count() const
	{
		return th1.members_count();
	}
	
	complex operator()(const complex &z) const
	{
		const auto z_L = a * z + b;
		
		const auto z_c = std::conj(z_L);
		const auto z_t = (z_L - tau * z_c) / (1 - tau);
		
		const auto invP = inverse(P);
		
		const auto zeta_c = invP(z_c);
		const auto zeta_t = invP(z_t);
		
		const auto th1_c = th1(zeta_c);
		const auto th2_c = th2(zeta_c);
		const auto th1_t = th1(zeta_t);
		const auto th2_t = th2(zeta_t);
		
		return th1_t / th2_t - th1_c / th2_c;
	}
	
	template<class U>
	inline auto map(const U &zz) const
	{
		U ww = zz;
		for (auto &e : ww) e = operator()(e);
		return ww;
	}
	
	template<class U>
	inline auto parallel_map(const U &zz) const
	{
		const std::size_t len = zz.size();
		const std::size_t div = 4;//std::sqrt(len);
		const std::size_t quot = len / div;
		const std::size_t rem = len % div;
		
		U ww(len);
		auto it_z = zz.begin(), it_w = ww.begin();
		const auto invP = inverse(P);
		
		#define INNER_PARALLEL_MAP(count)															\
		{																							\
			U a1_c(count), a1_t(count);																\
			for (																					\
				auto it_zeta_c = a1_c.begin(), it_zeta_t = a1_t.begin();							\
				it_zeta_c != a1_c.end();															\
				++it_z, ++it_zeta_c, ++it_zeta_t													\
			) {																						\
				const auto z_L = a * *it_z + b;														\
				const auto z_c = std::conj(z_L);													\
				const auto z_t = (z_L - tau * z_c) / (1 - tau);										\
																									\
				*it_zeta_c = invP(z_c);																\
				*it_zeta_t = invP(z_t);																\
			}																						\
																									\
			U a2_c = a1_c, a2_t = a1_t;																\
			std::thread thr1_c(theta_series<complex>::template transform<U>, &th1, std::ref(a1_c));	\
			std::thread thr2_c(theta_series<complex>::template transform<U>, &th2, std::ref(a2_c));	\
			std::thread thr1_t(theta_series<complex>::template transform<U>, &th1, std::ref(a1_t));	\
			std::thread thr2_t(theta_series<complex>::template transform<U>, &th2, std::ref(a2_t));	\
			thr1_c.join();																			\
			thr2_c.join();																			\
			thr1_t.join();																			\
			thr2_t.join();																			\
																									\
			auto it_a1_c = a1_c.begin();															\
			auto it_a2_c = a2_c.begin();															\
			auto it_a1_t = a1_t.begin();															\
			auto it_a2_t = a2_t.begin();															\
			while (it_a1_c != a1_c.end())															\
				*it_w++ = *it_a1_t++ / *it_a2_t++ - *it_a1_c++ / *it_a2_c++;						\
		}
		
		for (std::size_t i = 0; i < div; ++i)
			INNER_PARALLEL_MAP(quot);
		
		INNER_PARALLEL_MAP(rem)
		
		return ww;
	}
};