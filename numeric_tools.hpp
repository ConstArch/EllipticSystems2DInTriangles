#pragma once

#include <cmath>
#include <exception>
#include <stdexcept>

#define SWAP(a, b) { auto t = a; a = b; b = t; }
//#define SWAP_INT(a, b) { (a) ^= (b); (b) ^= (a); (a) ^= (b); }

const float       PIf = 3.141593F;
const double      PI  = 3.1415926535897932;
const long double PIl = 3.14159265358979323846264338327950288L;

const float       Ef = 2.718282F;
const double      E  = 2.7182818284590452;
const long double El = 2.71828182845904523536028747135266250L;

double EPSILON = 1e-10;

template<typename T, typename U>
constexpr inline T round_to(const U arg)
{
	return static_cast<T>(std::round(arg));
}

template<typename T, typename U>
constexpr inline T floor_to(const U arg)
{
	const auto r = std::round(arg);
	return static_cast<T>(std::abs(arg - r) < EPSILON || arg > r ? r : r - 1);
}

template<typename T, typename U>
constexpr inline T ceil_to(const U arg)
{
	const auto r = std::round(arg);
	return static_cast<T>(std::abs(arg - r) < EPSILON || arg < r ? r : r + 1);
}

template<typename T, typename U>
inline bool equal_fp(const T x, const U y)
{
	return std::abs(x - y) <= EPSILON;
}

template<typename T, typename U>
inline bool unequal_fp(const T x, const U y)
{
	return std::abs(x - y) > EPSILON;
}

template<typename T, typename U>
inline bool less_fp(const T x, const U y)
{
	return y - x > EPSILON;
}

template<typename T, typename U>
inline bool greater_fp(const T x, const U y)
{
	return x - y > EPSILON;
}

template<typename T, typename U>
inline bool less_equal_fp(const T x, const U y)
{
	return x - y <= EPSILON;
}

template<typename T, typename U>
inline bool greater_equal_fp(const T x, const U y)
{
	return y - x <= EPSILON;
}

/*template<void *distance, typename T, typename U>
inline bool equal_vectors_fp(const T x, const U y)
{
	return distance(x, y) <= EPSILON;
}

template<void *distance, typename T, typename U>
inline bool unequal_vectors_fp(const T x, const U y)
{
	return distance(x, y) > EPSILON;
}*/

template<typename distance_t>
class equal_vectors_fp
{
private:
	distance_t distance;
	
public:
	equal_vectors_fp(const distance_t dist) : distance(dist) {}
	
	template<typename T, typename U>
	inline bool operator()(const T x, const U y) const
	{
		return distance(x, y) <= EPSILON;
	}
};

template<typename distance_t>
class unequal_vectors_fp
{
private:
	distance_t distance;
	
public:
	unequal_vectors_fp(const distance_t dist) : distance(dist) {}
	
	template<typename T, typename U>
	inline bool operator()(const T x, const U y) const
	{
		return distance(x, y) > EPSILON;
	}
};

template<typename T, typename U>
bool equal_collections(const T &xx, const U &yy)
{
	if (xx.size() != yy.size())
		throw std::invalid_argument("The collections must be the same length.");
	
	auto itx = std::begin(xx);
	auto ity = std::begin(yy);
	const auto endx = std::end(xx);
	while (itx != endx)
		if (*itx++ != *ity++)
			return false;
	return true;
}

template<typename T, typename U, typename comparator_t>
bool equal_collections(const T &xx, const U &yy, const comparator_t equal)
{
	if (xx.size() != yy.size())
		throw std::invalid_argument("The collections must be the same length.");
	
	auto itx = std::begin(xx);
	auto ity = std::begin(yy);
	const auto endx = std::end(xx);
	while (itx != endx)
		if (!equal(*itx++, *ity++))
			return false;
	return true;
}

template<typename T, typename U>
inline bool equal_collections_fp(const T &xx, const U &yy)
{
	return equal_collections(
		xx, yy,
		equal_fp<typename T::value_type, typename U::value_type>
	);
}

template<typename T>
void delete_dublicates(std::vector<T> &array)
{
    std::size_t new_size = 1;
    for (std::size_t i = 1; i < array.size(); ++i)
    {
        std::size_t j = 0;
        for (; j < new_size; ++j)
        {
            if (array[i] == array[j])
                break;
        }
        if (j == new_size)
            array[new_size++] = array[i];
    }
    array.resize(new_size);
}

template<typename T, typename comparator_t>
void delete_dublicates(std::vector<T> &array, const comparator_t equal)
{
    std::size_t new_size = 1;
    for (std::size_t i = 1; i < array.size(); ++i)
    {
        std::size_t j = 0;
        for (; j < new_size; ++j)
        {
            if (equal(array[i], array[j]))
                break;
        }
        if (j == new_size)
            array[new_size++] = array[i];
    }
    array.resize(new_size);
}

template<typename T>
constexpr T int_pow(T base, const int exp)
{
	T result = (T)1;
	if (exp >= 0)
		for (int i = 0; i < exp; ++i) result *= base;
	else
	{
		base = (T)1 / base;
		for (int i = 0; i > exp; --i) result *= base;
	}
	return result;
}

template<typename T>
constexpr T pow_sum(T base, int exp_1, int exp_2)
{
	if (exp_1 > exp_2) SWAP(exp_1, exp_2);

	T result, q;
	int abs_1 = std::abs(exp_1), abs_2 = std::abs(exp_2);
	if ((exp_1 ^ exp_2) & 0x80000000)
	{
		T inv_base = (T)1 / base;
		if (std::abs(base) < 1)
		{
			SWAP(base, inv_base);
			SWAP(abs_1, abs_2);
		}

		result = (T)1;

		q = inv_base;
		for (int i = 0; i < abs_1; ++i, q *= inv_base) result += q;
		
		q = base;
		for (int i = 0; i < abs_2; ++i, q *= base) result += q;
	}
	else
	{
		if (abs_1 > abs_2) SWAP(abs_1, abs_2);

		if (exp_1 < 0) base = (T)1 / base;
		
		result = int_pow(base, abs_1);
		q = result * base;
		for (int i = abs_1; i < abs_2; ++i, q *= base) result += q;
	}
	return result;
}