#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
//#include <complex>
//#include <typeinfo>
//#include <type_traits>
#include "linear_fractional_transformation.hpp"

/*template<class T>
std::istream &operator>>(std::istream &in, std::complex<T> &z)
{
	T x, y;
	in >> x >> y;
	z = std::complex(x, y);
	return in;
}*/

template<class T>
std::istream &operator>>(
	std::istream &in,
	linear_fractional_transformation<T> &f
) {
	return in >> f.a >> f.b >> f.c >> f.d;
}

template<class T>
std::ostream &operator<<(
	std::ostream &out,
	const linear_fractional_transformation<T> &f
) {
	return out
		<< f.a << ' '
		<< f.b << ' '
		<< f.c << ' '
		<< f.d;
}

enum marker_t
{
	ROW = 1,
	COLUMN = 2,
	TABLE = ROW | COLUMN
};

template<class T>
struct default_formatter
{
	const T &data;
	marker_t marker;

	default_formatter(const T &d, const marker_t m = ROW) : data(d), marker(m) {}
};

template<class T>
struct formatter : public default_formatter<T>
{
	std::ios::fmtflags flags;
	std::streamsize width, precision;

	formatter(
		const T &d,
		const marker_t m,
		const std::ios::fmtflags f,
		const std::streamsize w,
		const std::streamsize p
	) : default_formatter<T>(d, m), flags(f), width(w), precision(p) {}
};

/*template<class T>
std::ostream &operator<<(
	std::ostream &out,
	const default_formatter<std::complex<T>> &form
) {
	if (form.marker & ROW)
		out << form.data.real() << ' ' << form.data.imag();
	else
		out << form.data.real() << std::endl << form.data.imag();
	
	return out;
}

template<class T>
std::ostream &operator<<(
	std::ostream &out,
	const formatter<std::complex<T>> &form
) {
	auto oldf = out.flags();
	auto oldp = out.precision();
	out.flags(form.flags);
	out.precision(form.precision);

	if (form.marker & ROW)
		out
			<< std::setw(form.width) << form.data.real() << ' '
			<< std::setw(form.width) << form.data.imag();
	else
		out
			<< std::setw(form.width) << form.data.real() << std::endl
			<< std::setw(form.width) << form.data.imag();
	
	out.flags(oldf);
	out.precision(oldp);
	
	return out;
}*/

template<class T>
std::ostream &operator<<(
	std::ostream &out,
	const default_formatter<linear_fractional_transformation<T>> &form
) {
	switch (form.marker)
	{
	case ROW:
		out
			<< form.data.a << ' '
			<< form.data.b << ' '
			<< form.data.c << ' '
			<< form.data.d;
		break;
	
	case COLUMN:
		out
			<< form.data.a << std::endl
			<< form.data.b << std::endl
			<< form.data.c << std::endl
			<< form.data.d;
		break;
	
	case TABLE:
		out
			<< form.data.a << ' '
			<< form.data.b << std::endl
			<< form.data.c << ' '
			<< form.data.d;
		break;
	}
	
	return out;
}

template<class T>
std::ostream &operator<<(
	std::ostream &out,
	const formatter<linear_fractional_transformation<T>> &form
) {
	auto oldf = out.flags();
	auto oldp = out.precision();
	out.flags(form.flags);
	out.precision(form.precision);

	switch (form.marker)
	{
	case ROW:
		out
			<< std::setw(form.width) << form.data.a << ' '
			<< std::setw(form.width) << form.data.b << ' '
			<< std::setw(form.width) << form.data.c << ' '
			<< std::setw(form.width) << form.data.d;
		break;
	
	case COLUMN:
		out
			<< std::setw(form.width) << form.data.a << std::endl
			<< std::setw(form.width) << form.data.b << std::endl
			<< std::setw(form.width) << form.data.c << std::endl
			<< std::setw(form.width) << form.data.d;
		break;
	
	case TABLE:
		out
			<< std::setw(form.width) << form.data.a << ' '
			<< std::setw(form.width) << form.data.b << std::endl
			<< std::setw(form.width) << form.data.c << ' '
			<< std::setw(form.width) << form.data.d;
		break;
	}
	
	out.flags(oldf);
	out.precision(oldp);
	
	return out;
}

template<class T>
class typed_ifstream : public std::ifstream
{
public:
	typed_ifstream() : std::ifstream() {}
	
	template<class... Args>
	typed_ifstream(Args... args) : std::ifstream(args...) {}

	auto &read(T &arg)
	{
		return std::ifstream::read((char*)(&arg), sizeof(T));
	}

	auto &read_vector(std::vector<T> &arg)
	{
		/*static_assert(
			std::is_same<T, typename U::value_type>::value,
			std::string("The types ") + typeid(T).name() +
			" and " + typeid(U).name() + " must be same."
		);//*/
		return std::ifstream::read((char*)arg.data(), sizeof(T) * arg.size());
	}
};

template<class T>
class typed_ofstream : public std::ofstream
{
public:
	typed_ofstream() : std::ofstream() {}
	
	template<class... Args>
	typed_ofstream(Args... args) : std::ofstream(args...) {}

	auto &write(const T &arg)
	{
		return std::ofstream::write((char*)(&arg), sizeof(T));
	}

	auto &write_vector(const std::vector<T> &arg)
	{
		/*static_assert(
			std::is_same<T, typename U::value_type>::value,
			std::string("The types ") + typeid(T).name() +
			" and " + typeid(U).name() + " must be same."
		);//*/
		return std::ofstream::write((char*)arg.data(), sizeof(T) * arg.size());
	}
};