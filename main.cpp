#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <string>
#include "linear_fractional_transformation.hpp"
#include "solution.hpp"
#include "io_tools.hpp"

using real = double;
using complex = std::complex<real>;
using transform = linear_fractional_transformation<complex>;

using namespace std::complex_literals;

int main(int argc, char **argv)
{
	const char *args_file_address, *mesh_file_address;
	if (argc >= 3)
	{
		args_file_address = argv[1];
		mesh_file_address = argv[2];
	}
	else
	{
		args_file_address = "args.txt";
		mesh_file_address = "mesh.txt";
	}

	real tau;
	triangle<real> tr;
	complex z_singular;

	transform P, H1, H2;
	unsigned int level, m;

	std::ifstream fin;
	
	fin.open(args_file_address);
	if (!fin.is_open())
	{
		std::cerr << "File \'" << args_file_address << "\' not found.\n";
		return 0;
	}
	fin >> tau >> tr >> z_singular >> level >> m;
	fin.close();

	P = transform(2. + 7.i, 9., 6.i, 11.);
	H1 = transform(0., 1., 1., -inverse(P)((z_singular - tr.A) / (tr.B - tr.A)));
	H2 = transform();

	/*std::cout
		<< "tau = " << tau << std::endl
		<< "zeta = " << zeta << std::endl
		<< "P = " << P << std::endl
		<< "level = " << level << std::endl
		<< "m = " << m << std::endl
		<< "H1 = " << H1 << std::endl
		<< "H2 = " << H2 << std::endl;//*/

	auto t = clock();
	solution<real> f(tau, tr, P, level, m, H1, H2);
	auto dt = (double)(clock() - t) / CLOCKS_PER_SEC;
	std::cout
		<< "Approximate solution with " << f.members_count()
		<< " members of the series is built in ";
	if (dt < 1)
		std::cout << dt * 1000 << " ms.\n";
	else
		std::cout << dt << " sec.\n";

	/*complex z;
	std::cin >> z;
	std::cout
		<< "z := " << z << std::endl
		<< "f(z) -> " << f(z) << std::endl;//*/

	real x_min, x_max, y_min, y_max;

	std::size_t x_count, y_count;

	fin.open(mesh_file_address);
	if (!fin.is_open())
	{
		std::cerr << "File \'" << mesh_file_address << "\' not found.\n";
		return 0;
	}
	fin >> x_min >> x_max >> x_count >> y_min >> y_max >> y_count;
	fin.close();

	/*std::cout
		<< "x_min = " << x_min << ", x_max = " << x_max
		<< ", x_count = " << x_count << std::endl
		<< "y_min = " << y_min << ", y_max = " << y_max
		<< ", y_count = " << y_count << std::endl;//*/

	const real
		x_step = (x_max - x_min) / (x_count - 1),
		y_step = (y_max - y_min) / (y_count - 1);

	std::vector<real> xx(x_count), yy(y_count);
	for (std::size_t i = 0; i < x_count; ++i)
		xx.at(i) = x_min + i * x_step;
	for (std::size_t i = 0; i < y_count; ++i)
		yy.at(i) = y_min + i * y_step;

	std::vector<complex> mesh;
	mesh.reserve(x_count * y_count);
	for (const auto x : xx)
		for (const auto y : yy)
			mesh.emplace_back(x, y);
	
	t = clock();
	std::vector<complex> values = f.parallel_map(mesh);
	dt = (double)(clock() - t) / CLOCKS_PER_SEC;
	std::cout
		<< "The values are calculated on "
		<< values.size() << " points in ";
	if (dt < 1)
		std::cout << dt * 1000 << " ms.\n";
	else
		std::cout << dt << " sec.\n";
	
	typed_ofstream<complex> fout;
	fout.open(
		std::string("mesh[") + mesh_file_address + "].dat",
		std::ios::out | std::ios::binary | std::ios::trunc
	);
	fout.write_vector(mesh);
	fout.close();
	fout.open(
		std::string("values[") + args_file_address + "][" + mesh_file_address + "].dat",
		std::ios::out | std::ios::binary | std::ios::trunc
	);
	fout.write_vector(values);
	fout.close();//*/

	return 0;
}