#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <random>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <FftComplex.h>

using std::cout;
using std::uintmax_t;
using std::vector;

std::default_random_engine randGen((std::random_device())());

// The size of the array for FFT must be in n^2 order

static size_t bitReversal(size_t value, int width)
{
    size_t result = 0;
    for (int i = 0; i < width; i++, value >>= 1)
    {
        result = (result << 1) | (value & 1U);
    }
    return result;
}
static size_t powerOfTwo(size_t length)
{

    size_t rounding = (int)ceil((log(length) / log(2)));
    return pow(2, rounding);
}

static bool checkpowerofTwo(size_t length)
{
    return length != 0 && (length & (length - 1)) == 0;
}

vector<std::complex<double>> createComplex(vector<double> input)
{

    size_t length = checkpowerofTwo(length) ? input.size() : powerOfTwo(input.size());

    std::complex<double> c1(0.0, 0.0);

    vector<std::complex<double>> complexVector(length, c1);

    // cout<<length<<"ye lengtj"<<std::endl;
    // cout<<powerOfTwo(input.size())<<"ye lengtj"<<std::endl;

    for (size_t i = 0; i < length; i++)
    {

        // if the size of array is not in n^2 order then fill
        // the remaining space with zeroes.

        double value = i >= input.size() ? 0.0 : input[i];

        complexVector[i] = value;
    }

    return complexVector;
}

// data -> array that represents the vector fo complex samples
// reverse -> 1 to calculate FFT and -1 to calculate reverse FFT

void Fft::transform(vector<std::complex<double>> &data, bool inverse)
{

    size_t n = data.size();

    int sign = inverse ? -1 : 1;
    int levels = 0;
    for (size_t temp = n; temp > 1U; temp >>= 1)
        levels++;
    if (static_cast<size_t>(1U) << levels != n)
        throw std::domain_error("Length is not a power of 2");

    for (size_t j = 0; j < n; j++)
    {

        size_t x = bitReversal(j, levels);
        if (x > j)
        {

            std::swap(data[j], data[x]);
        }
    }

    for (size_t size = 2; size <= n; size <<= 1)
    {
        float expFactor = (-2.0 * M_PI / size) * sign;
        std::complex<double> angle(cos(expFactor), sin(expFactor));

        for (size_t s = 0; s < n; s += size)
        {
            std::complex<double> w(1);
            for (size_t j = 0; j < size / 2; j++)
            {
                std::complex<double> u = data[s + j], v = data[s + j + size / 2] * w;
                data[s + j] = u + v;
                data[s + j + size / 2] = u - v;
                w *= angle;
            }
        }
    }
}

void demoTester()
{

    vector<double> inputs;
    for (int i = 0; i < 20; i++)
    {
        std::uniform_real_distribution<double> range(-1.0, 1.0);
        inputs.push_back(range(randGen));
    }

    vector<std::complex<double>> v1 = createComplex(inputs);

    Fft::transform(v1, false);

    for (int i = 0; i < v1.size(); i++)
    {
        std::cout << v1[i];
    }
}

// int main()
// {

//     // std::cout << powerOfTwo(20);
//     // std::cout << checkpowerofTwo(3);

//     demoTester();
// }