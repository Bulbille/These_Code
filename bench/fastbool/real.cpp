#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>

#define LIKELY(x) __builtin_expect((x), 1)
#define UNLIKELY(x) __builtin_expect((x), 0)
#define NO_INLINE __attribute__((noinline))

// original source : 
// https://martin.ankerl.com/2018/12/08/fast-random-bool/?fbclid=IwAR2r38PZymhkLk_s4JCnXBy1Gi4Fa3nPTnFsekrYnhjSdX2qpDDc_HiCWL8

// extremely fast random number generator that also produces very high quality random.
// see PractRand: http://pracrand.sourceforge.net/PractRand.txt
class sfc64 {
    public:
        using result_type = uint64_t;

        static constexpr uint64_t(min)() { return 0; }
        static constexpr uint64_t(max)() { return UINT64_C(-1); }

        sfc64() : sfc64(std::random_device{}()) {}

        explicit sfc64(uint64_t seed) : m_a(seed), m_b(seed), m_c(seed), m_counter(1) {
            for (int i = 0; i < 12; ++i) {
                operator()();
            }
        }

        uint64_t operator()() noexcept {
            auto const tmp = m_a + m_b + m_counter++;
            m_a = m_b ^ (m_b >> right_shift);
            m_b = m_c + (m_c << left_shift);
            m_c = rotl(m_c, rotation) + tmp;
            return tmp;
        }

    private:
        template <typename T> T rotl(T const x, int k) { return (x << k) | (x >> (8 * sizeof(T) - k)); }

        static constexpr int rotation = 24;
        static constexpr int right_shift = 11;
        static constexpr int left_shift = 3;
        uint64_t m_a;
        uint64_t m_b;
        uint64_t m_c;
        uint64_t m_counter;
};


int main()
{
    const int nrolls=1e8;  // number of experiments
    const int nstars=1e4;     // maximum number of stars to distribute
    const int nintervals=1000; // number of intervals

    sfc64 generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    int p[nintervals]={};

    for (int i=0; i<nrolls; ++i) {
        double number = distribution(generator);
        ++p[int(nintervals*number)];
    }

    std::cout << "uniform_real_distribution (0.0,1.0):" << std::endl;
    std::cout << std::fixed; std::cout.precision(1);

    for (int i=0; i<nintervals; ++i) {
        std::cout << float(i)/nintervals << "-" << float(i+1)/nintervals << ": ";
        std::cout << std::string(p[i]*nstars/nrolls,'*') << std::endl;
    }

    return 0;
}

