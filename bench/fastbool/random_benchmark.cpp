#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <bitset>

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

// Unbiased, fastest variant. Gets rid of the counter by sacrificing 1 bit of randomness.
template <typename U = uint64_t> class RandomizerWithShiftT {
    public:
        template <typename Rng> bool operator()(Rng &rng) {
            if (UNLIKELY(1 == m_rand)) {
                m_rand = std::uniform_int_distribution<U>{}(rng) | s_mask_left1;
            }
            bool const ret = m_rand & 1;
            m_rand >>= 1;
            return ret;
        }

    private:
        static constexpr const U s_mask_left1 = U(1) << (sizeof(U) * 8 - 1);
        U m_rand = 1;
};
//*
template <typename U = uint64_t> class RandomizerInt {
    public :
        RandomizerInt(int size) : m_size(size) {};

        template <typename rand, typename Rng> U operator()(rand &boolrand, Rng &rng){
            int res = 0;
            int a ;
            for(int i = 0 ; i<m_size;i++){
                a = boolrand(rng);
                res += a << i ;
                std::cout << i << " " << a << " " << std::bitset<8>(res) <<"\n";
            }

            std::cout << "Final : " << std::bitset<8>(res) << "\n";
            return res;
        }
    private :
        int m_size;
};
//*/

int main(int, char **) {
    size_t iters = UINT64_C(100000000);
    sfc64 generator;
    RandomizerWithShiftT<uint64_t> random;
    RandomizerInt<int> rand4(5);

    std::cout <<  std::bitset<8>(random(generator)) << "\n";
    std::cout <<  std::bitset<8>(rand4(random,generator)) << " " << rand4(random,generator) << "\n";
    
    std::cout << (2<<-1) << " " << (2<<1) << " " << (2<<2) << " " << (2<<3) << " " << (2<<4) << "\n";
}

