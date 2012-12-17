/*
 *  general_utils.h
 *  opt-fast
 *
 *  john.mu@ieee.org
 *
 */

/* The MIT License

   Copyright (c) 2012 John C. Mu.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
 */


#ifndef GENERAL_UTILS_H

#define GENERAL_UTILS_H

#include "stl.h"


inline double logsumexp(double a, double b){
    double output = 0.0;

    if(a>b){
        output = a + log(1 + exp(b-a));
    }else{
        output = b + log(1 + exp(a-b));
    }

    return output;
}

inline double logsumexp(vector<double> a){
    double output = 0.0;

    double max_val = a[0];
    for(int i = 0;i<(int)a.size();i++){
        if(a[i]>max_val){
            max_val = a[i];
        }
    }

    output = max_val;


    double val = 0.0;

    for(int i = 0;i<(int)a.size();i++){
        val += exp(a[i]-max_val);
    }

    output += log(val);

    return output;
}

inline bool cut_region2(vector<vector<double> > &data,vector<vector<double> > &out0,
        vector<vector<double> > &out1, int dim, double lim) {

    double first = data[0][dim];
    bool all_same = true;

    for (int i = 0; i < (int) data.size(); i++) {

        if(all_same){
            if(fabs(first - data[i][dim]) > 1E-20) all_same = false;
        }

        if (data[i][dim] < lim) {
            out0.push_back(data[i]);
        }else if (data[i][dim] >= lim) {
            out1.push_back(data[i]);

        }



    }

    return !all_same;

}

inline bool cut_region(vector<vector<double> > &data,
        vector<vector<double> > &out ,int dim, int cut, double lim){

    double first = data[0][dim];
    bool all_same = true;

    for (int i = 0;i<(int)data.size();i++){

        if(all_same){
            if(fabs(first - data[i][dim]) > 1E-9) all_same = false;
        }

        if(cut == 0){
            if(data[i][dim] < lim) {
                out.push_back(data[i]);
            }

        }else if(cut == 1){
            if(data[i][dim] >= lim) {
                out.push_back(data[i]);
            }

        }else{
            cerr << "Wrong cut! " << cut << '\n';
            exit(2);
        }
    }


    return !all_same;

}


inline int count_region(vector<vector<double> > &data,
        int dim, int cut, double lim){
    int output = 0;

    for (int i = 0;i<(int)data.size();i++){


        if(cut == 0){
            if(data[i][dim] < lim) {
                output++;

            }

        }else if(cut == 1){
            if(data[i][dim] >= lim) {
                output++;

            }


        }else{
            cerr << "Wrong cut! " << cut << '\n';
            exit(2);
        }
    }

    return output;
}


inline void print_data(vector<vector<double> > &data){
    for(int i = 0;i<(int)data.size();i++){
        for(int  j=0;j<(int)data[i].size();j++){
            cerr << data[i][j] << ',';
        }
        cerr << '\n';
    }
}


//http://graphics.stanford.edu/~seander/bithacks.html
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
static const char LogTable256[256] =
    {
        -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
        LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
        LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
    };



class bit_str{
public:
    uint32_t data;
    static const int num_bits = sizeof(uint64_t)*8;



    bit_str(){
        data = 1u;
    }

    void clear(){
        data = 1u;
    }

    int size() const {

        int r; // r will be lg(v)
        register unsigned int t, tt; // temporaries

        tt = data >> 16;
        if (tt != 0) {
            r = (t = tt >> 8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
        } else {
            r = (t = data >> 8) ? 8 + LogTable256[t] : LogTable256[data];
        }

        return r;

    }

    bool push_back(int i){

        int s = size();

        if(s == 31){
            cerr << "Bit Str full!\n";
            return false;
        }

        if(i == 0){
            data = data & (~(((uint32_t)1)<<s));
        }else{
            data = data | (((uint32_t)1)<<s);
        }
        data = data | (((uint32_t)1)<<(s+1));

        return true;
    }

    void pop_back(){
        int s = size();


        if(s > 0){
            data = data & ((((uint32_t)1)<<s)-1);
            data = data | (((uint32_t)1)<<(s-1));
        }
    }

    bool operator==(const bit_str &b) const{
        return data == b.data;
    }

    bool operator!=(const bit_str &b) const{
        return data != b.data;
    }

    bool operator<(const bit_str &b) const{
        return data < b.data;
    }

    bool operator[](int i) const{
        return data & (1u<<i);
    }

    bool is_child(const bit_str &a) const{
        int self_size = size();
        int a_size = a.size();

        if (a_size < self_size){
            return false;
        }else if(a_size == self_size){
            return a.data == data;
        }else{

            uint32_t mask = (1u << self_size)-1;

            return ((data & mask) == (a.data & mask));
        }

    }

};



// Source: http://mlawire.blogspot.com/2009/07/c-whitespace-trimming-functions.html
inline void ltrim(string& str)
{
	string::size_type pos = 0;
	while (pos < str.size() && (isspace(str[pos]))) pos++;
	str.erase(0, pos);
}
inline void rtrim(string& str)
{
	string::size_type pos = str.size();
	while (pos > 0 && (isspace(str[pos - 1]))) pos--;
	str.erase(pos);
}
inline void trim2(string& str)
{
	ltrim(str);
	rtrim(str);
}

inline vector<string> split(const string &s, char delim)
{
    vector<string> elems;
	stringstream ss(s);
    string item;
    while(getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

// Split using all white space as delims
inline vector<string> split(const string &s)
{
    vector<string> elems;
	stringstream ss(s);
    string item;
    while(ss >> item) {
        elems.push_back(item);
    }
    return elems;
}


template <typename T> string toStr(T tmp) {
    ostringstream out;
    out << tmp;
    return out.str();
}

template <typename T> T strTo(string tmp) {
    T output;
    istringstream in(tmp);
    in >> output;
    return output;
}


// http://www.johndcook.com/cpp_phi.html

inline double phi(double x) {
    // constants
    double a1 = 0.254829592;
    double a2 = -0.284496736;
    double a3 = 1.421413741;
    double a4 = -1.453152027;
    double a5 = 1.061405429;
    double p = 0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x) / sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

    return 0.5 * (1.0 + sign * y);
}


// from knuth

inline uint64_t combination(uint32_t n, uint32_t k) {
    if (k > n) {
        return 0;
    }

    uint64_t r = 1;

    for (uint32_t d = 1; d <= k; ++d) {
        r *= n--;
        r /= d;
    }

    return r;
}

inline uint64_t power(unsigned int n, unsigned int m) {
    uint64_t output = 1;

    if (n == 0 || n == 1) {
        return 1u;
    }

    for (; m > 0; m--) {
        output *= n;
    }

    return output;
}



// http://www.johndcook.com/blog/2008/04/24/how-to-calculate-binomial-probabilities/

inline double binom_p(double p, double q, int m, int n) {
    double temp = lgamma(m + n + 1.0);
    temp -= lgamma(n + 1.0) + lgamma(m + 1.0);
    temp += m * log(p) + n * log(q);
    return exp(temp);
}

inline double binom_prob(uint32_t x, uint32_t n, double p) {
    return binom_p(p, 1 - p, x, n - x);
}


// MT code taken from
// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c
// Ported to C++ by John Mu

class MT_random {
private:
    static const uint64_t NN = 312;
    static const uint64_t MM = 156;
    static const uint64_t MATRIX_A = 0xB5026F5AA96619E9ULL;
    static const uint64_t UM = 0xFFFFFFFF80000000ULL; /* Most significant 33 bits */
    static const uint64_t LM = 0x7FFFFFFFULL; /* Least significant 31 bits */




    /* The array for the state vector */
    uint64_t mt[NN];
    /* mti==NN+1 means mt[NN] is not initialized */
    uint64_t mti;

    void init_all() {
        mti = NN + 1;
    }

public:

    MT_random() {

        init_all();
        //init_genrand64(time(NULL));
        init_genrand64(0);
    }

    MT_random(unsigned long long seed) {

        init_all();
        init_genrand64(seed);
    }

    /* initializes mt[NN] with a seed */
    void init_genrand64(unsigned long long seed) {
        mt[0] = seed;
        for (mti = 1; mti < NN; mti++)
            mt[mti] = (6364136223846793005ULL * (mt[mti - 1] ^ (mt[mti - 1] >> 62)) + mti);
    }

    /* initialize by an array with array-length */
    /* init_key is the array for initializing keys */

    /* key_length is its length */
    void init_by_array64(unsigned long long init_key[],
            unsigned long long key_length) {
        unsigned long long i, j, k;
        init_genrand64(19650218ULL);
        i = 1;
        j = 0;
        k = (NN > key_length ? NN : key_length);
        for (; k; k--) {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 62)) * 3935559000370003845ULL))
                    + init_key[j] + j; /* non linear */
            i++;
            j++;
            if (i >= NN) {
                mt[0] = mt[NN - 1];
                i = 1;
            }
            if (j >= key_length) j = 0;
        }
        for (k = NN - 1; k; k--) {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 62)) * 2862933555777941757ULL))
                    - i; /* non linear */
            i++;
            if (i >= NN) {
                mt[0] = mt[NN - 1];
                i = 1;
            }
        }

        mt[0] = 1ULL << 63; /* MSB is 1; assuring non-zero initial array */
    }

    /* generates a random number on [0, 2^64-1]-interval */
    unsigned long long genrand64_int64(void) {
        uint64_t i;
        uint64_t x;
        static uint64_t mag01[2] = {0ULL, MATRIX_A};

        if (mti >= NN) { /* generate NN words at one time */

            /* if init_genrand64() has not been called, */
            /* a default initial seed is used     */
            if (mti == NN + 1)
                init_genrand64(5489ULL);

            for (i = 0; i < NN - MM; i++) {
                x = (mt[i] & UM) | (mt[i + 1] & LM);
                mt[i] = mt[i + MM] ^ (x >> 1) ^ mag01[(int) (x & 1ULL)];
            }
            for (; i < NN - 1; i++) {
                x = (mt[i] & UM) | (mt[i + 1] & LM);
                mt[i] = mt[i + (MM - NN)] ^ (x >> 1) ^ mag01[(int) (x & 1ULL)];
            }
            x = (mt[NN - 1] & UM) | (mt[0] & LM);
            mt[NN - 1] = mt[MM - 1] ^ (x >> 1) ^ mag01[(int) (x & 1ULL)];

            mti = 0;
        }

        x = mt[mti++];

        x ^= (x >> 29) & 0x5555555555555555ULL;
        x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
        x ^= (x << 37) & 0xFFF7EEE000000000ULL;
        x ^= (x >> 43);

        return x;
    }

    /* generates a random number on [0, 2^63-1]-interval */
    long long genrand64_int63(void) {
        return (long long) (genrand64_int64() >> 1);
    }

    /* generates a random number on [0,1]-real-interval */
    double genrand64_real1(void) {
        return (genrand64_int64() >> 11) * (1.0 / 9007199254740991.0);
    }

    /* generates a random number on [0,1)-real-interval */
    double genrand64_real2(void) {
        return (genrand64_int64() >> 11) * (1.0 / 9007199254740992.0);
    }

    /* generates a random number on (0,1)-real-interval */
    double genrand64_real3(void) {
        return ((genrand64_int64() >> 12) + 0.5) * (1.0 / 4503599627370496.0);
    }

    // this is not the "correct" way... but good enough

    double genrand_norm(double mean, double sd) {
        return (genrand_norm() * sd) +mean;
    }

    double genrand_norm() {
        double sum = 0;

        for (int i = 0; i < 48; i++) {
            sum = sum + (double) (rand() + 1.0) / (RAND_MAX + 1.0);
        }

        // only for 48
        return (sum - 24) / 2;
    }

    double genrand_exp(double mean) {
        return -mean * log(1 - genrand64_real2());
    }

    int64_t genrand_binom(uint64_t n, double p) {

        if (p < 0 || p > 1) {
            return -1;
        }

        int64_t total = 0;
        for (uint64_t i = 0; i < n; i++) {
            total += (genrand64_real1() < p ? 1 : 0);
        }

        return total;
    }

    // generate random int inclusive of range [start,end]

    uint64_t genrand_int_range(uint64_t start, uint64_t end) {

        if(start == end){
            return start;
        }

        return (uint64_t) ((genrand64_real2()*(end - start + 1)) + start);
    }

    // bernouli p is true

    bool genrand_bern(double p) {
        return genrand64_real2() <= p;
    }



};



// This class is specific to linux :(

class mu_timer {
private:
    struct timeval start_tv;
    bool enabled;
public:

    mu_timer() {
        gettimeofday(&start_tv, NULL);
        enabled = true;
    }

    void set_enabled(bool e){
        enabled = e;
    }

    void reset() {
        if (enabled) {
            gettimeofday(&start_tv, NULL);
        }
    }

    double elapsed_time() {
        if (enabled) {
            struct timeval tv;
            gettimeofday(&tv, NULL);

            return (double) (tv.tv_sec - start_tv.tv_sec) + (tv.tv_usec - start_tv.tv_usec) / 1000000.0;
        } else {
            return 0;
        }
    }

    void print_elapsed_time(ostream &o, string s) {
        if (enabled) {
            o << s << " time: " << elapsed_time() << " s." << '\n';
        }
    }
};







static __inline__ unsigned long GetCC(void) {
    unsigned a, d;
    asm volatile("rdtsc" : "=a" (a), "=d" (d));
    return ((unsigned long) a) | (((unsigned long) d) << 32);
}



class seen_t {
private:
    uint32_t* locs; // hash table
    vector<uint32_t> locs_added;

public:
    static const uint32_t hash_size = 1048576;
    static const uint32_t shift_bits = 32 - 20;
    static const uint32_t magic = 2654435761u;

    //static const uint32_t mask = 0xFFFFF;

    seen_t() {
        locs = new uint32_t[hash_size];
        memset(locs, hash_size, sizeof (uint32_t));
        locs_added.reserve(5000);
    }

    ~seen_t() {
        clear();
        delete [] locs;
    }

    // true if inserted

    bool insert(uint32_t a) {
        uint32_t hash_val = (a * magic) >> shift_bits;


        if (locs[hash_val] == 0) {
            locs[hash_val] = a;
            locs_added.push_back(hash_val);
            return true;
        } else if (locs[hash_val] == a) {
            return false;
        } else {
            hash_val++;

            uint32_t num_seen = 1;
            while (num_seen < hash_size) {
                uint32_t hash_val_temp = hash_val;
                if (hash_val >= hash_size) {
                    hash_val_temp = hash_val_temp - hash_size;
                }


                if (locs[hash_val_temp] == 0) {
                    locs[hash_val_temp] = a;
                    locs_added.push_back(hash_val_temp);
                    return true;
                } else if (locs[hash_val_temp] == a) {
                    return false;
                }

                hash_val++;
                num_seen++;

            }

            return true; // hash is full :( ... should not happen
        }
    }

    // true if found

    bool find(uint32_t a) {
        uint32_t hash_val = (a * magic) >> shift_bits;


        if (locs[hash_val] == 0) {
            return false;
        } else if (locs[hash_val] == a) {
            return true;
        } else {
            hash_val++;

            uint32_t num_seen = 1;
            while (num_seen < hash_size) {
                uint32_t hash_val_temp = hash_val;
                if (hash_val >= hash_size) {
                    hash_val_temp = hash_val_temp - hash_size;
                }


                if (locs[hash_val_temp] == 0) {
                    return false;
                } else if (locs[hash_val_temp] == a) {
                    return true;
                }

                hash_val++;
                num_seen++;

            }

            return true; // hash is full :( ... should not happen
        }
    }

    void clear() {
        vector<uint32_t>::iterator added_it = locs_added.begin();
        while (added_it != locs_added.end()) {

            locs[*added_it] = 0;

            added_it++;
        }

        locs_added.clear();
    }

};




#endif




