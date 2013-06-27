/*
 *  stl.h
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



#ifndef STL_H

#define STL_H

#include <cstddef>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <ctime>
#include <map>
#include <sstream>
#include <set>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <assert.h>
#include <ctime>
#include <cctype>
#include <climits>
#include <stdint.h>
#include <bitset>
#include <queue>
#include <deque>
#include <list>
#include <queue>
#include <zlib.h>



// linux specific stuff
#include <sys/time.h>
#include <pthread.h>

// oh well
using namespace std;

namespace c{
    static const int cuts = 2; // A cut splits it into 2, no triary cuts for now
    static const double l2  = log(2);
    static const double pi  = 3.141592653589793238462643383279502884;
    static const double lpi = log(3.141592653589793238462643383279502884);
    static const double inf = INFINITY;
    static const uint32_t ra_null_val = 4294967295u;
    static const string PROG_NAME = "opt-fast";
    static const string BUILD = "0.5-r60";
    static const int MAX_COORD_DEPTH = 30;
    
    static const int magic = 1741324111; // for files created by this program
}


#endif
