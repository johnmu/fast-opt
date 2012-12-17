/*
 *  gamma_table.h
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



#ifndef GAMMA_TABLE_H
#define	GAMMA_TABLE_H

#include "stl.h"

class gamma_table{
private:
    vector<double> lgamma_vals_half;
    vector<double> lgamma_vals_one;
    int max;
public:
    gamma_table(int max){
        this->max = max;
        lgamma_vals_half.resize(max+1);
        lgamma_vals_one.resize(max+1);

        for(int i = 0;i<=max;i++){
            lgamma_vals_half[i] = lgamma(i+0.5);
            lgamma_vals_one[i]  = lgamma(i+1);
        }
    }

    // this is only for a special case of D
    double compute_lD2(int n, int n_j1, int n_j2){
        return lgamma_vals_half[n_j1] + lgamma_vals_half[n_j2] - lgamma_vals_one[n];
    }





};


#endif	/* GAMMA_TABLE_H */

