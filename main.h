/*
 *  main.h
 *  fast-opt
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



#ifndef MAIN_H

#define MAIN_H



#include "stl.h"
#include "general_utils.h"
#include "opt_tree.h"
#include "llopt_tree.h"
#include "disopt_tree.h"
#include "lsopt_tree.h"
#include "dfopt.h"
#include "map_tree.h"
#include "density_store.h"

inline vector<vector<double> > read_data(string filename, bool end_line) {
    int dim = 0;

    ifstream infile(filename.c_str());

    if (!infile.is_open()) {
        cerr << "ERROR: Could not open " << filename << '\n';
        exit(1);
    }

    string line;
    getline(infile, line);

    trim2(line);

    if(line.length() == 0){
        cerr << "ERROR: Empty file: " << filename << '\n';
        exit(1);
    }

    vector<string> line_list = split(line);

    dim = line_list.size();

    infile.close();

    infile.open(filename.c_str());

    vector<vector<double> > data;

    while (!infile.eof()) {

        getline(infile, line);
        trim2(line);

        if (line.length() == 0) continue;

        vector<string> ll = split(line);

        if ((int) ll.size() != dim) {
            cerr << "ERROR: Bad line dim: " << line << '\n';
            exit(2);
        }

        vector<double> d;
        d.resize(dim);

        bool good_data = true;
        for (int i = 0; i < dim; i++) {
            d[i] = strTo<double>(ll[i]);
            if (end_line) {
                if ((i != (dim - 1)) && ((d[i] > 1.0) || (d[i] < 0))) {
                    good_data = false;
                    break;
                }
            } else {
                if (d[i] > 1.0 || d[i] < 0 ) {
                    good_data = false;
                    break;
                }
            }
        }
        if(good_data){
            data.push_back(d);
        }else{
            cerr << "Warning: Data out of range("<< dim <<"): " << line << '\n';
        }

    }


    return data;
}



int llopt(vector<string> params);
int lsopt(vector<string> params);
int sopt(vector<string> params);
int opt(vector<string> params);
int opt_comp(vector<string> params);
int dfopt(vector<string> params);
int disopt(vector<string> params);
int copula(vector<string> params);
int hell_dist(vector<string> params);
int classify(vector<string> params);
int density(vector<string> params);
int print_partitions(vector<string> params);

int density_old(vector<string> params);
int bench(vector<string> params);

void print_usage_and_exit();

#endif



