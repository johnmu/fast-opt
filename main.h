/*
 *  main.h
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

// num_den is number of densities stored in this file
void init_file_out(ofstream &den_file, string file_name, int num_den) {

    string temp = file_name;
    den_file.open(temp.c_str(), ios::out | ios::binary);

    // write den file headers
    den_file.write((char*) &c::magic, sizeof (c::magic));
    den_file.write((char*) &num_den, sizeof (num_den));

}

// num_den is number of densities stored in this file
void init_file_in(ifstream &den_file, string file_name, int &num_den) {

    string temp = file_name;
    den_file.open(temp.c_str(), ios::in | ios::binary);

    // write den file headers
    int magic_val = 0;
    den_file.read((char*) &magic_val, sizeof (magic_val));
    
    if(magic_val != c::magic){
        cerr << "ERROR: Filetype mismatch\n";
        den_file.close();
        return;
    }
    
    den_file.read((char*) &num_den, sizeof (num_den));
    
    cerr<< "num_den: " << num_den << '\n';

}

double compute_density(vector<double> &point, map_tree *map_region_tree, 
        map_tree** m_map_region_tree, cdf** marginal, bool copula) {

    int num_dim = point.size();
    
    double log_density = 0;
    // product of the marginal densities
    if (copula) {
        for (int i = 0; i < num_dim; i++) {
            vector<double> single_point;
            single_point.push_back(point[i]);

            log_density += log(m_map_region_tree[i]->get_density(single_point));
        }
    }

    // copula transform the data point
    vector<double> trans_data = point;
    if (copula) {
        for (int i = 0; i < num_dim; i++) {
            trans_data[i] = marginal[i]->transform(trans_data[i]);
        }
    }

    log_density += map_region_tree->get_density(trans_data);

    return exp(log_density);
    
}

int load_densities(string joint_filename, string marginal_filename,
        map_tree *map_region_tree, opt_region_hash<uint32_t> *map_regions,
        map_tree** m_map_region_tree, opt_region_hash<uint32_t>** m_map_regions,
        cdf** marginal, bool copula){
    
    int num_dim = 0;
    {
        ifstream den_file;
        
        init_file_in(den_file, joint_filename, num_dim);

        if (num_dim != 1) {
            cerr << "More than one density in joint file\n";
            return 1;
        }

        map_region_tree->load(den_file);
        
        cerr << "DONE map_region_tree\n";
        
        map_regions->load(den_file);
        
        cerr << "DONE map_regions\n";

        den_file.close();
    }
    
    // load the marginal densities
    if(copula){
        cerr << "Loading Marginals...\n";
        
        ifstream den_file;
        
        init_file_in(den_file, marginal_filename, num_dim);

        if (num_dim != map_region_tree->get_num_children()) {
            cerr << "marginal densities num not consistent with joint\n";
            return 1;
        }
        
        m_map_region_tree = new map_tree*[num_dim];
        m_map_regions = new opt_region_hash<uint32_t>*[num_dim];
        
        for(int i = 0;i<num_dim;i++){
            // load each dimension
            m_map_region_tree[i] = new map_tree(0);
            m_map_regions[i] = new opt_region_hash<uint32_t>(2);
            
            m_map_region_tree[i]->load(den_file);
            m_map_regions[i]->load(den_file);
        }
        
        den_file.close();
        
    }
    
    // create the CDFs
    marginal = new cdf*[num_dim];
    if (copula) {
        for (int i = 0; i < num_dim; i++) {
            marginal[i] = new cdf(*(m_map_region_tree[i]), *(m_map_regions[i]));
        }
    }
    
    return 0;
}


int llopt(vector<string> params);
int lsopt(vector<string> params);
int sopt(vector<string> params);
int opt(vector<string> params);
int dfopt(vector<string> params);
int disopt(vector<string> params);
int copula(vector<string> params);
int hell_dist(vector<string> params);
int classify(vector<string> params);
int density(vector<string> params);
int density_old(vector<string> params);
int bench(vector<string> params);

void print_usage_and_exit();

#endif



