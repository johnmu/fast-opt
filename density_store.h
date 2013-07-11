/* 
 * File:   density_store.h
 * Author: johnmu
 *
 * Created on June 28, 2013, 3:33 PM
 */

#ifndef DENSITY_STORE_H
#define	DENSITY_STORE_H

#include "stl.h"
#include "general_utils.h"
#include "opt_tree.h"
#include "llopt_tree.h"
#include "disopt_tree.h"
#include "lsopt_tree.h"
#include "dfopt.h"
#include "map_tree.h"


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

    if (magic_val != c::magic) {
        cerr << "ERROR: Filetype mismatch\n";
        den_file.close();
        return;
    }

    den_file.read((char*) &num_den, sizeof (num_den));

    //cerr<< "num_den: " << num_den << '\n';

}



// This class stores densities that could be copulas

class density_store {
private:
    bool copula;

    map_tree map_region_tree;
    opt_region_hash<uint32_t> map_regions;

    map_tree** m_map_region_tree;
    opt_region_hash<uint32_t>** m_map_regions;
    cdf** marginal;

    int num_dim;
    int num_regions;
    vector<int> m_num_regions;

    void init() {
        m_map_region_tree = NULL;
        m_map_regions = NULL;
        marginal = NULL;
        num_dim = 0;
        num_regions = 0;
        copula = false;
    }

    int load_densities(string joint_filename, string marginal_filename) {

        { // block
            int num_dim = 0;
            ifstream den_file;

            init_file_in(den_file, joint_filename, num_dim);

            if (num_dim != 1) {
                cerr << "More than one density in joint file\n";
                return 1;
            }

            map_region_tree.load(den_file);
            map_regions.load(den_file);

            den_file.close();
            
            num_regions = map_regions.get_num_regions();
        }

        // load the marginal densities
        int num_dim = 0;
        if (copula) {
            cerr << "Loading Marginals...\n";

            ifstream den_file;

            init_file_in(den_file, marginal_filename, num_dim);

            if (num_dim != map_region_tree.get_num_children()) {
                cerr << "marginal densities num not consistent with joint\n";
                return 1;
            }

            m_map_region_tree = new map_tree*[num_dim];
            m_map_regions = new opt_region_hash<uint32_t>*[num_dim];

            for (int i = 0; i < num_dim; i++) {
                // load each dimension
                m_map_region_tree[i] = new map_tree(0, num_dim);
                m_map_regions[i] = new opt_region_hash<uint32_t>(2);

                m_map_region_tree[i]->load(den_file);
                m_map_regions[i]->load(den_file);
                
                m_num_regions.push_back(m_map_regions[i]->get_num_regions());
            }

            den_file.close();

        }
        
        

        // create the CDFs
        if (copula) {
            marginal = new cdf*[num_dim];
            for (int i = 0; i < num_dim; i++) {
                marginal[i] = new cdf(*(m_map_region_tree[i]), *(m_map_regions[i]));
            }
        }

        return 0;
    }

    void destroy() {
        for (int i = 0; i < num_dim; i++) {
            if (m_map_region_tree != NULL) delete m_map_region_tree[i];
            if (m_map_regions != NULL) delete m_map_regions[i];
            if (marginal != NULL) delete marginal[i];
        }

        if (m_map_region_tree != NULL) delete [] m_map_region_tree;
        if (m_map_regions != NULL) delete [] m_map_regions;
        if (marginal != NULL) delete [] marginal;
    }

public:

    density_store() {
        init();
    }

    density_store(string joint_filename, string marginal_filename, bool copula) {
        init();
        this->copula = copula;

        int status = load_densities(joint_filename, marginal_filename);

        if (status != 0) {
            cerr << "ERROR: Failed loading of file: " << joint_filename << '\n';
        }

        num_dim = map_region_tree.get_num_children();


    }

    int load_files(string joint_filename, string marginal_filename, bool copula) {
        destroy();
        init();
        this->copula = copula;

        int status = load_densities(joint_filename, marginal_filename);

        if (status != 0) {
            cerr << "ERROR: Failed loading of file: " << joint_filename << '\n';
            return status;
        }

        num_dim = map_region_tree.get_num_children();
        return 0;
    }

    double compute_density(const vector<double> &point, double pseduo_count) {

        if (num_dim != (int) point.size()) {
            cerr << "Error: Compute density dimension mismatch\n";
            return -c::inf;
        }

        int rep = 1;
        double total_density = 0.0;

        for (int x = 0; x < rep; x++) {
            
            vector<double> pet_point = point;
            
            double dist = 0.0;
            
            if(x>0){
                for (int i = 0; i < num_dim; i++) {
                    pet_point[i] += (0.1*(rand()/(double)RAND_MAX))-0.05;
                    
                    double val = (pet_point[i]- point[i]);
                    dist += val*val;
                }
                
                dist /= num_dim;
            }
            
            double log_density = 0;
            // product of the marginal densities
            if (copula) {
                for (int i = 0; i < num_dim; i++) {
                    vector<double> single_point;
                    single_point.push_back(pet_point[i]);

                    double curr_density = m_map_region_tree[i]->get_density(single_point,pseduo_count,m_num_regions[i]);

                    log_density += log(curr_density);
                }
            }

            // copula transform the data point
            vector<double> trans_data = pet_point;
            if (copula) {
                for (int i = 0; i < num_dim; i++) {
                    trans_data[i] = marginal[i]->transform(trans_data[i]);
                }
            }

            double joint_density = map_region_tree.get_density(trans_data,pseduo_count,num_regions);

            log_density += log(joint_density);

            total_density += exp(log_density)*exp(-dist);
        }

        return total_density/rep;
    }

    ~density_store() {
        destroy();
    }

};

#endif	/* DENSITY_STORE_H */

