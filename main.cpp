/*
 *  main.cpp
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


#include "main.h"


int main(int argc, char** argv) {

    cerr << "**" + c::PROG_NAME << '\n';
    cerr << "Version: " << c::BUILD << "\n";
    cerr << "John C. Mu\n";
    cerr << "http://www.stanford.edu/group/wonglab" << "\n";
    cerr << "\n";



    string mode = "";
    string temp = "";
    int error_num = 0;
    vector<string> params;

    if (argc < 2) {
        print_usage_and_exit();
    }

    mode = argv[1];

    for (int i = 2; i < argc; i++) {
        temp = argv[i];
        params.push_back(temp);
    }


    if (mode == "opt") {
        error_num = opt(params);
    } else if (mode == "llopt") {
        error_num = llopt(params);
    } else if (mode == "lsopt") {
        error_num = lsopt(params);
    } else if (mode == "dfopt") {
        error_num = dfopt(params);
    } else if (mode == "disopt") {
        error_num = disopt(params);
    } else if (mode == "hell_dist") {
        error_num = hell_dist(params);
    } else if (mode == "classify") {
        error_num = classify(params);
    } else if (mode == "density") {
        error_num = density(params);
    } else if (mode == "bench") {
        error_num = bench(params);
    } else {
        print_usage_and_exit();
    }

    return error_num;

}





int opt(vector<string> params) {

    string usage_text = "Usage: " + c::PROG_NAME + " opt <percent_points> <data_file> <output_name>\n"
            + "       percent_points -- Ratio of total data to stop at (0.01 = 1%)\n"
            + "            data_file -- One sample each row (Restricted to [0,1] cube)\n"
            + "MAP partitions output to STDOUT and <output_name>.den .Log to STDERR \n"
            + "Run full-OPT, very fast but uses a lot of memory. If dimension \n"
            + "greater than 5 use the other methods\n";

    if (params.size() != 3) {
        cerr << usage_text << endl;
        return 3;
    }
    
    string out_filename = params[2];



    vector<vector<double> > data = read_data(params[1],false);

    int dim = (int)data[0].size();
    int N   = data.size();

    cerr << "Running Full-OPT" << '\n';
    cerr << N << " points in " << dim << " dimensions.\n";

    double stop_ratio = strTo<double>(params[0]);
    int stop_points = 0;
    
    if(stop_ratio<1){
        stop_points = (int)(N*stop_ratio);
    }else{
        stop_points = (int)stop_ratio;
    }
    if(stop_points < 1) stop_points = 1;

    cerr << "Stopping at " << stop_points << " points.\n";

    double total_time = 0.0;
    
    mu_timer mt;

    opt_tree opt_slow(dim,stop_points,1000);

    mt.reset();
    opt_slow.construct_full_tree(data);
    total_time += mt.elapsed_time();
    mt.print_elapsed_time(cerr, "OPT tree");

    mt.reset();
    map_tree map_region_tree(data.size());
    opt_region_hash<uint32_t> map_regions(24);
    opt_slow.construct_MAP_tree(map_region_tree, map_regions, data.size());
    total_time += mt.elapsed_time();
    mt.print_elapsed_time(cerr, "MAP tree");

    cerr << "OPT construction time: " << total_time << " s.\n";
    
    cerr << "lPhi: " << opt_slow.get_lphi() << endl;
    print_MAP_density(cout, map_regions.print_density(),map_region_tree.get_ra(),data.size());

    // write out the density to file
    ofstream den_file;
    init_file_out(den_file,out_filename+".den",1);
    
    map_region_tree.save(den_file);
    map_regions.save(den_file);
    
    den_file.close();

    return 0;
}


int llopt(vector<string> params) {


    string usage_text = "Usage: " + c::PROG_NAME + " llopt [-np] <percent_points> <levels> <top_percent_points> <data_file>\n"
            + "       -np            -- Don't prune the tree, a lot faster but more memory usage\n"
            + "       percent_points -- Ratio of total region data to stop at for each look-ahead (0.01 = 1%)\n"
            + "   top_percent_points -- Ratio of total data to stop at (0.01 = 1% or 2 = 2 points)\n"
            + "               levels -- Maximum levels for each look-ahead\n "
            + "            data_file -- One sample each row (Restricted to [0,1] cube)\n"
            + "MAP partitions output to STDOUT and <output_name>.den .Log to STDERR \n"
            + "Run limited look-ahead OPT. Good for moderate dimensions (6-15)\n";

    if (params.size() < 4 || params.size() > 5) {
        cerr << usage_text << endl;
        return 3;
    }

    bool prune_tree = true;
    int params_offset = 0;
    if (params.size() == 5) {
        if (params[0] == "-np") {
            prune_tree = false;
            params_offset = 1;
        } else {
            cerr << usage_text << endl;
            return 3;
        }
    }


    vector<vector<double> > data = read_data(params[params_offset + 3],false);

    int dim = (int)data[0].size();
    int N   = data.size();

    cerr << "Running limited look-ahead OPT" << '\n';
    cerr << N << " points in " << dim << " dimensions.\n";

    int ll_levels = strTo<int>(params[params_offset + 1]);
    double stop_ratio = strTo<double>(params[params_offset + 0]);
    double top_stop_ratio = strTo<double>(params[params_offset + 2]);

    cerr << "Each small OPT stopping at " << ll_levels << " levels.\n";
    cerr << "Each level stopping at " << stop_ratio*100 << "% of points.\n";
    cerr << "Total stopping at " << top_stop_ratio*100 << "% of points.\n";


    mu_timer mt;

    llopt_tree llopt(dim,stop_ratio,top_stop_ratio,ll_levels,20*dim);
    map_tree map_region_tree(data.size());
    opt_region_hash<uint32_t> map_regions(24);

    mt.reset();
    llopt.construct_llopt_tree(&data, map_region_tree, map_regions,prune_tree);
    mt.print_elapsed_time(cerr, "LLOPT("+ toStr<int>(ll_levels) +") construction");

    print_MAP_density(cout, map_regions.print_density(),map_region_tree.get_ra(),data.size());

    return 0;
}


int lsopt(vector<string> params) {


    string usage_text = "Usage: " + c::PROG_NAME + " lsopt <percent_points> <lookahead_depth> <top_percent_points> <iterations> <convergence_iterations> <data_file>\n"
            + "         percent_points -- Ratio of total region data to stop at for each look-ahead (0.01 = 1%)\n"
            + "        lookahead_depth -- Depth to stop at per look-ahead\n"
            + "     top_percent_points -- Ratio of total data to stop at (0.01 = 1%)\n"
            + "             iterations -- Maximum number of iterations for each look-ahead\n "
            + " convergence_iterations -- Number of iterations the same before convergence\n "
            + "              data_file -- One sample each row\n"
            + "MAP partitions output to STDOUT. Log to STDERR \n"
            + "Run limited look-ahead sampled-OPT. Good for higher dimensions (>10).\n"
            + "iterations = [50-1000], convergence_iterations=[20-100]";


    if (params.size() != 6) {
        cerr << usage_text << endl;
        return 3;
    }

    int iterations = strTo<int>(params[3]);
    int convergence_iterations = strTo<int>(params[4]);
    int lookahead_depth = strTo<int>(params[1]);

    vector<vector<double> > data = read_data(params[5],false);

    int dim = (int)data[0].size();
    int N   = data.size();

    cerr << "Running limited look-ahead sampled-OPT" << '\n';
    cerr << "Max iterations: "<< iterations << '\n';
    cerr << "Iterations to converge: "<< convergence_iterations << '\n';
    cerr << N << " points in " << dim << " dimensions.\n";

    double stop_ratio = strTo<double>(params[0]);
    double top_stop_ratio = strTo<double>(params[2]);

    cerr << "Each level stopping at " << stop_ratio*100 << "% of points.\n";
    cerr << "Total stopping at " << top_stop_ratio*100 << "% of points.\n";


    mu_timer mt;

    lsopt_tree lsopt(dim,stop_ratio,top_stop_ratio,lookahead_depth,1000);
    map_tree map_region_tree(data.size());
    opt_region_hash<uint32_t> map_regions(24);

    mt.reset();
    lsopt.construct_lsopt_tree(data, iterations,convergence_iterations, map_region_tree, map_regions);
    mt.print_elapsed_time(cerr, "LSOPT construction");

    print_MAP_density(cout, map_regions.print_density(),map_region_tree.get_ra(),data.size());

    return 0;
}

int dfopt(vector<string> params) {

    string usage_text = "Usage: " + c::PROG_NAME + " dfopt <top_percent_points>" 
            + " <percent_points> <levels> <data_file> <output_name>\n"
            + "     top_percent_points -- Points to stop at (0.01 = 1%, 1 = 1point)\n"
            + "     percent_points     -- Points in each level to stop at (0.01 = 1%, 1 = 1point))\n"
            + "     levels             -- look ahead levels\n"
            + "     data_file          -- One sample per row (Restricted to [0,1] cube)\n"
            + "MAP partitions output to STDOUT and <output_name>.den .Log to STDERR \n"
            + "Run depth-first-OPT, very slow but technically can do any number of dimensions\n "
            + "if you wait several years. Run in a similar fashion to limited look-ahead. \n";

    if (params.size() != 5) {
        cerr << usage_text << endl;
        return 3;
    }
    
    string out_filename = params[4];

    vector<vector<double> > data = read_data(params[3],false);

    int dim = (int)data[0].size();
    int N   = data.size();

    cerr << "Running depth-first-OPT" << '\n';
    cerr << N << " points in " << dim << " dimensions.\n";

    double top_stop_ratio = strTo<double>(params[0]);
    if(top_stop_ratio < 1){
        cerr << "Stopping at " << top_stop_ratio*100.0 << "% points.\n";
    }else{
        cerr << "Stopping at " << (int)top_stop_ratio << " points.\n";
    }
    double stop_ratio = strTo<double>(params[1]);
    if(stop_ratio < 1){
        cerr << "Each level stopping at " << stop_ratio*100.0 << "% points.\n";
    }else{
        cerr << "Each level stopping at " << (int)stop_ratio << " points.\n";
    }
    int levels = strTo<int>(params[2]);
    cerr << "Look-ahead " << levels << " levels.\n";

    mu_timer mt;

    dfopt_tree dfopt(dim,top_stop_ratio,stop_ratio,100,levels);
    map_tree map_region_tree(data.size());
    opt_region_hash<uint32_t> map_regions(24);

    mt.reset();
    dfopt.construct_dfopt_tree(data, map_region_tree, map_regions);
    mt.print_elapsed_time(cerr, "DFOPT construction");

    print_MAP_density(cout, map_regions.print_density(),map_region_tree.get_ra(),data.size());

    // write out the density to file
    ofstream den_file;
    init_file_out(den_file,out_filename+".den",1);
    
    map_region_tree.save(den_file);
    map_regions.save(den_file);
    
    den_file.close();
    
    return 0;
}




int disopt(vector<string> params) {


    string usage_text = "Usage: " + c::PROG_NAME + " disopt [-np] <percent_points> <levels> <top_percent_points> <data_file>\n"
            + "       -np            -- Don't prune the tree, a bit faster but more memory usage\n"
            + "       percent_points -- Ratio of total region data to stop at for each look-ahead (0.01 = 1%)\n"
            + "   top_percent_points -- Ratio of total data to stop at (0.01 = 1% or 2 = 2 points)\n"
            + "               levels -- Maximum levels for each look-ahead\n "
            + "            data_file -- One sample each row\n"
            + "MAP partitions output to STDOUT. Log to STDERR \n"
            + "Run discrepancy limited look-ahead OPT. Good for moderate dimensions (6-15)\n";

    if (params.size() < 4 || params.size() > 5) {
        cerr << usage_text << endl;
        return 3;
    }

    bool prune_tree = true;
    int params_offset = 0;
    if (params.size() == 5) {
        if (params[0] == "-np") {
            prune_tree = false;
            params_offset = 1;
        } else {
            cerr << usage_text << endl;
            return 3;
        }
    }


    vector<vector<double> > data = read_data(params[params_offset + 3],false);

    int dim = (int)data[0].size();
    int N   = data.size();

    cerr << "Running discrepancy limited look-ahead OPT" << '\n';
    cerr << N << " points in " << dim << " dimensions.\n";

    int ll_levels = strTo<int>(params[params_offset + 1]);
    double stop_ratio = strTo<double>(params[params_offset + 0]);
    double top_stop_ratio = strTo<double>(params[params_offset + 2]);

    cerr << "Each small OPT stopping at " << ll_levels << " levels.\n";
    cerr << "Each level stopping at " << stop_ratio*100 << "% of points.\n";
    cerr << "Total stopping at " << top_stop_ratio*100 << "% of points.\n";


    mu_timer mt;

    disopt_tree disopt(dim,stop_ratio,top_stop_ratio,ll_levels,20*dim);
    map_tree map_region_tree(data.size());
    opt_region_hash<uint32_t> map_regions(24);

    mt.reset();
    disopt.construct_disopt_tree(data, map_region_tree, map_regions,prune_tree);
    mt.print_elapsed_time(cerr, "disOPT("+ toStr<int>(ll_levels) +") construction");

    print_MAP_density(cout, map_regions.print_density(),map_region_tree.get_ra(),data.size());

    return 0;
}



int hell_dist(vector<string> params){
    string usage_text = "Usage: " + c::PROG_NAME + " hell_dist <true_samples> <MAP_partitions> \n"
            + "     true_samples   -- Each row one data point followed by true density\n"
            + "     MAP_partitions -- Learned MAP partitions\n"
            + "Compute sampled Hellinger distance\n";

    if (params.size() != 2) {
        cerr << usage_text << endl;
        return 3;
    }

    vector<vector<double> > true_samples = read_data(params[0],true);

    int true_N = (int)true_samples.size();
    int dim    = (int)true_samples[0].size()-1;

    cerr << true_N << " samples in " << dim << " dimensions.\n";

    vector<vector<double> > map_partitions = read_data(params[1],true);

    int N_part = (int)map_partitions.size();
    if((2*(dim)) != ((int)map_partitions[0].size()-1)){
        cerr << "ERROR: Dimension mismatch!\n";
        exit(2);

    }


    double val = 0.0;
    int not_found = 0;
    for(int i = 0;i<true_N;i++){
        double true_den = true_samples[i][dim];
        double map_den  = 0.0;
        bool found = false;
        for(int j = 0;!found&&j<N_part;j++){
            found = true;
            for(int k = 0;found&&k<dim;k++){

                double low  = map_partitions[j][2*k];
                double high = map_partitions[j][(2*k)+1];

                if(true_samples[i][k] <= low || true_samples[i][k] > high){
                    found = false;
                }
            }

            if(found){
                map_den = map_partitions[j][2*dim];
            }
        }

        if(found){
                val += sqrt(map_den/true_den);

                //cerr << val << "," << map_den << "," << true_den << '\n';
        }else{
            cerr << "Sample not found: " << i  << ":" << true_den << '\n';
            not_found++;
        }
    }

    if (val/(double)(true_N-not_found) > 1){
        cerr << "Not enough sample points in true distribution\n";
        cerr << "or density does not integrate to 1.\n";
        exit(2);
    }

    cout << sqrt(1-(val/(double)(true_N-not_found))) << '\n';

    return 0;
}


int classify(vector<string> params){
    string usage_text = "Usage: " + c::PROG_NAME
            + " classify [-c] <learned_dist_0> ... <learned_dist_n-1> <test_data> \n"
            + "     learned_dist_0...n-1 -- MAP partitions from each class\n"
            + "     test_data            -- Each row one data point\n "
            + "Do classification. Doesn't deal with the case of equal densities yet\n";

    if (params.size() < 3) {
        cerr << usage_text << endl;
        return 3;
    }
    
    int param_offset = 0;
    bool confusion = false;
    if(params[0] == "-c"){
        param_offset = 1;
        confusion = true;
    }
    
    vector<vector<double> > test_data = read_data(params[0+param_offset],confusion);

    int num_class = params.size()-1;

    cerr << "Total classes: " << num_class << '\n';

    vector<vector<double> > test_data = read_data(params[num_class],false);

    vector<string> joint_filenames = split(params[1+param_offset],',');
    
    vector<string> marginal_filenames;
    bool copula = false;
    if((params.size()-param_offset) == 3){
        copula = true;
        marginal_filenames = split(params[2+param_offset],',');
    }
    
    if(copula && joint_filenames.size() != marginal_filenames.size()){
        cerr << "ERROR: Must if marginals are given, must be same number as copula densities.\n";
        return 1;
    }
    
    int num_class = (int)joint_filenames.size();
>>>>>>> copula

    int test_N = (int)test_data.size();
    int dim    = (int)test_data[0].size();

    cerr << test_N << " test data points in " << dim << " dimensions.\n";

    vector<vector<vector<double> > > class_dist(num_class);
    for(int i = 0;i<num_class;i++){
        class_dist[i] = read_data(params[i], true);

        if ((2 * (dim)) != ((int) class_dist[i][0].size() - 1)) {
            cerr << "ERROR: Dimension mismatch!\n";
            exit(2);

        }
    }



    // Loop through the data and classify!
    // This should be changed to binary tree
    for(int i = 0;i<test_N;i++){

        vector<double> class_density(num_class,0.0);
        int max_class = -1;
        double max_den = -c::inf;

        for(int c = 0;c < num_class; c++) {
            double map_den = 0.0;
            bool found = false;
            for (int j = 0; !found && j < (int)class_dist[c].size(); j++) {
                found = true;
                for (int k = 0; found && k < dim; k++) {

                    double low = class_dist[c][j][2 * k];
                    double high = class_dist[c][j][(2 * k) + 1];

                    if (test_data[i][k] <= low || test_data[i][k] > high) {
                        found = false;
                    }
                }

                if (found) {
                    map_den = class_dist[c][j][2 * dim];
                }
            }

            if (!found) {
                cerr << "Sample not found: class:" << c << ", idx:" << i << '\n';
            }else{
                class_density[c] = map_den;

                if(map_den > max_den){
                    max_class = c;
                    max_den = map_den;
                }
            }
        }

        // Print it out
        cout << max_class;
        for(int c = 0;c < num_class; c++) {
            cout << '\t' << class_density[c] ;
        }
        cout << '\n';


    }

    return 0;
}



int density(vector<string> params){
    string usage_text = "Usage: " + c::PROG_NAME
            + " classify <MAP_partitions> <sample_data> \n"
            + "     MAP_partitions -- MAP partitions of a distribution\n"
            + "     sample_data    -- Each row one data point\n "
            + "Do classification.\n";

    if (params.size() != 2) {
        cerr << usage_text << endl;
        return 3;
    }


    vector<vector<double> > test_data = read_data(params[1],false);

    int test_N = (int)test_data.size();
    int dim    = (int)test_data[0].size();

    cerr << test_N << " data points in " << dim << " dimensions.\n";

    vector<vector<double> > MAP_dist = read_data(params[0], true);;


    // Loop through the data and classify!
    // This should be changed to binary tree
    for (int i = 0; i < test_N; i++) {


        double map_den = 0.0;
        bool found = false;
        for (int j = 0; !found && j < (int) MAP_dist.size(); j++) {
            found = true;
            for (int k = 0; found && k < dim; k++) {

                double low = MAP_dist[j][2 * k];
                double high = MAP_dist[j][(2 * k) + 1];

                if (test_data[i][k] <= low || test_data[i][k] > high) {
                    found = false;
                }
            }

            if (found) {
                map_den = MAP_dist[j][2 * dim];
            }
        }

        if (!found) {
            cerr << "Sample not found, idx:" << i << '\n';
        } else {
            cout << map_den << '\n';
        }

    }

    return 0;
}


int bench(vector<string> params){
    string usage_text = "Usage: " + c::PROG_NAME
            + " bench <num_points> <dim> \n"
            + "     num_points -- Number of random points to generate\n"
            + "     dim        -- Total dimension, we only count one though\n "
            + "Do benchmarking!.\n";

    if (params.size() != 2) {
        cerr << usage_text << endl;
        return 3;
    }

    int N   = strTo<int>(params[0]);
    int dim = strTo<int>(params[1]);
    
    MT_random randgen(0);
    
    vector<vector<double> > data(N,vector<double>(dim,0.0));
    vector<int> data_one(N,0);

    for(int i = 0; i < N; i++){
        data_one[i] = i;
        for (int j = 0;j < dim;j++){
            data[i][j] = randgen.genrand64_real1();
        }
    }

    mu_timer mt;
    
    mt.reset();
    {
        int cut = 0;
        double lim = 0.5;
        vector<vector<double> > out;
        double first = data[0][dim];
        bool all_same = true;

        for (int i = 0; i < (int) data.size(); i++) {

            if (all_same) {
                if (fabs(first - data[i][dim]) > 1E-19) all_same = false;
            }

            if (cut == 0) {
                if (data[i][dim] < lim) {
                    out.push_back(data[i]);
                }

            } else if (cut == 1) {
                if (data[i][dim] >= lim) {
                    out.push_back(data[i]);
                }

            } 
        }

    }
    mt.print_elapsed_time(cerr, "Full vector try 1");
    
    
    mt.reset();
    {
        int cut = 0;
        double lim = 0.5;
        vector<vector<double> > out;
        double first = data[0][dim];
        bool all_same = true;

        for (int i = 0; i < (int) data.size(); i++) {

            if (all_same) {
                if (fabs(first - data[i][dim]) > 1E-19) all_same = false;
            }

            if (cut == 0) {
                if (data[i][dim] < lim) {
                    out.push_back(data[i]);
                }

            } else if (cut == 1) {
                if (data[i][dim] >= lim) {
                    out.push_back(data[i]);
                }

            } 
        }

    }
    mt.print_elapsed_time(cerr, "Full vector try 2");
    
    
    mt.reset();
    {
        int cut = 0;
        double lim = 0.5;
        vector<int> out;
        double first = data[data_one[0]][dim];
        bool all_same = true;

        for (int i = 0; i < (int) data_one.size(); i++) {

            if (all_same) {
                if (fabs(first - data[data_one[i]][dim]) > 1E-19) all_same = false;
            }

            if (cut == 0) {
                if (data[data_one[i]][dim] < lim) {
                    out.push_back(data_one[i]);
                }

            } else if (cut == 1) {
                if (data[data_one[i]][dim] >= lim) {
                    out.push_back(data_one[i]);
                }

            } 
        }

    }
    mt.print_elapsed_time(cerr, "single vector try 1");
    
    
    
    mt.reset();
    {
        int cut = 0;
        double lim = 0.5;
        vector<int> out;
        double first = data[data_one[0]][dim];
        bool all_same = true;

        for (int i = 0; i < (int) data_one.size(); i++) {

            if (all_same) {
                if (fabs(first - data[data_one[i]][dim]) > 1E-19) all_same = false;
            }

            if (cut == 0) {
                if (data[data_one[i]][dim] < lim) {
                    out.push_back(data_one[i]);
                }

            } else if (cut == 1) {
                if (data[data_one[i]][dim] >= lim) {
                    out.push_back(data_one[i]);
                }

            } 
        }

    }
    mt.print_elapsed_time(cerr, "single vector try 2");
    
    
    return 0;
}


void print_usage_and_exit() {
    cerr << "Usage: " + c::PROG_NAME + " <option>" << "\n";
    cerr << "Options:" << "\n";
    cerr << "-== Density Estimation ==-" << '\n';
    cerr << "  opt        -- MAP partitions from full OPT" << "\n";
    cerr << "  llopt      -- MAP partitions from LL-OPT" << "\n";
    //cerr << "  lsopt      -- MAP partitions from LL-sampled OPT [experimental]" << "\n";
    cerr << "  dfopt      -- MAP partitions from depth-first OPT" << "\n";
    cerr << "  disopt     -- MAP partitions from discrepancy-LL-OPT [experimental]" << "\n";
    cerr << "\n";
    cerr << "-== Other tools ==-" << '\n';
    cerr << "  hell_dist   -- Compute sample Hellinger distance from a known density" << "\n";
    cerr << "  classify    -- Do classification with MAP partitions" << "\n";
    cerr << "  density     -- Get the density at particular points" << "\n";
    //cerr << "  bench       -- Benchmark counting speed [temporary]" << "\n";
    exit(2);
}




