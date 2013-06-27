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
    } else if (mode == "copula") {
        error_num = copula(params);
    } else if (mode == "hell_dist") {
        error_num = hell_dist(params);
    } else if (mode == "classify") {
        error_num = classify(params);
    } else if (mode == "density") {
        error_num = density(params);
    } else if (mode == "density_old") {
        error_num = density_old(params);
    } else if (mode == "bench") {
        error_num = bench(params);
    } else {
        print_usage_and_exit();
    }

    return error_num;

}





int opt(vector<string> params) {

    string usage_text = "Usage: " + c::PROG_NAME + " opt <percent_points> <data_file>\n"
            + "       percent_points -- Ratio of total data to stop at (0.01 = 1%)\n"
            + "            data_file -- One sample each row (Restricted to [0,1] cube)\n"
            + "MAP partitions output to STDOUT. Log to STDERR \n"
            + "Run full-OPT, very fast but uses a lot of memory. If dimension \n"
            + "greater than 5 use the other methods\n";

    if (params.size() != 2) {
        cerr << usage_text << endl;
        return 3;
    }

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
    map_tree map_region_tree(N,dim);
    opt_region_hash<uint32_t> map_regions(20);
    opt_slow.construct_MAP_tree(map_region_tree, map_regions, N);
    total_time += mt.elapsed_time();
    mt.print_elapsed_time(cerr, "MAP tree");

    cerr << "OPT construction time: " << total_time << " s.\n";
    
    cerr << "lPhi: " << opt_slow.get_lphi() << endl;
    print_MAP_density(cout, map_regions.get_regions(),map_region_tree.get_ra(),data.size());


    return 0;
}


int llopt(vector<string> params) {


    string usage_text = "Usage: " + c::PROG_NAME + " llopt [-np] <percent_points> <levels> <top_percent_points> <data_file> <output_name>\n"
            + "       -np            -- Don't prune the tree, a lot faster but more memory usage\n"
            + "       percent_points -- Ratio of total region data to stop at for each look-ahead (0.01 = 1%)\n"
            + "   top_percent_points -- Ratio of total data to stop at (0.01 = 1% or 2 = 2 points)\n"
            + "               levels -- Maximum levels for each look-ahead\n "
            + "            data_file -- One sample each row (Restricted to [0,1] cube)\n"
            + "MAP partitions output to STDOUT and <output_name>.den . Log to STDERR \n"
            + "Run limited look-ahead OPT. Good for moderate dimensions (6-15)\n";

    if (params.size() < 5 || params.size() > 6) {
        cerr << usage_text << endl;
        return 3;
    }

    bool prune_tree = true;
    int params_offset = 0;
    if (params.size() == 6) {
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
    
    string out_filename = params[params_offset + 4];

    cerr << "Each small OPT stopping at " << ll_levels << " levels.\n";
    cerr << "Each level stopping at " << stop_ratio*100 << "% of points.\n";
    cerr << "Total stopping at " << top_stop_ratio*100 << "% of points.\n";


    mu_timer mt;

    llopt_tree llopt(dim,stop_ratio,top_stop_ratio,ll_levels,20*dim);
    map_tree map_region_tree(N,dim);
    opt_region_hash<uint32_t> map_regions(20);

    mt.reset();
    llopt.construct_llopt_tree(&data, map_region_tree, map_regions,prune_tree);
    mt.print_elapsed_time(cerr, "LLOPT("+ toStr<int>(ll_levels) +") construction");

    // write to cout
    print_MAP_density(cout, map_regions.get_regions(),map_region_tree.get_ra(),N);

    // write out the density to file
    ofstream den_file;
    init_file_out(den_file,out_filename+".den",1);
    
    map_region_tree.save(den_file);
    map_regions.save(den_file);
    
    den_file.close();
    
    return 0;
}


int lsopt(vector<string> params) {


    string usage_text = "Usage: " + c::PROG_NAME + " lsopt <percent_points> <lookahead_depth> <top_percent_points> <iterations> <convergence_iterations> <data_file>\n"
            + "         percent_points -- Ratio of total region data to stop at for each look-ahead (0.01 = 1%)\n"
            + "        lookahead_depth -- Depth to stop at per look-ahead\n"
            + "     top_percent_points -- Ratio of total data to stop at (0.01 = 1%)\n"
            + "             iterations -- Maximum number of iterations for each look-ahead\n "
            + " convergence_iterations -- Number of iterations the same before convergence\n "
            + "              data_file -- One sample each row (Restricted to [0,1] cube)\n"
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
    map_tree map_region_tree(N,dim);
    opt_region_hash<uint32_t> map_regions(20);

    mt.reset();
    lsopt.construct_lsopt_tree(data, iterations,convergence_iterations, map_region_tree, map_regions);
    mt.print_elapsed_time(cerr, "LSOPT construction");

    print_MAP_density(cout, map_regions.get_regions(),map_region_tree.get_ra(),N);

    
    return 0;
}

int dfopt(vector<string> params) {

    string usage_text = "Usage: " + c::PROG_NAME + " dfopt <top_percent_points>" 
            + " <percent_points> <levels> <data_file>\n"
            + "     top_percent_points -- Points to stop at (0.01 = 1%, 1 = 1point)\n"
            + "     percent_points     -- Points in each level to stop at (0.01 = 1%, 1 = 1point))\n"
            + "     levels             -- look ahead levels\n"
            + "     data_file          -- One sample per row (Restricted to [0,1] cube)\n"
            + "MAP partitions output to STDOUT. Log to STDERR\n "
            + "Run depth-first-OPT, very slow but technically can do any number of dimensions\n "
            + "if you wait several years. Run in a similar fashion to limited look-ahead. \n";

    if (params.size() != 4) {
        cerr << usage_text << endl;
        return 3;
    }

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
    map_tree map_region_tree(N,dim);
    opt_region_hash<uint32_t> map_regions(20);

    mt.reset();
    dfopt.construct_dfopt_tree(data, map_region_tree, map_regions);
    mt.print_elapsed_time(cerr, "DFOPT construction");

    print_MAP_density(cout, map_regions.get_regions(),map_region_tree.get_ra(),N);

    return 0;
}




int disopt(vector<string> params) {


    string usage_text = "Usage: " + c::PROG_NAME + " disopt [-np] <percent_points> <levels> <top_percent_points> <data_file>\n"
            + "       -np            -- Don't prune the tree, a bit faster but more memory usage\n"
            + "       percent_points -- Ratio of total region data to stop at for each look-ahead (0.01 = 1%)\n"
            + "   top_percent_points -- Ratio of total data to stop at (0.01 = 1% or 2 = 2 points)\n"
            + "               levels -- Maximum levels for each look-ahead\n "
            + "            data_file -- One sample each row (Restricted to [0,1] cube)\n"
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
    map_tree map_region_tree(N,dim);
    opt_region_hash<uint32_t> map_regions(20);

    mt.reset();
    disopt.construct_disopt_tree(data, map_region_tree, map_regions,prune_tree);
    mt.print_elapsed_time(cerr, "disOPT("+ toStr<int>(ll_levels) +") construction");

    print_MAP_density(cout, map_regions.get_regions(),map_region_tree.get_ra(),N);

    return 0;
}


int copula(vector<string> params){
    string usage_text = "Usage: " + c::PROG_NAME + " copula <percent_points> <data_file> <output_name>\n"
            + "       percent_points -- Data count to stop at (0.01 = 1% or 2 = 2 points)\n"
            + "            data_file -- One sample each row (Restricted to [0,1] cube)\n"
            + "Transformed data to <output_name>.copula.txt\n"
            + "Marginal densities to <output_name>.copula.den\n"
            + "Log to STDERR \n"
            + "Do copula transform by estimating marginal densities with full OPT \n";

    if (params.size() != 3) {
        cerr << usage_text << endl;
        return 3;
    }

    string data_file = params[1];
    string out_filename = params[2];
    
    vector<vector<double> > data = read_data(data_file,false);

    int dim = (int)data[0].size();
    uint32_t N   = (uint32_t)data.size();

    cerr << "Running copula via Full-OPT" << '\n';
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

    // create file to save all the densities
    
    ofstream den_file;
    init_file_out(den_file,out_filename+".copula.den",dim);

    
    for (int i = 0; i < dim; i++) {
        
        cerr << "Computing marginal for dim: " << i+1 << '\n';
        // Extract one dimension of data
        
        vector<vector<double> > one_data(N,vector<double>(1,0.0));
        
        for(uint32_t j = 0;j<N;j++){
            one_data[j][0] = data[i][j];
        }
        
        // run full OPT on data

        double total_time = 0.0;
        mu_timer mt;
        opt_tree opt_slow(1, stop_points, 1000);

        mt.reset();
        opt_slow.construct_full_tree(one_data);
        total_time += mt.elapsed_time();
        mt.print_elapsed_time(cerr, "OPT tree");

        mt.reset();
        map_tree map_region_tree(N,dim);
        opt_region_hash<uint32_t> map_regions(20);
        opt_slow.construct_MAP_tree(map_region_tree, map_regions, N);
        total_time += mt.elapsed_time();
        mt.print_elapsed_time(cerr, "MAP tree");

        cerr << "Total construction time: " << total_time << " s.\n";
        cerr << "lPhi: " << opt_slow.get_lphi() << endl;
        //print_MAP_density(cout, map_regions.get_regions(), map_region_tree.get_ra(), data.size());

        // Construct CDF from the regions
        cdf marginal(map_region_tree,map_regions);
        
        // Transform the data
        for(uint32_t j = 0;j<N;j++){
            one_data[j][0] = marginal.transform(one_data[j][0]);
        }
        
        // Replace the data that was read in
        
        for(uint32_t j = 0;j<N;j++){
            data[i][j] = one_data[j][0];
        }
        
        // Write out the marginal into the den file
        map_region_tree.save(den_file);
        map_regions.save(den_file);
        
    }
    
    den_file.close();
    
    // Write out the transformed data
    ofstream out_file;
    {
        string temp = out_filename + ".copula.txt";
        out_file.open(temp.c_str(), ios::out);
    }
    
    out_file << std::scientific;
    for(uint32_t i = 0;i<N;i++){
        out_file << data[i][0];
        for(int j = 1;j<dim;j++){
            out_file << ' ' << data[i][j];
        }
        out_file << '\n';
    }
    
    out_file.close();
    
    
    return 0;
}



int hell_dist(vector<string> params){
    string usage_text = "Usage: " + c::PROG_NAME + " hell_dist <true_samples> <joint/copula_den> [marginal_den] \n"
            + "     true_samples     -- Each row one data point followed by true density\n"
            + "     joint/copula_den -- Joint or Copula density, Copula needs marginals\n"
            + "     marginal_den     -- Marginal densities, if given then copula is used\n"
            + "Compute sampled Hellinger distance\n";

    if (params.size() < 2 || params.size() > 3) {
        cerr << usage_text << endl;
        return 3;
    }

    vector<vector<double> > true_samples = read_data(params[0],true);

    int true_N = (int)true_samples.size();
    int dim    = (int)true_samples[0].size()-1;

    cerr << true_N << " samples in " << dim << " dimensions.\n";

    
    string joint_filename = params[1];
    
    string marginal_filename = "";
    bool copula = false;
    if(params.size() == 3){
        copula = true;
        marginal_filename = params[2];
    }

    // load the joint/copula densities
    map_tree map_region_tree(1,1);
    opt_region_hash<uint32_t> map_regions(20);
    
    map_tree** m_map_region_tree;
    opt_region_hash<uint32_t>** m_map_regions;
    cdf** marginal;
    
    int num_dim = 0;
    
    int status = load_densities(joint_filename, marginal_filename,
        &map_region_tree, &map_regions, m_map_region_tree, m_map_regions,
        marginal, copula);
    
    if(status != 0){
        return status;
    }
    
    num_dim = map_region_tree.get_num_children();
    
    
    double val = 0.0;
    int not_found = 0;
    for(int i = 0;i<true_N;i++){
        double true_den = true_samples[i][dim];
        double map_den  = 0.0;
        
        vector<double> point(true_samples[i].begin(),true_samples[i].end()-1);
        
        map_den = compute_density(point, &map_region_tree,m_map_region_tree,marginal,copula);

        val += sqrt(map_den / true_den); // importance sample
    }

    if (val/(double)(true_N-not_found) > 1){
        cerr << "Not enough sample points in true distribution\n";
        cerr << "or density does not integrate to 1.\n";
        exit(2);
    }

    cout << sqrt(1-(val/(double)(true_N-not_found))) << '\n';

    
    for(int i = 0;i<num_dim;i++){
        delete m_map_region_tree[i];
        delete m_map_regions[i];
        delete marginal[i];
    }
    
    delete [] m_map_region_tree;
    delete [] m_map_regions;
    delete [] marginal;
    
    return 0;
}


int classify(vector<string> params){
    string usage_text = "Usage: " + c::PROG_NAME
            + " classify <test_data> <learned_dist_0,...,learned_dist_n-1> [marginal_dist_0,...,marginal_dist_n-1] \n"
            + "     learned_dist_0...n-1  -- MAP partitions from each class (or copula partition)\n"
            + "     marginal_dist_0...n-1 -- Marginals generated by copula transform\n"
            + "     test_data             -- Each row one data point\n "
            + "Do classification. Doesn't deal with the case of equal densities yet\n";
    
    if (params.size() < 3) {
        cerr << usage_text << endl;
        return 3;
    }

    int num_class = params.size()-1;

    cerr << "Total classes: " << num_class << '\n';

    vector<vector<double> > test_data = read_data(params[num_class],false);

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



int density_old(vector<string> params){
    string usage_text = "Usage: " + c::PROG_NAME
            + " density_old <MAP_partitions> <sample_data> \n"
            + "     MAP_partitions -- MAP partitions of a distribution\n"
            + "     sample_data    -- Each row one data point\n "
            + "Get the density at locations.\n";

    if (params.size() != 2) {
        cerr << usage_text << endl;
        return 3;
    }


    vector<vector<double> > test_data = read_data(params[1],false);

    int test_N = (int)test_data.size();
    int dim    = (int)test_data[0].size();

    cerr << test_N << " data points in " << dim << " dimensions.\n";

    vector<vector<double> > MAP_dist = read_data(params[0], true);;


    // Loop through the data
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

int density(vector<string> params){
    string usage_text = "Usage: " + c::PROG_NAME
            + " density <sample_data> <joint_den/copula_den> [marginal_den]  \n"
            + "     sample_data    -- Each row one data point\n "
            + "     joint_den      -- Joint density (of copula) from OPT\n "
            + "     marginal_den   -- Marginal densities from OPT, for copula [Optional]\n "
            + "Gets the density at each sample point. If copula density is given\n"
            + "then a copula transform is performed with the marginals\n";

    if (params.size() < 2 || params.size() > 3) {
        cerr << usage_text << endl;
        return 3;
    }


    vector<vector<double> > test_data = read_data(params[0],false);

    int test_N = (int)test_data.size();
    int dim    = (int)test_data[0].size();

    cerr << test_N << " data points in " << dim << " dimensions.\n";

    string joint_filename = params[1];
    
    string marginal_filename = "";
    bool copula = false;
    if(params.size() == 3){
        copula = true;
        marginal_filename = params[2];
    }

    // load the joint/copula densities
    map_tree map_region_tree(1,1);
    opt_region_hash<uint32_t> map_regions(20);
    
    map_tree** m_map_region_tree = NULL;
    opt_region_hash<uint32_t>** m_map_regions = NULL;
    cdf** marginal = NULL;
    
    int num_dim = 0;
    
    int status = load_densities(joint_filename, marginal_filename,
        &map_region_tree, &map_regions, m_map_region_tree, m_map_regions,
        marginal, copula);
    
    if(status != 0){
        return status;
    }
    
    num_dim = map_region_tree.get_num_children();
    
    // Check things were loaded sensibly
    //cerr << "num_dim: " << num_dim << '\n';
    //cerr << "num points: " << map_region_tree.get_num_points() << '\n';
    //cerr << "---------\n";
    //print_MAP_density(cerr, map_regions.get_regions(),map_region_tree.get_ra(),map_region_tree.get_num_points());
    
    // Loop through the data
    for (vector<vector<double> >::iterator it = test_data.begin();
            it != test_data.end(); it++) {
        // for each data point
        
        double den = compute_density(*it, &map_region_tree,m_map_region_tree,marginal,copula);
        
        cerr << "den: " << den << '\n';
        
        cout << std::scientific
                << den
                << '\n';
    }
    
    for(int i = 0;i<num_dim;i++){
        if(m_map_region_tree != NULL) delete m_map_region_tree[i];
        if(m_map_regions != NULL) delete m_map_regions[i];
        if(marginal != NULL) delete marginal[i];
    }
    
    if(m_map_region_tree != NULL) delete [] m_map_region_tree;
    if(m_map_regions != NULL) delete [] m_map_regions;
    if(marginal != NULL) delete [] marginal;

    return 0;
}

int bench(vector<string> params){
    string usage_text = "Usage: " + c::PROG_NAME
            + " bench <num_points> <dim> \n"
            + "     num_points -- Number of random points to generate\n"
            + "     dim        -- Total dimension, we only count one though\n "
            + "Do benchmarking of counting!.\n";

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
    cerr << "  lsopt      -- MAP partitions from LL-sampled OPT [experimental]" << "\n";
    cerr << "  dfopt      -- MAP partitions from depth-first OPT" << "\n";
    //cerr << "  disopt     -- MAP partitions from discrepancy-LL-OPT [experimental]" << "\n";
    cerr << "  copula     -- Perform copula transform with full OPT" << "\n";
    cerr << "\n";
    cerr << "-== Other tools ==-" << '\n';
    cerr << "  hell_dist   -- Compute sample Hellinger distance from a known density" << "\n";
    cerr << "  classify    -- Do classification with MAP partitions [experimental]" << "\n";
    cerr << "  density     -- Get the density at particular points [experimental]" << "\n";
    cerr << "  bench       -- Benchmark counting speed [temporary]" << "\n";
    exit(2);
}




