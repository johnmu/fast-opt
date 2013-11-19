/*
 *  main.cpp
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
    } else if (mode == "opt_comp") {
        error_num = opt_comp(params);
    } else if (mode == "coopt") {
        error_num = coopt(params);
    } else if (mode == "coopt_scan_old") {
        error_num = coopt_scan_old(params);
    } else if (mode == "coopt_scan") {
        error_num = coopt_scan(params);
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
    } else if (mode == "print") {
        error_num = print_partitions(params);
    } else if (mode == "density_old") {
        error_num = density_old(params);
    } else if (mode == "bench") {
        error_num = bench(params);
    } else if (mode == "normalize") {
        error_num = normalize(params);
    } else if (mode == "demean") {
        error_num = demean(params);
    } else if (mode == "vec_quant_sam_quals") {
        error_num = vec_quant_sam_quals(params);
    } else {
        print_usage_and_exit();
    }

    return error_num;

}

int opt(vector<string> params) {

    string usage_text = "Usage: " + c::PROG_NAME + " opt <percent_points> <data_file> <output_name>\n"
            + "       percent_points -- Ratio of total data to stop at (0.01 = 1%, or 2 = 2 points)\n"
            + "            data_file -- One sample each row (Restricted to [0,1] cube)\n"
            + "          output_name -- Result output to <output_name>.den\n"
            + "Log output to STDERR \n"
            + "Run full-OPT, very fast but uses a lot of memory. If dimension \n"
            + "greater than 4 use the other methods. Best choice for 1 or 2 dimensional data.\n";

    if (params.size() != 3) {
        cerr << usage_text << endl;
        return 3;
    }

    string out_filename = params[2];

    vector<vector<double> > data;
    vector<uint32_t> skipped; 
    read_data(params[1],data,skipped);

    int dim = (int) data[0].size();
    int N = data.size();

    cerr << "Running Full-OPT" << '\n';
    cerr << N << " points in " << dim << " dimensions.\n";

    double stop_ratio = strTo<double>(params[0]);
    int stop_points = 0;

    if (stop_ratio < 1) {
        stop_points = (int) (N * stop_ratio);
    } else {
        stop_points = (int) stop_ratio;
    }
    if (stop_points < 1) stop_points = 1;

    cerr << "Stopping at " << stop_points << " points.\n";

    double total_time = 0.0;

    mu_timer mt;


    opt_tree opt_slow(dim, stop_points, 1000);

    mt.reset();
    opt_slow.construct_full_tree(data);
    total_time += mt.elapsed_time();
    mt.print_elapsed_time(cerr, "OPT tree");

    mt.reset();
    map_tree map_region_tree(N, dim);
    opt_region_hash<uint32_t> map_regions(20);
    opt_slow.construct_MAP_tree(map_region_tree, map_regions, N);
    total_time += mt.elapsed_time();
    mt.print_elapsed_time(cerr, "MAP tree");

    cerr << "OPT construction time: " << total_time << " s.\n";

    cerr << "Log Phi: " << opt_slow.get_lphi() << endl;

    // write out the density to file
    ofstream den_file;
    init_file_out(den_file, out_filename + ".den", 1);

    map_region_tree.save(den_file);
    map_regions.save(den_file);

    den_file.close();

    return 0;
}

int llopt(vector<string> params) {


    string usage_text = "Usage: " + c::PROG_NAME + " llopt [-np] <percent_points> <levels> <top_percent_points> <data_file> <output_name>\n"
            + "                  -np -- Don't prune the tree, a lot faster but more memory usage\n"
            + "       percent_points -- Ratio of total region data to stop at for each look-ahead (0.01 = 1%)\n"
            + "   top_percent_points -- Ratio of total data to stop at (0.01 = 1% or 2 = 2 points)\n"
            + "               levels -- Maximum levels for each look-ahead\n "
            + "            data_file -- One sample each row (Restricted to [0,1] cube)\n"
            + "          output_name -- MAP partitions output to <output_name>.den\n"
            + "Log output to STDERR \n"
            + "Run limited look-ahead OPT. Good for moderate dimensions (5-50).\n"
            + "For higher numbers of dimensions (20-50) use levels = 2.\n"
            + "Memory usage and running time scales exponentially with look-ahead levels.\n";

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

    vector<vector<double> > data;
    vector<uint32_t> skipped; 
    read_data(params[params_offset + 3],data,skipped);

    int dim = (int) data[0].size();
    int N = data.size();

    cerr << "Running limited look-ahead OPT" << '\n';
    cerr << N << " points in " << dim << " dimensions.\n";

    int ll_levels = strTo<int>(params[params_offset + 1]);
    double stop_ratio = strTo<double>(params[params_offset + 0]);
    double top_stop_ratio = strTo<double>(params[params_offset + 2]);

    string out_filename = params[params_offset + 4];

    cerr << "Each small OPT stopping at " << ll_levels << " levels.\n";
    cerr << "Each level stopping at " << stop_ratio * 100 << "% of points.\n";

    if (top_stop_ratio >= 1) {
        cerr << "Total stopping at " << top_stop_ratio << " points.\n";
    } else {
        cerr << "Total stopping at " << top_stop_ratio * 100 << "% of points.\n";
    }

    mu_timer mt;

    llopt_tree llopt(dim, stop_ratio, top_stop_ratio, ll_levels, 20 * dim);
    map_tree map_region_tree(N, dim);
    opt_region_hash<uint32_t> map_regions(20);

    mt.reset();
    llopt.construct_llopt_tree(&data, map_region_tree, map_regions, prune_tree);
    mt.print_elapsed_time(cerr, "LLOPT(" + toStr<int>(ll_levels) + ") construction");

    // write out the density to file
    ofstream den_file;
    init_file_out(den_file, out_filename + ".den", 1);

    map_region_tree.save(den_file);
    map_regions.save(den_file);

    den_file.close();

    return 0;
}

// not used for now... need to convert to new format

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

    
    vector<vector<double> > data;
    vector<uint32_t> skipped;
    read_data(params[5],data,skipped);

    int dim = (int) data[0].size();
    int N = data.size();

    cerr << "Running limited look-ahead sampled-OPT" << '\n';
    cerr << "Max iterations: " << iterations << '\n';
    cerr << "Iterations to converge: " << convergence_iterations << '\n';
    cerr << N << " points in " << dim << " dimensions.\n";

    double stop_ratio = strTo<double>(params[0]);
    double top_stop_ratio = strTo<double>(params[2]);

    cerr << "Each level stopping at " << stop_ratio * 100 << "% of points.\n";
    cerr << "Total stopping at " << top_stop_ratio * 100 << "% of points.\n";


    mu_timer mt;

    lsopt_tree lsopt(dim, stop_ratio, top_stop_ratio, lookahead_depth, 1000);
    map_tree map_region_tree(N, dim);
    opt_region_hash<uint32_t> map_regions(20);

    mt.reset();
    lsopt.construct_lsopt_tree(data, iterations, convergence_iterations, map_region_tree, map_regions);
    mt.print_elapsed_time(cerr, "LSOPT construction");

    print_MAP_density(cout, map_regions.get_regions(), map_region_tree.get_ra(), N);


    return 0;
}

int dfopt(vector<string> params) {

    string usage_text = "Usage: " + c::PROG_NAME + " dfopt <percent_points>"
            + " <levels> <top_percent_points> <data_file> <output_name>\n"
            + "     percent_points     -- Points in each level to stop at (0.01 = 1%, 1 = 1point))\n"
            + "     levels             -- look ahead levels\n"
            + "     top_percent_points -- Points to stop at (0.01 = 1%, 1 = 1point)\n"
            + "     data_file          -- One sample per row (Restricted to [0,1] cube)\n"
            + "     output_name        -- One sample per row (Restricted to [0,1] cube)\n"
            + "MAP partitions output to STDOUT and <output_name>.den .Log to STDERR \n"
            + "Run depth-first-OPT, very slow but technically can do any number of dimensions\n "
            + "if you wait several years. Run in a similar fashion to limited look-ahead. \n"
            + "Not recommended for any dataset, unless computer has limited memory. \n";

    if (params.size() != 5) {
        cerr << usage_text << endl;
        return 3;
    }

    string out_filename = params[4];

    vector<vector<double> > data;
    vector<uint32_t> skipped;
    read_data(params[3],data,skipped);

    int dim = (int) data[0].size();
    int N = data.size();

    cerr << "Running depth-first-OPT" << '\n';
    cerr << N << " points in " << dim << " dimensions.\n";

    double top_stop_ratio = strTo<double>(params[2]);
    if (top_stop_ratio < 1) {
        cerr << "Stopping at " << top_stop_ratio * 100.0 << "% points.\n";
    } else {
        cerr << "Stopping at " << (int) top_stop_ratio << " points.\n";
    }
    double stop_ratio = strTo<double>(params[0]);
    if (stop_ratio < 1) {
        cerr << "Each level stopping at " << stop_ratio * 100.0 << "% points.\n";
    } else {
        cerr << "Each level stopping at " << (int) stop_ratio << " points.\n";
    }
    int levels = strTo<int>(params[1]);
    cerr << "Look-ahead " << levels << " levels.\n";

    mu_timer mt;

    dfopt_tree dfopt(dim, top_stop_ratio, stop_ratio, 100, levels);
    map_tree map_region_tree(N, dim);
    opt_region_hash<uint32_t> map_regions(20);

    mt.reset();
    dfopt.construct_dfopt_tree(data, map_region_tree, map_regions);
    mt.print_elapsed_time(cerr, "DFOPT construction");

    // write out the density to file
    ofstream den_file;
    init_file_out(den_file, out_filename + ".den", 1);

    map_region_tree.save(den_file);
    map_regions.save(den_file);

    den_file.close();

    return 0;
}

// this wasn't really a good idea anyways

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

    vector<vector<double> > data;
    vector<uint32_t> skipped;
    read_data(params[params_offset + 3],data,skipped);

    int dim = (int) data[0].size();
    int N = data.size();

    cerr << "Running discrepancy limited look-ahead OPT" << '\n';
    cerr << N << " points in " << dim << " dimensions.\n";

    int ll_levels = strTo<int>(params[params_offset + 1]);
    double stop_ratio = strTo<double>(params[params_offset + 0]);
    double top_stop_ratio = strTo<double>(params[params_offset + 2]);

    cerr << "Each small OPT stopping at " << ll_levels << " levels.\n";
    cerr << "Each level stopping at " << stop_ratio * 100 << "% of points.\n";
    cerr << "Total stopping at " << top_stop_ratio * 100 << "% of points.\n";


    mu_timer mt;

    disopt_tree disopt(dim, stop_ratio, top_stop_ratio, ll_levels, 20 * dim);
    map_tree map_region_tree(N, dim);
    opt_region_hash<uint32_t> map_regions(20);

    mt.reset();
    disopt.construct_disopt_tree(data, map_region_tree, map_regions, prune_tree);
    mt.print_elapsed_time(cerr, "disOPT(" + toStr<int>(ll_levels) + ") construction");

    print_MAP_density(cout, map_regions.get_regions(), map_region_tree.get_ra(), N);

    return 0;
}

int copula(vector<string> params) {
    string usage_text = "Usage: " + c::PROG_NAME + " copula <percent_points> <data_file> <output_name>\n"
            + "       percent_points -- Data count to stop at (0.01 = 1% or 2 = 2 points)\n"
            + "            data_file -- One sample each row (Restricted to [0,1] cube)\n"
            + "Transformed data to <output_name>.copula.txt\n"
            + "Marginal densities to <output_name>.marginal.den\n"
            + "Log to STDERR \n"
            + "Do copula transform by estimating marginal densities with full OPT \n";

    if (params.size() != 3) {
        cerr << usage_text << endl;
        return 3;
    }

    string data_file = params[1];
    string out_filename = params[2];

    vector<vector<double> > data;
    vector<uint32_t> skipped;
    read_data(data_file,data,skipped);

    int dim = (int) data[0].size();
    uint32_t N = (uint32_t) data.size();

    cerr << "Running copula via Full-OPT" << '\n';
    cerr << N << " points in " << dim << " dimensions.\n";

    double stop_ratio = strTo<double>(params[0]);
    int stop_points = 0;

    if (stop_ratio < 1) {
        stop_points = (int) (N * stop_ratio);
    } else {
        stop_points = (int) stop_ratio;
    }
    if (stop_points < 1) stop_points = 1;

    cerr << "Stopping at " << stop_points << " points.\n";

    // create file to save all the densities

    ofstream den_file;
    init_file_out(den_file, out_filename + ".marginal.den", dim);


    for (int i = 0; i < dim; i++) {

        cerr << "Computing marginal for dim: " << i + 1 << '\n';
        // Extract one dimension of data

        vector<vector<double> > one_data(N, vector<double>(1, 0.0));

        for (uint32_t j = 0; j < N; j++) {
            one_data[j][0] = data[j][i];
        }


        // run full OPT on data

        double total_time = 0.0;
        mu_timer mt;
        opt_tree opt_slow(1, stop_points, 1000);

        mt.reset();
        opt_slow.construct_full_tree(one_data);
        total_time += mt.elapsed_time();
        //mt.print_elapsed_time(cerr, "OPT tree");

        mt.reset();
        map_tree map_region_tree(N, 1);
        opt_region_hash<uint32_t> map_regions(20);
        opt_slow.construct_MAP_tree(map_region_tree, map_regions, N);
        total_time += mt.elapsed_time();
        //mt.print_elapsed_time(cerr, "MAP tree");

        cerr << "Total construction time: " << total_time << " s.\n";
        cerr << "Log Phi: " << opt_slow.get_lphi() << endl;
        //print_MAP_density(cout, map_regions.get_regions(), map_region_tree.get_ra(), data.size());

        // Construct CDF from the regions
        cdf marginal(map_region_tree, map_regions);

        // Transform the data
        for (uint32_t j = 0; j < N; j++) {
            one_data[j][0] = marginal.transform(one_data[j][0]);
        }

        // Replace the data that was read in

        for (uint32_t j = 0; j < N; j++) {
            data[j][i] = one_data[j][0];
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
    for (uint32_t i = 0; i < N; i++) {
        out_file << data[i][0];
        for (int j = 1; j < dim; j++) {
            out_file << ' ' << data[i][j];
        }
        out_file << '\n';
    }

    out_file.close();


    return 0;
}

int hell_dist(vector<string> params) {
    string usage_text = "Usage: " + c::PROG_NAME + " hell_dist <true_samples> <joint/copula_den/density_list> [marginal_den] \n"
            + "     true_samples     -- Each row one data point, last column is true density (Restricted to [0,1] cube)\n"
            + "     joint/copula_den -- Joint or Copula density, Copula needs marginals, use - for uniform\n"
            + "     density_list     -- If not .den file, will try to read as a list of densities\n"
            + "     marginal_den     -- Marginal densities, if given then copula is used\n"
            + "Compute sampled Hellinger distance via importance sampling. "
            + "The true samples must be sampled from the true distribution.\n"
            + "Computed as H(f,g) = \\sqrt(1 - \\int{\\sqrt(f(x)g(x))dx})";

    if (params.size() < 2 || params.size() > 3) {
        cerr << usage_text << endl;
        return 3;
    }

    vector<vector<double> > true_samples;
    vector<uint32_t> skipped;
    read_data(params[0], true_samples, skipped);

    int true_N = (int) true_samples.size();
    int dim = (int) true_samples[0].size() - 1;

    cerr << true_N << " samples in " << dim << " dimensions.\n";

    string joint_filename = params[1];

    if (joint_filename == "-") {
        joint_filename = "";
    }

    string marginal_filename = "";
    bool copula = false;
    if (params.size() == 3) {
        copula = true;
        marginal_filename = params[2];
    }

    // load the joint/copula densities
    bool den_list = false;
    vector<double> density_list;
    density_store dens(joint_filename, marginal_filename);

    if(!dens.is_loaded() && marginal_filename.length() == 0){
        // try to load as density list
        density_list.reserve(true_samples.size());
        ifstream infile(joint_filename.c_str());
        
        if(!infile.is_open()){
            cerr << "Cannot open density file: " << joint_filename << '\n';
            return 1;
        }
        
        string temp = "";
        uint32_t skip_idx = 0;
        uint32_t idx = 0;
        while(!infile.eof()){
            temp = "";
            getline(infile,temp);
            trim2(temp);
            if(temp.length() == 0){
                continue;
            }
            bool skip =  false;
            if(skip_idx < skipped.size()){
                if(idx == skipped[skip_idx]){
                    skip_idx++;
                    skip = true;
                }
            }
            if(!skip){
                density_list.push_back(strTo<double>(temp));
            }
            idx++;
        }
        
        infile.close();
        
        cerr << "Read " << density_list.size() << " density values\n";
        
        if(density_list.size() != true_samples.size()){
            cerr << "Incorrect density file length or format...\n";
            return 1;
        }
        den_list = true;
    }else if(!dens.is_loaded()){
        cerr << "Error loading files...\n";
        return 1;
    }
    
    double val = 0.0;

    for (int i = 0; i < true_N; i++) {
        double true_den = true_samples[i][dim];
        double map_den = 0.0;

        if(!den_list) {
            vector<double> point(true_samples[i].begin(), true_samples[i].end() - 1);
            map_den = dens.compute_density(point, 0);
        } else {
            map_den = density_list[i];
        }

        val += sqrt(map_den / true_den); // importance sample
    }

    if (val / (double) (true_N) > 1) {
        cerr << "Not enough sample points in true distribution\n";
        cerr << "or density does not integrate to 1.\n";
        exit(2);
    }

    cout << sqrt(1 - (val / (double) (true_N))) << '\n';

    return 0;
}

int classify(vector<string> params) {
    string usage_text = "Usage: " + c::PROG_NAME
            + " classify [-c] <test_data> <prior_0,...,prior_n-1> <learned_dist_0,...,learned_dist_n-1> [marginal_dist_0,...,marginal_dist_n-1] \n"
            + "                       -c -- Classes given in test data, confusion matrix will be printed [optional]"
            + "                test_data -- Each row one data point, last column can be class (Restricted to [0,1] cube)\n "
            + "            prior_0...n-1 -- The class prior, automatically normalized\n"
            + "     learned_dist_0...n-1 -- MAP partitions from each class (or copula partition), use - for uniform\n"
            + "    marginal_dist_0...n-1 -- Marginals generated by copula transform\n"
            + "Do classification. If densities equal, random one is selected\n"
            + "Prediction output to STDOUT, Confusion matrix to STDERR\n";

    if (params.size() < 3 || params.size() > 5) {
        cerr << usage_text << endl;
        return 3;
    }

    int param_offset = 0;
    bool confusion = false;
    if (params[0] == "-c") {
        param_offset = 1;
        confusion = true;
    }

    vector<double> prior;
    {
        vector<string> prior_str = split(params[1 + param_offset], ',');

        double prior_sum = 0;
        for (int i = 0; i < (int) prior_str.size(); i++) {
            double val = strTo<double>(prior_str[i]);
            prior.push_back(val);
            prior_sum += val;
        }

        // do normalization
        for (int i = 0; i < (int) prior_str.size(); i++) {
            prior[i] = prior[i] / prior_sum;
        }
    }

    vector<vector<double> > test_data;
    vector<uint32_t> skipped;
    read_data(params[0 + param_offset],test_data,skipped,confusion);

    int test_N = (int) test_data.size();
    int dim = (int) test_data[0].size() - param_offset;

    cerr << test_N << " data points in " << dim << " dimensions.\n";

    vector<string> joint_filenames = split(params[2 + param_offset], ',');

    vector<string> marginal_filenames;
    bool copula = false;
    if ((params.size() - param_offset) == 4) {
        copula = true;
        marginal_filenames = split(params[3 + param_offset], ',');
    }

    if (copula && joint_filenames.size() != marginal_filenames.size()) {
        cerr << "ERROR: Must if marginals are given, must be same number as copula densities.\n";
        return 1;
    }

    if (joint_filenames.size() != prior.size()) {
        cerr << "ERROR: Prior must be same length as classes\n";
        return 1;
    }

    int num_class = (int) joint_filenames.size();

    cerr << "Total classes: " << num_class << '\n';

    if (confusion) {
        // check the classes are right
        for (int i = 0; i < test_N; i++) {
            double class_val = test_data[i].back();
            int class_val_i = (int) class_val;
            if (fabs(class_val - (double) class_val_i) > 0.0001) {
                cerr << "ERROR: Classes given are not integer\n";
                return 1;
            }

            if (class_val_i > num_class - 1) {
                cerr << "ERROR: Classes given outside range [0," << (num_class - 1) << "]\n";
                return 1;
            }
        }
    }

    // load the joint/copula densities

    vector<density_store> class_dist(num_class);

    for (int i = 0; i < num_class; i++) {

        if (joint_filenames[i] == "-") {
            joint_filenames[i] = "";
        }

        if (copula) {
            class_dist[i].load_files(joint_filenames[i], marginal_filenames[i]);
        } else {
            class_dist[i].load_files(joint_filenames[i], "");
        }

    }

    MT_random rand_gen(0);

    // Loop through the data and classify!
    // This should be changed to binary tree
    vector<int> result;
    for (int i = 0; i < test_N; i++) {

        vector<double> class_density(num_class, 0.0);
        int max_class = -1;
        double max_den = -c::inf;

        vector<double> point;
        if (confusion) {
            point = vector<double>(test_data[i].begin(), test_data[i].end() - 1);
        } else {
            point = test_data[i];
        }

        for (int c = 0; c < num_class; c++) {
            double map_den = 0.0;

            // bayes rule
            map_den = prior[c] * class_dist[c].compute_density(point, 0.5);

            class_density[c] = map_den;

            if (map_den > max_den) {
                max_class = c;
                max_den = map_den;
            }
        }

        // check for equality
        vector<int> equal_class;
        for (int c = 0; c < num_class; c++) {
            if (class_density[c] == max_den) {
                equal_class.push_back(c);
            }
        }

        if (equal_class.size() > 1) {
            // choose random one
            // this should be based on the prior [FIXME]
            max_class = equal_class[rand_gen.genrand_int_range(0, equal_class.size() - 1)];
        }

        // Print it out
        cout << max_class;
        for (int c = 0; c < num_class; c++) {
            cout << '\t' << class_density[c];
        }
        cout << '\n';

        result.push_back(max_class);
    }

    if (confusion) {
        // generate confusion matrix
        vector<vector<int> > conf_mat(num_class, vector<int>(num_class, 0));
        int num_correct = 0;
        for (int i = 0; i < test_N; i++) {
            int correct = (int) test_data[i].back();
            int predict = result[i];

            if (correct == predict) {
                num_correct++;
            }

            conf_mat[correct][predict]++;
        }

        // output classification rate
        cerr << "Classification Rate: " << ((double) num_correct / (double) test_N) << '\n';

        // output confusion matrix
        cerr << " ,";
        for (int i = 0; i < num_class; i++) {
            cerr << i << ',';
        }
        cerr << '\n';
        for (int i = 0; i < num_class; i++) {
            cerr << i << ',';
            for (int j = 0; j < num_class; j++) {
                cerr << conf_mat[i][j] << ',';
            }
            cerr << '\n';
        }
    }

    return 0;
}

// the old way of getting density from the partition... has some numerical issues

int density_old(vector<string> params) {
    string usage_text = "Usage: " + c::PROG_NAME
            + " density_old <MAP_partitions> <sample_data> \n"
            + "     MAP_partitions -- MAP partitions of a distribution\n"
            + "        sample_data -- Each row one data point (Restricted to [0,1] cube)\n "
            + "Get the density at locations.\n";

    if (params.size() != 2) {
        cerr << usage_text << endl;
        return 3;
    }


    vector<vector<double> > test_data;
    vector<uint32_t> test_skipped;
    read_data(params[1],test_data,test_skipped);

    int test_N = (int) test_data.size();
    int dim = (int) test_data[0].size();

    cerr << test_N << " data points in " << dim << " dimensions.\n";

    vector<vector<double> > MAP_dist;
    vector<uint32_t> MAP_skipped;
    read_data(params[0],MAP_dist,MAP_skipped,true);

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

int density(vector<string> params) {
    string usage_text = "Usage: " + c::PROG_NAME
            + " density <sample_data> <joint_den/copula_den> [marginal_den]  \n"
            + "      sample_data -- Each row one data point (Restricted to [0,1] cube)\n "
            + "        joint_den -- Joint density (of copula) from OPT, use - for uniform\n "
            + "     marginal_den -- Marginal densities from OPT, for copula [Optional]\n "
            + "Gets the density at each sample point. If marginal density is given\n"
            + "then a copula transform is performed with the marginals\n";

    if (params.size() < 2 || params.size() > 3) {
        cerr << usage_text << endl;
        return 3;
    }
    
    
    vector<vector<double> > test_data;
    vector<uint32_t> skipped;
    read_data(params[0],test_data,skipped);
    
    int test_N = (int) test_data.size();
    int dim = (int) test_data[0].size();

    cerr << test_N << " data points in " << dim << " dimensions.\n";

    string joint_filename = params[1];

    if (joint_filename == "-") {
        joint_filename = "";
    }

    string marginal_filename = "";
    bool copula = false;
    if (params.size() == 3) {
        copula = true;
        marginal_filename = params[2];
    }

    // load the joint/copula densities

    density_store dens(joint_filename, marginal_filename);

    // Loop through the data
    for (vector<vector<double> >::iterator it = test_data.begin();
            it != test_data.end(); it++) {
        // for each data point

        double den = dens.compute_density(*it, 0);
        cout << std::scientific << den << '\n';
    }

    return 0;
}

        
int normalize(vector<string> params) {
    string usage_text = "Usage: " + c::PROG_NAME
            + " normalize <-n/-u> <sample_data> <prefix>\n"
            + "            -n/-u -- Normalize/un-Normalize"
            + "      sample_data -- Each row one data point\n "
            + "           prefix -- [for -n] create <prefix>_norm.txt and create <prefix>_limits.txt\n "
            + "                  -- [for -u] create <prefix>_orig.txt and read <prefix>_limits.txt\n "
            + "Normalize data points to [0,1] cube\n";

    if (params.size() != 3) {
        cerr << usage_text << endl;
        return 3;
    }

    bool run_norm = true; 
    if(params[0] == "-n"){
        run_norm = true;
    }else if(params[0] == "-u"){
        run_norm = false;
    }else{
        cerr << "Error: Unknown option\n";
        cerr << usage_text << endl;
        return 3;
    }
    
    string prefix = params[2];
    string temp = "";
    
    if(!run_norm){
        // check for existence of limits file
        temp = prefix+"_limits.txt";
        ifstream infile(temp.c_str(),ios::in);
        if(!infile.is_open()){
            infile.close();
            cerr << "Error: Cannot open limits file: " << temp <<"\n";
            return 1;
        }
    }
    
    vector<vector<double> > data;
    vector<uint32_t> skipped;
    read_data(params[1],data,skipped, false, true);

    int N = (int) data.size();
    int dim = (int) data[0].size();

    cerr << N << " data points in " << dim << " dimensions.\n";

    if (run_norm) {
        // if normalize
        // compute limits
        vector<double> min_vals(dim,c::inf);
        vector<double> max_vals(dim,-c::inf);
        for (int i = 0;i < N;i++){
            for (int j = 0;j < dim;j++){
                if(min_vals[j] > data[i][j]){
                    min_vals[j] = data[i][j];
                }
                if(max_vals[j] < data[i][j]){
                    max_vals[j] = data[i][j];
                }
            }
        }

        // output limits
        temp = prefix+"_limits.txt";
        ofstream outfile;
        outfile.open(temp.c_str(),ios::out);
        
        outfile << dim << '\n';
        for(int i = 0;i<dim;i++){
            outfile << std::scientific <<  min_vals[i] << " " << max_vals[i] << '\n';
        }
        
        outfile.close();

        // output normalized values
        temp = prefix+"_norm.txt";
        outfile.open(temp.c_str(), ios::out);
        outfile << std::scientific;
        
        for (int i = 0; i < N; i++) {

            for (int j = 0; j < dim; j++) {
                double norm_val = 0.0;
                double diff = (max_vals[j] - min_vals[j]);
                if (diff != 0) {
                    norm_val = (data[i][j] - min_vals[j]) / diff;
                } else {
                    norm_val = (data[i][j] - min_vals[j]);
                }
                outfile << std::scientific << norm_val;
                if(j != dim-1) outfile  << ' ';
            }
            if(i != N-1)outfile << '\n';
        }
        
        outfile.close();

    } else {
        // if unnormalize
        // read in limits
        temp = prefix+"_limits.txt";
        ifstream infile;
        infile.open(temp.c_str(),ios::in);
        int num_vals = 0;
        infile >> num_vals;
        
        if(num_vals != dim){
            cerr << "Error: Dimension mismatch: " << num_vals << ',' << dim <<  '\n';
        }
        
        vector<double> min_vals(dim,c::inf);
        vector<double> max_vals(dim,-c::inf);
        
        for(int i = 0;i<dim;i++){
            infile >> min_vals[i] >> max_vals[i];
        }
        infile.close();

        // output un-normalized values
        temp = prefix + "_orig.txt";
        ofstream outfile;
        outfile.open(temp.c_str(), ios::out);
        outfile << std::scientific;

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < dim; j++) {
                double orig_val = 0.0;
                double diff = (max_vals[j] - min_vals[j]);
                if (diff != 0) {
                    orig_val = (data[i][j]*diff) + min_vals[j];
                } else {
                    orig_val = (data[i][j] + min_vals[j]);
                }
                outfile << std::scientific << orig_val;
                if (j != dim - 1) outfile << ' ';
            }
            if (i != N - 1)outfile << '\n';
        }

        outfile.close();
    }
    return 0;
}

int demean(vector<string> params) {
    string usage_text = "Usage: " + c::PROG_NAME
            + " demean <sample_data>\n"
            + "      sample_data -- Each row one data point\n "
            + "Remove the mean from the data\n";

    if (params.size() != 1) {
        cerr << usage_text << endl;
        return 3;
    }
    
    vector<vector<double> > data;
    vector<uint32_t> skipped;
    read_data(params[0],data,skipped, false, true);

    vector<string> ll = split(params[0],'.');
    string prefix = ll[0];
    
    int N = (int) data.size();
    int dim = (int) data[0].size();

    cerr << N << " data points in " << dim << " dimensions.\n";

    vector<double> means(dim, 0);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            means[j] += data[i][j]/N;
        }
    }

    string temp = prefix + "_demean.txt";
    ofstream outfile;
    outfile.open(temp.c_str(), ios::out);
    outfile << ios::scientific;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            outfile << data[i][j]-means[j]; 
            if (j != dim - 1) outfile  << ' ';
        }
        if (i != N - 1)outfile << '\n';
    }

    outfile.close();
    
    return 0;
}


int opt_comp(vector<string> params) {

    string usage_text = "Usage: " + c::PROG_NAME + " opt_comp <percent_points> <data_file_1> <data_file_2> <output_name>\n"
            + "       percent_points -- Ratio of total data to stop at (0.01 = 1%, or 2 = 2 points)\n"
            + "            data_file -- One sample each row (Restricted to [0,1] cube)\n"
            + "          output_name -- Result output to <output_name>.den\n"
            + "Log output to STDERR \n"
            + "Run full-OPT, very fast but uses a lot of memory. If dimension \n"
            + "greater than 4 use the other methods. Best choice for 1 or 2 dimensional data.\n"
            + "Recommend always stopping at 2 points.\n";

    if (params.size() != 4) {
        cerr << usage_text << endl;
        return 3;
    }

    string out_filename = params[3];

    vector<vector<double> > data_1;
    vector<uint32_t> skipped_1;
    read_data(params[1],data_1,skipped_1, false);
    
    vector<vector<double> > data_2;
    vector<uint32_t> skipped_2;
    read_data(params[2],data_2,skipped_2, false);

    int dim = (int) data_1[0].size();

    if (dim != (int) data_2[0].size()) {
        cerr << "Error: dimension mismatch between data\n";
    }

    uint32_t N[2] = {(uint32_t) data_1.size(), (uint32_t) data_2.size()};

    cerr << "Running Two Sample Comparison Full-OPT" << '\n';
    cerr << "Data_1: " << N[0] << " points in " << dim << " dimensions.\n";
    cerr << "Data_2: " << N[1] << " points in " << dim << " dimensions.\n";

    double stop_ratio = strTo<double>(params[0]);
    int stop_points = 0;

    if (stop_ratio < 1) {
        stop_points = (int) (min(N[0], N[1]) * stop_ratio);
    } else {
        stop_points = (int) stop_ratio;
    }
    if (stop_points < 1) stop_points = 1;

    cerr << "Stopping at " << stop_points << " points.\n";

    double total_time = 0.0;

    mu_timer mt;

    {
        cerr << "Constructing OPT from data_2\n";
        opt_tree opt_slow(dim, stop_points, 1000);
        mt.reset();
        opt_slow.construct_full_tree(data_2);
        total_time += mt.elapsed_time();
        mt.print_elapsed_time(cerr, "data_2 OPT tree");

        double orig_lphi = opt_slow.get_lphi();
        cerr << "Log Phi: " << opt_slow.get_lphi() << endl;

        mt.reset();
        map_tree map_region_tree(N[0], dim);
        opt_region_hash<uint32_t> map_regions(20);
        opt_slow.construct_MAP_tree(map_region_tree, map_regions, N[0]);
        total_time += mt.elapsed_time();
        mt.print_elapsed_time(cerr, "MAP tree");

        //print_MAP_density(cerr, map_regions.get_regions(),
        //            map_region_tree.get_ra(), map_region_tree.get_num_points());

        cerr << "Comparing to data_1\n";
        opt_tree opt_slow_comp(dim, stop_points, 1000, &opt_slow);
        mt.reset();
        opt_slow_comp.construct_full_tree(data_1);
        total_time += mt.elapsed_time();
        mt.print_elapsed_time(cerr, "data_1 comp OPT tree");

        double comp_lphi = opt_slow_comp.get_lphi();
        cerr << "comp Log Phi: " << opt_slow_comp.get_lphi() << endl;

        double lphi_ratio = orig_lphi - comp_lphi;
        cerr << "lphi_ratio: " << lphi_ratio - c::l2 << " = " << exp(lphi_ratio) / 2 << '\n';

        mt.reset();
        map_tree comp_map_region_tree(N[0], dim);
        opt_region_hash<uint32_t> comp_map_regions(20);
        opt_slow_comp.construct_MAP_tree(comp_map_region_tree, comp_map_regions, N[0]);
        total_time += mt.elapsed_time();
        mt.print_elapsed_time(cerr, "comp MAP tree");

        //print_MAP_density(cerr, comp_map_regions.get_regions(),
        //            comp_map_region_tree.get_ra(), comp_map_region_tree.get_num_points());

    }

    {
        cerr << "Constructing OPT from data_1\n";
        opt_tree opt_slow(dim, stop_points, 1000);
        mt.reset();
        opt_slow.construct_full_tree(data_1);
        total_time += mt.elapsed_time();
        mt.print_elapsed_time(cerr, "data_1 OPT tree");

        double orig_lphi = opt_slow.get_lphi();
        cerr << "Log Phi: " << opt_slow.get_lphi() << endl;

        cerr << "Comparing to data_2\n";
        opt_tree opt_slow_comp(dim, stop_points, 1000, &opt_slow);
        mt.reset();
        opt_slow_comp.construct_full_tree(data_2);
        total_time += mt.elapsed_time();
        mt.print_elapsed_time(cerr, "data_2 comp OPT tree");

        double comp_lphi = opt_slow_comp.get_lphi();
        cerr << "comp Log Phi: " << opt_slow_comp.get_lphi() << endl;

        double lphi_ratio = orig_lphi - comp_lphi;
        cerr << "lphi_ratio: " << lphi_ratio - c::l2 << " = " << exp(lphi_ratio) / 2 << '\n';

        mt.reset();
        map_tree map_region_tree(N[1], dim);
        opt_region_hash<uint32_t> map_regions(20);
        opt_slow_comp.construct_MAP_tree(map_region_tree, map_regions, N[1]);
        total_time += mt.elapsed_time();
        mt.print_elapsed_time(cerr, "comp MAP tree");

        //print_MAP_density(cerr, map_regions.get_regions(),
        //            map_region_tree.get_ra(), map_region_tree.get_num_points());
    }


    return 0;
}

int coopt(vector<string> params) {

    string usage_text = "Usage: " + c::PROG_NAME + " coopt <percent_points> <data_file_1> <data_file_2> <output_name>\n"
            + "       percent_points -- Ratio of total data to stop at (0.01 = 1%, or 2 = 2 points)\n"
            + "            data_file -- One sample each row (Restricted to [0,1] cube)\n"
            + "          output_name -- Result output to <output_name>.den\n"
            + "Log output to STDERR \n"
            + "Run full co-OPT, very fast but uses a lot of memory. If dimension \n"
            + "greater than 4 use the other methods. Best choice for 1-3 dimensional data.\n"
            + "Recommend always stopping at 2 points.\n";

    if (params.size() != 4) {
        cerr << usage_text << endl;
        return 3;
    }

    string out_filename = params[3];

    vector<vector<double> > data[2] = {vector<vector<double> >(),vector<vector<double> >()};
    vector<uint32_t> skipped[2] = {vector<uint32_t>() ,vector<uint32_t>() };
    
    read_data(params[1], data[0],skipped[0],false);
    read_data(params[2], data[1],skipped[1],false);
            
    int dim = (int) data[0][0].size();

    if (dim != (int) data[1][0].size()) {
        cerr << "Error: dimension mismatch between data\n";
    }

    uint32_t N[2] = {(uint32_t) data[0].size(), (uint32_t) data[1].size()};

    cerr << "Running Two Sample Comparison Full-OPT" << '\n';
    cerr << "Data_1: " << N[0] << " points in " << dim << " dimensions.\n";
    cerr << "Data_2: " << N[1] << " points in " << dim << " dimensions.\n";

    double stop_ratio = strTo<double>(params[0]);
    int stop_points = 0;

    if (stop_ratio < 1) {
        stop_points = (int) (min(N[0], N[1]) * stop_ratio);
    } else {
        stop_points = (int) stop_ratio;
    }
    if (stop_points < 1) stop_points = 1;

    cerr << "Stopping at " << stop_points << " points.\n";

    double total_time = 0.0;

    mu_timer mt;

    {
        cerr << "Constructing coupling OPT tree\n";
        copt_tree opt_slow_comp(dim, stop_points, 1000);
        mt.reset();
        opt_slow_comp.construct_full_tree(data);
        total_time += mt.elapsed_time();
        mt.print_elapsed_time(cerr, "coupling OPT tree");

        cerr << "Log coupling prob: " << opt_slow_comp.get_log_coupling_prob() << endl;

        mt.reset();
        map_tree comp_map_region_tree(N[0] + N[1], dim);
        opt_region_hash<uint32_t> comp_map_regions(20);
        opt_slow_comp.construct_MAP_tree(comp_map_region_tree, comp_map_regions, N[0] + N[1]);
        total_time += mt.elapsed_time();
        mt.print_elapsed_time(cerr, "comp MAP tree");

        print_MAP_density(cerr, comp_map_regions.get_regions(),
                comp_map_region_tree.get_ra(), comp_map_region_tree.get_num_points());

    }

    return 0;
}

int coopt_scan_old(vector<string> params) {

    string usage_text = "Usage: " + c::PROG_NAME + " coopt_scan_old <percent_points> <max_depth> <window_size> <scan_resolution> <data_file>\n"
            + "       percent_points -- Ratio of total data to stop at (0.01 = 1%, or 2 = 2 points)\n"
            + "            max_depth -- Maximum depth in the tree, smallest region area will be 2^(-max_depth)\n"
            + "          window_size -- Size of window to scan for changes, break at window_size/2\n"
            + "      scan_resolution -- Number of data points to move each time\n"
            + "            data_file -- One sample each row (Restricted to [0,1] cube)\n"
            + "Log posterior coupling probability to STDOUT\n"
            + "Scans through data_file and runs a coupling OPT for each window with the split at window_size/2"
            + "Outputs one partition for each coupling OPT \n"
            + "Runs full-OPT, uses a lot of memory for high dimensions. Not recommended for dimension \n"
            + "greater than 5. Best choice for 1 - 3 dimensional data.\n";

    if (params.size() != 5) {
        cerr << usage_text << endl;
        return 3;
    }
    
    int max_depth = strTo<int>(params[1]);

    int window_size = strTo<int>(params[2]);
    int scan_res = strTo<int>(params[3]);
    
    vector<vector<double> > data;
    vector<uint32_t> skipped;
    read_data(params[4],data,skipped, false);

    int dim = (int) data[0].size();
    int N = data.size();

    cerr << "Running coupling OPT scan" << '\n';
    cerr << "Data: " << N << " points in " << dim << " dimensions.\n";

    double stop_ratio = strTo<double>(params[0]);
    int stop_points = 0;

    if (stop_ratio < 1) {
        stop_points = (int) (N * stop_ratio);
    } else {
        stop_points = (int) stop_ratio;
    }
    if (stop_points < 1) stop_points = 1;

    cerr << "Stopping at " << stop_points << " points.\n";

    double total_time = 0.0;

    mu_timer mt;

    // loop to split the dataset
    if (N < window_size) {
        cerr << "N too small\n";
        return 1;
    }

    int idx = 0;
    while ((idx + window_size) < N) {
        // this splitting is really inefficient
        
        mt.reset();
        vector<vector<double> > split_data[2] = {vector<vector<double> >(), vector<vector<double> >()};
        split_data[0].reserve(window_size / 2 + 1);
        split_data[1].reserve(window_size / 2 + 1);

        for (int i = 0; i < window_size; i++) {
            if (i < window_size / 2) {
                split_data[0].push_back(data[idx + i]);
            } else {
                split_data[1].push_back(data[idx + i]);
            }
        }
        total_time += mt.elapsed_time();
        mt.print_elapsed_time(cerr, "Split data");

        cerr << "idx: " << idx << " - " << idx + window_size << '\n';

        mt.reset();
        copt_tree opt_slow_comp(dim, stop_points, max_depth);
        opt_slow_comp.construct_full_tree(split_data);
        total_time += mt.elapsed_time();
        mt.print_elapsed_time(cerr, "coupling OPT tree");

        double coup_prob = opt_slow_comp.get_log_coupling_prob();
        cerr << "Log coupling prob: " << coup_prob << endl;

        /*
        cerr << "partitions!!~!@~!\n";
        vector<pair<opt_region, ctree_node*> > print_stuff = opt_slow_comp.region_cache.get_regions();
        for (int i = 0; i < (int) print_stuff.size(); i++) {
            print_stuff[i].first.print_region_limits();
            ctree_node* node = print_stuff[i].second;
            int count[2];
            node->get_count(count);
            cerr << " ; " << print_stuff[i].second;
            cerr << " : " << count[0] << "," << count[1];
            cerr << " : " << node->get_lP() << " : " << node->get_lphi();
            cerr << '\n';
        }
        */
        
        cout << idx + window_size / 2 << " " << coup_prob << '\n';

        mt.reset();
        map_tree comp_map_region_tree(N, dim);
        opt_region_hash<uint32_t> comp_map_regions(20);
        opt_slow_comp.construct_MAP_tree(comp_map_region_tree, comp_map_regions, N);
        total_time += mt.elapsed_time();
        mt.print_elapsed_time(cerr, "comp MAP tree");

        string temp = "old_partition_" + toStr<int>(idx + window_size / 2) + ".txt";
        ofstream out_file(temp.c_str());
        print_MAP_density(out_file, comp_map_regions.get_regions(),
                comp_map_region_tree.get_ra(), comp_map_region_tree.get_num_points());
        print_MAP_density(cerr, comp_map_regions.get_regions(),
                comp_map_region_tree.get_ra(), comp_map_region_tree.get_num_points());
        out_file.close();


        idx = idx + scan_res;
        if (((idx + window_size) >= N)&& (idx + window_size < N + scan_res - 1)) {
            idx = N - window_size - 1;
        }
    }

    cerr.unsetf(ios::scientific);
    cerr << "Total time: " << total_time << '\n';


    return 0;
}

int coopt_scan(vector<string> params) {

    string usage_text = "Usage: " + c::PROG_NAME + " coopt_scan <percent_points> <max_depth> <window_size> <output_interval> <data_file>\n"
            + "       percent_points -- Ratio of total data to stop at (0.01 = 1%, or 2 = 2 points)\n"
            + "            max_depth -- Maximum depth in the tree, smallest region area will be 2^(-max_depth)\n"
            + "          window_size -- Size of window to scan for changes, break at window_size/2\n"
            + "      scan_resolution -- Number of data points to move each time\n"
            + "            data_file -- One sample each row (Restricted to [0,1] cube)\n"
            + "Log posterior coupling probability to STDOUT\n"
            + "Scans through data_file and runs a coupling OPT for each window with the split at window_size/2"
            + "Outputs one partition for each coupling OPT \n"
            + "Runs full-OPT, uses a lot of memory for high dimensions. Not recommended for dimension \n"
            + "greater than 5. Best choice for 1 - 3 dimensional data.\n";

    if (params.size() != 5) {
        cerr << usage_text << endl;
        return 3;
    }

    int window_size = strTo<int>(params[2]);
    int half_wind = (window_size / 2);
    int output_interval = strTo<int>(params[3]);
    
    int max_depth = strTo<int>(params[1]);
    
    if (half_wind <= 0) {
        cerr << "Error: window size too small\n";
        return 1;
    }
    
    vector<vector<double> > data;
    vector<uint32_t> skipped;
    read_data(params[4],data,skipped, false);

    int dim = (int) data[0].size();
    int N = data.size();

    cerr << "Running coupling OPT scan" << '\n';
    cerr << "Maximum depth: " << max_depth << '\n';
    cerr << "Window size: " << window_size << '\n';
    cerr << "Half Window size: " << half_wind << '\n';
    cerr << "Data: " << N << " points in " << dim << " dimensions.\n";

    double stop_ratio = strTo<double>(params[0]);
    int stop_points = 0;

    if (stop_ratio < 1) {
        stop_points = (int) (window_size * stop_ratio);
    } else {
        stop_points = (int) stop_ratio;
    }
    if (stop_points < 1) stop_points = 1;

    cerr << "Stopping at " << stop_points << " points.\n";

    mu_timer mt;

    if (N < window_size) {
        cerr << "N too small\n";
        return 1;
    }

    // run copt for the first window
    cerr << "idx: " << half_wind << '\n';
    online_copt_tree* online_comp = new online_copt_tree(dim, stop_points, max_depth, window_size);
    uint32_t start_idx[2] = {0, half_wind};
    uint32_t end_idx[2] = {half_wind - 1, window_size - 1};

    online_comp->construct_full_tree(data, start_idx, end_idx);

    /*
    cerr << "partitions!!~!@~!\n";
    vector<pair<opt_region, uint32_t> > print_stuff = online_comp->get_region_cache()->get_regions();
    region_allocator<online_ctree_node>* print_ra = online_comp->get_ra();
    for (int i = 0;i<(int)print_stuff.size();i++){
        print_stuff[i].first.print_region_limits();
        online_ctree_node* node = (*print_ra)[print_stuff[i].second];
        int count[2];
        node->get_count(count);
        cerr << " ; " << print_stuff[i].second;
        cerr << " : " << node->get_sequence_id();
        cerr << " : " << count[0] << "," << count[1];
        cerr << " : " << node->get_lP() << " : " << node->get_lphi();
        cerr << '\n';
    }*/
     

    // output the partition
    // need to change to output immediately rather than store in hash table
    map_tree comp_map_region_tree(N, dim);
    opt_region_hash<uint32_t> comp_map_regions(15);
    online_comp->construct_MAP_tree(comp_map_region_tree, comp_map_regions, window_size);

    string temp = "partition_" + toStr<int>(half_wind) + ".txt";
    ofstream out_file(temp.c_str());
    print_MAP_density(out_file, comp_map_regions.get_regions(),
            comp_map_region_tree.get_ra(), comp_map_region_tree.get_num_points());

    print_MAP_density(cerr, comp_map_regions.get_regions(),
            comp_map_region_tree.get_ra(), comp_map_region_tree.get_num_points());


    cerr << "lP: " << online_comp->get_lP() << '\n';
    cerr << "lphi: " << online_comp->get_lphi() << '\n';
    cerr << "log couping: " << online_comp->get_log_coupling_prob() << '\n';
    cout << half_wind << " " << online_comp->get_log_coupling_prob()  << '\n';
    out_file.close();

    int idx = 1;
    //bool any_diff = true; // need to keep a dequeue for this
    uint32_t del_nodes[2] = {0,0};
    while ((idx + window_size) < N) {
        // swap out the middle data point

        cerr << "idx: " << half_wind + idx << '\n';

        //if(any_diff){
            del_nodes[0] =  uint32_t(idx - 1);
            del_nodes[1] = half_wind + idx - 1;
        //}
        
        // [ add0,  del0,  add1,  del1]
        uint32_t pts[4] = {half_wind + idx - 1, del_nodes[0], window_size + idx - 1, del_nodes[1]};
        
        //any_diff = false;
        //for(int i = 0;(!any_diff) && i<2;i++){
        //    for(int j = 0;(!any_diff) && j<dim;j++){
        //        if(data[pts[(2*i)+0]][j] != data[pts[(2*i)+1]][j]){
        //            any_diff = true;
        //        }
        //    }
       // }

        //if (any_diff) {
            online_comp->update_points(data, pts, idx);
            online_comp->prune_tree(data, pts, idx); // this must be run every iteration (not yet tested for other case)
        //}
        /*
        cerr << "partitions!!~!@~!\n";
        vector<pair<opt_region, uint32_t> > print_stuff = online_comp->get_region_cache()->get_regions();
        region_allocator<online_ctree_node>* print_ra = online_comp->get_ra();
        for (int i = 0; i < (int) print_stuff.size(); i++) {
            print_stuff[i].first.print_region_limits();
            online_ctree_node* node = (*print_ra)[print_stuff[i].second];
            int count[2];
            node->get_count(count);
            cerr << " ; " << print_stuff[i].second;
            cerr << " : " << node->get_sequence_id();
            cerr << " : " << count[0] << "," << count[1];
            cerr << " : " << node->get_lP() << " : " << node->get_lphi();
            cerr << '\n';
        }*/
        
        
        cout << half_wind + idx << " " << online_comp->get_log_coupling_prob()  << '\n';
        
        // output the partition
        if (((idx + window_size) == (N - 1)) || ((idx % output_interval) == 0)) {
            map_tree comp_map_region_tree(N, dim);
            opt_region_hash<uint32_t> comp_map_regions(15);
            online_comp->construct_MAP_tree(comp_map_region_tree, comp_map_regions, window_size);

            string temp = "partition_" + toStr<int>(half_wind + idx) + ".txt";
            ofstream out_file(temp.c_str());
            print_MAP_density(out_file, comp_map_regions.get_regions(),
                    comp_map_region_tree.get_ra(), comp_map_region_tree.get_num_points());
            print_MAP_density(cerr, comp_map_regions.get_regions(),
                    comp_map_region_tree.get_ra(), comp_map_region_tree.get_num_points());
            cerr << "lP: " << online_comp->get_lP() << '\n';
            cerr << "lphi: " << online_comp->get_lphi() << '\n';
            cerr << "log couping: " << online_comp->get_log_coupling_prob() << '\n';
            out_file.close();
        }

        idx = idx + 1;
    }

    cerr.unsetf(ios::scientific);
    cerr << "Total time: " << mt.elapsed_time() << "s\n";
    delete online_comp;


    return 0;
}

int print_partitions(vector<string> params) {
    string usage_text = "Usage: " + c::PROG_NAME
            + " print <density_file> [output_prefix]\n"
            + "      density_file -- The .den file, can have multiple densities\n"
            + "     output_prefix -- output_prefix.txt, one per density otherwise to STDOUT\n"
            + "Outputs the partitions in standard format\n"
            + "start_0 end_0 ... start_n-1 end_n-1 density count\n";

    if (params.size() < 1 || params.size() > 2) {
        cerr << usage_text << endl;
        return 3;
    }

    bool output_file = false;
    string out_prefix = "";
    if (params.size() == 2) {
        out_prefix = params[1];
        output_file = true;
    }

    // load the joint/copula densities
    density_store dens(params[0]);

    // print the densities
    if (!output_file) {
        dens.print_density(cout);
    } else {
        dens.print_density(out_prefix);
    }

    return 0;
}

// urgh, this benchmarking code can be improved

int bench(vector<string> params) {
    string usage_text = "Usage: " + c::PROG_NAME
            + " bench <num_points> <dim> \n"
            + "     num_points -- Number of random points to generate\n"
            + "     dim        -- Total dimension, we only count one though\n "
            + "Do benchmarking of counting!.\n";

    if (params.size() != 2) {
        cerr << usage_text << endl;
        return 3;
    }

    int N = strTo<int>(params[0]);
    int dim = strTo<int>(params[1]);

    MT_random randgen(0);

    vector<vector<double> > data(N, vector<double>(dim, 0.0));
    vector<int> data_one(N, 0);

    for (int i = 0; i < N; i++) {
        data_one[i] = i;
        for (int j = 0; j < dim; j++) {
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


int vec_quant_sam_quals(vector<string> params) {
    string usage_text = "Usage: " + c::PROG_NAME + " vec_quant_sam_quals <SAM/qual_file> <offset> <transfomred_quals> <den_file>\n"
            + "Offset typically 33 (illumina 1.8+,sanger) or 64 (illumina 1.3+)\n"
            + "replace quals in SAM file with new quantized ones!";


    if (params.size() != 4) {
        cerr << usage_text << endl;
        return 3;
    }

    string filename = params[0];
    int offset = strTo<int>(params[1]);
    
    string quals_file = params[2];
    string partitions = params[3];
    
    ifstream infile;

    int file_type = 2; // 1 is sam file, 2 is quals file
    
    infile.open(filename.c_str(), ios::in);
    if (!infile.is_open()) {
        cerr << "Error: cannot open input file: " << filename << '\n';
        return 3;
    }else{
        // determine file type
        // read a few lines from the file and if there are no '@' we will assume it is a quals file
        string line = "";
        int count = 1;
        while (!infile.eof()) {
            getline(infile, line);

            trim2(line);

            if ((int) line.length() < 1) {
                continue;
            }

            if (line[0] == '@') {
                file_type = 1;
                break;
            }
            if(count > 10){
                break;
            }
            count++;
        }
    }
    infile.close();
    
    infile.open(quals_file.c_str(), ios::in);
    if (!infile.is_open()) {
        cerr << "Error: cannot open quals file: " << quals_file << '\n';
        return 3;
    }
    infile.close();
    
    infile.open(partitions.c_str(), ios::in);
    if (!infile.is_open()) {
        cerr << "Error: cannot open partition file: " << partitions << '\n';
        return 3;
    }
    infile.close();
    
    // Read in transformed quals file
    vector<vector<double> > quals;
    vector<uint32_t> skipped_quals;
    read_data(quals_file, quals,skipped_quals,false);
    
    // Read in partitions
    density_store dens(partitions, "");
    
    
    // Read in quals from SAM file
    cerr << "Reading quals...\n";
    vector<vector<double> > sam_quals;
    
    string line = "";

    infile.open(filename.c_str(), ios::in);
    int idx = 0;

    while (!infile.eof()) {
        getline(infile, line);

        trim2(line);

        if ((int) line.length() < 1) {
            continue;
        }

        if (file_type == 1 && line[0] == '@') {
            continue;
        }

        vector<string> line_list = split(line);

        if (file_type == 1) {
            if (line_list.size() < 10) {
                cerr << "Bad line: " << line << '\n';
                continue;
            }

            // read in quals, need to invert

            // get flag
            int flag = strTo<int>(line_list[1]);

            string new_quals = line_list[10];

            // invert if necessary
            bool invert = false;

            if (!((flag & 4) && (flag & 8))) {
                // at least on read is mapped
                if (flag & 64) {
                    // first in pair
                    invert = flag & 16;
                } else if (flag & 128) {
                    // second in pair
                    invert = (!(flag & 16));
                }

            }

            if (invert) {
                reverse(new_quals.begin(), new_quals.end());
            }

            sam_quals.push_back(vector<double>(new_quals.length()));

            int len = (int) new_quals.length();
            for (int i = 0; i < len; i++) {
                sam_quals[idx][i] = (double) (new_quals[i] - offset);
            }
        } else {
            sam_quals.push_back(vector<double>(line_list.size()));

            int len = (int) line_list.size();
            for (int i = 0; i < len; i++) {
                sam_quals[idx][i] = strTo<double>(line_list[i]);
            }
        }
        idx++;
    }

    infile.close();

    cerr << "Assigning quals to regions...\n";
    vector<uint32_t> indexes(sam_quals.size());
    vector<vector<vector<double> > > regions;
    opt_region_hash<uint32_t> region_cache(12);
    
    uint32_t curr_idx = 0;
    for(size_t i = 0;i < quals.size();i++){
        opt_region reg = dens.find_region(quals[i]);
        pair<uint32_t,bool> out = region_cache.find(reg);
        if(!out.second){
            regions.push_back(vector<vector<double> >());
            regions[curr_idx].push_back(sam_quals[i]);
            region_cache.insert(reg,curr_idx);
            indexes[i] = curr_idx;
            curr_idx++;
        }else{
            regions[out.first].push_back(sam_quals[i]);
            indexes[i] = out.first;
        }
    }
    
    cerr << "Number of regions: " << regions.size() << '\n';
    
    
    cerr << "Computing means...\n";
    // compute medians or means or etc...
    vector<vector<double> > quantizers(regions.size());
    int dim = (int)regions[0][0].size();
    for(size_t i = 0;i<regions.size();i++){
        quantizers[i] = vector<double>(dim);
        int num_values = (int) regions[i].size();
        for(int j = 0;j<num_values;j++){
            for(int k=0;k<dim;k++){
                quantizers[i][k] += ((double)regions[i][j][k]/ num_values);
            }
        }
    }
    
    cerr << "output the quantizers and indexes...\n";
    {
        string temp = quals_file + "_quant.txt";
        ofstream outfile(temp.c_str(), ios::out);

        for (size_t i = 0; i < quantizers.size(); i++) {
            for (int k = 0; k < dim; k++) {
                outfile << (char)((int)quantizers[i][k] + 33);
            }
            outfile << '\n';
        }

        outfile.close();
        
        
        temp = quals_file + "_idx.txt";
        outfile.open(temp.c_str(), ios::out);

        for (size_t i = 0; i < indexes.size(); i++) {
            outfile.put((char)(indexes[i]));
            outfile.put((char)(indexes[i]>>8));
            outfile.put((char)(indexes[i]>>16));
        }

        outfile.close();
    }
    
    cerr << "output the original quals in 2-bit encoding...\n";
    {
        string temp = filename + "_2bit.txt";
        ofstream outfile(temp.c_str(), ios::out);
        
        for(int i = 0;i<(int)sam_quals.size();i++){
            char val = 0;
            for (int k = 0; k < dim; k++) {
                int qual = (int) sam_quals[i][k];
                int shift = 2 * (k % 4);
                if (qual < 10) {
                    val += 0 << shift;
                } else if (qual < 20) {
                    val += 1 << shift;
                } else if (qual < 30) {
                    val += 2 << shift;
                } else {
                    val += 3 << shift;
                }

                if (k % 4 == 3) {
                    outfile.put(val);
                    val = 0;
                }
            }
            if(dim%4 != 0){
                outfile.put(val);
            }
        }
        
        outfile.close();
    }
    
    // 4 9 14 19 24 29 34 41
    cerr << "output the original quals in 3-bit encoding...\n";
    {
        string temp = filename + "_3bit.txt";
        ofstream outfile(temp.c_str(), ios::out);
        
        for(int i = 0;i<(int)sam_quals.size();i++){
            char val = 0;
            for (int k = 0; k < dim; k++) {
                int qual = (int) sam_quals[i][k];
                int shift = 3 * (k % 2);
                if (qual < 5) {
                    val += 0 << shift;
                } else if (qual < 10) {
                    val += 1 << shift;
                } else if (qual < 15) {
                    val += 2 << shift;
                } else if (qual < 20) {
                    val += 3 << shift;
                } else if (qual < 25) {
                    val += 4 << shift;
                } else if (qual < 30) {
                    val += 5 << shift;
                } else if (qual < 35) {
                    val += 6 << shift;
                } else {
                    val += 7 << shift;
                }

                if (k % 2 == 1) {
                    outfile.put(val);
                    val = 0;
                }
            }
            if(dim%2 != 0){
                outfile.put(val);
            }
        }
        
        outfile.close();
    }
    

    if (file_type == 1) {
        line = "";

        infile.open(filename.c_str(), ios::in);
        idx = 0;
        cerr << "outputting samfile with new quals...\n";
        while (!infile.eof()) {
            getline(infile, line);

            trim2(line);

            if ((int) line.length() < 1) {
                continue;
            }

            if (line[0] == '@') {
                cout << line << '\n';
                continue;
            }

            vector<string> line_list = split(line);

            if (line_list.size() < 10) {
                cerr << "Bad line: " << line << '\n';
                cout << line << '\n';
                continue;
            }


            // get flag
            int flag = strTo<int>(line_list[1]);

            string new_quals = line_list[10];

            // Replace the quals
            int len = (int) new_quals.length();
            for (int i = 0; i < len; i++) {
                new_quals[i] = (char) (quantizers[indexes[idx]][i] + offset);
            }
            idx++;

            // invert if necessary
            bool invert = false;

            if (!((flag & 4) && (flag & 8))) {
                // at least on read is mapped
                if (flag & 64) {
                    // first in pair
                    invert = flag & 16;
                } else if (flag & 128) {
                    // second in pair
                    invert = (!(flag & 16));
                }

            }

            if (invert) {
                reverse(new_quals.begin(), new_quals.end());
            }

            // print it out

            line_list[10] = new_quals;

            cout << line_list[0];

            for (int i = 1; i < (int) line_list.size(); i++) {
                cout << '\t' << line_list[i];
            }

            cout << "\n";
        }

        infile.close();
    }else{
        for(size_t i = 0;i<indexes.size();i++){
            for(int j = 0;j<dim;j++){
                cout << floor(quantizers[indexes[i]][j]);
                if(j!=dim) cout << ' ';
            }
            cout << "\n";
        }
    }

    return 0;
}

void print_usage_and_exit() {
    cerr << "Usage: " + c::PROG_NAME + " <option>" << "\n";
    cerr << "Options:" << "\n";
    cerr << "-== Density Estimation and Related ==-" << '\n';
    cerr << "  opt        -- MAP partitions from full OPT" << "\n";
    cerr << "  llopt      -- MAP partitions from LL-OPT" << "\n";
    //cerr << "  lsopt      -- MAP partitions from LL-sampled OPT [experimental]" << "\n";
    cerr << "  dfopt      -- MAP partitions from depth-first OPT" << "\n";
    //cerr << "  disopt     -- MAP partitions from discrepancy-LL-OPT [experimental]" << "\n";
    cerr << "  copula     -- Copula transform with full OPT" << "\n";
    cerr << "\n";
    cerr << "-== Two sample comparison ==-" << '\n';
    cerr << "  coopt       -- Two sample comparison with the full co-OPT [experimental]" << "\n";
    cerr << "  coopt_scan  -- Change-point detection with online full co-OPT [experimental]" << "\n";
    cerr << "\n";
    cerr << "-== Other tools ==-" << '\n';
    cerr << "  hell_dist  -- Compute sample Hellinger distance from a known density" << "\n";
    cerr << "  classify   -- Do classification with MAP partitions" << "\n";
    cerr << "  density    -- Get the density at particular points" << "\n";
    cerr << "  print      -- Print partitions from a .den file" << "\n";
    cerr << "  normalize  -- Normalize data points to [0,1]" << "\n";
    cerr << "  demean     -- Subtract mean from each column" << "\n";
    //cerr << "  bench       -- Benchmark counting speed [temporary]" << "\n";
    exit(2);
}




