/*
 * dfopt.h
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



#ifndef DFOPT_H
#define	DFOPT_H

#include "general_utils.h"
#include "opt_utils.h"

struct df_working_unit_t {
    vector<vector<double> > data;
    current_region curr_reg;
    opt_region working_reg;

    df_working_unit_t(){

    }

    df_working_unit_t(int num_children) {
        curr_reg.init(num_children);
        working_reg.init(num_children);
    }

    df_working_unit_t(vector<vector<double> > &data, current_region &curr_reg
            , opt_region &working_reg) {


        this->data = data;
        this->curr_reg = curr_reg;
        this->working_reg = working_reg;
    }
};

struct node_t {
    double lphi[2];
    int count[2];

    node_t() {
        lphi[0] = -c::inf;
        lphi[1] = -c::inf;
        count[0] = -1;
        count[1] = -1;
    }
};

struct dfpile_t {
    vector<vector<double> > data;
    vector<node_t> val; // lphi and count of children
    int dim;
    int cut;

    dfpile_t() {
        dim = -1;
        cut = -1;
    }

    dfpile_t(int num_children) {
        dim = -1;
        cut = -1;
        val.resize(num_children);
    }

    dfpile_t(vector<vector< double> > &data) {
        this->data = data;
        dim = -1;
        cut = -1;
        val.resize(data.size());
    }
};

class dfopt_tree {
private:

    int num_children;

    
    double top_count_ratio;
    int top_max_depth;
    double count_ratio;
    int max_depth;

    void init(int num_children, double top_count_ratio,double count_ratio,
                int top_levels,int levels) {
        this->num_children = num_children;
        this->count_ratio = count_ratio;
        this->max_depth = levels;
        this->top_count_ratio = top_count_ratio;
        this->top_max_depth = top_levels;
    }

public:

    dfopt_tree(int num_children, double top_count_ratio,double count_ratio,
                int top_levels,int levels) {
        init(num_children, top_count_ratio, count_ratio, top_levels,levels);
    }

    dfopt_tree(int num_children) {

        init(num_children, 0.001, 0.002, 100,4);

    }

    void store_region(uint32_t map_node, opt_region_hash<uint32_t> &map_regions,
            region_allocator<map_tree_node> *map_ra, opt_region &working_reg, int count, int depth) {

        map_regions.insert(working_reg, map_node);

        (*map_ra)[map_node]->set_area(-depth);
        (*map_ra)[map_node]->set_count(count);


    }

    double compute_lphi(dfpile_t &curr_pile, int depth, gamma_table &gt) {
        
        //cerr << "lphi depth:" << depth << '\n';
        
        vector<double> lphi_list;
        lphi_list.reserve(num_children + 1);

        int count = curr_pile.data.size();

        //cerr << "count: " << count << '\n';
        
        double max_val = (count * depth * c::l2) - c::l2;
        lphi_list.push_back(max_val);

        //cerr << "1max_val: " << max_val << '\n';
        
        double ld = -log(num_children) - c::lpi - c::l2;

        for (int i = 0; i < num_children; i++) {
            // check for null
            if (isinf(curr_pile.val[i].lphi[0]) || isinf(curr_pile.val[i].lphi[1])) {
                cerr << "NULL child!!! " << i << ',' << depth << '\n';
                exit(2);
            }

            int child_1_count = curr_pile.val[i].count[0];
            int child_2_count = curr_pile.val[i].count[1];
            
            //cerr << "child_1_count(" << i <<"): " << child_1_count << '\n';
            //cerr << "child_2_count(" << i <<"): " << child_2_count << '\n';

            //cerr << "lphi0: " << curr_pile.val[i].lphi[0] << '\n';
            //cerr << "lphi1: " << curr_pile.val[i].lphi[1] << '\n';
            
            if (child_1_count < 0 || child_2_count < 0) {
                cerr << "neg count!!! " << i << ',' << depth << '\n';
                exit(2);
            }

            double val = ld;
            val += curr_pile.val[i].lphi[0] + curr_pile.val[i].lphi[1];
            val += gt.compute_lD2(count, child_1_count, child_2_count);

            //cerr << "1val: " << val << '\n';
            
            lphi_list.push_back(val);

            if (val > max_val) {
                max_val = val;
            }

        }

        double lphi = max_val;
        double sum = 0;
        for (int i = 0; i < (num_children + 1); i++) {
            sum += exp(lphi_list[i] - max_val);
        }

        if (sum > 0) lphi += log(sum);
        
        //cerr << "1lphi: " << lphi << '\n';
        

        return lphi;
    }



    double compute_dflphi(df_working_unit_t& wu, int top_depth, gamma_table &gt,
            int64_t &num_nodes, int64_t &num_zero_nodes) {


        int N = wu.data.size();

        int count_lim = 1;

        if (count_ratio < 1) {
            count_lim = (int) floor(N * count_ratio);
            if (count_lim < 1) {
                count_lim = 1;
            }
        }else{
            count_lim = (int)count_ratio;
        }

        //cerr << "count_lim: " << count_lim << '\n';

        vector<dfpile_t> pile;
        pile.push_back(dfpile_t(wu.data));
        pile[0].dim = 0;
        pile[0].cut = 0;

        current_region curr_reg = wu.curr_reg;
        opt_region working_reg = wu.working_reg;

        int depth = 0;

        bool done = false;
        double lphi_out = -c::inf;

        while (!done) {

            if (pile.size() == 0) {
                done = true;
                continue;
            }

            // print the pile
#ifdef DEBUG
            cerr << "---- " << depth << " , " << pile.size() << " ---\n";
            for (int i = 0; i < (int) pile.size(); i++) {
                cerr << "Data size: " << pile[i].data.size() << '\n';
                cerr << "Dim:       " << pile[i].dim << '\n';
                cerr << "Cut:       " << pile[i].cut << '\n';
                if(pile[i].data.size()<=3){
                    print_data(pile[i].data);
                }
            }
            cerr << "=====\n";
#endif
            int curr_dim = pile[depth].dim;
            int curr_cut = pile[depth].cut;

            if (!isinf(pile[depth].val[curr_dim].lphi[curr_cut])) {

                //cerr << "Increment cut: " << pile[depth].cut << " dim: "<< pile[depth].dim <<"\n";

                if (pile[depth].cut < c::cuts - 1) {
                    pile[depth].cut++;
                } else if (pile[depth].dim < num_children - 1) {
                    pile[depth].dim++;
                    pile[depth].cut = 0;
                } else {
                    // reached end of node!! back up

                    //cerr << "FULL depth:" << depth << " top_depth:" << top_depth << "\n";

                  
                    lphi_out = compute_lphi(pile[depth], depth + top_depth, gt);

                    if (depth > 0) {
                        int prev_dim = pile[depth - 1].dim;
                        int prev_cut = pile[depth - 1].cut;
                        pile[depth - 1].val[prev_dim].lphi[prev_cut] = lphi_out;
                    }

                    //cerr << "lphi_out: " << lphi_out << "\n";
                    
                    //cerr << "BackUP\n";

                    depth--;
                    pile.pop_back();
                    if (depth < 0) continue;

                    curr_reg.uncut(pile[depth].dim, pile[depth].cut);
                    working_reg.uncut(pile[depth].dim);

                }
                
                //cerr << "Now at:: cut: " << pile[depth].cut << " dim: "<< pile[depth].dim <<"\n";
                
            }else{


                vector<vector<double> > temp_data;
                bool is_diff = cut_region(pile[depth].data,temp_data,
                        curr_dim, curr_cut, curr_reg.get_lim(curr_dim));


                int curr_count = (int) temp_data.size();

                //cerr << "curr_dim: " << curr_dim << " curr_cut: " << curr_cut << '\n';
                //cerr << "curr_reg.get_lim(curr_dim): " << curr_reg.get_lim(curr_dim) << '\n';
                //cerr << "curr_reg.get_resolution(curr_dim): " << curr_reg.get_resolution(curr_dim) << '\n';

                //cerr << "curr_count: " << curr_count << '\n';

                pile[depth].val[curr_dim].count[curr_cut] = curr_count;


                if (curr_count <= count_lim || depth >= max_depth || !is_diff) {

                    //cerr << "leaf\n";

                    num_nodes++;
                    
                    if (curr_count <= count_lim) {
                        num_zero_nodes++;
                    }

                    if (num_nodes % 1000000 == 0) {
                        cerr << "Nodes(" << pile[0].dim << "):"
                                << num_nodes << " : " << num_zero_nodes
                                << " : " << (num_nodes - num_zero_nodes) << '\n';
                    }

                    //cerr << "lphi curr_count: " << curr_count << '\n';
                    //cerr << "depth: " << depth << '\n';
                    //cerr << "top_depth: " << top_depth << '\n';
                    //cerr << "leaf lphi: " << curr_count * (depth + top_depth + 1) * c::l2 << '\n';
                    
                    pile[depth].val[curr_dim].lphi[curr_cut] = curr_count * (depth + top_depth + 1) * c::l2;
                }else{


                    //cerr << "GO down\n";

                    curr_dim = pile[depth].dim;
                    curr_cut = pile[depth].cut;

                    //cerr << "curr_dim2: " << curr_dim << " curr_cut2: " << curr_cut << '\n';

                    // go down
                    depth++;

                    working_reg.cut(curr_dim, curr_cut);
                    curr_reg.cut(curr_dim, curr_cut);

                    pile.push_back(dfpile_t(temp_data));
                    pile[depth].dim = 0;
                    pile[depth].cut = 0;

                    num_nodes++;


                    if (num_nodes % 1000000 == 0) {
                        cerr << "Nodes(" << pile[0].dim << "):"
                                << num_nodes << " : " << num_zero_nodes
                                << " : " << (num_nodes - num_zero_nodes) << '\n';
                    }

                }
            }


        } // end while(done)



        return lphi_out;
    }

    void construct_dfopt_tree(vector<vector<double> > &data,
            map_tree &map_region_tree, opt_region_hash<uint32_t> &map_regions) {

        MT_random rand_gen;

        // need to treat the empty tree/data case

        int64_t num_nodes = 0;
        int64_t num_zero_nodes = 0;

        int N = data.size();

        gamma_table gt(N);

        region_allocator<map_tree_node> *map_ra = map_region_tree.get_ra();


        current_region start_region(num_children);
        opt_region start_working(num_children);

        vector<df_working_unit_t> curr_good_regions;
        curr_good_regions.push_back(df_working_unit_t(data, start_region,
                start_working));

        vector<uint32_t> curr_map_nodes;
        curr_map_nodes.push_back(map_region_tree.get_full_tree());

        
        int top_count_lim = 1;

        if (top_count_ratio < 1) {
            top_count_lim = (int) floor(N * top_count_ratio);
            if (top_count_lim < 1) {
                top_count_lim = 1;
            }
        }else{
            top_count_lim = (int)top_count_ratio;
        }
        
        //cerr << "top_count_lim: " << top_count_lim << '\n';
        
        int map_depth = 0;


        // Now create the MAP tree!
        // Do a breadth first search on the MAP tree
        bool done = false;

        while (!done) {

            // if no good regions we are done!
            if (curr_good_regions.size() == 0) {
                done = true;
                continue;
            }

            //cerr << "======MAP_DEPTH: " << map_depth << '\n';
            //cerr << "Good Size: " << curr_good_regions.size() << '\n';

            // if too deep we stop all the regions
            if (map_depth >= max_depth) {
                done = true;
                
                //cerr << "STOP ALL REGIONS!" << '\n';

                // stop all the regions...
                vector<uint32_t>::iterator mit = curr_map_nodes.begin();
                for (vector<df_working_unit_t>::iterator wit = curr_good_regions.begin();
                        wit != curr_good_regions.end(); wit++) {

                    store_region(*mit, map_regions, map_ra, wit->working_reg,
                            wit->data.size(), map_depth);
                    mit++;
                }

                continue;

            }
            
            


            // for each good region
            vector<df_working_unit_t> next_good_regions;
            vector<uint32_t> next_map_nodes;

            vector<uint32_t>::iterator mit = curr_map_nodes.begin();
            int count = 1;
            for (vector<df_working_unit_t>::iterator wit = curr_good_regions.begin();
                    wit != curr_good_regions.end(); wit++) {


                
                // if region has few points we stop
                if (map_depth >= top_max_depth ||
                        (int) wit->data.size() <= top_count_lim ||
                        wit->working_reg.full()) {
                    // STOP
                    
                    //cerr << "STOP region!" << '\n';

                    store_region(*mit, map_regions, map_ra, wit->working_reg,
                            wit->data.size(), map_depth);

                    mit++; //grr
                    continue;

                }

                // Compute the lphi of the children of this good region
                vector<df_working_unit_t> children_wu0(num_children);
                vector<df_working_unit_t> children_wu1(num_children);

                vector<double> children_lphi0(num_children);
                vector<double> children_lphi1(num_children);

                for(int d = 0;d<num_children;d++){

                    bool is_diff = cut_region2(wit->data, children_wu0[d].data, children_wu1[d].data,
                            d, wit->curr_reg.get_lim(d));

                    //cerr << "d: " << d << " is_diff: " << is_diff << '\n';

                    children_wu0[d].curr_reg = wit->curr_reg;
                    children_wu0[d].working_reg = wit->working_reg;
                    children_wu0[d].curr_reg.cut(d, 0);
                    children_wu0[d].working_reg.cut(d, 0);

                    children_wu1[d].curr_reg = wit->curr_reg;
                    children_wu1[d].working_reg = wit->working_reg;
                    children_wu1[d].curr_reg.cut(d, 1);
                    children_wu1[d].working_reg.cut(d, 1);

                    bool count_ok_0 = true;
                    bool count_ok_1 = true;
                    if(count_ratio < 0){
                        count_ok_0 = ((int) children_wu0[d].data.size()
                            <= (int) max((int) floor(children_wu0[d].data.size() * count_ratio), 1));
                        count_ok_1 = ((int) children_wu1[d].data.size()
                            <= (int) max((int) floor(children_wu1[d].data.size() * count_ratio), 1));
                    }else{
                        count_ok_0 = (int) children_wu0[d].data.size() <= count_ratio;
                        count_ok_1 = (int) children_wu1[d].data.size() <= count_ratio;
                    }
                    
                    
                    if (!is_diff || count_ok_0) {
                        children_lphi0[d] = children_wu0[d].data.size() * (map_depth + 1) * c::l2;
                    } else {
                        children_lphi0[d] = compute_dflphi(children_wu0[d], (map_depth + 1), gt,
                                num_nodes, num_zero_nodes);
                    }

                    if (!is_diff || count_ok_1) {
                        children_lphi1[d] = children_wu1[d].data.size() * (map_depth + 1) * c::l2;
                    } else {
                        children_lphi1[d] = compute_dflphi(children_wu1[d], (map_depth + 1), gt,
                                num_nodes, num_zero_nodes);
                    }

                }

                //cerr << "Nodes{g:" << count << "}:"
                //            << num_nodes << " : " << num_zero_nodes
                //            << " : " << (num_nodes - num_zero_nodes) << '\n';

                // Compute lphi of current node from children
                vector<double> lphi_list;
                lphi_list.reserve(num_children + 1);

                int count_temp = wit->data.size();

                double max_val = (count_temp * map_depth * c::l2) - c::l2;
                lphi_list.push_back(max_val);
                
                //cerr << "max_val: " << max_val << '\n';

                double ld = -log(num_children) - c::lpi - c::l2;

                for (int i = 0; i < num_children; i++) {
                    // check for null
                    if (isinf(children_lphi0[i]) || isinf(children_lphi1[i])) {
                        cerr << "mapNULL child!!! " << i << ',' << map_depth << '\n';
                        exit(2);
                    }

                    int child_1_count = children_wu0[i].data.size();
                    int child_2_count = children_wu1[i].data.size();

                    if (child_1_count < 0 || child_2_count < 0) {
                        cerr << "map neg count!!! " << i << ',' << map_depth << '\n';
                        exit(2);
                    }

                    double val = ld;
                    val += children_lphi0[i] + children_lphi1[i];
                    val += gt.compute_lD2(count_temp, child_1_count, child_2_count);

                    //cerr << "val: " << val << '\n';
                    
                    lphi_list.push_back(val);

                    if (val > max_val) {
                        max_val = val;
                    }

                }

                double lphi = max_val;
                double sum = 0;
                for (int i = 0; i < (num_children + 1); i++) {
                    sum += exp(lphi_list[i] - max_val);
                }

                if (sum > 0) lphi += log(sum);


                //cerr << "lphi: " << lphi << '\n';


                // compute the posterior pho and decide if we stop
                double post_rho = -c::l2;
                post_rho += map_depth * c::l2 * wit->data.size(); // phi_0
                post_rho -= lphi;

                //cerr << "post_rho: " << post_rho << '\n';

                // if not stop we split and add the two regions to the good region list
                if (post_rho>-c::l2) {
                    // STOP

                    store_region(*mit, map_regions, map_ra, wit->working_reg,
                            wit->data.size(), map_depth);


                } else {
                    // Choose MAP dimension

                    int map_dim = 0;

                    double max_post_prob = -c::inf;

                    //cerr << "curr_p_count: " <<  wit->data.size() << '\n';
                    
                    for (int i = 0; i < num_children; i++) {
                        double post_prob = gt.compute_lD2(wit->data.size(),children_wu0[i].data.size()
                                , children_wu1[i].data.size());
                        post_prob += children_lphi0[i];
                        post_prob += children_lphi1[i];


                        //cerr << "post_prob[" << children_wu0[i].data.size() << "|" << children_wu1[i].data.size() << "](0) = " << post_prob << '\n';


                        if (post_prob > max_post_prob) {
                            map_dim = i;
                            max_post_prob = post_prob;
                        }
                    }


                    //cerr << "MAP dimension = " << map_dim << '\n';


                    // Split


                    next_good_regions.push_back(children_wu0[map_dim]);

                    next_good_regions.push_back(children_wu1[map_dim]);


                    pair<uint32_t, map_tree_node*> new_map_node0 = map_ra->create_node();
                    pair<uint32_t, map_tree_node*> new_map_node1 = map_ra->create_node();

                    (*map_ra)[(*mit)]->set_dim(map_dim);
                    (*map_ra)[(*mit)]->set_child(0, new_map_node0.first);
                    (*map_ra)[(*mit)]->set_child(1, new_map_node1.first);


                    next_map_nodes.push_back(new_map_node0.first);
                    next_map_nodes.push_back(new_map_node1.first);


                }

                mit++;
                count++;
            }

            // move next into current
            map_depth++;
            curr_good_regions = next_good_regions;
            curr_map_nodes = next_map_nodes;

        }


    }


};





#endif	/* DFOPT_H */

