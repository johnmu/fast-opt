/*
 *  llcopt_tree.h
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



#ifndef LLCOPT_TREE_H
#define	LLCOPT_TREE_H



#include "stl.h"
#include "general_utils.h"
#include "opt_utils.h"
#include "gamma_table.h"
#include "opt_tree.h"




struct llc_tree_node_sparse{
private:
    int count[2]; // number of points in this region for each type
    
    float lP;   // conditional marginal likelihood of not coupling
    float lphi; // conditional marginal for the OPT (combined sample) 

public:

    llc_tree_node_sparse(){
        count[0] = -1;
        count[1] = -1;
        lP = -c::inf;
        lphi = -c::inf;
    }
    
    llc_tree_node_sparse(const llc_tree_node_sparse& a){
        count[0] = a.count[0];
        count[1] = a.count[1];
        
        lP = a.lP;
        lphi = a.lphi;
    }

    void set_uniform(int depth){
        lP = abs(count[0]+count[1]) * depth * c::l2;
        lphi = lP;
    }

    float get_lP(){
        return lP;
    }
    
    float get_lphi(){
        return lphi;
    }
    
    void set_lP(float lP){
        this->lP = lP;
    }
    
    void set_lphi(float lphi){
        this->lphi = lphi;
    }

    void set_lP(double lP){
        this->lP = lP;
    }

    void set_count(int count[2]){
        this->count[0] = count[0];
        this->count[1] = count[1];
    }
    
    void get_count(int count_out[2]){
        count_out[0] = this->count[0];
        count_out[1] = this->count[1];
    }
    
    int get_count(){
        return count[0]+count[1];
    }
    
};


struct llc_node_t{
    llc_tree_node_sparse data[2];
    
    llc_node_t(){
        data[0] = llc_tree_node_sparse();
        data[1] = llc_tree_node_sparse();
    }
    
    llc_node_t(const llc_node_t &a){
        data[0] = a.data[0];
        data[1] = a.data[1];
    }
};

struct cdfpile_t{
    vector<uint32_t> data[2];
    vector<llc_node_t> nodes;
    int dim;
    int cut;
    
    cdfpile_t(const cdfpile_t &a){
        data[0] = a.data[0];
        data[1] = a.data[1];
        nodes = a.nodes;
        dim = a.dim;
        cut = a.cut;
    }
    
    cdfpile_t(){
        data[0] = vector<uint32_t>();
        data[1] = vector<uint32_t>();
        dim = 0;
        cut = 0;
    }
    
    cdfpile_t(int num_children){
        data[0] = vector<uint32_t>();
        data[1] = vector<uint32_t>();
        nodes.resize(num_children);
        dim = 0;
        cut = 0;
    }
};

struct llc_working_unit_t{
    vector<uint32_t> data[2];
    current_region curr_reg;
    opt_region working_reg;

    llc_working_unit_t(int num_children){
        data[0] = vector<uint32_t>();
        data[1] = vector<uint32_t>();
        curr_reg.init(num_children);
        working_reg.init(num_children);
    }
    
    llc_working_unit_t(const llc_working_unit_t &a){
        data[0] = a.data[0];
        data[1] = a.data[1];;
        curr_reg = a.curr_reg;
        working_reg = a.working_reg;
    }

    llc_working_unit_t(vector<uint32_t > data[2], current_region &curr_reg, opt_region &working_reg){

        this->data[0] = data[0];
        this->data[1] = data[1];
        this->curr_reg = curr_reg;
        this->working_reg = working_reg;
    }
};


struct llc_pile_t{
    vector<llc_working_unit_t> good_regions;
    vector<uint32_t> map_nodes;
};


struct llc_mouse_params_t{
    int start_dimension;
    int num_children;
    int count_lim;
    int top_depth;  // the depth of the parent
    int top_count_lim;
    int max_depth;
    int top_max_depth;
    llc_working_unit_t* wup;
    
    vector<vector<double> > *all_data;
    
    llc_node_t *output;
};



class llcopt_tree{
private:
    region_allocator<llc_tree_node_sparse> ra;
    int num_children;

    double count_lim;
    double top_count_lim;
    int max_depth;
    int top_max_depth;
    
    gamma_table gt;

    void init(int num_children, double count_lim, double top_count_lim,
                int max_depth, int top_max_depth) {
        this->num_children = num_children;
        
        if(top_count_lim < 2) top_count_lim = 2;
        if(count_lim > top_count_lim) count_lim = top_count_lim;
        else if(count_lim < 2) count_lim = 2;
        
        this->count_lim = count_lim;
        this->top_count_lim = top_count_lim;
        this->max_depth = max_depth;
        this->top_max_depth = top_max_depth;
    }

public:

    llcopt_tree(int num_children, double count_lim,double top_count_lim,int max_depth, int top_max_depth) {
        init(num_children, count_lim, top_count_lim, max_depth,top_max_depth);
    }

    llcopt_tree(int num_children) {
        init(num_children, 2, 2 , 2,num_children*20);
    }


    void store_region(uint32_t map_node, opt_region_hash<uint32_t> &map_regions,
            region_allocator<map_tree_node> *map_ra,opt_region &working_reg, int count, int depth){

        map_regions.insert(working_reg, map_node);

        (*map_ra)[map_node]->set_area(-depth);
        (*map_ra)[map_node]->set_count(count);
    }


    // compute lPhi and lP
    void compute_lPs(llc_tree_node_sparse &self, vector<llc_node_t> &children,int depth) {

        int num_children = children.size();
        
        vector<double> lphi_list;
        lphi_list.reserve(num_children+1);
        
        // Base measure
        int count[2];
        self.get_count(count);
        int total_count = abs(count[0]+count[1]);
        double max_val = (total_count * depth * c::l2) - (c::l2);
        lphi_list.push_back (max_val);
        
        // The random constants
        // lambda, D([1/2,1/2]) and 1/2
        double ld = -log(num_children) - c::lpi - c::l2;

        for(int i = 0;i<num_children;i++){
            int child_1_count = abs(children[i].data[0].get_count());
            int child_2_count = abs(children[i].data[1].get_count());

            double val = ld;
            val += children[i].data[0].get_lphi();
            val += children[i].data[1].get_lphi();
            val += gt.compute_lD2(total_count,child_1_count,child_2_count);

            lphi_list.push_back(val);

            if(val > max_val){
                max_val = val;
            }
        }

        float lphi = max_val;
        double sum = 0;
        for(int i = 0;i<(num_children+1);i++){
            sum += exp(lphi_list[i] - max_val);
        }
        if(sum > 0) lphi += log(sum);
        
        self.set_lphi(lphi);
        
        // this is for the coupling case
        
        max_val = lphi - c::l2;
        lphi_list[0] = max_val;
        
        // The random constants
        // lambda, (D([1/2,1/2]))^2 and 1/2
        ld = -log(num_children) - c::lpi - c::lpi - c::l2;

        for(int i = 0;i<num_children;i++){
            int child_1_count[2];
            children[i].data[0].get_count(child_1_count);
            child_1_count[0] = abs(child_1_count[0]);
            child_1_count[1] = abs(child_1_count[1]);
            int child_2_count[2];
            children[i].data[1].get_count(child_2_count);
            child_2_count[0] = abs(child_2_count[0]);
            child_2_count[1] = abs(child_2_count[1]);

            double val = ld;
            val += children[i].data[0].get_lP();
            val += children[i].data[1].get_lP();
            val += gt.compute_lD2(count[0],child_1_count[0],child_2_count[0]);
            val += gt.compute_lD2(count[1],child_1_count[1],child_2_count[1]);

            lphi_list[i+1] = val;

            if(val > max_val){
                max_val = val;
            }
        }
        
        float lP = max_val;
        sum = 0;
        for(int i = 0;i<(num_children+1);i++){
            sum += exp(lphi_list[i] - max_val);
        }

        if(sum > 0) lP += log(sum);
        
        self.set_lP(lP);
    }


    void* small_opt_thread(void* params) {

        llc_mouse_params_t* p = (llc_mouse_params_t*)params;

        int start_dimension = p->start_dimension;
        int top_depth = p->top_depth;
        int num_children = p->num_children;
        int count_lim = p->count_lim;
        int max_depth = p->max_depth;
        int top_max_depth = p->top_max_depth;
        llc_working_unit_t* wup = p->wup;
        llc_node_t* output = p->output;
        
        vector<vector<double> > *all_data = p->all_data;

        vector<cdfpile_t> pile;
        pile.push_back(cdfpile_t(num_children));
        pile[0].dim = start_dimension;
        pile[0].cut = 0;
        pile[0].data[0] = wup->data[0];
        pile[0].data[1] = wup->data[1];

        current_region curr_reg = wup->curr_reg;
        opt_region working_reg = wup->working_reg;

        int depth = 0;

        bool done = false;
        
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
            
            if(pile[0].dim != start_dimension){
                *output = pile[0].nodes[start_dimension];
                pile.pop_back();
                continue;
            }

            if (!isinf(pile[depth].nodes[curr_dim].data[curr_cut].get_lP())) {

                //cerr << "Increment cut: " << pile[depth].cut << " dim: "<< pile[depth].dim <<"\n";

                if (pile[depth].cut < c::cuts - 1) {
                    pile[depth].cut++;
                } else if (pile[depth].dim < num_children - 1) {
                    pile[depth].dim++;
                    pile[depth].cut = 0;
                } else {
                    
                    // reached end of node!! back up

                    //cerr << "FULL depth:" << depth << " top_depth:" << top_depth << "\n";
                  
                    llc_tree_node_sparse self;
                    int prev_count[2] = {pile[depth].data[0].size(),pile[depth].data[1].size()};
                    self.set_count(prev_count);
                    
                    compute_lPs(self,pile[depth].nodes,depth + top_depth);

                    if (depth > 0) {
                        int prev_dim = pile[depth - 1].dim;
                        int prev_cut = pile[depth - 1].cut;
                        pile[depth - 1].nodes[prev_dim].data[prev_cut] = self;
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

                vector<uint32_t> new_data[2];

                for (int k = 0; k < 2; k++) {
                    new_data[k] = vector<uint32_t>();
                    if (pile[depth].data[k].size() > 0) {
                        cut_region_one(*all_data, pile[depth].data[k],
                                new_data[k], curr_dim, curr_cut, curr_reg.get_lim(curr_dim));
                    } 
                }

                //cerr << "curr_dim: " << curr_dim << " curr_cut: " << curr_cut << '\n';
                //cerr << "curr_reg.get_lim(curr_dim): " << curr_reg.get_lim(curr_dim) << '\n';
                //cerr << "curr_reg.get_resolution(curr_dim): " << curr_reg.get_resolution(curr_dim) << '\n';

                //cerr << "curr_count: " << curr_count << '\n';

                int curr_count[2];
                curr_count[0] = new_data[0].size();
                curr_count[1] = new_data[1].size();
                pile[depth].nodes[curr_dim].data[curr_cut].set_count(curr_count);


                if ((abs(curr_count[0]) + abs(curr_count[1])) <= count_lim
                        || (curr_count[0] <=0 && curr_count[1] <=0)
                        || (depth +top_depth) >= top_max_depth
                        || (depth) >= max_depth
                        || working_reg.full()) {

                    //cerr << "leaf\n";


                    //cerr << "lphi curr_count: " << curr_count << '\n';
                    //cerr << "depth: " << depth << '\n';
                    //cerr << "top_depth: " << top_depth << '\n';
                    //cerr << "leaf lphi: " << curr_count * (depth + top_depth + 1) * c::l2 << '\n';
                    
                    pile[depth].nodes[curr_dim].data[curr_cut].set_uniform(depth + top_depth + 1);
                }else{
                    //cerr << "GO down\n";

                    curr_dim = pile[depth].dim;
                    curr_cut = pile[depth].cut;

                    //cerr << "curr_dim2: " << curr_dim << " curr_cut2: " << curr_cut << '\n';

                    // go down
                    depth++;

                    working_reg.cut(curr_dim, curr_cut);
                    curr_reg.cut(curr_dim, curr_cut);

                    pile.push_back(cdfpile_t(num_children));
                    pile[depth].dim = 0;
                    pile[depth].cut = 0;
                    pile[depth].data[0] = new_data[0];
                    pile[depth].data[1] = new_data[1];
                }
            }


        } // end while(done)



        return NULL;
    }


    void do_small_opt(vector<vector<double> > *all_data, llc_working_unit_t &w,
            vector<llc_node_t> &nodes, MT_random &rand_gen,int start_depth){


        llc_mouse_params_t* params = new llc_mouse_params_t[num_children];

        for (int d = 0; d < num_children; d++) {
            params[d].count_lim = count_lim;
            params[d].top_count_lim = top_count_lim;
            params[d].max_depth = max_depth;
            params[d].top_max_depth = top_max_depth;
            params[d].num_children = num_children;
            params[d].start_dimension = d;
            params[d].top_depth = start_depth;
            params[d].wup = &w;
            params[d].all_data = all_data;
            params[d].output = &nodes[d];
        }

        for (int d = 0; d < num_children; d++) {
            //cerr << "d=" << d << '\n';
            small_opt_thread((void*) &(params[d]));
            
        }
       
        delete [] params;
    }

int get_map_dim(llc_working_unit_t &w,vector<llc_node_t> &nodes,int map_depth) {
        // Choose MAP dimension

        int map_dim = 0;
        double max_post_prob = -c::inf;

        int count0[2];
        int count1[2];
        
        int curr_count[2];
        curr_count[0] = w.data[0].size();
        curr_count[1] = w.data[1].size();

        for (int i = 0; i < num_children; i++) {
            nodes[i].data[0].get_count(count0);
            nodes[i].data[1].get_count(count1);

            double post_prob = 0;
            post_prob += gt.compute_lD2(curr_count[0], count0[0], count1[0]);
            post_prob += gt.compute_lD2(curr_count[1], count0[1], count1[1]);
            post_prob += nodes[i].data[0].get_lP();
            post_prob += nodes[i].data[1].get_lP();

            if (post_prob > max_post_prob) {
                map_dim = i;
                max_post_prob = post_prob;
            }
        }

        return map_dim;
    }

    void construct_llcopt_tree(vector<vector<double> > *all_data, uint32_t first_len,
            map_tree &map_region_tree, opt_region_hash<uint32_t> &map_regions) {

        MT_random rand_gen;

        // need to treat the empty tree/data case

        uint32_t N_all = all_data->size();
        uint32_t N[2] = {first_len,N_all-first_len};
        
        vector<uint32_t> data[2];
        uint32_t idx = 0;
        for(int k = 0;k<2;k++){
            data[k].resize(N[k]);
            for (uint32_t i = 0; i < N[k]; i++) {
                data[k][i] = idx;
                idx++;
            }
        }
        
        cerr << "Full tree stopping at " << count_lim << " points\n";
        cerr << "Each look-ahead stopping at " << max_depth << " levels\n";
        cerr << "Full tree stopping at " << top_max_depth << " levels\n";

        gt.init(N_all);

        region_allocator<map_tree_node> *map_ra = map_region_tree.get_ra();

        current_region start_region(num_children);
        opt_region     start_working(num_children);

        vector<llc_pile_t> llpile;
        llpile.push_back(llc_pile_t());

        llpile[0].good_regions.push_back(llc_working_unit_t(data,start_region,start_working));
        llpile[0].map_nodes.push_back(map_region_tree.get_full_tree());

        int map_depth = 0;
        double total_area = 0.0;
        int prev_area_level = 0;

        // Now create the MAP tree!
        // Do a breadth first search on the MAP tree
        bool done = false;
        int count = 1;
        while(!done){

            // if no good regions we are done!
            if(llpile.size() == 0){
                done = true;
                continue;
            }

            if(llpile[map_depth].good_regions.size() == 0){
                llpile.pop_back();
                map_depth--;
                continue;
            }

            vector<llc_working_unit_t>::iterator wu_it = (llpile[map_depth].good_regions.end() - 1);
            vector<uint32_t>::iterator map_node_it = (llpile[map_depth].map_nodes.end() - 1);

            // if too deep we stop all the regions
            int curr_count[2] = {wu_it->data[0].size(),wu_it->data[1].size()};
            
            if (map_depth >= top_max_depth || (curr_count[0]+curr_count[1]) <= count_lim
                    || wu_it->working_reg.full()) {

                //cerr << "mini-STOP count: " << ra[wu_it->node_idx]->count << '\n';

                total_area += exp(wu_it->working_reg.get_area() * c::l2);
                if (floor(total_area * 100) >= prev_area_level) {
                    cerr << "Depth(" << map_depth << "):Area(" << 100 * total_area << "%)\n";
                    prev_area_level = floor(total_area * 100) + 1;
                }
                store_region(*map_node_it, map_regions, map_ra, wu_it->working_reg,
                        curr_count[0]-curr_count[1], map_depth);

                llpile[map_depth].good_regions.pop_back();
                llpile[map_depth].map_nodes.pop_back();
                count++;

                continue;
            }
            
            // for each good region

            //cerr << "Good [" << map_depth << "](" << count << "/"
            //        << llpile[map_depth].good_regions.size()
            //        << "), Data size: " << curr_data_size << '\n';

            // do a smaller OPT for the good region to determine which
            // dimension to cut
            vector<llc_node_t> nodes;
            do_small_opt(all_data, *wu_it,nodes,rand_gen, map_depth);

            //cerr << "\nNodes:" << num_nodes
            //        << " : " << num_zero_nodes
            //        << " : " << (num_nodes - num_zero_nodes) << '\n';

            llc_tree_node_sparse self;
            self.set_count(curr_count);
            compute_lPs(self, nodes,map_depth);

            // compute the posterior pho and decide if we stop
            double post_rho = -c::l2;

            // base measure
            post_rho += self.get_lphi();
            post_rho -= self.get_lP();

            // if not stop we split and add the two regions to the good region list
            if (post_rho>-c::l2) {
                // STOP

                total_area += exp(wu_it->working_reg.get_area() * c::l2);
                if (floor(total_area * 100) >= prev_area_level) {
                    cerr << "Depth(" << map_depth << "):Area(" << 100 * total_area << "%)\n";
                    prev_area_level = floor(total_area * 100) + 1;
                }

                store_region(*map_node_it, map_regions, map_ra, wu_it->working_reg,
                        curr_count[0]-curr_count[1], map_depth);

                llpile[map_depth].good_regions.pop_back();
                llpile[map_depth].map_nodes.pop_back();
                count++;

            } else {
                // Choose MAP dimension
                int map_dim = get_map_dim(*wu_it, nodes, map_depth);

                //cerr << "children: " << map_child_idx_0 << "," << map_child_idx_1 << '\n';

                // Split
                vector<uint32_t > new_data_0[2];
                vector<uint32_t > new_data_1[2];

                // too much data copying going on here
                for (int k = 0; k < 2; k++) {
                    cut_region2_one(*all_data, wu_it->data[k], new_data_0[k], new_data_1[k],
                            map_dim, wu_it->curr_reg.get_lim(map_dim));
                }
                
                current_region next_curr_reg0 = wu_it->curr_reg;
                opt_region next_working_reg0 = wu_it->working_reg;
                next_curr_reg0.cut(map_dim, 0);
                next_working_reg0.cut(map_dim, 0);

                current_region next_curr_reg1 = wu_it->curr_reg;
                opt_region next_working_reg1 = wu_it->working_reg;
                next_curr_reg1.cut(map_dim, 1);
                next_working_reg1.cut(map_dim, 1);

                // after this point wu_it is invalid
                llpile.push_back(llc_pile_t());
                map_depth++;

                map_node_it = (llpile[map_depth - 1].map_nodes.end() - 1);

                //cerr << "ADD PILE " << '\n';

                llpile[map_depth].good_regions.push_back(llc_working_unit_t(new_data_0,
                        next_curr_reg0, next_working_reg0));

                llpile[map_depth].good_regions.push_back(llc_working_unit_t(new_data_1,
                        next_curr_reg1, next_working_reg1));

                pair<uint32_t, map_tree_node*> new_map_node0 = map_ra->create_node();
                pair<uint32_t, map_tree_node*> new_map_node1 = map_ra->create_node();

                //cerr << "children: " << new_map_node0.first
                //        << "," << new_map_node1.first << '\n';

                //cerr << "*map_node_it: " << *map_node_it << '\n';

                (*map_ra)[(*map_node_it)]->set_dim(map_dim);
                (*map_ra)[(*map_node_it)]->set_child(0, new_map_node0.first);
                (*map_ra)[(*map_node_it)]->set_child(1, new_map_node1.first);

                llpile[map_depth].map_nodes.push_back(new_map_node0.first);
                llpile[map_depth].map_nodes.push_back(new_map_node1.first);

                llpile[map_depth - 1].good_regions.pop_back();
                llpile[map_depth - 1].map_nodes.pop_back();


            }

        }

    }


};

#endif	/* LLCOPT_TREE_H */
