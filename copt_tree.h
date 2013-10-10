/*
 *  copt_tree.h
 *  fast-opt
 *
 *  john.mu@ieee.org
 *
 */

/* The MIT License

   Copyright (c) 2013 John C. Mu.

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


#ifndef COPT_TREE_H
#define	COPT_TREE_H

//#define DEBUG
//#define DEBUG_MAP

#include "general_utils.h"
#include "gamma_table.h"
#include "opt_utils.h"
#include "map_tree.h"


class ctree_node
{
private:
    int count[2]; // number of points in this region for each type
    
    // this tree is designed so that it is a DAG
    // It is possible for branches to merge if the partitions are the same
    // pointer -> dimensions -> cuts (always 2 cuts for now)
    ctree_node*** children;

    // essentially the number of dimensions
    int num_children;
    
    float lP;   // conditional marginal likelihood of not coupling
    float lphi; // conditional marginal for the OPT (combined sample) 
    
    void init(){
        count[0] = -1;
        count[1] = -1;
        lP = -c::inf;
        lphi = -c::inf;
        children = NULL;
        num_children = 0;
    }

public:
    ctree_node(){
        init();
    }

    ctree_node(int num_children){
        init();
        this->num_children = num_children;

        children = new ctree_node**[num_children];

        for(int i = 0;i<num_children;i++){
            children[i] = new ctree_node*[c::cuts];
            for(int j = 0;j<c::cuts;j++){
                children[i][j] = NULL;
            }
        }
    }

    ~ctree_node(){
        if(children == NULL) return;

        for(int i = 0;i<num_children;i++){
            delete [] children[i];
        }
        delete [] children;
    }

    bool is_leaf(){
        return children == NULL;
    }

    void set_uniform(int depth){
        lP = (count[0]+count[1]) * depth * c::l2;
        lphi = lP;
    }

    // compute lPhi and lP
    void compute_lPs(int depth, gamma_table& gt){
        
        vector<double> lphi_list;
        lphi_list.reserve(num_children+1);
        
        // Base measure
        int total_count = (count[0]+count[1]);
        double max_val = (total_count * depth * c::l2) - (c::l2);
        lphi_list.push_back (max_val);
        
        // The random constants
        // lambda, D([1/2,1/2]) and 1/2
        double ld = -log(num_children) - c::lpi - c::l2;

        for(int i = 0;i<num_children;i++){
            // check for null
            if(children[i][0] == NULL || children[i][1] == NULL){
                cerr << "lphi NULL child!!! " << i << ',' << depth << '\n';
                exit(2);
            }

            int child_1_count = children[i][0]->get_count();
            int child_2_count = children[i][1]->get_count();

            if(child_1_count < 0 || child_2_count < 0){
                cerr << "lphi neg count!!! " << i << ',' << depth << '\n';
                exit(2);
            }

            double val = ld;
            val += children[i][0]->get_lphi();
            val += children[i][1]->get_lphi();
            val += gt.compute_lD2(total_count,child_1_count,child_2_count);

            lphi_list.push_back(val);

            if(val > max_val){
                max_val = val;
            }
        }

        lphi = max_val;
        double sum = 0;
        for(int i = 0;i<(num_children+1);i++){
            sum += exp(lphi_list[i] - max_val);
        }
        if(sum > 0) lphi += log(sum);
        
        // this is for the coupling case
        
        max_val = lphi - c::l2;
        lphi_list[0] = max_val;
        
        // The random constants
        // lambda, (D([1/2,1/2]))^2 and 1/2
        ld = -log(num_children) - c::lpi - c::lpi - c::l2;

        for(int i = 0;i<num_children;i++){
            // check for null
            if(children[i][0] == NULL || children[i][1] == NULL){
                cerr << "clphi NULL child!!! " << i << ',' << depth << '\n';
                exit(2);
            }

            int child_1_count[2];
            children[i][0]->get_count(child_1_count);
            int child_2_count[2];
            children[i][1]->get_count(child_2_count);

            if(child_1_count[0] < 0 || child_2_count[0] < 0
                    || child_1_count[1] < 0 || child_2_count[1] < 0){
                cerr << "cneg count!!! " << i << ',' << depth << '\n';
                exit(2);
            }

            double val = ld;
            val += children[i][0]->get_lP();
            val += children[i][1]->get_lP();
            val += gt.compute_lD2(count[0],child_1_count[0],child_2_count[0]);
            val += gt.compute_lD2(count[1],child_1_count[1],child_2_count[1]);

            lphi_list[i+1] = val;

            if(val > max_val){
                max_val = val;
            }

        }
        
        lP = max_val;
        sum = 0;
        for(int i = 0;i<(num_children+1);i++){
            sum += exp(lphi_list[i] - max_val);
        }

        if(sum > 0) lP += log(sum);
    }

    double get_lP(){
        return lP;
    }
    
    double get_lphi(){
        return lphi;
    }

    int get_num_children(){
        return num_children;
    }

    int get_cuts(){
        return c::cuts;
    }

    void set_lP(double lP){
        this->lP = lP;
    }
    
    void set_lPhi(double lphi){
        this->lphi = lphi;
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

    void set_child(int dim,int cut,ctree_node* node){
        if(children != NULL)children[dim][cut] = node;
    }

    ctree_node* get_child(int dim,int cut){
        if(children == NULL) return NULL;
        return children[dim][cut];
    }
};



class copt_tree
{
private:
    ctree_node* root;
    int num_children;

    int count_lim;
    int max_depth;

    void init(int num_children, int count_lim, int max_depth){
        this->num_children = num_children;
        root = new ctree_node(num_children);
        this->count_lim = count_lim;
        
        if(this->count_lim < 2){
            this->count_lim = 2;
            cerr << "Warning: Minimum count limit is 2, may be fixed in future version\n";
            // this is due to counting separately rather than together
        }
        
        this->max_depth = max_depth;
        region_cache.init_table(27);
    }

public:
    
    opt_region_hash<ctree_node*> region_cache;

    copt_tree(int num_children, int count_lim, int max_depth){
        init(num_children,count_lim,max_depth);
    }

    copt_tree(int num_children){
        init(num_children,5,1000);
    }

    ~copt_tree(){
        delete root;
    }

    void construct_full_tree(vector<vector<double> > all_data[2]){
        int64_t num_nodes = 0;
        int64_t num_zero_nodes = 0;

        int N[2] = {(int)all_data[0].size(),(int)all_data[1].size()};

        gamma_table gt(N[0]+N[1]);
        root->set_count(N);

        vector<cpile_t<ctree_node*,uint32_t > > pile;
        pile.push_back(cpile_t<ctree_node*,uint32_t >());
        
        for (int k = 0; k < 2; k++) {
            pile[0].data[k] = vector<uint32_t>(N[k],0);
            for (int i = 0; i < N[k]; i++) {
                pile[0].data[k][i] = i;
            }
        }
        
        pile[0].node = root;
        pile[0].dim  = 0;
        pile[0].cut  = 0;

        current_region curr_reg(num_children);
        opt_region working_reg(num_children);

        

        int depth = 0;
        
        bool done = false;
        while (!done){

            if(pile.size() == 0){
                done = true;
                continue;
            }

            int curr_dim = pile[depth].dim;
            int curr_cut = pile[depth].cut;
            ctree_node* curr_node = pile[depth].node;

            int curr_count[2];
            curr_node->get_count(curr_count);
            // work out what to count
            
            cerr << "== depth = " << depth << '\n';
            cerr << "curr_dim = " << curr_dim << '\n';
            cerr << "curr_cut = " << curr_cut << '\n';
            cerr << "curr_count[0] = " << curr_count[0] << '\n';
            cerr << "curr_count[1] = " << curr_count[1] << '\n';
            for(int k = 0;k<2;k++){
                cerr << "data[" << k << "]:\n";
                for(int i = 0;i<(int)pile[depth].data[k].size();i++){
                    for(int j = 0;j<num_children;j++){
                        cerr << all_data[k][pile[depth].data[k][i]][j] << ',';
                    }
                    cerr << "\n";
                }
            }

            bool back_up = false;
            // check if current node is leaf or at end
            if(curr_node->is_leaf()
                    || (curr_count[0]+curr_count[1]) <= count_lim
                    || depth >= max_depth
                    || working_reg.full()){
                // back up
                back_up = true;

                cerr << "LEAF!\n";
                
                // assume cuts are the same, so don't nee to search for lP0
                curr_node->set_uniform(depth);
                
            }else if(curr_node->get_child(curr_dim,curr_cut) != NULL){
                // move to next node
                
                if(pile[depth].cut < c::cuts - 1){
                    pile[depth].cut++;
                }else if(pile[depth].dim < num_children - 1){
                    pile[depth].dim++;
                    pile[depth].cut = 0;
                }else{
                    // reached end of node!! back up
                    cerr << "BACKUP\n";
                    
                    back_up = true;
                    curr_node->compute_lPs(depth,gt);
                    
                    cerr << "lphi = " << curr_node->get_lphi() << '\n';
                    cerr << "lP =   " << curr_node->get_lP() << '\n';
                }
            }

            if (back_up) {

                depth--;
                pile.pop_back();
                if(depth < 0) continue;
                
                curr_reg.uncut(pile[depth].dim,pile[depth].cut);
                working_reg.uncut(pile[depth].dim);

                continue;
            }

            curr_dim = pile[depth].dim;
            curr_cut = pile[depth].cut;

            // do the counting
            
            working_reg.cut(curr_dim,curr_cut);
            uint32_t working_hash = region_cache.hash(working_reg);

            pair<ctree_node*,bool> new_node = region_cache.find(working_reg,working_hash);

            if (!new_node.second) {

                pile.push_back(cpile_t<ctree_node*,uint32_t >());
                depth++;

                bool is_diff_sep[2] = {true, true};
                for (int k = 0; k < 2; k++) {
                    if(pile[depth-1].data[k].size() > 0){
                        is_diff_sep[k] = cut_region_one(all_data[k], pile[depth - 1].data[k], pile[depth].data[k],
                            curr_dim, curr_cut, curr_reg.get_lim(curr_dim));
                    }else{
                        is_diff_sep[k] = false;
                    }
                }
                
                bool is_diff = is_diff_sep[0]||is_diff_sep[1];

                curr_reg.cut(curr_dim, curr_cut);

                pile[depth].dim = 0;
                pile[depth].cut = 0;

                int curr_count[2];
                curr_count[0] = pile[depth].data[0].size();
                curr_count[1] = pile[depth].data[1].size();
                
                cerr << "##curr_count[0] = " << curr_count[0] << '\n';
                cerr << "##curr_count[1] = " << curr_count[1] << '\n';
               
                // must match the backup criteria
                // kind of un-elegant that we need this...
                if (!is_diff || (curr_count[0]+curr_count[1]) <= count_lim 
                        || depth >= max_depth || working_reg.full()){
                    new_node.first = new ctree_node();
                    num_zero_nodes++;
                }else {
                    new_node.first = new ctree_node(num_children);
                }

                new_node.first->set_count(curr_count);

                num_nodes++;

                curr_node->set_child(curr_dim, curr_cut, new_node.first);
                pile[depth].node = new_node.first;
                region_cache.insert(working_reg,new_node.first,working_hash);

                if (num_nodes % 1000000 == 0) {
                    cerr << "Nodes(" <<pile[0].dim << "):"
                            << num_nodes << " : " << num_zero_nodes
                            << " : " << (num_nodes-num_zero_nodes) <<'\n';
                }

            } else {
#ifdef DEBUG
                cerr << "found node: " << curr_dim << "," << curr_cut << '\n';
#endif
                curr_node->set_child(curr_dim, curr_cut, new_node.first);
                working_reg.uncut(curr_dim);
            }
        }

        cerr << "Nodes:" << num_nodes
                << ", Zero nodes:" << num_zero_nodes
                << ", Non-Zero nodes:" << (num_nodes-num_zero_nodes) <<'\n';


    }

    /////////////////////////
    // the tree and the regions
    // N is number of data points
    void construct_MAP_tree(map_tree &map_region_tree,opt_region_hash<uint32_t> &map_regions,int N){

        // Check if tree and region is empty
        
        // copy the root over

        // the second part is actually not used :/
        vector<pile_t<uint32_t,char> > map_pile;
        map_pile.push_back(pile_t<uint32_t,char>());
        map_pile[0].node = map_region_tree.get_full_tree();
        map_pile[0].dim  = -1;
        map_pile[0].cut  = -1;

        region_allocator<map_tree_node> *map_ra = map_region_tree.get_ra();

        // initialise state variables
        int depth = 0;
        gamma_table gt(N);

        vector<pile_t<ctree_node*,char > > pile;
        pile.push_back(pile_t<ctree_node*,char >());

        pile[0].data = vector<char>();
        pile[0].node = root;
        pile[0].dim  = -1;
        pile[0].cut  = -1;

        opt_region working_reg(num_children);

        bool done = false;
        while(!done) {

            if (pile.size() == 0) {
                done = true;
                continue;
            }
            // this is only kind of temporary

            int curr_dim = pile[depth].dim;
            int curr_cut = pile[depth].cut;
            ctree_node* curr_node = pile[depth].node;

            uint32_t curr_map_node = map_pile[depth].node;

            bool back_up = false;
            bool add_region = false;
            int map_dim = -1;

            // work out whether we stop at this node
            int curr_count[2];
            curr_node->get_count(curr_count);
            
            if (curr_node->is_leaf()|| (curr_count[0]+curr_count[1]) <= count_lim 
                    || depth >= max_depth || working_reg.full()) {
                // we are already at a uniform node

                // add to regions
                add_region = true;
                back_up = true;
                
            } else if(curr_dim == -1){
                double post_rho = -c::l2;
                
                // base measure

                post_rho += curr_node->get_lphi();
                post_rho -= curr_node->get_lP();
                
                if(post_rho>-c::l2){
                    // add to regions
                    add_region = true;
                    back_up = true;
                }else{
                    // choose a dimension

                    map_dim = 0;
                    ctree_node* child_0 = curr_node->get_child(0,0);
                    ctree_node* child_1 = curr_node->get_child(0,1);

                    int count0[2];
                    child_0->get_count(count0);
                    int count1[2];
                    child_1->get_count(count1);
                    int curr_count[2];
                    curr_node->get_count(curr_count);
                    
                    double max_post_prob = 0;
                    max_post_prob += gt.compute_lD2(curr_count[0],count0[0],count1[0]);
                    max_post_prob += gt.compute_lD2(curr_count[1],count0[1],count1[1]);
                    max_post_prob += child_0->get_lP();
                    max_post_prob += child_1->get_lP();

                    for(int i = 1;i<num_children; i++) {
                        child_0 = curr_node->get_child(i, 0);
                        child_1 = curr_node->get_child(i, 1);
                        
                        child_0->get_count(count0);
                        child_1->get_count(count1);

                        double post_prob = 0;
                        post_prob += gt.compute_lD2(curr_count[0],count0[0],count1[0]);
                        post_prob += gt.compute_lD2(curr_count[1],count0[1],count1[1]);
                        post_prob += child_0->get_lP();
                        post_prob += child_1->get_lP();

                        if(post_prob > max_post_prob){
                            map_dim = i;
                            max_post_prob = post_prob;
                        }
                    }

                    curr_dim = map_dim;
                    pile[depth].dim = map_dim;

                    curr_cut = 0;
                    pile[depth].cut = 0;

                }

            }else{
                // we have already chosen a dimension to cut...

                if(curr_cut < c::cuts - 1) {
                    pile[depth].cut++;
                    curr_cut++;
                }else{
                    back_up = true;
                }
            }


            if(add_region){

                // add region to the hash
                map_regions.insert(working_reg, curr_map_node);

                // compute the density
                (*map_ra)[curr_map_node]->set_area(-depth);
                int curr_count[2];
                curr_node->get_count(curr_count);
                (*map_ra)[curr_map_node]->set_count((curr_count[0]-curr_count[1]));
            }

            if(back_up){
                depth--;
                pile.pop_back();
                if (depth < 0) continue;
                working_reg.uncut(pile[depth].dim);
                continue;
            }

            // go down
            depth++;

            pair<uint32_t,map_tree_node*> new_map_node = map_ra->create_node();

            (*map_ra)[curr_map_node]->set_dim(curr_dim);
            (*map_ra)[curr_map_node]->set_child(curr_cut,new_map_node.first);

            map_pile.push_back(pile_t<uint32_t,char >());
            map_pile[depth].node = new_map_node.first;

            working_reg.cut(curr_dim, curr_cut);

            pile.push_back(pile_t<ctree_node*,char >());

            pile[depth].dim  = -1;
            pile[depth].cut  = -1;
            pile[depth].node = curr_node->get_child(curr_dim,curr_cut);
        }
    }


    int get_num_children(){
        return num_children;
    }

    ctree_node* get_full_tree(){
        return root;
    }

    double get_lP(){
        return root->get_lP();
    }
    
    double get_log_coupling_prob(){
        return root->get_lphi() - root->get_lP() - c::l2;
    }
};



#endif	/* COPT_TREE_H */

