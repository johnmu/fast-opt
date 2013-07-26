/*
 *  opt_tree.h
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


#ifndef OPT_TREE_H
#define	OPT_TREE_H

//#define DEBUG
//#define DEBUG_MAP

#include "general_utils.h"
#include "gamma_table.h"
#include "opt_utils.h"
#include "map_tree.h"


class tree_node
{
private:
    int count; // number of points in this region, this is also the dimension
                // to cut it is MAP tree.
    
    // this tree is designed so that it is a DAG
    // It is possible for branches to merge if the partitions are the same
    // pointer -> dimensions -> cuts (always 2 cuts for now)
    tree_node*** children;

    // essentially the number of dimensions
    int num_children;
    
    double lphi;  // this is also the density if MAP tree

    void init(){
        count = -1;
        lphi = -c::inf;
    }

public:
    tree_node(){
        children = NULL;
        num_children = 0;

        init();
    }

    tree_node(int num_children){
        this->num_children = num_children;

        children = new tree_node**[num_children];

        for(int i = 0;i<num_children;i++){
            children[i] = new tree_node*[c::cuts];
            for(int j = 0;j<c::cuts;j++){
                children[i][j] = NULL;
            }
        }

        init();
    }

    ~tree_node(){

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
        lphi = count * depth * c::l2;
    }

    void compute_lphi(int depth, gamma_table& gt, double lphi0){
        
        vector<double> lphi_list;
        lphi_list.reserve(num_children+1);
        
        // Base measure
        double max_val = 0.0;
        if (lphi0 == -c::inf) {
            max_val = (count * depth * c::l2) - c::l2;
        }else{
            max_val = lphi0;
        }
        
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

#ifdef DEBUG
            cerr << "child_1_count(" << i <<"): " << child_1_count << '\n';
            cerr << "child_2_count(" << i <<"): " << child_2_count << '\n';
            
            cerr << "lphi0: " << children[i][0]->get_lphi() << '\n';
            cerr << "lphi1: " << children[i][1]->get_lphi() << '\n';
#endif
            
            if(child_1_count < 0 || child_2_count < 0){
                cerr << "neg count!!! " << i << ',' << depth << '\n';
                exit(2);
            }

            double val = ld;
            val += children[i][0]->get_lphi();
            val += children[i][1]->get_lphi();
            val += gt.compute_lD2(count,child_1_count,child_2_count);

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

    void set_lphi(double lphi){
        this->lphi = lphi;
    }

    void set_count(int count){
        this->count = count;
    }

    int get_count(){
        return count;
    }

    void set_child(int dim,int cut,tree_node* node){
        if(children != NULL)children[dim][cut] = node;
    }

    tree_node* get_child(int dim,int cut){
        if(children == NULL) return NULL;
        return children[dim][cut];
    }

};



class opt_tree
{
private:
    tree_node* root;
    int num_children;

    int count_lim;
    int max_depth;
    
    opt_tree* base_measure;

    void init(int num_children, int count_lim, int max_depth,opt_tree* base_measure){
        this->num_children = num_children;
        root = new tree_node(num_children);
        this->count_lim = count_lim;
        this->max_depth = max_depth;
        this->base_measure = base_measure;
    }

public:

    opt_tree(int num_children, int count_lim, int max_depth,opt_tree* base_measure){
        init(num_children,count_lim,max_depth,base_measure);
    }
    
    opt_tree(int num_children, int count_lim, int max_depth){
        init(num_children,count_lim,max_depth,NULL);
    }

    opt_tree(int num_children){
        init(num_children,5,1000,NULL);
    }

    ~opt_tree(){
        delete root;
    }

    void construct_full_tree(vector<vector<double> > &all_data){
        int64_t num_nodes = 0;
        int64_t num_zero_nodes = 0;

        uint32_t N = (uint32_t)all_data.size();

        gamma_table gt(N);
        root->set_count(N);

        vector<pile_t<tree_node*,uint32_t > > pile;
        pile.push_back(pile_t<tree_node*,uint32_t >());
        
        pile[0].data = vector<uint32_t>(N,0);
        
        for(uint32_t i = 0;i<N;i++){
            pile[0].data[i] = i;
        }
        
        pile[0].node = root;
        pile[0].dim  = 0;
        pile[0].cut  = 0;

        current_region curr_reg(num_children);
        opt_region working_reg(num_children);

        opt_region_hash<tree_node*> region_cache(27);

        int depth = 0;
        
        bool done = false;
        while (!done){

            if(pile.size() == 0){
                done = true;
                continue;
            }

            // print the pile
#ifdef DEBUG
            cerr << "===---- " << depth << " , " << pile.size() << " ---\n";
            cerr << num_nodes << " : " << num_zero_nodes << '\n';
            for(int i = 0;i<(int)pile.size();i++){
                cerr << "Data size: " << pile[i].data.size() << '\n';
                cerr << "Dim:       " << pile[i].dim << '\n';
                cerr << "Cut:       " << pile[i].cut << '\n';
                //if(pile[i].data.size()<20){
                //    print_data(pile[i].data);
                //}
            }
            cerr << "*******\n";
#endif
            // this is only kind of temporary

            int curr_dim = pile[depth].dim;
            int curr_cut = pile[depth].cut;
            tree_node* curr_node = pile[depth].node;

            // work out what to count

            bool back_up = false;

            // check if current node is leaf or at end
            if(curr_node->is_leaf()
                    || curr_node->get_count() <= count_lim
                    || depth >= max_depth
                    || working_reg.full()){
                // back up
                back_up = true;
                curr_node->set_uniform(depth);
#ifdef DEBUG
                cerr << "BACKUP ZERO\n";
#endif

            }else if(curr_node->get_child(curr_dim,curr_cut) != NULL){
                // move to next node
#ifdef DEBUG
                cerr << "CHILD FOUND\n";
                cerr << "pile[depth].cut: " << pile[depth].cut << '\n';
                cerr << "pile[depth].dim: " << pile[depth].dim << '\n';
#endif

                if(pile[depth].cut < c::cuts - 1){
                    pile[depth].cut++;
                }else if(pile[depth].dim < num_children - 1){
                    pile[depth].dim++;
                    pile[depth].cut = 0;
                }else{
                    // reached end of node!! back up
                    back_up = true;
                    
                    double lphi0 = -c::inf;
                    if(!(base_measure == NULL)){
                        lphi0 = base_measure->get_reg_lphi(working_reg);
                        if(lphi0 == -c::inf){
                            lphi0 = 0;
                        }
                    }
                    
                    curr_node->compute_lphi(depth,gt,lphi0);
#ifdef DEBUG
                    cerr << "BACKUP CHILD\n";
#endif
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
#ifdef DEBUG
            cerr << "CUT: " << curr_dim << "," << curr_cut << '\n';
#endif
            
            working_reg.cut(curr_dim,curr_cut);
            uint32_t working_hash = region_cache.hash(working_reg);
#ifdef DEBUG
            cerr << "working_hash: " << working_hash << '\n';
#endif
            pair<tree_node*,bool> new_node = region_cache.find(working_reg,working_hash);

            if (!new_node.second) {

#ifdef DEBUG
                cerr << "Add node: " << curr_dim << "," << curr_cut << '\n';
#endif
                pile.push_back(pile_t<tree_node*,uint32_t >());
                depth++;

                bool is_diff = cut_region_one(all_data,pile[depth - 1].data,pile[depth].data,
                        curr_dim, curr_cut, curr_reg.get_lim(curr_dim));

                curr_reg.cut(curr_dim, curr_cut);

                pile[depth].dim = 0;
                pile[depth].cut = 0;

                int curr_count = pile[depth].data.size();

                // must match the backup criteria
                // kind of un-elegant that we need this...
                if (!is_diff || curr_count <= count_lim 
                        || depth >= max_depth || working_reg.full()){
                    new_node.first = new tree_node();
                }else {
                    new_node.first = new tree_node(num_children);
                }

                new_node.first->set_count(curr_count);

                num_nodes++;
                if (!is_diff || new_node.first->get_count() <= count_lim)num_zero_nodes++;

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

        vector<pile_t<tree_node*,char > > pile;
        pile.push_back(pile_t<tree_node*,char >());

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
            tree_node* curr_node = pile[depth].node;

            uint32_t curr_map_node = map_pile[depth].node;

#ifdef DEBUG_MAP
            cerr << "===-----------------------" << '\n';
            cerr << "===--- DEPTH = " << depth << '\n';
            cerr << "working_reg:";
            working_reg.print_region_limits(cerr);
            cerr << '\n';
#endif

            bool back_up = false;
            bool add_region = false;
            int map_dim = -1;

            // work out whether we stop at this node
            if (curr_node->is_leaf()|| curr_node->get_count() <= count_lim 
                    || depth >= max_depth || working_reg.full()) {
                // we are already at a uniform node

                // add to regions
                add_region = true;
                back_up = true;
#ifdef DEBUG_MAP
                cerr << "BOTTOM " << '\n';
#endif
            } else if(curr_dim == -1){
                double post_rho = -c::l2;
                
                // base measure
                if(base_measure == NULL){
                    post_rho += depth * c::l2 * curr_node->get_count(); // phi_0
                }else{
                    double lphi0 = base_measure->get_reg_lphi(working_reg);
                    if(!lphi0 == -c::inf){
                        post_rho += lphi0;
                    }
                }
                
                post_rho -= curr_node->get_lphi();
#ifdef DEBUG_MAP
                cerr << "get_lphi = " << curr_node->get_lphi() << '\n';
                cerr << "post_rho = " << post_rho << '\n';
#endif
                if(post_rho>-c::l2){
                    // add to regions
                    add_region = true;
                    back_up = true;
                }else{
                    // choose a dimension

                    map_dim = 0;
                    tree_node* child_0 = curr_node->get_child(0,0);
                    tree_node* child_1 = curr_node->get_child(0,1);

                    double max_post_prob = gt.compute_lD2(curr_node->get_count()
                        ,child_0->get_count(),child_1->get_count());
                    max_post_prob += child_0->get_lphi();
                    max_post_prob += child_1->get_lphi();
#ifdef DEBUG_MAP
                    cerr << "post_prob["<< child_0->get_count() << "|" << child_1->get_count() <<"](0) = " << max_post_prob << '\n';
#endif
                    for(int i = 1;i<num_children; i++) {
                        child_0 = curr_node->get_child(i, 0);
                        child_1 = curr_node->get_child(i, 1);

                        double post_prob = gt.compute_lD2(curr_node->get_count()
                                , child_0->get_count(), child_1->get_count());
                        post_prob += child_0->get_lphi();
                        post_prob += child_1->get_lphi();
#ifdef DEBUG_MAP
                        cerr << "post_prob["<< child_0->get_count() << "|" << child_1->get_count() <<"]("<< i << ") = " << post_prob << '\n';
#endif
                        if(post_prob > max_post_prob){
                            map_dim = i;
                            max_post_prob = post_prob;
                        }
                    }
#ifdef DEBUG_MAP
                    cerr << "map_dim = " << map_dim << '\n';
#endif
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
#ifdef DEBUG_MAP
            cerr << "curr_dim = " << curr_dim << '\n';
            cerr << "curr_cut = " << curr_cut << '\n';
#endif

            if(add_region){
#ifdef DEBUG_MAP
                cerr << "ADD REGION" << '\n';
#endif
                // add region to the hash
                map_regions.insert(working_reg, curr_map_node);

                // compute the density

                (*map_ra)[curr_map_node]->set_area(-depth);
                (*map_ra)[curr_map_node]->set_count(curr_node->get_count());
            }

            if(back_up){
#ifdef DEBUG_MAP
                cerr << "BACKUP" << '\n';
#endif
                depth--;
                pile.pop_back();
                if (depth < 0) continue;
                working_reg.uncut(pile[depth].dim);

                continue;
            }
#ifdef DEBUG_MAP
            cerr << "GODOWN" << '\n';
#endif
            // go down
            depth++;

            pair<uint32_t,map_tree_node*> new_map_node = map_ra->create_node();

            (*map_ra)[curr_map_node]->set_dim(curr_dim);
            (*map_ra)[curr_map_node]->set_child(curr_cut,new_map_node.first);

            map_pile.push_back(pile_t<uint32_t,char >());
            map_pile[depth].node = new_map_node.first;

            working_reg.cut(curr_dim, curr_cut);

            pile.push_back(pile_t<tree_node*,char >());

            pile[depth].dim  = -1;
            pile[depth].cut  = -1;
            pile[depth].node = curr_node->get_child(curr_dim,curr_cut);
        }
    }

    // search through tree to get the lphi of a particular region
    // unless we reach a leaf.. then return -inf
    // can only be used after the tree has been completed
    double get_reg_lphi(opt_region &reg){
        // start from root and work way down one dimension at a time
        // it doesn't matter what order cause the multi-tree is already a DAG
        
        if(reg.num_children() != num_children){
            cerr << "Error: get_lphi: dimension mismatch\n";
            return -c::inf;
        }
        
        tree_node* curr_node = root;
        bool done = false;

        // check if we reached the end
        if (curr_node == NULL) {
            done = true;
        }

        int num_cuts = 0;
        bool reach_end_cut = false;
        bool reach_end_dim = false;
        for (int d = 0; !done && d < num_children; d++) {
            // for each dimension
            reach_end_cut = false;
            num_cuts = reg[d].size();
            for (int i = 0;!done && i<num_cuts;i++){
                int cut = (int)reg[d][i];
                curr_node = curr_node->get_child(0,cut);
                if(curr_node->is_leaf()){
                    done = true;
                }
                if(i == num_cuts-1){
                    reach_end_cut = true;
                }
            }
            if(d == num_children -1){
                reach_end_dim = true;
            }
        }
        
        if(reach_end_dim && reach_end_cut){
            return curr_node->get_lphi();
        }else{
            return -c::inf;
        }
        
    }

    int get_num_children(){
        return num_children;
    }

    tree_node* get_full_tree(){
        return root;
    }

    double get_lphi(){
        return root->get_lphi();
    }
};



#endif	/* OPT_TREE_H */

