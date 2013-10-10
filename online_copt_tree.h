/*
 *  online_copt_tree.h
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


#ifndef ONLINE_COPT_TREE_H
#define	ONLINE_COPT_TREE_H

//#define DEBUG
//#define DEBUG_MAP

#include "general_utils.h"
#include "gamma_table.h"
#include "opt_utils.h"
#include "map_tree.h"


class online_ctree_node
{
private:
    int count[2]; // number of points in this region for each type
    
    float lP;   // conditional marginal likelihood of not coupling
    float lphi; // conditional marginal for the OPT (combined sample) 
    
    // used to keep track of which nodes have been updated
    int sequence_id;              
    

public:
    
    // this stores the data points, but only used at the leaf nodes. 
    // it should be NULL for any non-leaf node
    vector<uint32_t>* data[2]; 
    
    online_ctree_node(){
        count[0] = -1;
        count[1] = -1;
        lP = -c::inf;
        lphi = -c::inf;
        sequence_id = 0;
        data[0] = NULL;
        data[1] = NULL;
    }
    
    online_ctree_node(const online_ctree_node& a){
        count[0] = a.count[0];
        count[1] = a.count[1];
        
        lP = a.lP;
        lphi = a.lphi;
        sequence_id = a.sequence_id;
        
        if(a.data[0] != NULL){
            data[0] = new vector<uint32_t>(a.data[0]->begin(),a.data[0]->end());
        }else{
            data[0] = NULL;
        }
        if(a.data[1] != NULL){
            data[1] = new vector<uint32_t>(a.data[1]->begin(),a.data[1]->end());
        }else{
            data[1] = NULL;
        }
                
    }

    ~online_ctree_node(){        
        if(data[0]!=NULL) delete data[0];
        if(data[1]!=NULL) delete data[1];
    }

    bool is_leaf(){
        return (data[0] != NULL);
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

    void set_sequence_id(int sequence_id){
        this->sequence_id = sequence_id;
    }
    
    void get_count(int count_out[2]){
        count_out[0] = this->count[0];
        count_out[1] = this->count[1];
    }
    
    int get_count(){
        return count[0]+count[1];
    }
    
    int get_sequence_id(){
        return sequence_id;
    }
};



class online_copt_tree
{
private:
    uint32_t root;
    int num_children;

    int count_lim;
    int max_depth;
    
    region_allocator<online_ctree_node> ra;
    opt_region_hash<uint32_t> region_cache;
    gamma_table gt;

    void init(int num_children, int count_lim, int max_depth,int max_N){
        this->num_children = num_children;
        root = c::ra_null_val;
        this->count_lim = count_lim;
        if(this->count_lim < 2){
            this->count_lim = 2;
            cerr << "Warning: Minimum count limit is 2, may be fixed in future version\n";
            // this is due to counting separately rather than together
        }
        this->max_depth = max_depth;
        region_cache.init_table(27);
        gt.init(max_N);
    }
    
    static uint32_t get_child(opt_region &working_reg, opt_region_hash<uint32_t> &region_cache,int dim, int cut){

        if(!working_reg.cut(dim,cut)){
            cerr << "CANNOT CUT: ";
            working_reg.print_region();
            cerr << '\n';
            exit(2);
        }

        pair<uint32_t,bool> out = region_cache.find(working_reg);
        working_reg.uncut(dim);
        
        if(out.second){
            return out.first;
        }else{
            return c::ra_null_val;
        }
    }
    
    static uint32_t get_parent(opt_region &working_reg, opt_region_hash<uint32_t> &region_cache,int dim){

        if(working_reg[dim].size() == 0){
            // no parent on this dimension
            return c::ra_null_val;
        }
        
        int prev_cut = working_reg.back(dim);
        working_reg.uncut(dim);

        pair<uint32_t,bool> out = region_cache.find(working_reg);
        working_reg.cut(dim,prev_cut);
        if(out.second){
            return out.first;
        }else{
            return c::ra_null_val;
        }
    }
    
    // compute lPhi and lP
    // modify so that we don't access the hash twice
    void compute_lPs(opt_region &working_reg, uint32_t curr_node, int depth, int seq_idx){
        
        vector<double> lphi_list;
        lphi_list.reserve(num_children+1);
        
        // Base measure
        int count[2];
        ra[curr_node]->get_count(count);
        int total_count = abs(count[0]+count[1]);
        double max_val = (total_count * depth * c::l2) - (c::l2);
        lphi_list.push_back (max_val);
        
        // The random constants
        // lambda, D([1/2,1/2]) and 1/2
        double ld = -log(num_children) - c::lpi - c::l2;

        for(int i = 0;i<num_children;i++){
            // check for null
            uint32_t child_id[2];
            child_id[0] = get_child(working_reg, region_cache,i,0);
            child_id[1] = get_child(working_reg, region_cache,i,1);
            
            ra[child_id[0]]->set_sequence_id(seq_idx);
            ra[child_id[1]]->set_sequence_id(seq_idx);
            
            int child_1_count = abs(ra[child_id[0]]->get_count());
            int child_2_count = abs(ra[child_id[1]]->get_count());

            if(child_1_count < 0 || child_2_count < 0){
                cerr << "lphi neg count!!! " << i << ',' << depth << '\n';
                exit(2);
            }

            double val = ld;
            val += ra[child_id[0]]->get_lphi();
            val += ra[child_id[1]]->get_lphi();
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
        
        ra[curr_node]->set_lphi(lphi);
        
        // this is for the coupling case
        
        max_val = lphi - c::l2;
        lphi_list[0] = max_val;
        
        // The random constants
        // lambda, (D([1/2,1/2]))^2 and 1/2
        ld = -log(num_children) - c::lpi - c::lpi - c::l2;

        for(int i = 0;i<num_children;i++){
            uint32_t child_id[2];
            child_id[0] = get_child(working_reg, region_cache,i,0);
            child_id[1] = get_child(working_reg, region_cache,i,1);

            int child_1_count[2];
            ra[child_id[0]]->get_count(child_1_count);
            child_1_count[0] = abs(child_1_count[0]);
            child_1_count[1] = abs(child_1_count[1]);
            int child_2_count[2];
            ra[child_id[1]]->get_count(child_2_count);
            child_2_count[0] = abs(child_2_count[0]);
            child_2_count[1] = abs(child_2_count[1]);

            if(child_1_count[0] < 0 || child_2_count[0] < 0
                    || child_1_count[1] < 0 || child_2_count[1] < 0){
                cerr << "cneg count!!! " << i << ',' << depth << '\n';
                exit(2);
            }

            double val = ld;
            val += ra[child_id[0]]->get_lP();
            val += ra[child_id[1]]->get_lP();
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
        
        ra[curr_node]->set_lP(lP);
    }

public:

    online_copt_tree(int num_children, int count_lim, int max_depth,int max_N){
        init(num_children,count_lim,max_depth, max_N);
    }

    online_copt_tree(int num_children,int max_N){
        init(num_children,5,1000, max_N);
    }

    ~online_copt_tree(){
    }

    // start and end are inclusive
    void construct_full_tree(vector<vector<double> > &all_data, uint32_t start[2], uint32_t end[2]){
        int64_t num_nodes = 0;

        for (int k = 0; k < 2; k++) {
            if(start[k] > end[k]){
                cerr << "Error: start after end: " << k << '\n';
            }
        }

        int N[2] = {(int)(end[0]-start[0]+1),(int)(end[1]-start[1]+1)};

        pair<uint32_t,online_ctree_node*> root_out = ra.create_node();
        root = root_out.first;
        root_out.second->set_count(N);

        vector<cpile_t<uint32_t,uint32_t > > pile(1);
        
        
        for (int k = 0; k < 2; k++) {
            pile[0].data[k].resize(N[k]);
            int idx = 0;
            for (uint32_t i = start[k]; i <= end[k]; i++) {
                pile[0].data[k][idx] = i;
                idx++;
            }
        }
        
        pile[0].node = root;
        pile[0].dim  = 0;
        pile[0].cut  = 0;

        current_region curr_reg(num_children);
        opt_region working_reg(num_children);
        
        region_cache.insert(working_reg,root);

        int depth = 0;
        
        bool done = false;
        while (!done){

            if(pile.size() == 0){
                done = true;
                continue;
            }

            int curr_dim = pile[depth].dim;
            int curr_cut = pile[depth].cut;
            online_ctree_node* curr_node = ra[pile[depth].node];

            int curr_count[2];
            curr_node->get_count(curr_count);
            // work out what to count

            bool back_up = false;
            // check if current node is leaf or at end
            if(curr_node->is_leaf()
                    || (curr_count[0]+curr_count[1]) <= count_lim
                    || depth >= max_depth
                    || working_reg.full()){
                // back up
                back_up = true;

                curr_node->set_uniform(depth);
                
                // save data in the leaf nodes
                curr_node->data[0] = new vector<uint32_t>(pile[depth].data[0]);
                curr_node->data[1] = new vector<uint32_t>(pile[depth].data[1]);

            } else if (curr_dim > num_children - 1) {
                // reached end of node!! back up
                back_up = true;

                compute_lPs(working_reg, pile[depth].node, depth,0);

            }

            if (back_up) {
                depth--;
                pile.pop_back();
                if (depth < 0) continue;
                
                curr_reg.uncut(pile[depth].dim,pile[depth].cut);
                working_reg.uncut(pile[depth].dim);

                if(pile[depth].cut < c::cuts - 1){
                    pile[depth].cut++;
                }else if(pile[depth].dim <= num_children - 1){
                    pile[depth].dim++;
                    pile[depth].cut = 0;
                }
                
                continue;
            }

            curr_dim = pile[depth].dim;
            curr_cut = pile[depth].cut;

            // do the counting

            working_reg.cut(curr_dim,curr_cut);
            uint32_t working_hash = region_cache.hash(working_reg);
            
            pair<uint32_t,bool> new_node = region_cache.find(working_reg,working_hash);

            if (!new_node.second) {

                pile.push_back(cpile_t<uint32_t,uint32_t>());
                depth++;

                bool is_diff_sep[2] = {true, true};
                for (int k = 0; k < 2; k++) {
                    if(pile[depth-1].data[k].size() > 0){
                        is_diff_sep[k] = cut_region_one(all_data, pile[depth - 1].data[k],
                                pile[depth].data[k],curr_dim, curr_cut, curr_reg.get_lim(curr_dim));
                    }else{
                        is_diff_sep[k] = false;
                    }
                }
                
                bool is_diff = is_diff_sep[0]||is_diff_sep[1];

                curr_reg.cut(curr_dim, curr_cut);

                pile[depth].dim = 0;
                pile[depth].cut = 0;

                int curr_count[2];
                if(is_diff){
                    curr_count[0] = pile[depth].data[0].size();
                    curr_count[1] = pile[depth].data[1].size();
                }else{
                    curr_count[0] = -pile[depth].data[0].size();
                    curr_count[1] = -pile[depth].data[1].size();
                }
                

                pair<uint32_t, online_ctree_node*> out = ra.create_node();
                new_node.first = out.first;

                out.second->set_count(curr_count);

                num_nodes++;
                
                pile[depth].node = new_node.first;
                region_cache.insert(working_reg,new_node.first,working_hash);

                if (num_nodes % 1000000 == 0) {
                    cerr << "Nodes(" <<pile[0].dim << "):"
                            << num_nodes <<'\n';
                }

            } else {

                working_reg.uncut(curr_dim);
                
                if(pile[depth].cut < c::cuts - 1){
                    pile[depth].cut++;
                }else if(pile[depth].dim <= num_children - 1){
                    pile[depth].dim++;
                    pile[depth].cut = 0;
                }
            }
        }

        cerr << "Nodes:" << num_nodes <<'\n';
    }

    // Remove the leaf nodes with children that are not needed anymore
    // [ add0,  del0,  add1,  del1]
    
    void prune_tree(vector<vector<double> > &all_data, uint32_t pts[4],int seq_idx){

        int N[2] = {0,0};
        ra[root]->get_count(N);

        vector<cpile_t<uint32_t,uint32_t> > pile;
        pile.push_back(cpile_t<uint32_t,uint32_t>());
        
        // here we use data to store which points are still active
        // 1 = active, 0 = inactive
        pile[0].data[0].push_back(1); //add
        pile[0].data[0].push_back(1); //del
        pile[0].data[1].push_back(1); //add
        pile[0].data[1].push_back(1); //del
        
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
            online_ctree_node* curr_node = ra[pile[depth].node];

            int curr_count[2];
            curr_node->get_count(curr_count);
            // work out what to count

            bool back_up = false;
            


            // check if current node is leaf or at end
            // do this check if the sequence id is reached

            // check leaf criteria
            if ((curr_count[0] + curr_count[1]) <= count_lim
                    || depth >= max_depth
                    || working_reg.full()) {
                // back up
                back_up = true;

                // check if this node is a leaf
                // If it is a leaf delete all children

                bool is_leaf = curr_node->is_leaf();
                if (is_leaf) {

                    // check if it has children
                    // probably need to check all children :/
                    // also return the data

                    vector<epile_t<uint32_t> > node_pile;
                    node_pile.push_back(epile_t<uint32_t>());
                    int edepth = 0;
                    node_pile[edepth].node = pile[depth].node;
                    node_pile[edepth].dim = 0;
                    node_pile[edepth].cut = 0;

                    while (node_pile.size() > 0) {
                        int edim = node_pile[edepth].dim;
                        int ecut = node_pile[edepth].cut;

                        if (edim < num_children) {

                            working_reg.cut(edim, ecut);

                            //uint32_t del_node = region_cache.erase(working_reg);
                            uint32_t reg_hash = region_cache.hash(working_reg);
                            pair<uint32_t, bool> out = region_cache.find(working_reg, reg_hash);
                            uint32_t del_node = out.first;
                            bool delete_node = true;
                            if (out.second) {
                                if (abs(ra[out.first]->get_sequence_id()) == seq_idx) {
                                    // don't delete

                                    delete_node = false;
                                } else {
                                    // check other criteria for not deleting
                                    // ie. any parent is not a leaf
                                    
                                    for(int k = 0;k<num_children;k++){
                                        uint32_t parent = get_parent(working_reg, region_cache,k);
                                        if(parent == c::ra_null_val){
                                            continue;
                                        }else{
                                            int parent_count[2];
                                            ra[parent]->get_count(parent_count);
                                            if (!((parent_count[0] + parent_count[1]) <= count_lim)) {
                                                // don't delete because parent is not a leaf
                                                delete_node = false;
                                            }
                                        }
                                    }

                                    if (delete_node) {

                                        region_cache.erase(working_reg, reg_hash);

                                    }
                                }
                            

                                if (ra[del_node]->is_leaf()) {
                                    if (delete_node) ra.delete_node(del_node);
                                } else {

                                    edepth++;
                                    node_pile.push_back(epile_t<uint32_t>());
                                    if (delete_node)node_pile[edepth].node = del_node;
                                    else node_pile[edepth].node = c::ra_null_val;
                                    node_pile[edepth].dim = 0;
                                    node_pile[edepth].cut = 0;
                                    node_pile[edepth].insert = node_pile[edepth - 1].insert;
                                    continue;
                                }
                            }
                            // uncut
                            working_reg.uncut(edim);

                        }


                        if (node_pile[edepth].dim == num_children) {
                            if (edepth != 0) {
                                if (node_pile[edepth].node != c::ra_null_val)ra.delete_node(node_pile[edepth].node);
                            }
                            node_pile.pop_back();
                            if(edepth>0)edepth--;
                            if(edepth>=0){
                                working_reg.uncut(node_pile[edepth].dim);
                            }
                        }
                        if (edepth >= 0) {
                            if (node_pile[edepth].cut < c::cuts - 1) {
                                node_pile[edepth].cut++;
                            } else if (node_pile[edepth].dim <= num_children - 1) {
                                node_pile[edepth].dim++;
                                node_pile[edepth].cut = 0;
                            }
                        }
                        
                    }
                }


            } else if (pile[depth].dim > num_children - 1) {
                // reached end of node!! back up

                // if it is not a leaf, delete the data points within
                if (curr_node->data[0] != NULL) {
                    delete curr_node->data[0];
                    delete curr_node->data[1];

                    curr_node->data[0] = NULL;
                    curr_node->data[1] = NULL;
                }
                
                back_up = true;

            }


            
            bool point_included = false;
            
            if (back_up) {
                depth--;
                pile.pop_back();
                if(depth < 0) continue;
                
                curr_reg.uncut(pile[depth].dim,pile[depth].cut);
                working_reg.uncut(pile[depth].dim);

            }else{
                // if sequence id is not reached check this
                // check if current node includes any of the points (in the current pile)
                // if doesn't include back up
                double lim = curr_reg.get_lim(curr_dim);

                for (int i = 0; i < 4 && !point_included; i++) {

                    if (pile[depth].data[i / 2][i % 2] == 1) {

                        if (curr_cut == 0) {
                            if (all_data[pts[i]][curr_dim] < lim) {
                                point_included = true;
                            }
                        } else if (curr_cut == 1) {
                            if (all_data[pts[i]][curr_dim] >= lim) {
                                point_included = true;
                            }
                        }
                    }
                }
            }

            if(back_up || (!point_included)){
                if(pile[depth].cut < c::cuts - 1){
                    pile[depth].cut++;
                }else if(pile[depth].dim <= num_children - 1){
                    pile[depth].dim++;
                    pile[depth].cut = 0;
                }
                
                continue;
            }
            
            curr_dim = pile[depth].dim;
            curr_cut = pile[depth].cut;

            // do the counting

            // determine if current node is a leaf
            // if it is a leaf we may need to delete the node or create new nodes
            // at most create one new level
            
            working_reg.cut(curr_dim,curr_cut);
            uint32_t working_hash = region_cache.hash(working_reg);

            pair<uint32_t,bool> new_node = region_cache.find(working_reg,working_hash);
            
            // In addition to checking if the new node exists we check the sequence number
            bool recount = false;
            if (!new_node.second){
                exit(1); // should always find
            }else{
                if(ra[new_node.first]->get_sequence_id()!=-seq_idx){
                    recount = true;
                }
            }
            
            if (recount){

                pile.push_back(cpile_t<uint32_t, uint32_t >());
                depth++;
                pile[depth].dim = 0;
                pile[depth].cut = 0;
                pile[depth].node = new_node.first;
                
                
                
                // determine which data points are active
                double lim = curr_reg.get_lim(curr_dim);
                for (int k = 0; k < 2; k++) {

                    pile[depth].data[k] = pile[depth-1].data[k];
                    for (int i = 0; i < 2; i++) {

                        if (pile[depth - 1].data[k][i] == 1) {
                            if (curr_cut == 0) {
                                if (!(all_data[pts[(2*k)+i]][curr_dim] < lim)) {
                                    pile[depth].data[k][i] = 0;

                                }
                            } else if (curr_cut == 1) {

                                if (!(all_data[pts[(2*k)+i]][curr_dim] >= lim)) {
                                    pile[depth].data[k][i] = 0;

                                }
                            }
                        }
                    }
                }

                online_ctree_node* new_ptr = ra[new_node.first];

                new_ptr->set_sequence_id(-seq_idx);

                // cut after the counting
                curr_reg.cut(curr_dim, curr_cut);

            } else {

                working_reg.uncut(curr_dim);
                
                if(pile[depth].cut < c::cuts - 1){
                    pile[depth].cut++;
                }else if(pile[depth].dim <= num_children - 1){
                    pile[depth].dim++;
                    pile[depth].cut = 0;
                }
            }
        }

    }
    
    
    
    // add and delete 1 point from each of the datasets
    // used to shift the window 
    // [ add0,  del0,  add1,  del1]
    
    void update_points(vector<vector<double> > &all_data, uint32_t pts[4],int seq_idx){

        int N[2] = {0,0};
        ra[root]->get_count(N);

        vector<cpile_t<uint32_t,uint32_t> > pile;
        pile.push_back(cpile_t<uint32_t,uint32_t>());
        
        // here we use data to store which points are still active
        // 1 = active, 0 = inactive
        pile[0].data[0].push_back(1); //add
        pile[0].data[0].push_back(1); //del
        pile[0].data[1].push_back(1); //add
        pile[0].data[1].push_back(1); //del
        
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
            online_ctree_node* curr_node = ra[pile[depth].node];

            int curr_count[2];
            curr_node->get_count(curr_count);
            // work out what to count
            bool back_up = false;
            


            // check if current node is leaf or at end
            // do this check if the sequence id is reached
            if (curr_node->get_sequence_id() == seq_idx || depth == 0) {

                // check leaf criteria
                if ((curr_count[0] + curr_count[1]) <= count_lim
                        || depth >= max_depth
                        || working_reg.full()) {
                    // back up
                    back_up = true;

                    // check if this node is a leaf
                    // If it is a leaf delete all children

                    bool is_leaf = curr_node->is_leaf();
                    if (!is_leaf) {
                        curr_node->data[0] = new vector<uint32_t>();
                        curr_node->data[1] = new vector<uint32_t>();
                        curr_node->data[0]->reserve(curr_count[0]+1);
                        curr_node->data[1]->reserve(curr_count[1]+1);
                        
                        // check if it has children
                        // probably need to check all children :/
                        // also return the data

                        vector<epile_t<uint32_t> > node_pile;
                        node_pile.push_back(epile_t<uint32_t>());
                        int edepth = 0;
                        node_pile[edepth].node = pile[depth].node;
                        node_pile[edepth].dim = 0;
                        node_pile[edepth].cut = 0;
                        node_pile[edepth].insert = true;

                        // first insert what is known
                        for (int k = 0; k < 2; k++) {
                            if (pile[depth].data[k][0] == 1) {
                                curr_node->data[k]->push_back(pts[2 * k]);
                            }
                        }
                        
                        while (node_pile.size() > 0) {
                            int edim = node_pile[edepth].dim;
                            int ecut = node_pile[edepth].cut;

                            if (edim < num_children) {

                                working_reg.print_region(cerr);
                                working_reg.cut(edim, ecut);

                                //uint32_t del_node = region_cache.erase(working_reg);
                                uint32_t reg_hash = region_cache.hash(working_reg);
                                pair<uint32_t,bool> out = region_cache.find(working_reg,reg_hash);
                                uint32_t del_node = out.first;
                                
                                if (out.second) {

                                    if (ra[del_node]->is_leaf()) {
                                        if (node_pile[edepth].insert) {
                                            // since this is a leaf, everything below must be a leaf

                                            for (int k = 0; k < 2; k++) {
 
                                                // if the sequence number is correct, simply insert it
                                                //if (false && ra[del_node]->get_sequence_id() == seq_idx) {
                                                if (false) {
                                                    // before deleting children need to store children data within
                                                    curr_node->data[k]->insert(curr_node->data[k]->end(),
                                                            ra[del_node]->data[k]->begin(), ra[del_node]->data[k]->end());

                                                } else {
                                                    // if it is incorrect, we need to re-count it, urgh

                                                    for (vector<uint32_t>::iterator it = ra[del_node]->data[k]->begin();
                                                            it != ra[del_node]->data[k]->end(); it++) {
                                                        // check each active pts to see which ones need to be skipped
                                                        bool skip = false;
                                                        for(int i = 0;i<2;i++){
                                                            if(pile[depth].data[k][i]==1 && *it == pts[2*k+i]){
                                                                skip = true;
                                                                break;
                                                            }
                                                        }
                                                        if(!skip){
                                                            curr_node->data[k]->push_back(*it);

                                                        }

                                                    }

                                                }
                                            }

                                        }
                                        //if(delete_node) ra.delete_node(del_node);
                                    } else {
                                        edepth++;
                                        node_pile.push_back(epile_t<uint32_t>());
                                        //if(delete_node)node_pile[edepth].node = del_node;
                                        //else node_pile[edepth].node = c::ra_null_val;
                                        node_pile[edepth].dim = 0;
                                        node_pile[edepth].cut = 0;
                                        node_pile[edepth].insert = node_pile[edepth - 1].insert;
                                        continue;
                                    }
                                }
                                // uncut
                                working_reg.uncut(edim);
      
                            }
                            
                            if (node_pile[edepth].dim == num_children) {

                                node_pile.pop_back();
                                edepth--;
                                if(edepth>=0)working_reg.uncut(node_pile[edepth].dim);
                            }

                            if (edepth >= 0) {
                                if (node_pile[edepth].cut < c::cuts - 1) {
                                    node_pile[edepth].cut++;
                                } else if (node_pile[edepth].dim <= num_children - 1) {
                                    node_pile[edepth].dim++;
                                    node_pile[edepth].cut = 0;
                                    node_pile[edepth].insert = false;
                                }
                            }

                        }
   
                        if(curr_node->data[0]->size() != curr_count[0]){
                            cerr << "COUNT_ERROR 0: " << curr_node->data[0]->size() << "," << curr_count[0] << '\n';
                        }
                        if(curr_node->data[1]->size() != curr_count[1]){
                            cerr << "COUNT_ERROR 1: " << curr_node->data[1]->size() << "," << curr_count[1] << '\n';
                        }
                        

                    }

                    // assume cuts are the same, so don't nee to search for lP0
                    curr_node->set_uniform(depth);

                } else if (pile[depth].dim > num_children - 1) {

                    // reached end of node!! back up

                    back_up = true;

                    compute_lPs(working_reg, pile[depth].node, depth,seq_idx);

                }
                
            }
            
            bool point_included = false;
            
            if (back_up) {

                curr_node->set_sequence_id(seq_idx);
                
                depth--;
                pile.pop_back();
                if(depth < 0) continue;
                
                curr_reg.uncut(pile[depth].dim,pile[depth].cut);
                working_reg.uncut(pile[depth].dim);


            }else{
                // if sequence id is not reached check this
                // check if current node includes any of the points (in the current pile)
                // if doesn't include back up
                double lim = curr_reg.get_lim(curr_dim);

                for (int i = 0; i < 4 && !point_included; i++) {

                    if (pile[depth].data[i / 2][i % 2] == 1) {

                        if (curr_cut == 0) {

                            if (all_data[pts[i]][curr_dim] < lim) {

                                point_included = true;
                            }
                        } else if (curr_cut == 1) {
                            if (all_data[pts[i]][curr_dim] >= lim) {
                                point_included = true;
                            }
                        }
                    }
                }

                if (!point_included) {

                    // test is the child exists
                    uint32_t child_id = get_child(working_reg, region_cache, curr_dim, curr_cut);
                    
                    if (child_id != c::ra_null_val) {
                        // if child exists, update sequence id

                        ra[child_id]->set_sequence_id(seq_idx);
                    } else {
                        // if child doesn't exist, treat as point included. 

                        point_included = true;
                    }
                }
            }

            if(back_up || (!point_included)){
                if(pile[depth].cut < c::cuts - 1){
                    pile[depth].cut++;
                }else if(pile[depth].dim <= num_children - 1){
                    pile[depth].dim++;
                    pile[depth].cut = 0;
                }
                
                continue;
            }
            
            curr_dim = pile[depth].dim;
            curr_cut = pile[depth].cut;

            // do the counting

            // determine if current node is a leaf
            // if it is a leaf we may need to delete the node or create new nodes
            // at most create one new level
            
            working_reg.cut(curr_dim,curr_cut);
            uint32_t working_hash = region_cache.hash(working_reg);

            pair<uint32_t,bool> new_node = region_cache.find(working_reg,working_hash);
            
            
            // In addition to checking if the new node exists we check the sequence number
            bool recount = false;
            if (!new_node.second){

                recount = true;
            }else{
                if(ra[new_node.first]->get_sequence_id()!=seq_idx){
                    recount = true;
                }
            }
            
            if (recount){
                pile.push_back(cpile_t<uint32_t, uint32_t >());
                depth++;
                pile[depth].dim = 0;
                pile[depth].cut = 0;
                
                
                
                // determine which data points are active
                double lim = curr_reg.get_lim(curr_dim);

                for (int k = 0; k < 2; k++) {
 
                    pile[depth].data[k] = pile[depth-1].data[k];
                    for (int i = 0; i < 2; i++) {

                        if (pile[depth - 1].data[k][i] == 1) {
                            if (curr_cut == 0) {

                                if (!(all_data[pts[(2*k)+i]][curr_dim] < lim)) {
                                    pile[depth].data[k][i] = 0;
 
                                }
                            } else if (curr_cut == 1) {
 
                                if (!(all_data[pts[(2*k)+i]][curr_dim] >= lim)) {
                                    pile[depth].data[k][i] = 0;

                                }
                            }
                        }
                    }
                }
                
                if (!new_node.second) {
                    // if the region does not exist
                    // we need to do counting and create a new node
                    vector<uint32_t>* new_data[2]; 
                    
                    bool is_diff_sep[2] = {true, true};
                    for (int k = 0; k < 2; k++) {
                        new_data[k] = new vector<uint32_t>();
                        if (curr_node->data[k]->size() > 0) {
                            is_diff_sep[k] = cut_region_one(all_data, *(curr_node->data[k]),
                                    *(new_data[k]), curr_dim, curr_cut, curr_reg.get_lim(curr_dim));
                        } else {
                            is_diff_sep[k] = false;
                        }
                    }
                    bool is_diff = is_diff_sep[0] || is_diff_sep[1];

                    int curr_count[2];
                    if (is_diff) {
                        curr_count[0] = new_data[0]->size();
                        curr_count[1] = new_data[1]->size();
                    } else {
                        curr_count[0] = -new_data[0]->size();
                        curr_count[1] = -new_data[1]->size();
                    }

                    // if we create a new node, remember to copy the leaf data into it. 

                    pair<uint32_t, online_ctree_node*> out = ra.create_node();
                    new_node.first = out.first;

                    out.second->set_count(curr_count);
                    out.second->data[0] = new_data[0];
                    out.second->data[1] = new_data[1];
                    out.second->set_sequence_id(seq_idx);
                    
                    region_cache.insert(working_reg,new_node.first,working_hash);

                } else {
                    
                    online_ctree_node* new_ptr = ra[new_node.first];
                    
                    // update counts
                    // modify counts here 
                    int count[2];
                    new_ptr->get_count(count);

                    for (int i = 0; i < 4; i++) {
                        int k = i/2;
                        int j = i % 2;
                        if (pile[depth].data[k][j] == 0) {
                            continue;
                        }

                        if (j == 0) {
                            if(count[k]>=0) count[k]++;
                            else count[k]--;
                        } else {
                            if(count[k]>0) count[k]--;
                            else if((count[k]<0)) count[k]++;
                            else {
                                cerr << "Error: Remove from zero?\n";
                                exit(1);
                            }
                        }
                    }

                    new_ptr->set_count(count);
                    new_ptr->set_sequence_id(seq_idx);

                    if (new_ptr->is_leaf()) {
                        // if the region is a leaf
                        // modify the data points in the leaf

                        for(int k = 0;k<2;k++){
                            if(pile[depth].data[k][1] == 1){
                                // del
                                vector<uint32_t>* data_ptr = new_ptr->data[k];
                                uint32_t del_pt = pts[(2*k)+1];
                                vector<uint32_t>::iterator del_it = data_ptr->end();
                                for (vector<uint32_t>::iterator it = data_ptr->begin();
                                        it != data_ptr->end();it++){
                                    if(*it == del_pt){
                                        del_it = it;
                                        break;
                                    }
                                }
                                if(del_it != data_ptr->end()){
                                    data_ptr->erase(del_it);
                                }else{
                                    // error!!!
                                    cerr << "Del error!\n";
                                }
                            }
                            
                            // add
                            if(pile[depth].data[k][0] == 1){
                                vector<uint32_t>* data_ptr = new_ptr->data[k];
                                data_ptr->push_back(pts[2*k]);
                            }
                        }
                        
                    } else {
                        // if the region is not a leaf
                        // do nothing?
                    }
                    
                    
                }
                
                // cut after the counting
                curr_reg.cut(curr_dim, curr_cut);
                
                pile[depth].node = new_node.first;
            } else {

                working_reg.uncut(curr_dim);
                
                if(pile[depth].cut < c::cuts - 1){
                    pile[depth].cut++;
                }else if(pile[depth].dim <= num_children - 1){
                    pile[depth].dim++;
                    pile[depth].cut = 0;
                }
            }
        }

    }
    
    /////////////////////////
    // the tree and the regions
    // N is number of data points
    void construct_MAP_tree(map_tree &map_region_tree,opt_region_hash<uint32_t> &map_regions,int N){

        // Check if tree and region is empty
        if(root == c::ra_null_val){
            cerr << "Error: online coupling OPT not constructed yet\n";
            exit(2);
        }

        // the second part is actually not used :/
        vector<pile_t<uint32_t,char> > map_pile;
        map_pile.push_back(pile_t<uint32_t,char>());
        map_pile[0].node = map_region_tree.get_full_tree();
        map_pile[0].dim  = -1;
        map_pile[0].cut  = -1;

        region_allocator<map_tree_node> *map_ra = map_region_tree.get_ra();

        // initialize state variables
        int depth = 0;

        vector<pile_t<uint32_t,char > > pile;
        pile.push_back(pile_t<uint32_t,char >());

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
            online_ctree_node* curr_node = ra[pile[depth].node];

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
                    uint32_t child_id[2];
                    child_id[0] = get_child(working_reg, region_cache, 0, 0);
                    child_id[1] = get_child(working_reg, region_cache, 0, 1);
                    
                    map_dim = 0;
                    online_ctree_node* child_0 = ra[child_id[0]];
                    online_ctree_node* child_1 = ra[child_id[1]];

                    int count0[2];
                    child_0->get_count(count0);
                    count0[0] = abs(count0[0]);
                    count0[1] = abs(count0[1]);
                    
                    int count1[2];
                    child_1->get_count(count1);
                    count1[0] = abs(count1[0]);
                    count1[1] = abs(count1[1]);
                    
                    int curr_count[2];
                    curr_node->get_count(curr_count);
                    curr_count[0] = abs(curr_count[0]);
                    curr_count[1] = abs(curr_count[1]);
                    
                    double max_post_prob = 0;
                    max_post_prob += gt.compute_lD2(curr_count[0],count0[0],count1[0]);
                    max_post_prob += gt.compute_lD2(curr_count[1],count0[1],count1[1]);
                    max_post_prob += child_0->get_lP();
                    max_post_prob += child_1->get_lP();

                    for(int i = 1;i<num_children; i++) {
                        child_id[0] = get_child(working_reg, region_cache, i, 0);
                        child_id[1] = get_child(working_reg, region_cache, i, 1);
                        
                        online_ctree_node* child_0 = ra[child_id[0]];
                        online_ctree_node* child_1 = ra[child_id[1]];
                        
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

            pile.push_back(pile_t<uint32_t,char >());

            pile[depth].dim  = -1;
            pile[depth].cut  = -1;
            pile[depth].node = get_child(working_reg, region_cache, curr_dim, curr_cut);
            
            working_reg.cut(curr_dim, curr_cut);
        }
    }


    int get_num_children(){
        return num_children;
    }

    online_ctree_node* get_full_tree(){
        return ra[root];
    }

    double get_lP(){
        return ra[root]->get_lP();
    }
    
    double get_lphi(){
        return ra[root]->get_lphi();
    }
    
    double get_log_coupling_prob(){
        online_ctree_node* curr_node = ra[root];
        return curr_node->get_lphi() - curr_node->get_lP() - c::l2;
    }
    
    opt_region_hash<uint32_t>* get_region_cache(){
        return &region_cache;
    }
    
    region_allocator<online_ctree_node>* get_ra(){
        return &ra;
    }
};



#endif	/* COPT_TREE_H */

