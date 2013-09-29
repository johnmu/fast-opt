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
    uint32_t sequence_id;              
    
    void init(){
        count[0] = -1;
        count[1] = -1;
        lP = -c::inf;
        lphi = -c::inf;
        //children = NULL;
        //num_children = 0;
        sequence_id = 0;
        data[0] = NULL;
        data[1] = NULL;
    }

public:
    
    // this stores the data points, but only used at the leaf nodes. 
    // it should be NULL for any non-leaf node
    vector<uint32_t>* data[2]; 
    
    online_ctree_node(){
        init();
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
        cerr << "delete: " <<  data[0] << ',' << data[1] << '\n';
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
    
    static uint32_t get_child(opt_region working_reg, opt_region_hash<uint32_t> &region_cache,
            int dim, int cut){

        if(!working_reg.cut(dim,cut)){
            cerr << "CANNOT CUT: ";
            working_reg.print_region();
            cerr << '\n';
            exit(2);
        }

        pair<uint32_t,bool> out = region_cache.find(working_reg);

        if(out.second){
            return out.first;
        }else{
            return c::ra_null_val;
        }
    }
    
    // compute lPhi and lP
    // modify so that we don't access the hash twice
    void compute_lPs(opt_region &working_reg, uint32_t curr_node, int depth){
        
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
            cerr << "N["<< k<<"]"<<N[k]<<'\n';
            for (uint32_t i = start[k]; i <= end[k]; i++) {
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
            online_ctree_node* curr_node = ra[pile[depth].node];

            int curr_count[2];
            curr_node->get_count(curr_count);
            // work out what to count

            cerr << "== depth = " << depth << '\n';
            cerr << "curr_dim = " << curr_dim << '\n';
            cerr << "curr_cut = " << curr_cut << '\n';
            cerr << "curr_count[0] = " << curr_count[0] << '\n';
            cerr << "curr_count[1] = " << curr_count[1] << '\n';
            
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
                cerr << "Is leaf! save data! "<< curr_node->is_leaf() <<"\n";
                curr_node->data[0] = new vector<uint32_t>(pile[depth].data[0]);
                curr_node->data[1] = new vector<uint32_t>(pile[depth].data[1]);

            } else if (curr_dim > num_children - 1) {
                // reached end of node!! back up
                cerr << "end of nodes backup!\n";

                back_up = true;
                //(opt_region &working_reg, uint32_t curr_node, int depth, int num_children
                compute_lPs(working_reg, pile[depth].node, depth);

            }

            if (back_up) {
                cerr << "DO BACKUP\n";
                
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
            cerr << "CUT!...\n";
            working_reg.cut(curr_dim,curr_cut);
            uint32_t working_hash = region_cache.hash(working_reg);

            pair<uint32_t,bool> new_node = region_cache.find(working_reg,working_hash);

            if (!new_node.second) {
                cerr << "new node not found...\n";
                
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
                cerr << "node was found...\n";
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
            
            // if sequence id is not reached check this
            // check if current node includes any of the 4 points
            // if doesn't include back up

            bool point_included = false;
            if (curr_node->get_sequence_id() != seq_idx) {
                double lim = curr_reg.get_lim(curr_dim);

                for (int i = 0; i < 4 && !point_included; i++) {
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

                if (!point_included) {
                    back_up = true;

                    // also set the sequence ID to current idx
                    curr_node->set_sequence_id(seq_idx);
                }
            }
            
            

            // check if current node is leaf or at end
            // do this check if the sequence id is reached
            if (curr_node->get_sequence_id() == seq_idx) {
                
                // check leaf criteria
                if ((curr_count[0] + curr_count[1]) <= count_lim
                        || depth >= max_depth
                        || working_reg.full()) {
                    // back up
                    back_up = true;

                    // check if this node is a leaf
                    // If it is a leaf delete all children

                    // check if it has children
                    // probably need to check all children :/
                    for (int i = 0; i < num_children; i++) {
                        // get the child
                        working_reg.cut(curr_dim, curr_cut);
                        uint32_t del_node = region_cache.erase(working_reg);
                        if (del_node != c::ra_null_val) {
                            ra.delete_node(del_node);
                        }
                        // uncut
                        working_reg.uncut(curr_dim);
                    }

                    // assume cuts are the same, so don't nee to search for lP0
                    curr_node->set_uniform(depth);

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
                    compute_lPs(working_reg, pile[depth].node, depth);

                }
                
            }

            if (back_up) {

                depth--;
                pile.pop_back();
                if(depth < 0) continue;
                
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

                curr_reg.cut(curr_dim, curr_cut);
                
                // determine which data points are active
                double lim = curr_reg.get_lim(curr_dim);
                for (int k = 0; k < 2; k++) {
                    for (int i = 0; i < 2; i++) {
                        pile[depth].data[k][i] = pile[depth-1].data[k][i];
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

                    pile[depth].dim = 0;
                    pile[depth].cut = 0;

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
                } else {
                    // update counts
                    // modify counts here 
                    int count[2];
                    curr_node->get_count(count);
                    for (int i = 0; i < 4; i++) {
                        int k = i/2;
                        int j = i % 2;
                        if (pile[depth].data[k][j] == 0) {
                            continue;
                        }

                        if (j == 0) {
                            count[k]++;
                        } else {
                            count[k]--;
                        }
                    }
                    
                    
                    
                    if (ra[new_node.first]->is_leaf()) {
                        // if the region is a leaf
                        // modify the data points in the leaf

                        for(int k = 0;k<2;k++){
                            if(pile[depth].data[k][1] == 1){
                                // del
                                vector<uint32_t>* data_ptr = curr_node->data[k];
                                uint32_t del_pt = pts[(2*k)+1];
                                vector<uint32_t>::iterator del_it = data_ptr->end();
                                for (vector<uint32_t>::iterator it = data_ptr->begin();
                                        it != data_ptr->end();it++){
                                    if(*it == del_pt){
                                        del_it = it;
                                        it = data_ptr->end();
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
                                vector<uint32_t>* data_ptr = curr_node->data[k];
                                data_ptr->push_back(pts[2*k]);
                            }
                        }
                        
                    } else {
                        // if the region is not a leaf
                        // do nothing?
                        
                    }
                }
                
                // do the counting of the data points here
                
                
                pile[depth].node = new_node.first;
                region_cache.insert(working_reg,new_node.first,working_hash);


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
            cerr << "Error: coupling OPT not constructed yet\n";
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
    
    double get_log_coupling_prob(){
        online_ctree_node* curr_node = ra[root];
        return curr_node->get_lphi() - curr_node->get_lP() - c::l2;
    }
};



#endif	/* COPT_TREE_H */

