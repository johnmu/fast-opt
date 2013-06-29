/*
 *  llopt_tree.h
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



#ifndef LLOPT_TREE_H
#define	LLOPT_TREE_H



#include "stl.h"
#include "general_utils.h"
#include "opt_utils.h"
#include "gamma_table.h"
#include "opt_tree.h"




struct ll_tree_node_sparse{

    float lphi; // this is also the density if MAP tree
    int count; // number of points in this region
               // the count is negative if we stop cutting

    void init(){
        count = -1;
        lphi = -c::inf;
    }

    ll_tree_node_sparse() {
        init();
    }

    ll_tree_node_sparse(int a) {
        init();
    }

    double get_lphi_unif(int depth){
        if(count >= 0){
            return (count * depth * c::l2) - c::l2;
        }else{
            return -c::inf;
        }
    }

    double get_lphi2(int depth) {
        if(!isinf(lphi))return lphi;
        else return get_lphi_unif(depth)+c::l2;
    }


};



struct ll_working_unit_t{
    vector<uint32_t > data;
    current_region curr_reg;
    opt_region working_reg;
    uint32_t node_idx;


    ll_working_unit_t(int num_children){
        curr_reg.init(num_children);
        working_reg.init(num_children);
        node_idx = c::ra_null_val;
    }

    ll_working_unit_t(vector<uint32_t > &data, current_region &curr_reg
        , opt_region &working_reg, uint32_t node_idx){


        this->data = data;
        this->curr_reg = curr_reg;
        this->working_reg = working_reg;
        this->node_idx = node_idx;
    }
};


struct ll_pile_t{
    vector<ll_working_unit_t> good_regions;
    vector<uint32_t> map_nodes;
};


struct ll_mouse_params_t{
    int start_dimension;
    int num_children;
    int count_lim;
    int top_depth;  // the depth of the parent
    int top_count_lim;
    int max_depth;
    int top_max_depth;
    ll_working_unit_t* wup;
    opt_region_hash<uint32_t>* region_cache;
    gamma_table *gt;
    region_allocator<ll_tree_node_sparse> *ra;
    int64_t *num_nodes;
    int64_t *num_zero_nodes;
    
    vector<vector<double> > *all_data;

    pthread_mutex_t* locker;
};



class llopt_tree{
private:
    region_allocator<ll_tree_node_sparse> ra;
    uint32_t root;
    int num_children;

    double count_ratio;
    double top_count_ratio;
    int max_depth;
    int top_max_depth;

    void init(int num_children, double count_ratio, double top_count_ratio,
                int max_depth, int top_max_depth) {
        this->num_children = num_children;
        pair<uint32_t,ll_tree_node_sparse*> out = ra.create_node(num_children);
        root = out.first;
        this->count_ratio = count_ratio;
        this->top_count_ratio = top_count_ratio;
        this->max_depth = max_depth;
        this->top_max_depth = top_max_depth;
    }

public:

    llopt_tree(int num_children, double count_ratio,double top_count_ratio,
                    int max_depth, int top_max_depth) {
        init(num_children, count_ratio, top_count_ratio, max_depth,top_max_depth);
    }

    llopt_tree(int num_children) {

        init(num_children, 0.01, 0.001 , 3,num_children*20);

    }


    void store_region(uint32_t map_node, opt_region_hash<uint32_t> &map_regions,
            region_allocator<map_tree_node> *map_ra,opt_region &working_reg, int count, int depth){

        map_regions.insert(working_reg, map_node);

        (*map_ra)[map_node]->set_area(-depth);
        (*map_ra)[map_node]->set_count(count);
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


    static void compute_lphi(region_allocator<ll_tree_node_sparse> &ra,opt_region &working_reg,
            opt_region_hash<uint32_t> &region_cache,
            uint32_t curr_node, int depth, gamma_table &gt, int calling_loc,int num_children) {

        //cerr << "compute_lphi(" << depth << "): ";
        //working_reg.print_region(cerr);
        //cerr << '\n';

        vector<double> lphi_list; // should pre-allocate
        double max_val = (ra[curr_node]->count * depth * c::l2) - c::l2;
        lphi_list.push_back(max_val);

        //if(calling_loc == 3) cerr << "max_val: " << max_val << '\n';

        double ld = -log(num_children) - c::lpi - c::l2;

        //if(calling_loc == 3) cerr << "ld: " << ld << '\n';

        for (int i = 0; i < num_children; i++) {

            uint32_t child_id[2];
            child_id[0] = get_child(working_reg, region_cache,i,0);
            child_id[1] = get_child(working_reg, region_cache,i,1);

            // check for null
            /*
            for (int f = 0; f < c::cuts; f++) {
                if (child_id[f] == c::ra_null_val) {
                    cerr << "(" << calling_loc << ")fully NULL child!!! " << f << "|" << i << ',' << depth << "|"
                            << child_id[0] << "," << child_id[1] << '\n';
                    cerr << "curr_count: " << ra[curr_node]->count << '\n';
                    for (int j = 0; j < num_children; j++) {
                        cerr << "child(" << j << "): "
                                << get_child(working_reg, region_cache, j, 0) << ","
                                << get_child(working_reg, region_cache, j, 1) << "\n";
                    }

                    cerr << "reg: ";
                    working_reg.print_region(cerr);
                    cerr << '\n';

                    exit(2);

                }

            }
            */
            
            
            

            int child_1_count = ra[child_id[0]]->count;
            int child_2_count = ra[child_id[1]]->count;

            
            if (child_1_count < 0 && child_2_count < 0) {
                cerr << "neg count!!! " << i << ',' << depth << '\n';
                exit(2);
            }


            //if(calling_loc == 3){
            //    cerr << i << ":" << child_1_count << ',' << child_2_count << '\n';
            //    cerr << "lphi: " << ra[child_id[0]]->get_lphi2(depth + 1) << ","
            //            << ra[child_id[1]]->get_lphi2(depth + 1) << '\n';
            //    cerr << "gt: " << gt.compute_lD2(ra[curr_node]->count, child_1_count, child_2_count)<<'\n';
            //}

            if(child_1_count < 0){
                child_1_count = -child_1_count;
            }
            
            if(child_2_count < 0){
                child_2_count = -child_2_count;
            }
            
            
            double val = ld;
            val += ra[child_id[0]]->get_lphi2(depth + 1);
            val += ra[child_id[1]]->get_lphi2(depth + 1);
            val += gt.compute_lD2(ra[curr_node]->count, child_1_count, child_2_count);

            //if(calling_loc == 3) cerr << "val: " << val << '\n';

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

        ra[curr_node]->lphi = lphi;

        //if(calling_loc == 3) cerr << "lphi: " << lphi << '\n';
    }


    static void* small_opt_thread(void* params) {

        ll_mouse_params_t* p = (ll_mouse_params_t*)params;

        int start_dimension = p->start_dimension;
        int top_depth = p->top_depth;
        int num_children = p->num_children;
        int count_lim = p->count_lim;
        int max_depth = p->max_depth;
        int top_max_depth = p->top_max_depth;
        ll_working_unit_t* wup = p->wup;
        opt_region_hash<uint32_t>* region_cache = p->region_cache;
        gamma_table *gt = p->gt;
        region_allocator<ll_tree_node_sparse> *ra = p->ra;
        int64_t *num_nodes = p->num_nodes;
        int64_t *num_zero_nodes = p->num_zero_nodes;
        //pthread_mutex_t* locker = p->locker;

        vector<vector<double> > *all_data = p->all_data;

        vector<pile_t<uint32_t,uint32_t> > pile;
        pile.push_back(pile_t<uint32_t,uint32_t> ());
        pile[0].node = wup->node_idx;
        pile[0].dim = start_dimension;
        pile[0].cut = 0;
        pile[0].data = wup->data;
        

        current_region curr_reg = wup->curr_reg;
        opt_region working_reg = wup->working_reg;

        bool done = false;
        int depth = 0;

        //cerr << "%%%%!!!----dimension" << start_dimension << "\n";

        while (!done) {

            if (pile.size() == 0) {
                done = true;
                continue;
            }


            uint32_t curr_node_idx = pile[depth].node;

            //pthread_mutex_lock(locker);
            ll_tree_node_sparse curr_node = *((*ra)[curr_node_idx]);
            //pthread_mutex_unlock(locker);

            bool backup = false;

            //cerr << "----\nDEPTH("<< depth <<"): " << (depth + top_depth)
            //        << " count: "<< curr_node.count << '\n';
            //cerr << "START dim:cut --- " << pile[depth].dim << ":" << pile[depth].cut << '\n';


            if (curr_node.count <= count_lim
                    || (depth +top_depth) >= top_max_depth
                    || (depth) >= max_depth
                    || working_reg.full()) {

                //cerr << "BOTTOM BACKUP\n";

                // make sure we at least cut once
                if(depth != 0) backup = true;
            } else {
                if ((depth == 0&&pile[depth].dim !=start_dimension)
                        || pile[depth].dim > num_children - 1){

                    //cerr << "CUTALL BACKUP\n";
                    if(depth != 0){
                        compute_lphi(*ra, working_reg, *region_cache,
                            curr_node_idx, depth + top_depth, *gt,
                            1, num_children);

                        //cerr << "COMPUTE LPHI: " << (*((*ra)[curr_node_idx])).lphi << '\n';

                    }
                    backup = true;
                }

            }

            bool is_diff = true;

            if (!backup) {
                int curr_dim = pile[depth].dim;
                int curr_cut = pile[depth].cut;

                //cerr << "NB dim:cut --- " << curr_dim << ":" << curr_cut << '\n';

                pile.push_back(pile_t<uint32_t,uint32_t > ());
                depth++;

                is_diff = cut_region_one(*all_data, pile[depth - 1].data, pile[depth].data,
                        curr_dim, curr_cut, curr_reg.get_lim(curr_dim));

                curr_reg.cut(curr_dim, curr_cut);
                working_reg.cut(curr_dim, curr_cut);

                pile[depth].dim = 0;
                pile[depth].cut = 0;

                int curr_count = pile[depth].data.size();

                uint32_t working_hash = region_cache->hash(working_reg);

                //pthread_mutex_lock(locker);
                pair<uint32_t, bool> new_node = region_cache->find(working_reg, working_hash);
                //pthread_mutex_unlock(locker);

                if (!new_node.second) {

                    //cerr << "not found: " << curr_count << "\n";

                    // MUTEX
                    //pthread_mutex_lock(locker);
                    pair<uint32_t, ll_tree_node_sparse*> out = ra->create_node(num_children);
                    new_node.first = out.first;
                    
                    if(is_diff){
                        out.second->count = curr_count;
                    }else{
                        out.second->count = -curr_count;
                    }

                    //cerr << "INSERT" << '\n';
                    region_cache->insert(working_reg, new_node.first, working_hash);


                    (*num_nodes)++;
                    if ((*ra)[new_node.first]->count <= count_lim)(*num_zero_nodes)++;

                    pile[depth].node = new_node.first;

                    // UNMUTEX
                    //pthread_mutex_unlock(locker);

                } else {

                    //cerr << "found node(" << (*ra)[new_node.first]->count << "): " << curr_dim << "," << curr_cut << '\n';

                    //ra[curr_node]->set_child(curr_dim, curr_cut, new_node.first);

                    pile[depth].node = new_node.first;
                    //cerr << "fUNCUT: " << curr_dim  << '\n';

                }
                


            }

            if (backup) {

                //cerr << "BACKUP\n";
                //pthread_mutex_lock(locker);
                // estimate phi

                //int curr_count = curr_node.count;
                //int curr_depth = (depth + top_depth);

                //cerr << "curr_count: " << curr_count << '\n';


                //pthread_mutex_unlock(locker);

                //cerr << "BEFORE depth: " << depth << '\n';

                depth--;
                pile.pop_back();
                if (depth < 0) continue;

                //cerr << "BEFORE dim:cut --- " << pile[depth].dim << ":" << pile[depth].cut << '\n';
                //working_reg.print_region();
                //cerr << '\n';
                //cerr << "AFTER depth: " << depth << '\n';

                curr_reg.uncut(pile[depth].dim, pile[depth].cut);
                working_reg.uncut(pile[depth].dim);

                //working_reg.print_region();
                //cerr << '\n';

                if (pile[depth].cut < c::cuts - 1) {
                    pile[depth].cut++;
                } else if (pile[depth].dim <= num_children - 1) {
                    pile[depth].dim++;
                    pile[depth].cut = 0;
                }

                //cerr << "AFTER dim:cut --- " << pile[depth].dim << ":" << pile[depth].cut << '\n';

                continue;
            }
        }
        return NULL;
    }


    void do_small_opt(vector<vector<double> > *all_data, ll_working_unit_t &w,
            opt_region_hash<uint32_t> &region_cache,gamma_table &gt,
            MT_random &rand_gen,
            int64_t &num_nodes, int64_t &num_zero_nodes, int start_depth,
            int top_count_lim) {

        int count_lim = (int)floor(w.data.size()*count_ratio);
        if(count_lim < top_count_lim) count_lim = top_count_lim;

        //pthread_t* t_group = new pthread_t[num_children]; // just do 3 threads for now :)
        //pthread_mutex_t* locker = new pthread_mutex_t();
        //pthread_mutex_init(locker, NULL);

        ll_mouse_params_t* params = new ll_mouse_params_t[num_children];

        for (int d = 0; d < num_children; d++) {
            params[d].count_lim = count_lim;
            params[d].top_count_lim = top_count_lim;
            params[d].gt = &gt;
            //params[d].locker = locker;
            params[d].max_depth = max_depth;
            params[d].top_max_depth = top_max_depth;
            params[d].num_children = num_children;
            params[d].num_nodes = &num_nodes;
            params[d].num_zero_nodes = &num_zero_nodes;
            params[d].ra = &ra;
            params[d].region_cache = &region_cache;
            params[d].start_dimension = d;
            params[d].top_depth = start_depth;
            params[d].wup = &w;
            params[d].all_data = all_data;
        }

        for (int d = 0; d < num_children; d++) {
            // Start threads

            small_opt_thread((void*) &(params[d]));
            
            //pthread_create(&(t_group[d]), NULL, small_opt_thread, (void*) &(params[d]));
        // no threading for now
        //}


        // wait for threads
        //for (int d = 0; d < num_children; d++) {
            //void* output;
            //pthread_join(t_group[d], &output);
        }
        
        //delete locker;
        delete [] params;
        //delete [] t_group;
    }

int get_map_dim(ll_working_unit_t &w,opt_region_hash<uint32_t> &region_cache,gamma_table &gt,
                    int map_depth) {
        // Choose MAP dimension

        int map_dim = 0;

        vector<uint32_t>child_idx_0(num_children);
        vector<uint32_t>child_idx_1(num_children);

        for (int i = 0; i < num_children; i++) {

            if(w.working_reg.full()){
                w.working_reg.print_region();
                cerr << "\n";
            }

            child_idx_0[i] = get_child(w.working_reg, region_cache, i, 0);
            child_idx_1[i] = get_child(w.working_reg, region_cache, i, 1);

        }




        double max_post_prob = -c::inf;

        for (int i = 0; i < num_children; i++) {
            double post_prob = gt.compute_lD2(ra[w.node_idx]->count
                    , ra[child_idx_0[i]]->count, ra[child_idx_1[i]]->count);
            post_prob += ra[child_idx_0[i]]->get_lphi2(map_depth + 1);
            post_prob += ra[child_idx_1[i]]->get_lphi2(map_depth + 1);
            //#ifdef DEBUG_MAP2
            //cerr << "gt: " << gt.compute_lD2(ra[wit->node_idx]->count
            //        , ra[child_idx_0[i]]->count, ra[child_idx_1[i]]->count) <<'\n';

            //cerr << "lphi_3_0: " <<  ra[child_idx_0[i]]->get_lphi3(map_depth + 1) << '\n';
            //cerr << "lphi_3_1: " <<  ra[child_idx_1[i]]->get_lphi3(map_depth + 1) << '\n';

            //cerr << "post_prob[" << ra[child_idx_0[i]]->count << "|" << ra[child_idx_1[i]]->count << "](" << i << ") = " << post_prob << '\n';
            //#endif
            if (post_prob > max_post_prob) {
                map_dim = i;
                max_post_prob = post_prob;
            }
        }

        //#ifdef DEBUG_MAP2
        //cerr << map_dim << ',';
        //#endif

        return map_dim;

    }

    void construct_llopt_tree(vector<vector<double> > *all_data,
            map_tree &map_region_tree, opt_region_hash<uint32_t> &map_regions,
            bool prune_tree) {

        MT_random rand_gen;

        // need to treat the empty tree/data case

        int64_t num_nodes = 0;
        int64_t num_zero_nodes = 0;

        uint32_t N = all_data->size();
        vector<uint32_t> data(N,0);
        for (uint32_t i = 0;i<N;i++){
            data[i] = i;
        }
        
        int count_lim = 0;

        if (top_count_ratio < 1) {
            count_lim = (int) floor((double) N * top_count_ratio); // this is for the whole tree
        } else {
            count_lim = (int)top_count_ratio;
        }
        
        if(count_lim < 1)count_lim = 1;

        cerr << "Full tree stopping at " << count_lim << " points\n";
        cerr << "Each look-ahead stopping at " << max_depth << " levels\n";
        cerr << "Full tree stopping at " << top_max_depth << " levels\n";

        gamma_table gt(N);
        ra[root]->count = N;

        region_allocator<map_tree_node> *map_ra = map_region_tree.get_ra();

        opt_region_hash<uint32_t> region_cache(27);

        current_region start_region(num_children);
        opt_region     start_working(num_children);


        vector<ll_pile_t> llpile;
        llpile.push_back(ll_pile_t());

        llpile[0].good_regions.push_back(ll_working_unit_t(data,start_region,
                start_working,root));

        llpile[0].map_nodes.push_back(map_region_tree.get_full_tree());

        int map_depth = 0;
        int max_depth = -1;

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

            vector<ll_working_unit_t>::iterator wu_it = (llpile[map_depth].good_regions.end() - 1);
            vector<uint32_t>::iterator map_node_it = (llpile[map_depth].map_nodes.end() - 1);

            // if too deep we stop all the regions
            if (map_depth >= top_max_depth || (int)wu_it->data.size() <= count_lim
                    || wu_it->working_reg.full()) {

                //cerr << "mini-STOP count: " << ra[wu_it->node_idx]->count << '\n';

                store_region(*map_node_it, map_regions, map_ra, wu_it->working_reg,
                        ra[wu_it->node_idx]->count, map_depth);

                llpile[map_depth].good_regions.pop_back();
                llpile[map_depth].map_nodes.pop_back();
                count++;

                continue;

            }
            
            // for each good region
            int curr_data_size = wu_it->data.size();
            if(map_depth>max_depth){
                max_depth = map_depth;
                cerr << "Depth("<< map_depth <<"):Size("<< curr_data_size <<")\n "; 
            }
            //cerr << "Good [" << map_depth << "](" << count << "/"
            //        << llpile[map_depth].good_regions.size()
            //        << "), Data size: " << curr_data_size << '\n';

            // do a smaller OPT for the good region to determine which
            // dimension to cut

            do_small_opt(all_data, *wu_it,
                    region_cache, gt, rand_gen,
                    num_nodes, num_zero_nodes, map_depth, count_lim);

            //cerr << "\nNodes:" << num_nodes
            //        << " : " << num_zero_nodes
            //        << " : " << (num_nodes - num_zero_nodes) << '\n';


            compute_lphi(ra,wu_it->working_reg, region_cache,
                    wu_it->node_idx, map_depth, gt, 3,num_children);

            //cerr << "lphi: " << ra[wu_it->node_idx]->get_lphi2(map_depth) << "\n";


            // compute the posterior pho and decide if we stop
            double post_rho = -c::l2;
            post_rho += map_depth * c::l2 * ra[wu_it->node_idx]->count; // phi_0
            post_rho -= ra[wu_it->node_idx]->get_lphi2(map_depth);

            // if not stop we split and add the two regions to the good region list
            if (post_rho>-c::l2) {
                // STOP

                //cerr << "STOP count: " << ra[wu_it->node_idx]->count << '\n';

                store_region(*map_node_it, map_regions, map_ra, wu_it->working_reg,
                        ra[wu_it->node_idx]->count, map_depth);

                llpile[map_depth].good_regions.pop_back();
                llpile[map_depth].map_nodes.pop_back();
                count++;

            } else {
                // Choose MAP dimension

                //cerr << "MAP dim: ";
                int map_dim = get_map_dim(*wu_it, region_cache, gt, map_depth);
                //cerr << '\n';

                //cerr << "DONE MAP\n";


                uint32_t map_child_idx_0 = get_child(wu_it->working_reg, region_cache, map_dim, 0);
                uint32_t map_child_idx_1 = get_child(wu_it->working_reg, region_cache, map_dim, 1);

                //cerr << "children: " << map_child_idx_0 << "," << map_child_idx_1 << '\n';

                // Split
                vector<uint32_t > new_data_0;
                vector<uint32_t > new_data_1;

                bool is_diff = cut_region2_one(*all_data,wu_it->data,new_data_0,new_data_1,
                        map_dim, wu_it->curr_reg.get_lim(map_dim));

                if(is_diff){

                    current_region next_curr_reg0 = wu_it->curr_reg;
                    opt_region next_working_reg0 = wu_it->working_reg;

                    next_curr_reg0.cut(map_dim, 0);
                    next_working_reg0.cut(map_dim, 0);

                    current_region next_curr_reg1 = wu_it->curr_reg;
                    opt_region next_working_reg1 = wu_it->working_reg;

                    next_curr_reg1.cut(map_dim, 1);
                    next_working_reg1.cut(map_dim, 1);


                    // after this point wu_it is invalid
                    llpile.push_back(ll_pile_t());
                    map_depth++;

                    map_node_it = (llpile[map_depth - 1].map_nodes.end() - 1);

                    //cerr << "ADD PILE " << '\n';

                    llpile[map_depth].good_regions.push_back(ll_working_unit_t(new_data_0,
                            next_curr_reg0, next_working_reg0, map_child_idx_0));

                    llpile[map_depth].good_regions.push_back(ll_working_unit_t(new_data_1,
                            next_curr_reg1, next_working_reg1, map_child_idx_1));

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

                    llpile[map_depth-1].good_regions.pop_back();
                    llpile[map_depth-1].map_nodes.pop_back();

                } else {
                    // store as if it is uniform
                    store_region(*map_node_it, map_regions, map_ra, wu_it->working_reg,
                        ra[wu_it->node_idx]->count, map_depth);

                    llpile[map_depth].good_regions.pop_back();
                    llpile[map_depth].map_nodes.pop_back();

                    count++;
                }


            }


            // remove all the nodes not consistent with current good regions

            if(prune_tree && ((curr_data_size > 5000 && ra.free_locs.size() < 1000000) || (curr_data_size > 10000))) {
                int64_t num_removed = 0;
                int64_t total_nodes = 0;

                for (uint32_t i = 0; i < region_cache.table_size; i++) {
                    if (region_cache.map_table[i] != NULL) {

                        map<opt_region, uint32_t> *new_map = new map<opt_region, uint32_t > ();

                        for (map<opt_region, uint32_t>::iterator it = region_cache.map_table[i]->begin();
                                it != region_cache.map_table[i]->end(); it++) {

                            bool is_child = false;

                            for (vector<ll_pile_t>::iterator ls_it = llpile.begin();
                                    ls_it != llpile.end() && !is_child; ls_it++) {
                                for (vector<ll_working_unit_t>::iterator reg_it = ls_it->good_regions.begin();
                                        reg_it != ls_it->good_regions.end() && !is_child; reg_it++) {

                                    if (reg_it->working_reg.is_child(it->first)) {
                                        is_child = true;
                                    }

                                }
                            }
                            if (is_child) {
                                new_map->insert(pair<opt_region, uint32_t > (it->first, it->second));
                            } else {
                                ra.delete_node(it->second);
                                num_removed++;
                            }

                            total_nodes++;

                        }

                        delete region_cache.map_table[i];
                        region_cache.map_table[i] = new_map;
                    }
                }

                //if(total_nodes>0)cerr << "Removed " << num_removed << "/" << total_nodes
                //        << "(" << (num_removed / (total_nodes / 100.0)) << "%)" << '\n';
            }

        }


    }


};

#endif	/* LLOPT_TREE_H */

