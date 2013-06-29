/*
 *  lsopt_tree.h
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



#ifndef LSOPT_TREE_H
#define	LSOPT_TREE_H

#include "stl.h"
#include "general_utils.h"
#include "opt_utils.h"
#include "gamma_table.h"
#include "opt_tree.h"




struct ls_tree_node_sparse{

    float lphi; // this is also the density if MAP tree
    float lphi_cut;
    int count; // number of points in this region
    int num_visits;

    void init(){
        num_visits = 0;
        count = -1;
        lphi = -c::inf;
        lphi_cut = -c::inf;
    }

    ls_tree_node_sparse() {
        init();
    }

    ls_tree_node_sparse(int a) {
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

    double get_lphi3(int depth) {
        if(!isinf(lphi))return lphi;
        else return logsumexp(lphi_cut,get_lphi_unif(depth));
    }



};



struct ls_working_unit_t{
    vector<vector<double> > data;
    current_region curr_reg;
    opt_region working_reg;
    uint32_t node_idx;


    ls_working_unit_t(int num_children){
        curr_reg.init(num_children);
        working_reg.init(num_children);
        node_idx = c::ra_null_val;
    }

    ls_working_unit_t(vector<vector<double> > &data, current_region &curr_reg
        , opt_region &working_reg, uint32_t node_idx){


        this->data = data;
        this->curr_reg = curr_reg;
        this->working_reg = working_reg;
        this->node_idx = node_idx;
    }
};


struct ls_pile_t{
    vector<ls_working_unit_t> good_regions;
    vector<uint32_t> map_nodes;
};


struct mouse_params_t{
    int start_dimension;
    int top_depth;
    int num_children;
    int count_lim;
    int top_count_lim;
    int max_depth;
    int lookahead_max_depth;
    ls_working_unit_t* wup;
    opt_region_hash<uint32_t>* region_cache;
    gamma_table *gt;
    region_allocator<ls_tree_node_sparse> *ra;
    int64_t *num_nodes;
    int64_t *num_zero_nodes;

    pthread_mutex_t* locker;
};



class lsopt_tree{
private:
    region_allocator<ls_tree_node_sparse> ra;
    uint32_t root;
    int num_children;

    double count_ratio;
    double top_count_ratio;
    int max_depth;
    int lookahead_max_depth;

    void init(int num_children, double count_ratio, double top_count_ratio, int lookahead_max_depth,int max_depth) {
        this->num_children = num_children;
        pair<uint32_t,ls_tree_node_sparse*> out = ra.create_node(num_children);
        root = out.first;
        this->count_ratio = count_ratio;
        this->top_count_ratio = top_count_ratio;
        this->max_depth = max_depth;
        this->lookahead_max_depth = lookahead_max_depth;
    }

public:

    lsopt_tree(int num_children, double count_ratio,double top_count_ratio, int lookahead_max_depth, int max_depth) {
        init(num_children, count_ratio, top_count_ratio, lookahead_max_depth, max_depth);
    }

    lsopt_tree(int num_children) {

        init(num_children, 0.01, 0.001 , 10, 100);

    }


    void store_region(uint32_t map_node, opt_region_hash<uint32_t> &map_regions,
            region_allocator<map_tree_node> *map_ra,opt_region &working_reg, int count, int depth){

        map_regions.insert(working_reg, map_node);

        (*map_ra)[map_node]->set_area(-depth);
        (*map_ra)[map_node]->set_count(count);
    }


    static uint32_t get_child(opt_region working_reg, opt_region_hash<uint32_t> &region_cache,
            int dim, int cut){
        return get_child(working_reg, region_cache,dim, cut, false);
    }

    static uint32_t get_child(opt_region working_reg, opt_region_hash<uint32_t> &region_cache,
            int dim, int cut, bool debug){

        if(debug){
            cerr << "before reg:";
            working_reg.print_region();
            cerr << '\n';
        }

        if(!working_reg.cut(dim,cut)){
            cerr << "CANNOT CUT: ";
            working_reg.print_region();
            cerr << '\n';

            exit(2);
        }

        if(debug){
            cerr << "after reg:";
            working_reg.print_region();
            cerr << '\n';
        }

        pair<uint32_t,bool> out = region_cache.find(working_reg);

        if(debug){
            cerr << out.first << ',' << out.second << '\n';
        }
        if(out.second){
            return out.first;
        }else{
            return c::ra_null_val;
        }
    }

    static void compute_approx_lphi_cut(region_allocator<ls_tree_node_sparse> &ra,opt_region &working_reg,
                opt_region_hash<uint32_t> &region_cache,
                uint32_t curr_node, int depth, gamma_table &gt, int calling_loc,
                vector<vector<double> > &data, current_region &curr_reg, int num_children){

        //cerr << "compute_approx_lphi_cut(" << depth << ")["<< ra[curr_node]->count << "]: ";
        //working_reg.print_region(cerr);
        //cerr << '\n';


        vector<double> lphi_list;
        double max_val = -c::inf;

        double ld = -log(num_children) - c::lpi - c::l2;

        //cerr << "ld: " << ld << '\n';

        for (int i = 0; i < num_children; i++) {
            // get the child ids
            uint32_t child_id[2];
            child_id[0] = get_child(working_reg, region_cache,i,0);
            child_id[1] = get_child(working_reg, region_cache,i,1);


            // check for null
            for (int f = 0; f < c::cuts; f++) {
                if (child_id[f] == c::ra_null_val) {
                    cerr << "a(" << calling_loc << ")NULL child!!! " << f << "|" << i << ',' << depth << "|"
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
            int child_1_count = ra[child_id[0]]->count;
            int child_2_count = ra[child_id[1]]->count;

            if (child_1_count < 0 || child_2_count < 0) {
                cerr << "aneg count!!! " << i << ',' << depth << '\n';
                exit(2);
            }


            //cerr << i << ":" << child_1_count << ',' << child_2_count << '\n';

            double val = ld;

            for (int x = 0; x < 2; x++) {
                if (isinf(ra[child_id[x]]->lphi)) {
                    if (isinf(ra[child_id[x]]->lphi_cut)) {

                        //cerr << "unif: "<< ra[child_id[x]]->get_lphi_unif(depth+1) + c::l2 << "\n";

                        val += ra[child_id[x]]->get_lphi_unif(depth+1) + c::l2;
                    } else {
                        // this is approximate

                        //cerr << "cut: " << logsumexp(ra[child_id[x]]->lphi_cut
                        //        , ra[child_id[x]]->get_lphi_unif(depth+1)) << "\n";

                        val += logsumexp(ra[child_id[x]]->lphi_cut
                                , ra[child_id[x]]->get_lphi_unif(depth+1));
                    }
                } else {

                    //cerr << "found: " << ra[child_id[x]]->lphi <<"\n";

                    val += ra[child_id[x]]->lphi;
                }
            }

            val += gt.compute_lD2(ra[curr_node]->count, child_1_count, child_2_count);

            //cerr << "gt: " << gt.compute_lD2(ra[curr_node]->count, child_1_count, child_2_count) << '\n';

            //cerr << "val: " << val << '\n';

            lphi_list.push_back(val);

            if (val > max_val) {
                max_val = val;
            }

        }

        double lphi_cut = max_val;
        double sum = 0;
        for (int i = 0; i < num_children; i++) {
            sum += exp(lphi_list[i] - max_val);
        }

        if (sum > 0) lphi_cut += log(sum);

        ra[curr_node]->lphi_cut = lphi_cut;

        //cerr << "lphi_cut: " << lphi_cut << '\n';

    }

    static void compute_lphi(region_allocator<ls_tree_node_sparse> &ra,opt_region &working_reg,
            opt_region_hash<uint32_t> &region_cache,
            uint32_t curr_node, int depth, gamma_table &gt, int calling_loc,
                vector<vector<double> > &data, current_region &curr_reg, int num_children) {

        //cerr << "compute_lphi(" << depth << "): ";
        //working_reg.print_region(cerr);
        //cerr << '\n';

        vector<double> lphi_list;
        double max_val = (ra[curr_node]->count * depth * c::l2) - c::l2;
        lphi_list.push_back(max_val);

        //cerr << "max_val: " << max_val << '\n';

        double ld = -log(num_children) - c::lpi - c::l2;

        //cerr << "ld: " << ld << '\n';

        for (int i = 0; i < num_children; i++) {

            uint32_t child_id[2];
            child_id[0] = get_child(working_reg, region_cache,i,0);
            child_id[1] = get_child(working_reg, region_cache,i,1);

            // check for null
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

            int child_1_count = ra[child_id[0]]->count;
            int child_2_count = ra[child_id[1]]->count;

            if (child_1_count < 0 || child_2_count < 0) {
                cerr << "neg count!!! " << i << ',' << depth << '\n';
                exit(2);
            }


            //cerr << i << ":" << child_1_count << ',' << child_2_count << '\n';

            double val = ld;
            val += ra[child_id[0]]->get_lphi2(depth + 1);
            val += ra[child_id[1]]->get_lphi2(depth + 1);
            val += gt.compute_lD2(ra[curr_node]->count, child_1_count, child_2_count);

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

        ra[curr_node]->lphi = lphi;

        //cerr << "lphi: " << lphi << '\n';

    }


/*
    void propagate_lphi(ls_working_unit_t &w,
            opt_region_hash<uint32_t> &region_cache,gamma_table &gt, int start_depth, int N) {

        vector<pile_t<uint32_t> > pile;
        pile.push_back(pile_t<uint32_t>());

        pile[0].data = vector<vector<double> >(); // should not need any data
        pile[0].node = w.node_idx;
        pile[0].dim = 0;
        pile[0].cut = 0;

        opt_region working_reg = w.working_reg;
        current_region curr_reg = w.curr_reg;

        int count_lim = (int)floor(N*count_ratio);
        if(count_lim < 1) count_lim = 1;

        int depth = 0;

        bool done = false;
        while (!done) {

            if (pile.size() == 0) {
                done = true;
                continue;
            }

#ifdef DEBUG3
            // print the pile
            cerr << "@@@---- dp:" << depth+start_depth << " , ps:" << pile.size() << " ---\n";
            //cerr << "REG: ";
            //working_reg.print_region(cerr);
            cerr << '\n';
            //for (int i = 0; i < (int) pile.size(); i++) {
                //cerr << "Data size: " << pile[i].data.size() << '\n';
            //    cerr << "Dim:       " << pile[i].dim << '\n';
            //    cerr << "Cut:       " << pile[i].cut << '\n';
                //if (ra[pile[i].node]->get_child(pile[i].dim, pile[i].cut) != c::ra_null_val) {
                //    cerr << "Child lphi: " << ra[ra[pile[i].node]->get_child(pile[i].dim, pile[i].cut)]->get_lphi3(depth+start_depth) << '\n';
                //}
                //if(pile[i].data.size()<20){
                //    print_data(pile[i].data);
                //}
            //}
#endif


            int curr_dim = pile[depth].dim;
            int curr_cut = pile[depth].cut;
            uint32_t curr_node = pile[depth].node;

            //cerr << "curr_dim: " << curr_dim << "| curr_cut: "<< curr_cut << '\n';

            // work out what to count

            bool back_up = false;
#ifdef DEBUG3
            if (ra[curr_node]->get_child(curr_dim, curr_cut) != NULL) {
                cerr << "@Child lphi: " << ra[curr_node]->get_child(curr_dim, curr_cut)->get_lphi() << '\n';
            } else {
                cerr << "@NULL node " << curr_dim << "," << curr_cut << " \n";
            }
#endif

            if(curr_dim >= num_children){
                // reached the end of node

                //cerr << "BACKUPEND\n";

                back_up = true;

                //cerr << "prev-lphicut: " << ra[curr_node]->lphi_cut << '\n';

                compute_approx_lphi_cut(working_reg, region_cache,curr_node,
                        depth+start_depth,gt,0,pile[depth].data,curr_reg); // last two params are wrong


#ifdef DEBUG3
                cerr << "LPHI---- " << depth+start_depth << " , " << pile.size() << " ---\n";
                cerr << "REG: ";
                working_reg.print_region(cerr);
                cerr << '\n';
                //cerr << num_nodes << " : " << num_zero_nodes << '\n';
                cerr << "lphicut: " << ra[curr_node]->lphi_cut << '\n';
                cerr << "lphi: " << ra[curr_node]->get_lphi3(depth+start_depth)
                        << "," << ra[curr_node]->get_lphi_unif(depth+start_depth)<< '\n';
                //for (int i = 0; i < (int) pile.size(); i++) {
                //    cerr << "Data size: " << pile[i].data.size() << '\n';
                //    cerr << "Dim:       " << pile[i].dim << '\n';
                //    cerr << "Cut:       " << pile[i].cut << '\n';
                //if(pile[i].data.size()<20){
                //    print_data(pile[i].data);
                //}
                //}


                cerr << "BACKUP CHILD\n";

#endif
            }

            // check if current node is leaf or at end
            if (get_child(working_reg, region_cache,0, 0) == c::ra_null_val||
                    !isinf(ra[curr_node]->lphi) ||
                    ra[curr_node]->count <= count_lim || depth+start_depth >= max_depth) {
                // back up
                back_up = true;

#ifdef DEBUG3
                cerr << "BACKUP ZERO\n";
#endif
            }


            if (back_up) {

                depth--;
                pile.pop_back();
                if (depth < 0) continue;

                //cerr << depth << " bUNCUT: " << curr_dim << '\n';
                working_reg.uncut(pile[depth].dim);

                if (pile[depth].cut < c::cuts - 1) {
                    pile[depth].cut++;
                } else{
                    pile[depth].dim++;
                    pile[depth].cut = 0;
                }

                continue;
            }


            curr_dim = pile[depth].dim;
            curr_cut = pile[depth].cut;

            // do the counting



            uint32_t new_node = get_child(working_reg, region_cache,curr_dim, curr_cut);
            //cerr << "CUT: " << curr_dim << "," << curr_cut << '\n';
            working_reg.cut(curr_dim, curr_cut);
            //uint32_t working_hash = region_cache.hash(working_reg);
#ifdef DEBUG
            cerr << "working_hash: " << working_hash << '\n';
#endif


#ifdef DEBUG3
            cerr << "FOUND node: " << curr_dim << "," << curr_cut << '\n';

            cerr << "new lphi:  " << new_node->get_lphi() << '\n';
            cerr << "new count: " << new_node->get_count() << '\n';
#endif


            //uint32_t new_node = ra[curr_node]->get_child(curr_dim, curr_cut);


            pile.push_back(pile_t<uint32_t > ());
            depth++;


            pile[depth].dim = 0;
            pile[depth].cut = 0;

            pile[depth].node = new_node;





        }


    }
 */

    static void* iteration_thread(void* params) {

        mouse_params_t* p = (mouse_params_t*)params;

        int start_dimension = p->start_dimension;
        int top_depth = p->top_depth;
        int num_children = p->num_children;
        int count_lim = p->count_lim;
        int top_count_lim = p->top_count_lim;
        int max_depth = p->max_depth;
        int lookahead_max_depth = p->lookahead_max_depth;
        ls_working_unit_t* wup = p->wup;
        opt_region_hash<uint32_t>* region_cache = p->region_cache;
        gamma_table *gt = p->gt;
        region_allocator<ls_tree_node_sparse> *ra = p->ra;
        int64_t *num_nodes = p->num_nodes;
        int64_t *num_zero_nodes = p->num_zero_nodes;
        pthread_mutex_t* locker = p->locker;

        MT_random rand_gen;

        vector<pile_t<uint32_t,vector<double> > > pile;
        pile.push_back(pile_t<uint32_t,vector<double>  > ());

        pile[0].node = wup->node_idx;
        pile[0].dim = start_dimension;
        pile[0].cut = 0;
        pile[0].data = wup->data;

        bool first_cut = false; // we choose the first dimension to cut.

        current_region curr_reg = wup->curr_reg;
        opt_region working_reg = wup->working_reg;

        vector<double> curr_dim_score(num_children, -1.0);

        bool done = false;
        int depth = 0;

        // iterate with a binary tree until bottom

        //cerr << "!!!----dimension" << start_d << "\n";

        while (!done) {

            if (pile.size() == 0) {
                done = true;
                continue;
            }


            int curr_dim = pile[depth].dim;
            int curr_cut = pile[depth].cut;
            uint32_t curr_node_idx = pile[depth].node;

            pthread_mutex_lock(locker);
            ls_tree_node_sparse curr_node = *((*ra)[curr_node_idx]);
            pthread_mutex_unlock(locker);

            bool backup = false;

            //cerr << "----\nDEPTH: " << (depth + top_depth)  << " count: "<< curr_node.count << '\n';


            if (curr_node.count <= count_lim
                    || (depth +top_depth) >= max_depth
                    || depth >= lookahead_max_depth
                    || working_reg.full()) {

                //cerr << "BOTTOM BACKUP\n";

                backup = true;
            } else if (first_cut && curr_dim >= 0) {
                if (curr_cut < c::cuts - 1) {
                    curr_cut++;
                    pile[depth].cut++;
                } else {

                    //cerr << "CUTALL BACKUP\n";

                    backup = true;
                }
            }


            if (!backup && !first_cut) {


                for (int f = 0; f < 2; f++) {
                    pthread_mutex_lock(locker);
                    uint32_t child = get_child(working_reg, *(region_cache), start_dimension, f);
                    pthread_mutex_unlock(locker);
                    //uint32_t child = ra[curr_node]->get_child(d, f);

                    if (child == c::ra_null_val) {
                        // count it
                        working_reg.cut(start_dimension, f);
                        uint32_t working_hash = region_cache->hash(working_reg);

                        pthread_mutex_lock(locker);
                        pair<uint32_t, bool> new_node = region_cache->find(working_reg, working_hash);
                        pthread_mutex_unlock(locker);

                        //int curr_count = 0;
                        if (!new_node.second) {

                            int curr_count = count_region(pile[depth].data, start_dimension
                                    , f, curr_reg.get_lim(start_dimension));


                            // MUTEX

                            pthread_mutex_lock(locker);
                            pair<uint32_t, ls_tree_node_sparse*> out = ra->create_node();

                            new_node.first = out.first;
                            out.second->count = curr_count;

                            (*num_nodes)++;
                            if (curr_count <= count_lim)(*num_zero_nodes)++;

                            region_cache->insert(working_reg, new_node.first, working_hash);

                            // END MUTEX
                            pthread_mutex_unlock(locker);

                        }

                        working_reg.uncut(start_dimension);

                    }

                }


            }


            // work out which dimension to cut

            if (!backup && curr_dim < 0) { // if we have not chosen a dimension yet.



                bool any_good = false;
                for (int d = 0; !backup && (d < num_children); d++) {
                    double child_lphi_unif[2] = {-1.0, -1.0};
                    double child_lphi_cut[2] = {-1.0, -1.0};
                    int child_lphi_count[2] = {-1, -1};
                    int child_visit_count[2] = {0, 0};

                    bool computed = true;

                    for (int f = 0; f < 2; f++) {

                        pthread_mutex_lock(locker);
                        uint32_t child = get_child(working_reg, *region_cache, d, f);

                        pthread_mutex_unlock(locker);
                        //uint32_t child = ra[curr_node]->get_child(d, f);

                        if (child == c::ra_null_val) {
                            // count it
                            working_reg.cut(d, f);
                            uint32_t working_hash = region_cache->hash(working_reg);

                            pthread_mutex_lock(locker);
                            pair<uint32_t, bool> new_node = region_cache->find(working_reg, working_hash);
                            pthread_mutex_unlock(locker);

                            ls_tree_node_sparse* childp;

                            int curr_count = 0;
                            if (!new_node.second) {

                                curr_count = count_region(pile[depth].data, d, f, curr_reg.get_lim(d));

                                //cerr << "curr_count: " << curr_count << '\n';

                                // MUTEX
                                pthread_mutex_lock(locker);
                                pair<uint32_t, ls_tree_node_sparse*> out = ra->create_node();


                                new_node.first = out.first;
                                out.second->count = curr_count;
                                computed = false;


                                (*num_nodes)++;
                                if (curr_count <= count_lim)(*num_zero_nodes)++;

                                region_cache->insert(working_reg, new_node.first, working_hash);

                                childp = out.second;

                                // END MUTEX
                                pthread_mutex_unlock(locker);


                                //cerr << "New ";

                            } else {
                                pthread_mutex_lock(locker);
                                childp = (*ra)[new_node.first];
                                pthread_mutex_unlock(locker);
                                //cerr << "Old ";
                            }

                            //cerr << "wr: ";
                            //working_reg.print_region(cerr);
                            //cerr << " = " << curr_count << ":" << ra[new_node.first]->count << '\n';
                            //ra[curr_node]->set_child(d, f, new_node.first);

                            curr_count = childp->count;
                            child_lphi_unif[f] = childp->get_lphi_unif((depth + top_depth) + 1);
                            if (isinf(childp->lphi_cut)) {
                                child_lphi_cut[f] = child_lphi_unif[f];
                            } else {
                                child_lphi_cut[f] = childp->lphi_cut;
                            }
                            child_lphi_count[f] = childp->count;


                            child_visit_count[f] = childp->num_visits;


                            // uncut for next loop
                            working_reg.uncut(d);


                        } else {

                            pthread_mutex_lock(locker);
                            ls_tree_node_sparse* childp = (*ra)[child];
                            pthread_mutex_unlock(locker);

                            if (isinf(childp->lphi)) {
                                if (childp->count > count_lim
                                        && (depth + top_depth) < max_depth
                                        && depth < lookahead_max_depth) {
                                    computed = false;
                                }
                            }

                            child_lphi_unif[f] = childp->get_lphi_unif((depth + top_depth) + 1);
                            if (isinf(childp->lphi_cut)) {
                                child_lphi_cut[f] = child_lphi_unif[f];
                            } else {
                                child_lphi_cut[f] = childp->lphi_cut;
                            }
                            child_lphi_count[f] = childp->count;


                            child_visit_count[f] = childp->num_visits;
                        }


                    }


                    if (computed)curr_dim_score[d] = -c::inf; // already done
                    else {
                        //curr_dim_score[d] = gt.compute_lD2(curr_node->get_count(),child_count[0],child_count[1]);
                        curr_dim_score[d] =
                                gt->compute_lD2(curr_node.count, child_lphi_count[0], child_lphi_count[1])
                                + logsumexp(child_lphi_cut[0], child_lphi_unif[0])
                                + logsumexp(child_lphi_cut[1], child_lphi_unif[1])
                                - (min(child_visit_count[0]+child_visit_count[1],2)
                                        *(1000.0 / ((depth + top_depth) + 1)));

                        //cerr << "curr_dim_score["<< curr_node.count
                        //        <<" ]("<<child_lphi_count[0] << "/"<<child_lphi_count[1] << ") "
                        //        <<"{" << child_visit_count[0] << "," << child_visit_count[1] << "}"
                        //        << d << ": " << curr_dim_score[d] << '\n';

                        any_good = true;
                    }

                    //cerr << "?curr_dim_score " << d << ": " << curr_dim_score[d] << '\n';
                }

                if (!any_good) {
                    //done = true;
                    //continue;

                    //cerr << "ALLGOOD BACKUP\n";

                    backup = true;
                } else {


                    double total = logsumexp(curr_dim_score);

                    //cerr << "total: " << total << '\n';
                    for (int d = 0; d < num_children; d++) {
                        //curr_dim_score[d] = 1 / (total - curr_dim_score[d]) + (1 / ((double) num_children * sqrt(i + 1)));
                        if (!isinf(curr_dim_score[d])){
                            curr_dim_score[d] = (1 / (1+log(total - curr_dim_score[d] + 1)))
                                + (1 / ((double) num_children));
                        }else{
                            curr_dim_score[d] = 0;
                        }
                        //curr_dim_score[d] = (1 / (total - curr_dim_score[d])) + (1 / ((double) ((depth + start_depth)+1) * sqrt(i + 1)));
                        //cerr << "#curr_dim_score " << d << ": " << curr_dim_score[d] << '\n';
                    }

                    curr_dim = choose_dim(curr_dim_score, rand_gen);
                    pile[depth].dim = curr_dim;

                    //cerr << "choice " << curr_dim << "\n";

                }
            }


            first_cut = true;


            //cerr << "curr_dim: " << curr_dim << ", curr_cut" << curr_cut << '\n';
            //cerr << "reg: ";
            //working_reg.print_region(cerr);
            //cerr << '\n';

            bool is_diff = true;

            if(!backup){

                pile.push_back(pile_t<uint32_t,vector<double> > ());
                depth++;

                is_diff = cut_region(pile[depth - 1].data, pile[depth].data,
                        curr_dim, curr_cut, curr_reg.get_lim(curr_dim));

                if(is_diff) {

                    curr_reg.cut(curr_dim, curr_cut);
                    working_reg.cut(curr_dim, curr_cut);

                    pile[depth].dim = -1;
                    pile[depth].cut = 0;

                    int curr_count = pile[depth].data.size();

                    uint32_t working_hash = region_cache->hash(working_reg);

                    pthread_mutex_lock(locker);
                    pair<uint32_t, bool> new_node = region_cache->find(working_reg, working_hash);



                    pthread_mutex_unlock(locker);

                    if (!new_node.second) {

                        //cerr << "not found: "<< curr_count <<"\n";

                        // MUTEX
                        pthread_mutex_lock(locker);
                        pair<uint32_t, ls_tree_node_sparse*> out = ra->create_node(num_children);
                        new_node.first = out.first;
                        out.second->count = curr_count;

                        region_cache->insert(working_reg, new_node.first, working_hash);


                        (*num_nodes)++;
                        if ((*ra)[new_node.first]->count <= count_lim)(*num_zero_nodes)++;

                        pile[depth].node = new_node.first;

                        // UNMUTEX
                        pthread_mutex_unlock(locker);


                    } else {

#ifdef DEBUG
                        cerr << "found node(" << (*ra)[new_node.first]->count << "): " << curr_dim << "," << curr_cut << '\n';
#endif
                        //ra[curr_node]->set_child(curr_dim, curr_cut, new_node.first);

                        pile[depth].node = new_node.first;
                        //cerr << "fUNCUT: " << curr_dim  << '\n';

                    }

                } else {
                    // all same

                    pile.pop_back();
                    depth--;
                }
            }


            if (backup) {

                //cerr << "BACKUP\n";


                pthread_mutex_lock(locker);
                if((*ra)[curr_node_idx]->count > count_lim)(*ra)[curr_node_idx]->num_visits++;

                // estimate phi


                int curr_count = curr_node.count;
                int curr_depth = (depth + top_depth);

                //cerr << "curr_count: " << curr_count << '\n';


                if (!is_diff
                        || curr_count <= count_lim
                        || depth >= lookahead_max_depth
                        || curr_depth >= max_depth) {

                    //cerr << "UNIFORM\n";

                    //curr_node->set_uniform(curr_depth);
                } else {
                    bool all_calc = true;


                    for (int d = 0; d < num_children; d++) {

                        for (int f = 0; f < 2; f++) {

                            uint32_t child = get_child(working_reg, *region_cache, d, f);
                            //uint32_t child = ra[curr_node]->get_child(d, f);

                            if (child == c::ra_null_val) {
                                // count it
                                working_reg.cut(d, f);
                                uint32_t working_hash = region_cache->hash(working_reg);
                                pair<uint32_t, bool> new_node = region_cache->find(working_reg, working_hash);

                                if (!new_node.second) {

                                    int curr_count = count_region(pile[depth].data, d, f, curr_reg.get_lim(d));

                                    // MUTEX

                                    pair<uint32_t, ls_tree_node_sparse*> out = ra->create_node();

                                    new_node.first = out.first;
                                    out.second->count = curr_count;
                                    (*num_nodes)++;
                                    if (curr_count <= count_lim)(*num_zero_nodes)++;

                                    region_cache->insert(working_reg, new_node.first, working_hash);

                                    //UNMUTEX


                                }
                                working_reg.uncut(d);

                                all_calc = false;

                            } else {

                                if (isinf((*ra)[child]->lphi) && (*ra)[child]->count > top_count_lim) {
                                    all_calc = false;
                                }
                            }

                        }


                    }



                    //cerr << "-?-\nDEPTH: " << (depth + start_depth) << '\n';
                    //cerr << "curr_dim: " << curr_dim << ", curr_cut" << curr_cut << '\n';
                    //cerr << "reg: ";
                    //working_reg.print_region(cerr);
                    //cerr << '\n';

                    if (all_calc) {
                        //ra[curr_node]->compute_lphi(curr_depth, gt,ra);

                        compute_lphi(*ra,working_reg, *region_cache, curr_node_idx, curr_depth, *gt,
                                1, pile[depth].data, curr_reg,num_children);


                        //cerr << "LPHI: " << ra[curr_node]->lphi <<"\n";

                    } else if (curr_depth > 0) {
                        //ra[curr_node]->compute_approx_lphi_cut(curr_depth, gt,ra);

                        compute_approx_lphi_cut(*ra,working_reg, *region_cache,
                                curr_node_idx, curr_depth, *gt, 1, pile[depth].data, curr_reg,num_children);

                        //cerr << "LPHI_CUT: " << ra[curr_node]->lphi_cut << " + "
                        //        << ra[curr_node]->get_lphi_unif(curr_depth) << " = "
                        //        << logsumexp(ra[curr_node]->lphi_cut,ra[curr_node]->get_lphi_unif(curr_depth)) <<"\n";

                    }
                }

                pthread_mutex_unlock(locker);

                depth--;
                pile.pop_back();
                if (depth < 0) continue;

                curr_reg.uncut(pile[depth].dim, pile[depth].cut);
                working_reg.uncut(pile[depth].dim);

                continue;
            }


        }

        return NULL;
    }


    void do_iteration(ls_working_unit_t &w,
            opt_region_hash<uint32_t> &region_cache,gamma_table &gt,
            MT_random &rand_gen,
            int64_t &num_nodes, int64_t &num_zero_nodes, int start_depth,
            int top_count_lim) {

        int count_lim = (int)floor(w.data.size()*count_ratio);
        if(count_lim < top_count_lim) count_lim = top_count_lim;


        pthread_t* t_group = new pthread_t[num_children]; // just do 3 threads for now :)
        pthread_mutex_t* locker = new pthread_mutex_t();
        pthread_mutex_init(locker, NULL);

        mouse_params_t* params = new mouse_params_t[num_children];

        for (int d = 0; d < num_children; d++) {
            params[d].count_lim = count_lim;
            params[d].top_count_lim = top_count_lim;
            params[d].gt = &gt;
            params[d].locker = locker;
            params[d].max_depth = max_depth;
            params[d].lookahead_max_depth = lookahead_max_depth;
            params[d].num_children = num_children;
            params[d].num_nodes = &num_nodes;
            params[d].num_zero_nodes = &num_zero_nodes;
            params[d].ra = &ra;
            params[d].region_cache = &region_cache;
            params[d].start_dimension = d;
            params[d].top_depth = start_depth;
            params[d].wup = &w;
        }

        for (int d = 0; d < num_children; d++) {
            // Start threads


            pthread_create(&(t_group[d]), NULL, iteration_thread, (void*) &(params[d]));
        }

        // wait for threads
        for (int d = 0; d < num_children; d++) {
            void* output;
            pthread_join(t_group[d], &output);
        }




        delete locker;
        delete [] params;
        delete [] t_group;
    }

int get_map_dim(ls_working_unit_t &w,opt_region_hash<uint32_t> &region_cache,gamma_table &gt,
                    int map_depth,bool debug) {
        // Choose MAP dimension

        int map_dim = 0;

        vector<uint32_t>child_idx_0(num_children);
        vector<uint32_t>child_idx_1(num_children);

        for (int i = 0; i < num_children; i++) {

            if(w.working_reg.full()){
                w.working_reg.print_region();
                cerr << "\n";
            }

            child_idx_0[i] = get_child(w.working_reg, region_cache, i, 0,debug);
            child_idx_1[i] = get_child(w.working_reg, region_cache, i, 1,debug);

            if(debug){
                w.working_reg.print_region();
                cerr << ":" << child_idx_0[i] << "," << child_idx_1[i] << '\n';

            }
        }

        if(debug){
            cerr << "w.node_idx: " << w.node_idx << '\n';
        }

        double max_post_prob = -c::inf;

        for (int i = 0; i < num_children; i++) {
            double post_prob = gt.compute_lD2(ra[w.node_idx]->count
                    , ra[child_idx_0[i]]->count, ra[child_idx_1[i]]->count);
            post_prob += ra[child_idx_0[i]]->get_lphi3(map_depth + 1);
            post_prob += ra[child_idx_1[i]]->get_lphi3(map_depth + 1);
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
        cerr << map_dim << ',';
        //#endif

        return map_dim;

    }

    void construct_lsopt_tree(vector<vector<double> > &data, int iterations,
            int convergence_iterations,
            map_tree &map_region_tree, opt_region_hash<uint32_t> &map_regions) {

        MT_random rand_gen;

        // need to treat the empty tree/data case

        int64_t num_nodes = 0;
        int64_t num_zero_nodes = 0;

        int N = data.size();

        int count_lim = (int)floor((double)N*top_count_ratio); // this is for the whole tree
        if(count_lim < 1)count_lim = 1;

        cerr << "Stopping at " << count_lim << " points\n";

        gamma_table gt(N);
        ra[root]->count = N;

        region_allocator<map_tree_node> *map_ra = map_region_tree.get_ra();

        opt_region_hash<uint32_t> region_cache(25);


        current_region start_region(num_children);
        opt_region     start_working(num_children);





        vector<ls_pile_t> lspile;
        lspile.push_back(ls_pile_t());

        lspile[0].good_regions.push_back(ls_working_unit_t(data,start_region,
                start_working,root));

        lspile[0].map_nodes.push_back(map_region_tree.get_full_tree());

        int map_depth = 0;



        // Now create the MAP tree!
        // Do a breadth first search on the MAP tree
        bool done = false;
        int count = 1;
        while(!done){


            // if no good regions we are done!
            if(lspile.size() == 0){
                done = true;
                continue;
            }

            if(lspile[map_depth].good_regions.size() == 0){
                lspile.pop_back();
                map_depth--;
                continue;
            }


            vector<ls_working_unit_t>::iterator wu_it = (lspile[map_depth].good_regions.end() - 1);
            vector<uint32_t>::iterator map_node_it = (lspile[map_depth].map_nodes.end() - 1);



            // if too deep we stop all the regions
            if (map_depth >= max_depth || (int)wu_it->data.size() <= count_lim
                    || wu_it->working_reg.full()) {

                cerr << "mini-STOP count: " << ra[wu_it->node_idx]->count << '\n';

                store_region(*map_node_it, map_regions, map_ra, wu_it->working_reg,
                        ra[wu_it->node_idx]->count, map_depth);

                lspile[map_depth].good_regions.pop_back();
                lspile[map_depth].map_nodes.pop_back();
                count++;

                continue;

            }


            // for each good region
            int curr_data_size = wu_it->data.size();

            cerr << "Good [" << map_depth << "](" << count << "/"
                    << lspile[map_depth].good_regions.size()
                    << "), Data size: " << curr_data_size << '\n';


            // iterate through each dimension
            // filling up the tree
            // Iterate until convergence

            int convergence_count = 0;
            int prev_map_dim = -1;

            for (int i = 0; i < iterations; i++) {
                do_iteration(*wu_it,
                        region_cache, gt, rand_gen,
                        num_nodes, num_zero_nodes, map_depth, count_lim);

                int temp_map_dim = get_map_dim(*wu_it, region_cache, gt, map_depth,false);

                if (i > 10) {
                    if (temp_map_dim != prev_map_dim) {
                        convergence_count = 0;
                        prev_map_dim = temp_map_dim;
                    } else {
                        convergence_count++;
                    }
                    if (convergence_count >= convergence_iterations) {
                        break;
                    }

                }


                if (i % 10 == 0 || i == iterations - 1) {
                    cerr << "\nIteration: " << i << ", ";
                    compute_approx_lphi_cut(ra,wu_it->working_reg, region_cache,
                            wu_it->node_idx, map_depth, gt, 2, wu_it->data, wu_it->curr_reg,num_children);

                    cerr << " lphi: " << ra[wu_it->node_idx]->get_lphi3(map_depth) << "\n";


                }

            }
            cerr << "\nNodes:" << num_nodes
                    << " : " << num_zero_nodes
                    << " : " << (num_nodes - num_zero_nodes) << '\n';

            // propagate lphis through current tree
            // compute lphi of current good node

            //cerr << "Propagate lphi\n";
            //propagate_lphi(*wit,region_cache,gt,map_depth);
            //ra[wit->node_idx]->compute_approx_lphi_cut(map_depth,gt,ra);



            compute_approx_lphi_cut(ra,wu_it->working_reg, region_cache,
                    wu_it->node_idx, map_depth, gt, 3, wu_it->data, wu_it->curr_reg,num_children);

            cerr << "lphi: " << ra[wu_it->node_idx]->get_lphi3(map_depth) << "\n";


            // compute the posterior pho and decide if we stop
            double post_rho = -c::l2;
            post_rho += map_depth * c::l2 * ra[wu_it->node_idx]->count; // phi_0
            post_rho -= ra[wu_it->node_idx]->get_lphi3(map_depth);

            //cerr << "---blah: " << ra[wit->node_idx]->get_lphi3(map_depth)
            //        << ',' << ra[wit->node_idx]->lphi_cut << '\n';
            //cerr << "phi_0: " << map_depth * c::l2 * ra[wit->node_idx]->count << '\n';
            //cerr << "post_rho: "<< post_rho << "|" << -c::l2 << "\n";


            // if not stop we split and add the two regions to the good region list
            if (post_rho>-c::l2) {
                // STOP

                cerr << "STOP count: " << ra[wu_it->node_idx]->count << '\n';

                store_region(*map_node_it, map_regions, map_ra, wu_it->working_reg,
                        ra[wu_it->node_idx]->count, map_depth);

                lspile[map_depth].good_regions.pop_back();
                lspile[map_depth].map_nodes.pop_back();
                count++;

            } else {
                // Choose MAP dimension

                cerr << "MAP dim: ";
                int map_dim = get_map_dim(*wu_it, region_cache, gt, map_depth,false);
                cerr << '\n';

                //cerr << "DONE MAP\n";


                uint32_t map_child_idx_0 = get_child(wu_it->working_reg, region_cache, map_dim, 0);
                uint32_t map_child_idx_1 = get_child(wu_it->working_reg, region_cache, map_dim, 1);

                //cerr << "children: " << map_child_idx_0 << "," << map_child_idx_1 << '\n';

                // Split
                vector<vector<double> > new_data_0;
                vector<vector<double> > new_data_1;

                bool is_diff = cut_region2(wu_it->data,new_data_0,new_data_1,
                        map_dim, wu_it->curr_reg.get_lim(map_dim));

                if (is_diff) {

                    current_region next_curr_reg0 = wu_it->curr_reg;
                    opt_region next_working_reg0 = wu_it->working_reg;

                    next_curr_reg0.cut(map_dim, 0);
                    next_working_reg0.cut(map_dim, 0);

                    current_region next_curr_reg1 = wu_it->curr_reg;
                    opt_region next_working_reg1 = wu_it->working_reg;

                    next_curr_reg1.cut(map_dim, 1);
                    next_working_reg1.cut(map_dim, 1);


                    // after this point wu_it is invalid
                    lspile.push_back(ls_pile_t());
                    map_depth++;

                    map_node_it = (lspile[map_depth - 1].map_nodes.end() - 1);

                    //cerr << "ADD PILE " << '\n';

                    lspile[map_depth].good_regions.push_back(ls_working_unit_t(new_data_0,
                            next_curr_reg0, next_working_reg0, map_child_idx_0));

                    lspile[map_depth].good_regions.push_back(ls_working_unit_t(new_data_1,
                            next_curr_reg1, next_working_reg1, map_child_idx_1));



                    pair<uint32_t, map_tree_node*> new_map_node0 = map_ra->create_node();
                    pair<uint32_t, map_tree_node*> new_map_node1 = map_ra->create_node();

                    //cerr << "children: " << new_map_node0.first
                    //        << "," << new_map_node1.first << '\n';

                    //cerr << "*map_node_it: " << *map_node_it << '\n';

                    (*map_ra)[(*map_node_it)]->set_dim(map_dim);
                    (*map_ra)[(*map_node_it)]->set_child(0, new_map_node0.first);
                    (*map_ra)[(*map_node_it)]->set_child(1, new_map_node1.first);

                    lspile[map_depth].map_nodes.push_back(new_map_node0.first);
                    lspile[map_depth].map_nodes.push_back(new_map_node1.first);

                    lspile[map_depth-1].good_regions.pop_back();
                    lspile[map_depth-1].map_nodes.pop_back();

                }else{
                    store_region(*map_node_it, map_regions, map_ra, wu_it->working_reg,
                            ra[wu_it->node_idx]->count, map_depth);

                    lspile[map_depth].good_regions.pop_back();
                    lspile[map_depth].map_nodes.pop_back();

                    count++;
                }


            }


            // remove all the nodes not consistent with current good regions

            if((curr_data_size > 100 && ra.free_locs.size() < 1000000) || (curr_data_size > 1000)) {
                int64_t num_removed = 0;
                int64_t total_nodes = 0;

                for (uint32_t i = 0; i < region_cache.table_size; i++) {
                    if (region_cache.map_table[i] != NULL) {

                        map<opt_region, uint32_t> *new_map = new map<opt_region, uint32_t > ();

                        for (map<opt_region, uint32_t>::iterator it = region_cache.map_table[i]->begin();
                                it != region_cache.map_table[i]->end(); it++) {

                            bool is_child = false;

                            for (vector<ls_pile_t>::iterator ls_it = lspile.begin();
                                    ls_it != lspile.end() && !is_child; ls_it++) {
                                for (vector<ls_working_unit_t>::iterator reg_it = ls_it->good_regions.begin();
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

                if(total_nodes>0)cerr << "Removed " << num_removed << "/" << total_nodes
                        << "(" << (num_removed / (total_nodes / 100.0)) << "%)" << '\n';
            }

        }


    }


};

#endif	/* LSOPT_TREE_H */

