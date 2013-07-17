/*
 *  opt_tree_sample.h
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



#ifndef OPT_TREE_SAMPLE_H
#define	OPT_TREE_SAMPLE_H


#include "general_utils.h"
#include "gamma_table.h"
#include "opt_utils.h"

class tree_node_sample {
private:
    int count; // number of points in this region, this is also the dimension
    // to cut in the MAP tree.
    tree_node_sample*** children;





    int num_children;
    double lphi; // this is also the density if MAP tree





public:

    //double* dim_scores; // these are the scores from the counts alone

    //double lphi_unif;
    double lphi_cut;

    int num_visits;

    void init() {
        num_visits = 0;
        count = -1;
        lphi = -c::inf;
        //lphi_unif = -c::inf;
        lphi_cut = -c::inf;
    }

    tree_node_sample() {
        children = NULL;
        num_children = 0;

        init();

    }

    tree_node_sample(int num_children) {
        this->num_children = num_children;

        children = new tree_node_sample**[num_children];
        //dim_scores = new double[num_children];

        for (int i = 0; i < num_children; i++) {

            //dim_scores[i] = -1.0;

            children[i] = new tree_node_sample*[c::cuts];
            for (int j = 0; j < c::cuts; j++) {
                children[i][j] = NULL;
            }
        }





        init();
    }

    ~tree_node_sample() {

        if (children == NULL) return;

        for (int i = 0; i < num_children; i++) {
            delete [] children[i];
        }
        delete [] children;

    }

    void set_uniform(int depth) {
        lphi = count * depth * c::l2;
        //lphi_unif = lphi;
        lphi_cut = -1.0e50; // a very very small number :)
    }

    void compute_lphi(int depth, gamma_table& gt) {
        vector<double> lphi_list;
        double max_val = (count * depth * c::l2) - c::l2;
        lphi_list.push_back(max_val);

        //cerr << "max_val: " << max_val << '\n';

        double ld = -log(num_children) - c::lpi - c::l2;

        //cerr << "ld: " << ld << '\n';

        for (int i = 0; i < num_children; i++) {
            // check for null
            if (children[i][0] == NULL || children[i][1] == NULL) {
                cerr << "NULL child!!! " << i << ',' << depth << '\n';
                exit(2);
            }

            int child_1_count = children[i][0]->get_count();
            int child_2_count = children[i][1]->get_count();

            if (child_1_count < 0 || child_2_count < 0) {
                cerr << "neg count!!! " << i << ',' << depth << '\n';
                exit(2);
            }


            //cerr << i << ":" << child_1_count << ',' << child_2_count << '\n';

            double val = ld;
            val += children[i][0]->get_lphi2(depth + 1);
            val += children[i][1]->get_lphi2(depth + 1);
            val += gt.compute_lD2(count, child_1_count, child_2_count);

            //cerr << "val: " << val << '\n';

            lphi_list.push_back(val);

            if (val > max_val) {
                max_val = val;
            }

        }

        lphi = max_val;
        double sum = 0;
        for (int i = 0; i < (num_children + 1); i++) {
            sum += exp(lphi_list[i] - max_val);
        }

        if (sum > 0) lphi += log(sum);



        //if(isnan(lphi)) cerr <<"nan sum:" << sum << "," << max_val << '\n';
        //if(isinf(lphi)) cerr <<"inf sum:" << sum << "," << max_val << '\n';


    }

    void compute_approx_lphi_cut(int depth, gamma_table& gt) {
        vector<double> lphi_list;
        double max_val = -c::inf;

        double ld = -log(num_children) - c::lpi - c::l2;

        //cerr << "ld: " << ld << '\n';

        for (int i = 0; i < num_children; i++) {
            // check for null
            if (children[i][0] == NULL || children[i][1] == NULL) {
                cerr << "NULL child!!! " << i << ',' << depth << '\n';
                exit(2);
            }

            int child_1_count = children[i][0]->get_count();
            int child_2_count = children[i][1]->get_count();

            if (child_1_count < 0 || child_2_count < 0) {
                cerr << "neg count!!! " << i << ',' << depth << '\n';
                exit(2);
            }


            //cerr << i << ":" << child_1_count << ',' << child_2_count << '\n';

            double val = ld;

            for (int x = 0; x < 2; x++) {
                if (isinf(children[i][x]->get_lphi())) {
                    if (isinf(children[i][x]->lphi_cut)) {
                        val += children[i][x]->get_lphi_unif(depth+1) + c::l2;
                    } else {
                        // this is approximate

                        val += logsumexp(children[i][x]->lphi_cut, children[i][x]->get_lphi_unif(depth+1));
                    }
                } else {
                    val += children[i][x]->get_lphi();
                }
            }

            val += gt.compute_lD2(count, child_1_count, child_2_count);

            //cerr << "val: " << val << '\n';

            lphi_list.push_back(val);

            if (val > max_val) {
                max_val = val;
            }

        }

        lphi_cut = max_val;
        double sum = 0;
        for (int i = 0; i < (num_children); i++) {
            sum += exp(lphi_list[i] - max_val);
        }

        if (sum > 0) lphi_cut += log(sum);



        //if(isnan(lphi)) cerr <<"nan sum:" << sum << "," << max_val << '\n';
        //if(isinf(lphi)) cerr <<"inf sum:" << sum << "," << max_val << '\n';


    }

    double get_lphi() {
        return lphi;
    }

    double get_lphi2(int depth) {
        if(!isinf(lphi))return lphi;
        else return get_lphi_unif(depth)+c::l2;
    }

    int get_num_children() {
        return num_children;
    }

    int get_cuts() {
        return c::cuts;
    }

    void set_lphi(double lphi) {
        this->lphi = lphi;
    }

    //void set_count(int count, int depth) {
    //    this->count = count;
        //this->lphi_unif = (count * depth * c::l2) - c::l2;
    //}

    double get_lphi_unif(int depth){
        if(count >= 0){
            return (count * depth * c::l2) - c::l2;
        }else{
            return -c::inf;
        }
    }

    void set_count(int count) {
        this->count = count;
    }

    int get_count() {
        return count;
    }

    void set_child(int dim, int cut, tree_node_sample* node) {
        if (children != NULL)children[dim][cut] = node;
    }

    tree_node_sample* get_child(int dim, int cut) {
        if (children == NULL) return NULL;
        return children[dim][cut];
    }

};



class opt_tree_sample {
private:
    tree_node_sample* root;
    int num_children;

    int count_lim;
    int max_depth;

    void init(int num_children, int count_lim, int max_depth) {
        this->num_children = num_children;
        root = new tree_node_sample(num_children);
        this->count_lim = count_lim;
        this->max_depth = max_depth;



    }

public:

    opt_tree_sample(int num_children, int count_lim, int max_depth) {
        init(num_children, count_lim, max_depth);
    }

    opt_tree_sample(int num_children) {

        init(num_children, 10, 100);

    }

    ~opt_tree_sample() {
        delete root;
    }

    void construct_sampled_tree(vector<vector<double> > &data, int iterations) {



        MT_random rand_gen;

        // need to treat the empty tree/data case


        int64_t num_nodes = 0;
        int64_t num_zero_nodes = 0;

        int N = data.size();

        gamma_table gt(N);
        root->set_count(N);

        opt_region_hash<tree_node_sample*> region_cache(24);


        // create a bunch of leaf nodes
        tree_node_sample** leaves = new tree_node_sample*[N+1];
        for(int i = 0;i<N+1;i++){
            leaves[i] = new tree_node_sample();
            leaves[i]->set_count(i);
        }

        // create paths in tree
        for (int i = 0; i < iterations; i++) {
            vector<pile_t<tree_node_sample*> > pile;
            pile.push_back(pile_t<tree_node_sample*>());

            pile[0].node = root;
            pile[0].dim = -1;
            pile[0].cut = 0;
            pile[0].data = data;

            current_region curr_reg(num_children);
            opt_region working_reg(num_children);

            vector<double> curr_dim_score(num_children, -1.0);

            bool done = false;
            int depth = 0;

            // iterate with a binary tree until bottom

            //cerr << "----" << i << "\n";

            while (!done) {

                if (pile.size() == 0) {
                    done = true;
                    continue;
                }

                int curr_dim = pile[depth].dim;
                int curr_cut = pile[depth].cut;
                tree_node_sample* curr_node = pile[depth].node;
                bool backup = false;



                //cerr << "----\nDEPTH: " << depth << '\n';



                if (curr_node->get_count() <= count_lim || depth >= max_depth) {

                    //cerr << "BOTTOM BACKUP\n";

                    backup = true;
                }else if(curr_dim >= 0){
                    if(curr_cut < c::cuts - 1){
                        curr_cut++;
                        pile[depth].cut++;
                    }else{

                        //cerr << "CUTALL BACKUP\n";

                        backup = true;
                    }
                }

                // work out which dimension to cut

                if (curr_dim < 0) { // if we have not chosen a dimension yet.

                    bool any_good = false;
                    for (int d = 0; !backup && (d < num_children); d++) {
                        double child_lphi_unif[2] = {-1.0, -1.0};
                        double child_lphi_cut[2] = {-1.0, -1.0};
                        int child_lphi_count[2] = {-1, -1};
                        int child_visit_count[2] = {0, 0};

                        bool computed = true;

                        for (int f = 0; f < 2; f++) {

                            tree_node_sample* child = curr_node->get_child(d, f);

                            if (child == NULL) {
                                // count it
                                working_reg.cut(d, f);
                                uint32_t working_hash = region_cache.hash(working_reg);
                                pair<tree_node_sample*,bool> new_node = region_cache.find(working_reg, working_hash);
                                working_reg.uncut(d);

                                int curr_count = 0;
                                if (!new_node.second) {

                                    curr_count = count_region(pile[depth].data, d, f, curr_reg.get_lim(d));

                                    if (curr_count <= count_lim || depth >= max_depth) {
                                        //new_node = new tree_node_sample();
                                        //new_node->set_count(curr_count, depth);

                                        new_node.first = leaves[curr_count];

                                    } else {
                                        new_node.first = new tree_node_sample(num_children);
                                        new_node.first->set_count(curr_count);
                                        computed = false;
                                    }
                                } else {
                                    curr_count = new_node.first->get_count();
                                }

                                curr_node->set_child(d, f, new_node.first);
                                child_lphi_unif[f] = new_node.first->get_lphi_unif(depth+1);
                                if (isinf(new_node.first->lphi_cut)) {
                                    child_lphi_cut[f] = child_lphi_unif[f];
                                } else {
                                    child_lphi_cut[f] = new_node.first->lphi_cut;
                                }
                                child_lphi_count[f] = new_node.first->get_count();


                                child_visit_count[f] = new_node.first->num_visits;

                            } else {

                                if (child->get_count() > count_lim && depth < max_depth && isinf(child->get_lphi())) {
                                    computed = false;
                                }

                                child_lphi_unif[f] = child->get_lphi_unif(depth+1);
                                if (isinf(child->lphi_cut)) {
                                    child_lphi_cut[f] = child_lphi_unif[f];
                                } else {
                                    child_lphi_cut[f] = child->lphi_cut;
                                }
                                child_lphi_count[f] = child->get_count();

                                child_visit_count[f] = child->num_visits;
                            }


                        }


                        if (computed)curr_dim_score[d] = 0.0; // already done
                        else {
                            //curr_dim_score[d] = gt.compute_lD2(curr_node->get_count(),child_count[0],child_count[1]);
                            curr_dim_score[d] =
                                    gt.compute_lD2(curr_node->get_count(), child_lphi_count[0], child_lphi_count[1])
                                    + logsumexp(child_lphi_cut[0], child_lphi_unif[0])
                                    + logsumexp(child_lphi_cut[1], child_lphi_unif[1])
                                    - (child_visit_count[0]*(10.0/(depth+1)))
                                    - (child_visit_count[1]*(10.0/(depth+1)));

                            //cerr << "curr_dim_score["<< curr_node->get_count()
                            //        <<" ]("<<child_lphi_count[0] << "/"<<child_lphi_count[1] << ") "
                            //        << d << ": " << curr_dim_score[d] << '\n';

                            any_good = true;
                        }
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
                            curr_dim_score[d] = (1 / sqrt(total - curr_dim_score[d] + 1)) + (1 / ((double) num_children * sqrt(i + 1)));
                            //curr_dim_score[d] = (1 / (total - curr_dim_score[d])) + (1 / ((double) (depth+1) * sqrt(i + 1)));
                            //cerr << "#curr_dim_score " << d << ": " << curr_dim_score[d] << '\n';
                        }

                        curr_dim = choose_dim(curr_dim_score, rand_gen);
                        pile[depth].dim = curr_dim;
                    }
                }



                //cerr << "curr_dim: " << curr_dim << ", curr_cut" << curr_cut << '\n';
                //cerr << "reg: ";
                //working_reg.print_region(cerr);
                //cerr << '\n';

                if(backup){

                    //cerr << "BACKUP\n";

                    curr_node->num_visits++;

                    // estimate phi

                    int curr_count = curr_node->get_count();
                    int curr_depth = depth;

                    if (curr_count <= count_lim || curr_depth >= max_depth) {

                        //cerr << "UNIFORM\n";

                        //curr_node->set_uniform(curr_depth);
                    } else {
                        bool all_calc = true;

                        for (int d = 0; d < num_children; d++) {
                            if (curr_node->get_child(d, 0) == NULL) {
                                all_calc = false;
                                break;
                            } else if (isinf(curr_node->get_child(d, 0)->get_lphi())) {
                                all_calc = false;
                                break;
                            }
                            if (curr_node->get_child(d, 1) == NULL) {
                                all_calc = false;
                                break;
                            } else if (isinf(curr_node->get_child(d, 1)->get_lphi())) {
                                all_calc = false;
                                break;
                            }
                        }

                        //cerr << "-??-\nDEPTH: " << depth << '\n';
                        //cerr << "curr_dim: " << curr_dim << ", curr_cut" << curr_cut << '\n';
                        //cerr << "reg: ";
                        //working_reg.print_region(cerr);
                        //cerr << '\n';

                        if (all_calc) {
                            curr_node->compute_lphi(curr_depth, gt);


                            //cerr << "LPHI: " << curr_node->get_lphi() <<"\n";

                        } else {
                            curr_node->compute_approx_lphi_cut(curr_depth, gt);

                            //cerr << "LPHI_CUT: " << curr_node->lphi_cut << " + "
                            //        << curr_node->get_lphi_unif(curr_depth) << " = "
                            //        << logsumexp(curr_node->lphi_cut,curr_node->get_lphi_unif(curr_depth)) <<"\n";

                        }
                    }

                    depth--;
                    pile.pop_back();
                    if (depth < 0) continue;

                    curr_reg.uncut(pile[depth].dim, pile[depth].cut);
                    working_reg.uncut(pile[depth].dim);

                    continue;
                }


                /*
                // work out which region of the dimension to explore
                tree_node_sample* child_0 = curr_node->get_child(curr_dim, 0);
                tree_node_sample* child_1 = curr_node->get_child(curr_dim, 1);

                if (child_0 != NULL && !isinf(child_0->get_lphi())) {
                    curr_cut = 1;
                } else if (child_1 != NULL && !isinf(child_1->get_lphi())) {
                    curr_cut = 0;
                } else {
                    //curr_cut = rand_gen.genrand_bern(child_1->get_count()/((double)child_0->get_count()+child_1->get_count()));
                    double rho_0 = 0.5;
                    if (!isinf(child_0->lphi_cut)) {

                        double cut_num = child_0->lphi_cut - (child_0->num_visits);



                        rho_0 = exp(cut_num - logsumexp(cut_num,child_0->lphi_unif));
                    }
                    double rho_1 = 0.5;
                    if (!isinf(child_1->lphi_cut)) {

                        double cut_num = child_1->lphi_cut - (child_1->num_visits);

                        rho_1 = exp(cut_num - logsumexp(cut_num,child_1->lphi_unif));
                    }

                    curr_cut = rand_gen.genrand_bern(rho_1 / (rho_0 + rho_1));
                }
                 */

                // do the cutting



                pile.push_back(pile_t<tree_node_sample*>());
                depth++;

                cut_region(pile[depth-1].data,pile[depth].data,
                        curr_dim, curr_cut, curr_reg.get_lim(curr_dim));

                curr_reg.cut(curr_dim, curr_cut);


                pile[depth].dim = -1;
                pile[depth].cut = 0;

                int curr_count = pile[depth].data.size();


                working_reg.cut(curr_dim, curr_cut);
                uint32_t working_hash = region_cache.hash(working_reg);
#ifdef DEBUG
                cerr << "working_hash: " << working_hash << '\n';
#endif
                pair<tree_node_sample*,bool> new_node = region_cache.find(working_reg, working_hash);

                if (!new_node.second) {

                    if (curr_count <= count_lim || depth >= max_depth) {
                        //new_node = new tree_node_sample();
                        //new_node->set_count(curr_count, depth);

                        new_node.first = leaves[curr_count];
                        //new_node->set_uniform(depth);

                        //done = true;
                    } else {


                        new_node.first = new tree_node_sample(num_children);
                        new_node.first->set_count(curr_count);

                        region_cache.insert(working_reg, new_node.first, working_hash);
                    }




                    num_nodes++;
                    if (new_node.first->get_count() <= count_lim)num_zero_nodes++;


                    curr_node->set_child(curr_dim, curr_cut, new_node.first);

                    pile[depth].node = new_node.first;



                    if (num_nodes % 100000 == 0) {
                        cerr << "Nodes(" << pile[0].dim << "):" << num_nodes
                                << " : " << num_zero_nodes
                                << " : " << (num_nodes-num_zero_nodes) << '\n';
                    }

                } else {

#ifdef DEBUG
                    cerr << "found node: " << curr_dim << "," << curr_cut << '\n';
#endif
                    curr_node->set_child(curr_dim, curr_cut, new_node.first);
                    pile[depth].node = new_node.first;
                    //cerr << "fUNCUT: " << curr_dim  << '\n';

                }





            }

            // go back up the pile and update region scores as well as compute
            // lphi when possible

            /*

            for (int x = (int) pile.size() - 1; x >= 0; x--) {
                tree_node_sample* curr_node = pile[x].node;
                int curr_count = curr_node->get_count();
                int curr_depth = x;

                curr_node->num_visits++;

                if (curr_count <= count_lim || curr_depth >= max_depth) {
                    curr_node->set_uniform(curr_depth);
                } else {
                    bool all_calc = true;

                    for (int d = 0; d < num_children; d++) {
                        if (curr_node->get_child(d, 0) == NULL) {
                            all_calc = false;
                            break;
                        } else if (isinf(curr_node->get_child(d, 0)->get_lphi())) {
                            all_calc = false;
                            break;
                        }
                        if (curr_node->get_child(d, 1) == NULL) {
                            all_calc = false;
                            break;
                        } else if (isinf(curr_node->get_child(d, 1)->get_lphi())) {
                            all_calc = false;
                            break;
                        }
                    }

                    if (all_calc) {
                        curr_node->compute_lphi(curr_depth, gt);
                    } else {
                        curr_node->compute_approx_lphi_cut(curr_depth, gt);
                    }
                }
            }
             */

            if (i % 100 == 0)cerr << "Iteration(" << root->lphi_cut << "): " << i << '\n';

        }


        cerr << "Done Iteration(" << root->lphi_cut << ") " << '\n';



        // fill out NULL nodes as uniform
        // compute lphi at same time





        vector<pile_t<tree_node_sample*> > pile;
        pile.push_back(pile_t<tree_node_sample*>());

        pile[0].data = data;
        pile[0].node = root;
        pile[0].dim = 0;
        pile[0].cut = 0;

        current_region curr_reg(num_children);
        opt_region working_reg(num_children);


        int depth = 0;
        //double lvolume = 0;

        bool done = false;


        while (!done) {

            if (pile.size() == 0) {
                done = true;
                continue;
            }

#ifdef DEBUG3
            // print the pile
            cerr << "@@@---- " << depth << " , " << pile.size() << " ---\n";
            cerr << "REG: ";
            working_reg.print_region(cerr);
            cerr << '\n';
            cerr << num_nodes << " : " << num_zero_nodes << '\n';
            for (int i = 0; i < (int) pile.size(); i++) {
                cerr << "Data size: " << pile[i].data.size() << '\n';
                cerr << "Dim:       " << pile[i].dim << '\n';
                cerr << "Cut:       " << pile[i].cut << '\n';
                if (pile[i].node->get_child(pile[i].dim, pile[i].cut) != NULL) {
                    cerr << "Child lphi: " << pile[i].node->get_child(pile[i].dim, pile[i].cut)->get_lphi() << '\n';
                }
                //if(pile[i].data.size()<20){
                //    print_data(pile[i].data);
                //}
            }
#endif



            // this is only kind of temporary

            int curr_dim = pile[depth].dim;
            int curr_cut = pile[depth].cut;
            tree_node_sample* curr_node = pile[depth].node;

            // work out what to count

            bool back_up = false;
#ifdef DEBUG3
            if (curr_node->get_child(curr_dim, curr_cut) != NULL) {
                cerr << "@Child lphi: " << curr_node->get_child(curr_dim, curr_cut)->get_lphi() << '\n';
            } else {
                cerr << "@NULL node " << curr_dim << "," << curr_cut << " \n";
            }
#endif

            // check if current node is leaf or at end
            if (curr_node->get_num_children() == 0 || curr_node->get_count() <= count_lim || depth >= max_depth) {
                // back up
                back_up = true;
                //curr_node->set_uniform(depth);
#ifdef DEBUG3
                cerr << "BACKUP ZERO\n";
#endif

            } else if (curr_node->get_child(curr_dim, curr_cut) != NULL
                    && (!isinf(curr_node->get_child(curr_dim, curr_cut)->get_lphi())
                    || curr_node->get_child(curr_dim, curr_cut)->get_num_children() == 0)) {
                // move to next node
#ifdef DEBUG3
                cerr << "CHILD FOUND\n";
                cerr << "pile[depth].cut: " << pile[depth].cut << '\n';
                cerr << "pile[depth].dim: " << pile[depth].dim << '\n';
#endif

                if (pile[depth].cut < c::cuts - 1) {
                    pile[depth].cut++;
                } else if (pile[depth].dim < num_children - 1) {
                    pile[depth].dim++;
                    pile[depth].cut = 0;
                } else {
                    // reached end of node!! back up
                    back_up = true;
                    curr_node->compute_lphi(depth, gt);
#ifdef DEBUG3
                    cerr << "LPHI---- " << depth << " , " << pile.size() << " ---\n";
                    cerr << "REG: ";
                    working_reg.print_region(cerr);
                    cerr << '\n';
                    cerr << num_nodes << " : " << num_zero_nodes << '\n';
                    cerr << "lphi: " << curr_node->get_lphi() << '\n';
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
            }

            if (back_up) {

                depth--;
                pile.pop_back();
                if (depth < 0) continue;



                //lvolume = lvolume + c::l2;


                curr_reg.uncut(pile[depth].dim, pile[depth].cut);
                //cerr << depth << " bUNCUT: " << curr_dim << '\n';
                working_reg.uncut(pile[depth].dim);

                continue;
            }


            curr_dim = pile[depth].dim;
            curr_cut = pile[depth].cut;

            // do the counting


            //cerr << "CUT: " << curr_dim << "," << curr_cut << '\n';
            working_reg.cut(curr_dim, curr_cut);
            uint32_t working_hash = region_cache.hash(working_reg);
#ifdef DEBUG
            cerr << "working_hash: " << working_hash << '\n';
#endif
            pair<tree_node_sample*,bool> new_node = region_cache.find(working_reg, working_hash);

            if (!new_node.second) {

#ifdef DEBUG3
                cerr << "Add node: " << curr_dim << "," << curr_cut << '\n';
#endif

                //pile.push_back(pile_sample_t());
                //depth++;

                int curr_count = count_region(pile[depth].data,
                        curr_dim, curr_cut, curr_reg.get_lim(curr_dim));

#ifdef DEBUG
                cerr << "ADD NODE curr_count: " << curr_count << '\n';
#endif

                //new_node = new tree_node_sample();

                //new_node->set_count(curr_count, depth);

                new_node.first = leaves[curr_count];

                //new_node->set_uniform(depth);

                num_nodes++;
                num_zero_nodes++;


                curr_node->set_child(curr_dim, curr_cut, new_node.first);

                //region_cache.insert(working_reg, new_node, working_hash);

                if (num_nodes % 100000 == 0) {
                    cerr << "Nodes(" << pile[0].dim << "):" << num_nodes
                            << " : " << num_zero_nodes
                            << " : " << (num_nodes-num_zero_nodes) << '\n';
                }

                working_reg.uncut(curr_dim);

            } else {

#ifdef DEBUG3
                cerr << "FOUND node: " << curr_dim << "," << curr_cut << '\n';

                cerr << "new lphi:  " << new_node->get_lphi() << '\n';
                cerr << "new count: " << new_node->get_count() << '\n';
#endif
                if (!isinf(new_node.first->get_lphi())) {
                    curr_node->set_child(curr_dim, curr_cut, new_node.first);

                    working_reg.uncut(curr_dim);

                } else {

                    pile.push_back(pile_t<tree_node_sample*>());
                    depth++;

                    cut_region(pile[depth - 1].data,pile[depth].data,
                            curr_dim, curr_cut, curr_reg.get_lim(curr_dim));


                    curr_reg.cut(curr_dim, curr_cut);

                    pile[depth].dim = 0;
                    pile[depth].cut = 0;

                    pile[depth].node = new_node.first;



                    curr_node->set_child(curr_dim, curr_cut, new_node.first);
                }



            }



        }





        cerr << "Nodes:" << num_nodes
                << ", Zero nodes:" << num_zero_nodes
                << ", Non-Zero nodes:" << (num_nodes-num_zero_nodes) <<'\n';

        // print the regions
        //region_cache.print_regions();
    }




    /////////////////////////
    // the tree and the regions
    // N is number of data points

    void construct_MAP_tree(map_tree &map_region_tree, opt_region_hash<uint32_t> &map_regions, int N) {


        // Check tree and region is empty

        // NEED TO DO

        // copy the root over



        vector<pile_t<uint32_t> > map_pile;
        map_pile.push_back(pile_t<uint32_t>());
        map_pile[0].node = map_region_tree.get_full_tree();
        map_pile[0].dim = -1;
        map_pile[0].cut = -1;

        region_allocator<map_tree_node> *map_ra = map_region_tree.get_ra();

        // initialise state variables
        int depth = 0;


        gamma_table gt(N);



        vector<pile_t<tree_node_sample*> > pile;
        pile.push_back(pile_t<tree_node_sample*>());

        pile[0].data = vector<vector<double> >();
        pile[0].node = root;
        pile[0].dim = -1;
        pile[0].cut = -1;

        //current_region curr_reg(num_children);
        opt_region working_reg(num_children);


        bool done = false;


        while (!done) {

            if (pile.size() == 0) {
                done = true;
                continue;
            }


            // this is only kind of temporary

            int curr_dim = pile[depth].dim;
            int curr_cut = pile[depth].cut;
            tree_node_sample* curr_node = pile[depth].node;

            uint32_t curr_map_node = map_pile[depth].node;

#ifdef DEBUG_MAP
            cerr << "===-----------------------" << '\n';
            cerr << "===--- DEPTH = " << depth << " count=" << curr_node->get_count() <<'\n';
            cerr << "working_reg:";
            working_reg.print_region_limits(cerr);
            cerr << '\n';
#endif

            bool back_up = false;
            bool add_region = false;
            int map_dim = -1;

            // work out whether we stop at this node
            if (curr_node->get_num_children() == 0 || curr_node->get_count() <= count_lim || depth >= max_depth) {
                // we are already at a uniform node

                // add to regions
                add_region = true;
                back_up = true;
#ifdef DEBUG_MAP
                cerr << "BOTTOM " << '\n';
#endif
            } else if (curr_dim == -1) {
                double post_rho = -c::l2;
                post_rho += depth * c::l2 * curr_node->get_count(); // phi_0
                post_rho -= curr_node->get_lphi2(depth);
#ifdef DEBUG_MAP
                cerr << "post_rho = " << post_rho << '\n';
#endif
                if (post_rho>-c::l2) {

                    // add to regions
                    add_region = true;
                    back_up = true;
                } else {
                    // choose a dimension

                    map_dim = 0;
                    tree_node_sample* child_0 = curr_node->get_child(0, 0);
                    tree_node_sample* child_1 = curr_node->get_child(0, 1);

                    double max_post_prob = gt.compute_lD2(curr_node->get_count()
                            , child_0->get_count(), child_1->get_count());
                    max_post_prob += child_0->get_lphi2(depth+1);
                    max_post_prob += child_1->get_lphi2(depth+1);
#ifdef DEBUG_MAP
                    cerr << "post_prob[" << child_0->get_count() << "|" << child_1->get_count() << "](0) = " << max_post_prob << '\n';
#endif
                    for (int i = 1; i < num_children; i++) {
                        child_0 = curr_node->get_child(i, 0);
                        child_1 = curr_node->get_child(i, 1);

                        double post_prob = gt.compute_lD2(curr_node->get_count()
                                , child_0->get_count(), child_1->get_count());
                        post_prob += child_0->get_lphi2(depth+1);
                        post_prob += child_1->get_lphi2(depth+1);
#ifdef DEBUG_MAP
                        cerr << "post_prob[" << child_0->get_count() << "|" << child_1->get_count() << "](" << i << ") = " << post_prob << '\n';
#endif
                        if (post_prob > max_post_prob) {
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

            } else {
                // we have already chosen a dimension to cut...

                if (curr_cut < c::cuts - 1) {
                    pile[depth].cut++;
                    curr_cut++;
                } else {
                    back_up = true;
                }
            }

#ifdef DEBUG_MAP
            cerr << "curr_dim = " << curr_dim << '\n';
            cerr << "curr_cut = " << curr_cut << '\n';
#endif
            if (add_region) {
#ifdef DEBUG_MAP
                cerr << "ADD REGION" << '\n';
#endif
                // add region to the hash
                map_regions.insert(working_reg, curr_map_node);

                // compute the density
                (*map_ra)[curr_map_node]->set_area(-depth);
                (*map_ra)[curr_map_node]->set_count(curr_node->get_count());

            }


            if (back_up) {
#ifdef DEBUG_MAP
                cerr << "BACKUP" << '\n';
#endif
                depth--;
                pile.pop_back();
                if (depth < 0) continue;

                //curr_reg.uncut(pile[depth].dim, pile[depth].cut);

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

            map_pile.push_back(pile_t<uint32_t>());
            map_pile[depth].node = new_map_node.first;



            working_reg.cut(curr_dim, curr_cut);

            pile.push_back(pile_t<tree_node_sample*>());



            //curr_reg.cut(curr_dim, curr_cut);


            pile[depth].dim = -1;
            pile[depth].cut = -1;
            pile[depth].node = curr_node->get_child(curr_dim, curr_cut);


        }


    }

    int get_num_children() {

        return num_children;
    }

    tree_node_sample* get_full_tree() {

        return root;
    }

    double get_lphi() {
        return root->get_lphi();
    }
};





#endif	/* OPT_TREE_SAMPLE_H */

