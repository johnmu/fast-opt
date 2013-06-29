/*
 *  map_tree.h
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



#ifndef MAP_TREE_H
#define	MAP_TREE_H


#include "general_utils.h"
#include "gamma_table.h"
#include "opt_utils.h"



class map_tree_node
{
private:
    int count;
    int area;

    int dim;
    child_t children;

    static const int num_children = 2;

    void init(){
        dim = -1;
        count = -1;
        area = 1;  // this means it is a leaf
    }

public:

    // stub :S
    map_tree_node(int a){
        init();
    }

    map_tree_node(){
        init();
    }


    void set_count(int count){
        this->count = count;
    }

    int get_count(){
        return count;
    }


    void set_area(int area){
        this->area = area;
    }

    int get_area(){
        return area;
    }

    void set_dim(int dim){
        this->dim = dim;
    }

    int get_dim(){
        return dim;
    }

    double get_area2(){
        return exp(area*c::l2);
    }


    // N is total number of points
    double get_density(int N){
        return exp(log(count/(double)N) - (area*c::l2));
    }

    void set_child(int cut,uint32_t node){
        children.val[cut] = node;
    }

    uint32_t get_child(int cut){
        return children.val[cut];
    }
    
    bool is_leaf(){
        return (dim==-1);
    }

};



inline void print_MAP_density(ostream &o,vector<pair<opt_region, uint32_t> > vec,region_allocator<map_tree_node> *ra, int N){

    double total_area = 0.0;
    double total_density = 0.0;
    int    total_count   = 0;

    for(int i = 0;i<(int)vec.size();i++){
        vec[i].first.print_region_limits(o);
        o << scientific << (*ra)[vec[i].second]->get_density(N);
        o << '\n';

        double area = exp(vec[i].first.get_area()*c::l2);

        total_area    += area;
        total_density += (*ra)[vec[i].second]->get_density(N) * area;
        total_count   += (*ra)[vec[i].second]->get_count();

    }

    cerr << "Regions: " << vec.size() <<   ", Total area: " << total_area
            << ", Total density: " << total_density
            << ", Total count: " << total_count << '\n';
}



class map_tree
{
private:
    uint32_t root;
    region_allocator<map_tree_node> ra;
    int num_points; // total number of data points, used to compute density

    void init(int num_points){
        pair<uint32_t,map_tree_node*> out = ra.create_node();
        root = out.first;
        this->num_points = num_points;
    }

public:

    map_tree(int num_points){
        init(num_points);
    }


    region_allocator<map_tree_node>* get_ra(){
        return &ra;
    }

    uint32_t get_full_tree(){
        return root;
    }
    
    double get_density(const vector<double> &data){
        
        if(ra[root]->get_dim() != (int)data.size()){
            cerr << "Warning: wrong dimension" << '\n';
            return -c::inf;
        }
        
        uint32_t curr_node = root;
        
        current_region curr_reg(data.size());
        
        while(!ra[curr_node]->is_leaf()){
            int curr_dim = ra[curr_node]->get_dim();
            int curr_cut = curr_reg.determine_cut(curr_dim,data[curr_dim]);
            curr_reg.cut(curr_dim,curr_cut);
            curr_node = ra[curr_node]->get_child(curr_cut);
        }
        
        return ra[curr_node]->get_density(num_points);
        
    }

};





#endif	/* MAP_TREE_H */

