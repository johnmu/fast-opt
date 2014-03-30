/*
 *  map_tree.h
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
        dim = -1; // this means it is a leaf
        count = -1;
        area = 0; // real_area = 2^(-area)  
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

    int get_count() const{
        return count;
    }

    void set_area(int area){
        this->area = area;
    }

    int get_area() const{
        return area;
    }

    void set_dim(int dim){
        this->dim = dim;
    }

    int get_dim() const{
        return dim;
    }

    double get_area2(){
        return exp(area*c::l2);
    }

    // N is total number of points
    double get_density(int N,double pseudo_count, int num_regions){
        return ((count+pseudo_count)/(double)(N+(num_regions*pseudo_count)))*exp(-(area*c::l2));
    }
    
    double get_prob(int N,double pseudo_count, int num_regions){
        return ((count+pseudo_count)/((double)N+(num_regions*pseudo_count)));
    }

    void set_child(int cut,uint32_t node){
        children.val[cut] = node;
    }

    uint32_t get_child(int cut){
        return children.val[cut];
    }
    
    bool is_leaf() const{
        return (dim==-1);
    }
    
    void save(ostream & out) const{
        out.write((char*)&count,sizeof(count));
        out.write((char*)&area,sizeof(area));
        out.write((char*)&dim,sizeof(dim));
        children.save(out);
    }
    
    void load(istream & in){
        in.read((char*)&count,sizeof(count));
        in.read((char*)&area,sizeof(area));
        in.read((char*)&dim,sizeof(dim));
        
        children.load(in);
    }
    
};


// need to use log_sum_exp for area!!
inline void print_MAP_density(ostream &o,vector<pair<opt_region, uint32_t> > vec,region_allocator<map_tree_node> *ra, int N){

    double total_area = 0.0;
    double total_prob = 0.0;
    int    total_count   = 0;

    for(int i = 0;i<(int)vec.size();i++){
        
         double area = exp(vec[i].first.get_area()*c::l2);
        
        vec[i].first.print_region_limits(o);
        o << scientific << (*ra)[vec[i].second]->get_density(N,0,0);
        o << ' ' << (*ra)[vec[i].second]->get_count();
        o << '\n';

        total_area    += area;
        total_prob += (*ra)[vec[i].second]->get_prob(N,0,0);
        total_count   += (*ra)[vec[i].second]->get_count();
    }

    cerr << "Regions: " << vec.size() <<   ", Total area: " << total_area
            << ", Total density: " << total_prob
            << ", Total count: " << total_count << '\n';
}

inline void print_den_influence(ostream &o,vector<pair<opt_region, uint32_t> > vec,region_allocator<map_tree_node> *ra, int N){

    int num_children = vec[0].first.num_children();
    vector<double> influence(num_children,0.0);

    for(int i = 0;i<(int)vec.size();i++){

        for (int j = 0; j < num_children; j++) {
            influence[j] += exp(vec[i].first.get_influence(j)*c::l2); // change to log_sum_exp
        }
    }

    for (int j = 0; j < num_children; j++) {
        o << log(influence[j]) << ',';
    }
    o << '\n';
}



class map_tree
{
private:
    uint32_t root;
    region_allocator<map_tree_node> ra;
    uint32_t num_points; // total number of data points, used to compute density
    int num_children; // The maximum dimension
    
    void init(uint32_t num_points,int num_children){
        pair<uint32_t,map_tree_node*> out = ra.create_node();
        root = out.first;
        this->num_points = num_points;
        this->num_children = num_children;
    }

public:

    map_tree(){
        // if this is used must load
    }
    
    map_tree(uint32_t num_points, int num_children){
        init(num_points,num_children);
    }

    region_allocator<map_tree_node>* get_ra(){
        return &ra;
    }

    uint32_t get_full_tree(){
        return root;
    }
    
    double get_density(const vector<double> &data,double pseudo_count,int num_regions){
        
        if(num_children != (int)data.size()){
            cerr << "Warning: wrong dimension(" << num_children << ")" << '\n';
            return -c::inf;
        }
        
        uint32_t curr_node = root;
        
        current_region curr_reg(num_children);
        
        while(!ra[curr_node]->is_leaf()){
            int curr_dim = ra[curr_node]->get_dim();
            int curr_cut = curr_reg.determine_cut(curr_dim,data[curr_dim]);
            
            curr_reg.cut(curr_dim,curr_cut);
            curr_node = ra[curr_node]->get_child(curr_cut);
        }

        return ra[curr_node]->get_density(num_points,pseudo_count,num_regions);
    }
    
    opt_region get_region(const vector<double> &data){
        
        if(num_children != (int)data.size()){
            cerr << "Warning: wrong dimension(" << num_children << ")" << '\n';
            return -c::inf;
        }
        
        uint32_t curr_node = root;
        
        current_region curr_reg(num_children);
        opt_region working_reg(num_children);
        
        while(!ra[curr_node]->is_leaf()){
            int curr_dim = ra[curr_node]->get_dim();
            int curr_cut = curr_reg.determine_cut(curr_dim,data[curr_dim]);
            
            curr_reg.cut(curr_dim,curr_cut);
            working_reg.cut(curr_dim,curr_cut);
            curr_node = ra[curr_node]->get_child(curr_cut);
        }

        return working_reg;
    }

    uint32_t get_num_points() const{
        return this->num_points;
    }
    
    int get_num_children() const{
        return this->num_children;
    }
    
    void save(ostream & out) const{
        out.write((char*)&root,sizeof(root));
        out.write((char*)&num_points,sizeof(num_points));
        out.write((char*)&num_children,sizeof(num_children));
        
        ra.save2(out);
                
    }
    
    void load(istream & in){
        in.read((char*)&root,sizeof(root));
        in.read((char*)&num_points,sizeof(num_points));
        in.read((char*)&num_children,sizeof(num_children));
        
        //cerr << "Root: " << root << '\n';
        //cerr << "num_points: " << num_points << '\n';
        
        ra.load2(in);
        
        
    }
};


// one-dimensional CDF constructed from OPT regions

struct cdf_t{
    double end;
    double den;
    
    cdf_t(){
        end = 1.0;
        den = 1.0;
    }
    
    cdf_t(const double &end,const double &den){
        this->end = end;
        this->den = den;
    }
    
    bool operator<(const cdf_t &a) const{
        return (end < a.end);
    }
    
};



class cdf{
private:
    vector<cdf_t> cdf_data;
    
public:
    
    // should add some sanity checks
    cdf(map_tree &map_region_tree, opt_region_hash<uint32_t> &map_regions){
        
        uint32_t N = map_region_tree.get_num_points();
        region_allocator<map_tree_node>* ra = map_region_tree.get_ra();
        
        // extract regions as (start, density) assume last one goes to 1
        vector<pair<opt_region, uint32_t> > regs = map_regions.get_regions();

        for (int i = 0; i < (int) regs.size(); i++) {
            // get the start location (assume only 1-D)
   
            pair<double, double> lims = regs[i].first.get_limits(0);
            double end = lims.second;
            
            // get the density
            double den = (*ra)[regs[i].second]->get_density(N,0,0);

            cdf_data.push_back(cdf_t(end,den));
        }
        
        // sort by end location
        sort(cdf_data.begin(),cdf_data.end());
        
        // sum up the densities
        int N_reg = (int)regs.size();
        double cum_sum = 0.0;
        for (int i = 0; i < N_reg; i++){
            
            double dist = 0;
            if(i == 0){
                dist = cdf_data[i].end;
            }else{
                dist = cdf_data[i].end - cdf_data[i-1].end;
            }

            cum_sum += cdf_data[i].den * dist;
            
            // fix rounding errors
            if(cum_sum>1)cum_sum = 1.0;
            
            cdf_data[i].den = cum_sum;
        }
        
    }
    
    double transform(double x) const{
        cdf_t temp;
        temp.end = x;
        temp.den = 0.0;
        
        vector<cdf_t>::const_iterator it = lower_bound(cdf_data.begin(),cdf_data.end(),temp);
        
        double start = 0.0;
        double start_den = 0.0;
        double end = it->end;
        double end_den = it->den;
        if(it != cdf_data.begin()){
            vector<cdf_t>::const_iterator it_prev = it - 1;
            start = it_prev->end;
            start_den = it_prev->den;
        }
        
        // interpolate to find the transformed value
        return ((x-start)/(end-start))*(end_den-start_den) + start_den;
    }
    
    void print_cdf(ostream &o) const{
        for(int i = 0;i<(int)cdf_data.size();i++){
            o << cdf_data[i].end << ":" << cdf_data[i].den << '\n';
        }
    }
};


#endif	/* MAP_TREE_H */

