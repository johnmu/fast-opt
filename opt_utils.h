/*
 *  opt_utils.h
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



#ifndef OPT_UTILS_H
#define	OPT_UTILS_H

#include "general_utils.h"

struct child_t{
    uint32_t val[c::cuts];

    child_t(){
        for(int i = 0;i<c::cuts;i++){
            val[i] = c::ra_null_val;
        }
    }
};



template <typename T>
struct pile_t{
    vector<vector< double> > data;
    T node;
    int dim;
    int cut;
};



//template <typename T>
class opt_region{
public:
    vector< bit_str > dim_cuts; // vector vector bool is SLOW

    void init(int num_children){
        dim_cuts.resize(num_children);
    }

    opt_region(){

    }

    opt_region(int num_children){
        init(num_children);
    }


    //opt_region(vector<pile_t<T> > &dim_pile, int num_children){
    //    init(num_children);
    //    for (int i = 0;i<(int) dim_pile.size();i++){
    //        dim_cuts[dim_pile[i].dim].push_back(dim_pile[i].cut);
    //    }
    //}

    ~opt_region(){

    }

    bool full(){
        for(vector< bit_str >::iterator it = dim_cuts.begin();
                it != dim_cuts.end();it++){
            if(it->size() >= c::MAX_COORD_DEPTH){
                cerr << "Full " << it->size() << ":";
                print_region();
                cerr << '\n';
                return true;
            }
        }
        return false;
    }

    bool cut(int dim, int cut){
        return dim_cuts[dim].push_back(cut);
    }

    void uncut(int dim){
        if(dim_cuts[dim].size() == 0){
            cerr << "Error, empty dim: " << dim << '\n';
        }
        //dim_cuts[dim].erase(dim_cuts[dim].size()-1);
        dim_cuts[dim].pop_back();
    }

    bool operator==(const opt_region& b)const {
        for(int i = 0;i<(int)dim_cuts.size();i++){
            if (this->dim_cuts[i] != b.dim_cuts[i]) return false;
        }
        return true;
    }

    bool operator<(const opt_region& b) const{
        //cerr << "Compare: ";
        //this->print_region();
        //cerr << " to ";
        //b.print_region();
        //cerr << "\n";

        bool out = dim_cuts < b.dim_cuts;

        //if(out){
        //    cerr << "ALL-LESS\n";
        //}else{
        //    cerr << "ALL-NOTLESS\n";
        //}

        return out;
    }


    void print_region() const{
        print_region(cerr);
    }

    void print_region(ostream &o) const{
        for(int i = 0;i<(int)dim_cuts.size();i++){
            for (int j = 0;j<(int)dim_cuts[i].size();j++){
                o << dim_cuts[i][j];
            }
            o << ',';
        }

    }


    void print_region_limits() const{
        print_region_limits(cout);
    }

    void print_region_limits(ostream &o) const{
        for(int i = 0;i<(int)dim_cuts.size();i++){

            double middle = 0.5;
            double length = 0.5;

            for (int j = 0;j<(int)dim_cuts[i].size();j++){
                length = length/2.0;
                if (dim_cuts[i][j] == 0){
                    middle = middle - length;
                }else{
                    middle = middle + length;
                }
            }
            o << scientific << (middle - length) << ' ' << (middle + length) << ' ';
        }

    }

    int get_area(){

        int area = 0;

        for(int i = 0;i<(int)dim_cuts.size();i++){
            area -= dim_cuts[i].size();
        }

        return area;
    }

    bool is_child(const opt_region &reg) const{

        if(reg.dim_cuts.size() != dim_cuts.size()){
            cerr << "opt_region: Bad dimension compare" << '\n';
            exit(2);
        }

        for(int i = 0;i<(int)dim_cuts.size();i++){
            if(!dim_cuts[i].is_child(reg.dim_cuts[i])){
                return false;
            }
        }

        return true;
    }

    int num_children() const{
        return (int)dim_cuts.size();
    }

};



template <typename T>
class region_allocator{
public:
    vector<T> store;
    vector<uint32_t> free_locs;

    pair<uint32_t,T*> create_node(int num_children){
        pair<uint32_t,T*> output;


        if (free_locs.size() == 0) {

            store.push_back(T(num_children));

            output.second = &store.back();
            output.first = (uint32_t) (store.size() - 1);

            if (output.first == c::ra_null_val) {
                cerr << "region_allocator: NULL create!\n";
                exit(2);
            }
        }else{
            output.first = free_locs.back();

            if (output.first == c::ra_null_val) {
                cerr << "region_allocator: NULL create free!\n";
                exit(2);
            }

            store[output.first] = T(num_children);
            output.second = &(store[output.first]);

            free_locs.pop_back();
        }

        return output;

    }


    pair<uint32_t,T*> create_node(){
        pair<uint32_t, T*> output;

        if (free_locs.size() == 0) {
            store.push_back(T());
            output.second = &store.back();
            output.first = (uint32_t) (store.size() - 1);

            if (output.first == c::ra_null_val) {
                cerr << "region_allocator: NULL create!\n";
                exit(2);
            }

        } else {
            output.first = free_locs.back();

            if (output.first == c::ra_null_val) {
                cerr << "region_allocator: NULL create freeee!\n";
                exit(2);
            }

            store[output.first] = T();
            output.second = &(store[output.first]);

            free_locs.pop_back();
        }

        return output;

    }

    void delete_node(uint32_t idx){
        free_locs.push_back(idx);
    }

    T* operator[](uint32_t idx){
        if(idx>=(uint32_t)store.size()){

            cerr << "region_allocator: Out of range: "<< idx<< "\n";

            return NULL;
        }else if(idx == c::ra_null_val){

            cerr << "region_allocator: NULL\n";

            return NULL;
        }

        return &store[idx];
    }
};



template <typename T>
class opt_region_hash {
public:
    static const uint32_t magic = 2654435761u;

    map<opt_region, T>** map_table;
    int table_bits;
    uint32_t table_size;
    uint32_t mask;

    opt_region_hash(int table_bits) {
        this->table_bits = table_bits;
        this->table_size = 1u << table_bits;
        this->mask = this->table_size - 1;

        map_table = new map<opt_region, T>*[table_size];
        for (uint32_t i = 0; i < table_size; i++) {
            map_table[i] = NULL;
        }
    }

    ~opt_region_hash() {
       //  iterate through and delete all the nodes
        for (uint32_t i = 0; i < table_size; i++) {
            if (map_table[i] != NULL) {
                delete map_table[i];
            }
        }
        delete [] map_table;
    }

    uint32_t hash(const opt_region& reg) const {
        uint32_t hash_val = 0;
        for (int i = 0; i < reg.num_children(); i++) {
            hash_val = hash_val ^ (reg.dim_cuts[i].data + i);
        }
        return (hash_val*magic) & mask;
    }

    pair<T,bool> find(opt_region& reg) {
        return find(reg, hash(reg));
    }

    pair<T,bool> find(opt_region& reg, uint32_t hash) {
        if (map_table[hash] == NULL) {
            return pair<T,bool>(T(),false);
        } else {

            typename map<opt_region, T>::iterator it = map_table[hash]->find(reg);
            if (it == map_table[hash]->end()) {
                //cerr << "NOTFOUND\n";
                return pair<T,bool>(T(),false);
            } else {
                //cerr << "FOUND\n";
                return pair<T,bool>(it->second,true);
            }
        }

        return pair<T,bool>(T(),false);
    }

    void insert(opt_region& reg, T node) {
        insert(reg, node, hash(reg));
    }

    void insert(opt_region& reg, T node, uint32_t hash) {
        if (map_table[hash] == NULL) {
            map_table[hash] = new map<opt_region, T>();
        }


        map_table[hash]->insert(pair<opt_region, T>(reg, node));
    }

    void print_regions() {
        for (uint32_t i = 0; i < table_size; i++) {
            if (map_table[i] != NULL) {
                cerr << i << "\n";
                for (typename map<opt_region, T>::iterator it = map_table[i]->begin();
                        it != map_table[i]->end(); it++) {
                    it->first.print_region();
                    cerr << '\n';
                }
            }
        }
    }

    vector<pair<opt_region, T> > print_density() {

        vector<pair<opt_region, T> > output;

        for (uint32_t i = 0; i < table_size; i++) {
            if (map_table[i] != NULL) {
                for (typename map<opt_region, T>::iterator it = map_table[i]->begin();
                        it != map_table[i]->end(); it++) {
                    output.push_back(pair<opt_region, T>(it->first,it->second));
                }
            }
        }

        return output;
    }


};



inline int choose_dim(const vector<double> &d, MT_random &rand_gen) {
    int N = d.size();

    double total = 0.0;

    for (int i = 0; i < N; i++) {
        total += d[i];
    }

    double val = rand_gen.genrand64_real1() * total;

    total = 0.0;

    for (int i = 0; i < N; i++) {
        total += d[i];
        if (total >= val) return i;
    }

    return N - 1;

}





class current_region{
private:
    vector<double> lim_pile;
    vector<double> resolution;
    int num_children;



public:

    void init(int num_children){
        this->num_children = num_children;
        this->lim_pile = vector<double>(num_children,0.5);
        this->resolution = vector<double>(num_children,0.5);
    }

    current_region(){

    }

    current_region(int num_children){
        init(num_children);
    }

    double get_lim(int dim){
        return lim_pile[dim];
    }

    double get_resolution(int dim){
        return resolution[dim];
    }

    void cut(int dim, int cut) {
        resolution[dim] = resolution[dim] / 2.0;
        if (cut == 0) {
            lim_pile[dim] = lim_pile[dim] - resolution[dim];
        } else {
            lim_pile[dim] = lim_pile[dim] + resolution[dim];
        }

    }

    void uncut(int dim, int cut) {
        if (cut == 0) {
            lim_pile[dim] = lim_pile[dim] + resolution[dim];
        } else {
            lim_pile[dim] = lim_pile[dim] - resolution[dim];
        }
        resolution[dim] = resolution[dim] * 2.0;
    }


    int count(const vector<vector<double> > &data){
        int output = 0;

        for(int i = 0;i<(int)data.size();i++){

            if((int)data[i].size() != num_children) return -1;

            bool add = true;

            for(int j = 0;j<(int)data[i].size();j++){
                if(data[i][j] < (lim_pile[j]-resolution[j])
                        || data[i][j] >= (lim_pile[j]+resolution[j])){

                    add = false;
                }
            }

            if(add) output++;
        }

        return output;
    }
    
    
    // determine which way the cut should be for MAP tree
    int determine_cut(int dim, double value){
        if(value >= lim_pile[dim]){
            return 1;
        }else{
            return 0;
        }
    }

};




#endif	/* OPT_UTILS_H */

