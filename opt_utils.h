/*
 *  opt_utils.h
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
    
    void save(ostream & out) const{
        for(int i = 0;i<c::cuts;i++){
            out.write((char*)&(val[i]),sizeof(val[i]));
        }
    }
    
    void load(istream & in){
        for(int i = 0;i<c::cuts;i++){
            in.read((char*)&(val[i]),sizeof(val[i]));
            //cerr << "Child[" << i << "]:" << val[i] << '\n';
        }
    }
};



template <typename T,typename U>
struct pile_t{
    vector<U> data;
    T node;
    int dim;
    int cut;
};

template <typename T,typename U>
struct cpile_t{
    vector<U> data[2];
    T node;
    int dim;
    int cut;
    
    cpile_t(const cpile_t &a){
        data[0] = a.data[0];
        data[1] = a.data[1];
        node = a.node;
        dim = a.dim;
        cut = a.cut;
    }
    
    cpile_t(){
        data[0] = vector<U>();
        data[1] = vector<U>();
    }
};

// stores the cuts that have been made
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

        bool out = dim_cuts < b.dim_cuts;

        return out;
    }
    
    bit_str & operator[](const int a){
        return dim_cuts[a];
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

            pair<double, double> lims = get_limits(i);
            o << scientific << lims.first << ' ' << lims.second << ' ';
        }

    }

    // get the limits of the d-th dimension

    pair<double, double> get_limits(int d) const{
        double middle = 0.5;
        double length = 0.5;

        for (int j = 0; j < (int) dim_cuts[d].size(); j++) {
            length = length / 2.0;
            if (dim_cuts[d][j] == 0) {
                middle = middle - length;
            } else {
                middle = middle + length;
            }
        }
        
        return pair<double, double>(middle - length,middle + length);
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
    
    
    void save(ostream & out) const{
        uint32_t len = dim_cuts.size();
        
        out.write((char*)&len,sizeof(len));
        
        for(uint32_t i = 0;i<len;i++){
            dim_cuts[i].save(out);
        }
    }
    
    void load(istream & in){
        uint32_t len = 0;
        
        in.read((char*)&len,sizeof(len));
        
        //cerr << "in len: " << len << '\n';
        
        dim_cuts.resize(len,bit_str());
        
        for(uint32_t i = 0;i<len;i++){
            //cerr << "len[i]= " << i << '\n';
            dim_cuts[i].load(in);
        }
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
    
    // common save and load operations (these should be private)
    void save_common(ostream & out) const{
        uint32_t free_len = (uint32_t)free_locs.size();
        uint32_t store_len = (uint32_t)store.size();
        
        out.write((char*)&free_len,sizeof(free_len));
        out.write((char*)&store_len,sizeof(store_len));
        
        for(uint32_t i = 0;i<free_len;i++){
            out.write((char*)&free_locs[i],sizeof(uint32_t));
        }
    }
    
    void load_common(istream & in){
        uint32_t free_len = 0;
        uint32_t store_len = 0;
        
        in.read((char*)&free_len,sizeof(free_len));
        in.read((char*)&store_len,sizeof(store_len));
        
        //cerr << "free_len: " << free_len << '\n';
        //cerr << "store_len: " << store_len << '\n';
        
        free_locs.clear();
        for(uint32_t i = 0;i<free_len;i++){
            free_locs.push_back(0);
        }
        
        for(uint32_t i = 0;i<free_len;i++){
            in.read((char*)&free_locs[i],sizeof(uint32_t));
        }
        
        store.clear();
        for(uint32_t i = 0;i<store_len;i++){
            store.push_back(T());
        }
    }
    
    // save and load for custom save/load
    void save2(ostream & out) const{
        save_common(out);
        uint32_t store_len = (uint32_t)store.size();
        
        for(uint32_t i = 0;i<store_len;i++){
            store[i].save(out);
        }
    }
    
    void load2(istream & in){
        load_common(in);
        uint32_t store_len = (uint32_t)store.size();
        
        //cerr << "store_len2: " << store_len << '\n';
        
        for(uint32_t i = 0;i<store_len;i++){
            //cerr << "store[i]: " << i << '\n';
            store[i].load(in);
        }
    }
    
    // default save/load
    void save(ostream & out) const{
        save_common(out);
        uint32_t store_len = (uint32_t)store.size();
        
        for(uint32_t i = 0;i<store_len;i++){
            out.write((char*)&store[i],sizeof(T));
        }
    }
    
    void load(istream & in){
        load_common(in);
        uint32_t store_len = (uint32_t)store.size();
        
        for(uint32_t i = 0;i<store_len;i++){
            in.read((char*)&store[i],sizeof(T));
        }
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
    
    void deinit(){
        //  iterate through and delete all the nodes
        
        if(map_table != NULL){
            for (uint32_t i = 0; i < table_size; i++) {
                if (map_table[i] != NULL) {
                    delete map_table[i];
                }
            }
            delete [] map_table;
        }
    }
    
    opt_region_hash(){
        // in this is used, must load
        map_table = NULL;
    }

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
       deinit();
    }
    
    void init_table(int table_bits){
        deinit();
        this->table_bits = table_bits;
        this->table_size = 1u << table_bits;
        this->mask = this->table_size - 1;

        map_table = new map<opt_region, T>*[table_size];
        for (uint32_t i = 0; i < table_size; i++) {
            map_table[i] = NULL;
        }
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
    
    // returns the node that was erased
    T erase(opt_region& reg) {
        return erase(reg, hash(reg));
    }

    T erase(opt_region& reg, uint32_t hash) {
        
        // BUG: This is bad fix this
        T del_node = (T)c::ra_null_val;
        
        if (map_table[hash] != NULL) {
            typename map<opt_region, T>::iterator it = map_table[hash]->find(reg);
            if (it != map_table[hash]->end()) {
                del_node = it->second;
                map_table[hash]->erase(it);
            } 
        }
        if(map_table[hash]->size()==0){
            delete map_table[hash];
            map_table[hash] = NULL;
        }
        return del_node;
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

    vector<pair<opt_region, T> > get_regions() {

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
    
    int get_num_regions(){
        int num_regions = 0;
        for (uint32_t i = 0; i < table_size; i++) {
            if (map_table[i] != NULL) {
                num_regions++;
            }
        }
        return num_regions;
    }

    void save(ostream & out) const{
        out.write((char*)&table_bits,sizeof(table_bits));
        out.write((char*)&table_size,sizeof(table_size));
        out.write((char*)&mask,sizeof(mask));
        
        // count the number of non-null entries
        uint32_t non_null_count = 0;
        for (uint32_t i = 0; i < table_size; i++) {
            if (map_table[i] != NULL) {
                non_null_count++;
            }
        }
        
        out.write((char*)&non_null_count,sizeof(non_null_count));
        
        // write out the non-null entries of the hash-table
        for (uint32_t i = 0; i < table_size; i++) {
            if (map_table[i] != NULL) {
                // write out the location
                out.write((char*)&i,sizeof(i));
                
                // write out the Map
                map<opt_region, T>* it = map_table[i];
                uint32_t map_len = it->size();
                
                out.write((char*)&map_len,sizeof(map_len));
                
                for(typename map<opt_region, T>::iterator map_it = it->begin();
                        map_it != it->end();map_it++){
                    map_it->first.save(out);
                    out.write((char*)&(map_it->second),sizeof(map_it->second));
                }
            }
        }
    }
    
    void load(istream & in){
        // delete the old table first
        if (map_table != NULL) {
            for (uint32_t i = 0; i < table_size; i++) {
                if (map_table[i] != NULL) {
                    delete map_table[i];
                }
            }
            delete [] map_table;
        }
        //cerr << "delete ok\n";
        
        in.read((char*)&table_bits,sizeof(table_bits));
        in.read((char*)&table_size,sizeof(table_size));
        in.read((char*)&mask,sizeof(mask));
        
        //cerr << "table_bits: " <<table_bits << '\n';
        //cerr << "table_size: " << table_size << '\n';
        //cerr << "mask: " << mask <<'\n';
        
        uint32_t non_null_count = 0;
        in.read((char*) &non_null_count, sizeof (non_null_count));

        //cerr << "non_null_count: " << non_null_count << '\n';

        // initialise the table
        map_table = new map<opt_region, T>*[table_size];
        for (uint32_t i = 0; i < table_size; i++) {
            map_table[i] = NULL;
        }
        
        // read in the values
        for (uint32_t i = 0; i < non_null_count; i++) {

            // read in the location
            uint32_t loc = 0;
            in.read((char*) &loc, sizeof (loc));
            
            //cerr << "loc: " << loc << '\n';
            
            // initialise the location
            map_table[loc] = new map<opt_region, T>();

            // read in the map
            map<opt_region, T>* it = map_table[loc];
            uint32_t map_len = 0;

            in.read((char*) &map_len, sizeof (map_len));
            
            //cerr << "map_len: " << map_len << '\n';

            for(uint32_t j = 0;j<map_len;j++){
                opt_region temp_region;
                T val;
                
                //cerr << "map[j]: " << j << '\n';
                
                temp_region.load(in);
                in.read((char*)&(val),sizeof(val));
                
                //cerr << "val: " << val << '\n';
                
                it->insert(pair<opt_region, T>(temp_region, val));
            }
            //cerr << "load map done" << '\n';

        }
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




// stores the boundaries of a region 
class current_region{
private:
    vector<double> lim_pile; // center of the region
    vector<double> resolution; // 1/2 the size of the region
    int num_children;

public:

    void init(int num_children){
        this->num_children = num_children;
        lim_pile = vector<double>(num_children,0.5);
        resolution = vector<double>(num_children,0.5);
    }

    current_region():lim_pile(),resolution(){
    }

    current_region(int num_children):lim_pile(num_children,0.5),resolution(num_children,0.5){
        this->num_children = num_children;
    }

    double get_lim(int dim){
        return lim_pile[dim];
    }

    double get_resolution(int dim){
        return resolution[dim];
    }
    
    vector<pair<double,double> > get_limits(){
        vector<pair<double,double> > result;
        
        for (int i = 0;i<num_children;i++){
            double res = resolution[i];
            double lim = lim_pile[i];
            result.push_back(pair<double,double>(lim-res,lim+res));
        }
        
        return result;
    }
    
    // this one returns min and range
    vector<pair<double,double> > get_limits2(){
        vector<pair<double,double> > result;
        
        for (int i = 0;i<num_children;i++){
            double res = resolution[i];
            double lim = lim_pile[i];
            result.push_back(pair<double,double>(lim-res,2*res));
        }
        
        return result;
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

