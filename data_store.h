/*
 *  data_store.h
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



#ifndef DATA_STORE_H
#define	DATA_STORE_H

#include "stl.h"

template <typename T>
class data_store {
private:
    vector<typename vector<T>::const_iterator> data;
public:
    data_store(){

    }

    data_store(vector<T> &d){
        data.reserve(d.size());
        for(typename vector<T>::const_iterator it = d.begin();it != d.end();it++ ){
            data.push_back(&(*it));
        }
    }

    data_store(data_store<T> &d){
        data.reserve(d.size());

        for(size_t i = 0;i<d.size();i++){
            data.push_back(d.get_ptr(i));
        }
    }

    T operator[](size_t i) const{
        return *(data[i]);
    }

    typename vector<T>::const_iterator get_ptr(size_t i){
        return data[i];
    }

    size_t size(){
        return data.size();
    }


    typename vector<typename vector<T>::const_iterator>::const_iterator begin(){
        return data.begin();
    }

    typename vector<typename vector<T>::const_iterator>::const_iterator end(){
        return data.end();
    }
};



#endif	/* DATA_STORE_H */

