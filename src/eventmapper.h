//
// Created by Vadim N. on 24/03/2015.
//

#ifndef _YMIR_EVENTMAPPER_H_
#define _YMIR_EVENTMAPPER_H_


#include "types.h"
#include "unordered_map"


namespace ymir {

    class EventMapper {
    public:


        eventind_t V_gene(segindex_t v_index) {

        }


        eventind_t V_deletions(segindex_t v_index, seq_len_t del_num) {

        }


        eventind_t J_gene(segindex_t j_index) {

        }


        eventind_t J_deletions(segindex_t j_index, seq_len_t del_num) {

        }


        eventind_t D_gene(segindex_t d_index) {

        }


        eventind_t D_deletions(segindex_t d_index, seq_len_t D5_del_num, seq_len_t D3_del_num) {

        }


        eventind_t insertion_length(segindex_t order, seq_len_t ins_num) {

        }


        eventind_t insertion_nuc(segindex_t order) {

        }

    protected:

        std::unordered_map<segindex_t, eventind_t> _map;
        
    };
}

#endif //_YMIR_EVENTMAPPER_H_
