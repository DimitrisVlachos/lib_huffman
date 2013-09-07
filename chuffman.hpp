#ifndef __chuffman_hpp__
#define __chuffman_hpp__
#include "bit_streams.hpp"
#include <stdint.h>
#include <vector>
#include <algorithm>

namespace chuffman_private {
	struct huff_node_t {
		uint64_t freq;
		uint64_t sym;
		uint32_t leaf;
		uint32_t _pad;

		huff_node_t* l;
		huff_node_t* r;
		 
		huff_node_t() : freq(0U),leaf(0U),sym(0U),l(0),r(0) {}
		huff_node_t(huff_node_t* in_l,huff_node_t* in_r) : freq(in_l->freq + in_r->freq),leaf(0U),sym(0U),l(in_l),r(in_r) {}
		huff_node_t(uint64_t in_sym,uint64_t in_freq) : freq(in_freq),leaf(1U),sym(in_sym),l(0),r(0) {}
		~huff_node_t() { delete l; delete r; }
	};

	struct huff_code_ref_t  {
		uint64_t code;
		uint64_t len;
		huff_code_ref_t() : code(0U),len(0U) {}
	};

	static inline uint64_t clc_blen(uint64_t in) { //TODO LUT
		register uint64_t bits = (uint64_t)1U,v = in;

		while ((v>>=(uint64_t)1U) != (uint64_t)0U)
			++bits;

		return bits;
	}

	static inline const bool cmp_node(const huff_node_t* l,const huff_node_t* r) { 
		return l->freq > r->freq; 
	}

	void gen_huff_code_list_recursive(huff_node_t* root,uint64_t code,uint32_t bit_len,std::vector<huff_code_ref_t>& code_map);
	bool store_recursive(bit_streams::bit_stream_writer_c* out,huff_node_t* root,const uint64_t sym_lo,const uint32_t sym_blen);
	huff_node_t* load_recursive(bit_streams::bit_stream_reader_c* in,const uint64_t sym_lo,const uint32_t sym_blen);
	huff_node_t* gen_huff_tree(std::vector<huff_node_t*>& nodes);
} 

class huffman_c {
	private:
	std::vector<chuffman_private::huff_code_ref_t> code_map;
	uint32_t min_code_len,max_code_len;
	uint32_t used_syms;
	chuffman_private::huff_node_t* tree_root;

	public:

	huffman_c() : min_code_len(0),max_code_len(0),used_syms(0),tree_root(0){}
	~huffman_c() { delete tree_root; }


	template <typename base_t>
	bool generate_tree_from_probabilities(const std::vector<base_t>& freq_table) {

		std::vector<chuffman_private::huff_node_t*> nodes;

		uint32_t used_syms = 0U;
		uint32_t hi_cl = 0U;
		uint32_t lo_cl = 0xffffffffU;

		code_map.clear();
		code_map.reserve(freq_table.size());
		nodes.reserve(256U);
 
		for (uint32_t i = 0U,j = freq_table.size();i < j;++i) {
			const base_t frq = freq_table[i];

			code_map.push_back(chuffman_private::huff_code_ref_t());
			if (0U == frq)
				continue;
			++used_syms;
			nodes.push_back( new chuffman_private::huff_node_t(i,frq) );
			if (i > hi_cl) 
				hi_cl = i;
			if (i < lo_cl) 
				lo_cl = i;
		}
		
		delete tree_root;
		tree_root = chuffman_private::gen_huff_tree(nodes);
		if (0==tree_root)
			return false;

		gen_huff_code_list_recursive(tree_root,0U,0U,code_map);

		if (lo_cl == 0xffffffffU)
			lo_cl = 0U;

		min_code_len = lo_cl;
		max_code_len = hi_cl;
		used_syms = used_syms;

		return true;
	}

	bool store(bit_streams::bit_stream_writer_c* out);
	bool load(bit_streams::bit_stream_reader_c* in);

	inline bool encode_sym(bit_streams::bit_stream_writer_c* out,const uint64_t sym) {
		out->write(code_map[sym].code,code_map[sym].len);
		return true;
	}

	inline uint64_t decode_sym(bit_streams::bit_stream_reader_c* in) {
		if (in->eof())
			return 0;

		register chuffman_private::huff_node_t* root = tree_root;

		while (root) {
			root = (in->read(1) == 1) ? root->r : root->l;
			if (root->leaf)
				return root->sym;	
		}

		return 0;
	}
};
#endif

