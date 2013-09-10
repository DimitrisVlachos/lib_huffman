
/*
	Author : Dimitris Vlachos (DimitrisV22@gmail.com @ github.com/DimitrisVlachos)
*/

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


	static void gen_huff_code_list_recursive(huff_node_t* root,uint64_t code,uint32_t bit_len,std::vector<huff_code_ref_t>& code_map) {
		if (root->leaf) {	//l==0
			code_map[root->sym].code = code;
			code_map[root->sym].len = bit_len;
			return;
		}

		gen_huff_code_list_recursive(root->l,code << (uint64_t)1U,bit_len  + 1U,code_map);
		gen_huff_code_list_recursive(root->r,(code << (uint64_t)1U) | (uint64_t)1U,bit_len + 1U,code_map);
	}

	template <class writer_type_c>
	static bool store_recursive(bit_streams::bit_stream_writer_c<writer_type_c>* out,huff_node_t* root,const uint64_t sym_lo,const uint32_t sym_blen) {
	
		//b0 : leaf , b1 : composite
		if (!root->leaf) {
			out->write(1,1);
			store_recursive<writer_type_c>(out,root->l,sym_lo,sym_blen);
			return store_recursive<writer_type_c>(out,root->r,sym_lo,sym_blen);
		}
		out->write(0,1);
		out->write(root->sym - sym_lo,sym_blen);
		return true;
	}

	template <class reader_type_c>
	static huff_node_t* load_recursive(bit_streams::bit_stream_reader_c<reader_type_c>* in,const uint64_t sym_lo,const uint32_t sym_blen) {
		if (1 == in->read(1)) {
			huff_node_t* l = load_recursive<reader_type_c>(in,sym_lo,sym_blen); if (!l) return 0;
			huff_node_t* r = load_recursive<reader_type_c>(in,sym_lo,sym_blen); if (!r) return 0;
			return new huff_node_t(l,r);
		}
		return new huff_node_t(sym_lo + in->read(sym_blen),0);
	}

	static huff_node_t* gen_huff_tree(std::vector<huff_node_t*>& nodes) {
		uint32_t len = nodes.size();

		if (0U == len)
			return 0;
		else if (1U == len)
			return nodes[0U];

		std::stable_sort(nodes.begin(),nodes.end(),cmp_node);

		huff_node_t* composite;

		do {
			composite = new huff_node_t(nodes[len - 1U],nodes[len - 2U]);
			if (!composite)
				return 0;

			nodes[len - 2U] = composite; //Replace dummy node with composite
			if (--len == 1)
				return nodes[0];
			std::stable_sort(nodes.begin(),nodes.begin() + len,cmp_node); 
		} while (len  > 1U);

		return nodes[0U];
	}
} 

template <class reader_type_c,class writer_type_c>
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

	inline bool store(bit_streams::bit_stream_writer_c<writer_type_c>* out) {
		const uint32_t m_bl = chuffman_private::clc_blen( max_code_len + (max_code_len==0) ); 
		out->write(max_code_len,16);
		out->write(min_code_len,m_bl);
		out->write(used_syms,m_bl);	
	
		const uint32_t blen = chuffman_private::clc_blen( (max_code_len - min_code_len) )  ;
		return chuffman_private::store_recursive(out,tree_root,min_code_len, blen + (0U == blen) );
	}

	inline bool load(bit_streams::bit_stream_reader_c<reader_type_c>* in) {
		max_code_len = in->read(16);
		const uint32_t m_bl = chuffman_private::clc_blen( max_code_len + (max_code_len==0) ); 
		min_code_len = in->read(m_bl);
		used_syms = in->read(m_bl);

		const uint32_t blen = chuffman_private::clc_blen( (max_code_len - min_code_len) ) ;
		delete tree_root;
		tree_root=chuffman_private::load_recursive(in,min_code_len,blen + (0U == blen));

		return tree_root != 0;
	}

	inline bool encode_sym(bit_streams::bit_stream_writer_c<writer_type_c>* out,const uint64_t sym) {
		out->write(code_map[sym].code,code_map[sym].len);
		return true;
	}

	inline uint64_t decode_sym(bit_streams::bit_stream_reader_c<reader_type_c>* in) {
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
