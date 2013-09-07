
/*
	Author : Dimitris Vlachos (DimitrisV22@gmail.com @ github.com/DimitrisVlachos)
*/

#include "chuffman.hpp"

namespace chuffman_private {
	void gen_huff_code_list_recursive(huff_node_t* root,uint64_t code,uint32_t bit_len,std::vector<huff_code_ref_t>& code_map) {
		if (root->leaf) {	//l==0
			code_map[root->sym].code = code;
			code_map[root->sym].len = bit_len;
			return;
		}

		gen_huff_code_list_recursive(root->l,code << (uint64_t)1U,bit_len  + 1U,code_map);
		gen_huff_code_list_recursive(root->r,(code << (uint64_t)1U) | (uint64_t)1U,bit_len + 1U,code_map);
	}

	bool store_recursive(bit_streams::bit_stream_writer_c* out,huff_node_t* root,const uint64_t sym_lo,const uint32_t sym_blen) {
	
		//b0 : leaf , b1 : composite
		if (!root->leaf) {
			out->write(1,1);
			store_recursive(out,root->l,sym_lo,sym_blen);
			return store_recursive(out,root->r,sym_lo,sym_blen);
		}
		out->write(0,1);
		out->write(root->sym - sym_lo,sym_blen);
		return true;
	}

	huff_node_t* load_recursive(bit_streams::bit_stream_reader_c* in,const uint64_t sym_lo,const uint32_t sym_blen) {
		if (1 == in->read(1)) {
			huff_node_t* l = load_recursive(in,sym_lo,sym_blen); if (!l) return 0;
			huff_node_t* r = load_recursive(in,sym_lo,sym_blen); if (!r) return 0;
			return new huff_node_t(l,r);
		}
		return new huff_node_t(sym_lo + in->read(sym_blen),0);
	}

	huff_node_t* gen_huff_tree(std::vector<huff_node_t*>& nodes) {
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

bool huffman_c::store(bit_streams::bit_stream_writer_c* out) {
	const uint32_t m_bl = chuffman_private::clc_blen( max_code_len + (max_code_len==0) ); 
	out->write(max_code_len,16);
	out->write(min_code_len,m_bl);
	out->write(used_syms,m_bl);	
	
	const uint32_t blen = chuffman_private::clc_blen( (max_code_len - min_code_len) )  ;
	return chuffman_private::store_recursive(out,tree_root,min_code_len, blen + (0U == blen) );
}

bool huffman_c::load(bit_streams::bit_stream_reader_c* in) {
	max_code_len = in->read(16);
	const uint32_t m_bl = chuffman_private::clc_blen( max_code_len + (max_code_len==0) ); 
	min_code_len = in->read(m_bl);
	used_syms = in->read(m_bl);

	const uint32_t blen = chuffman_private::clc_blen( (max_code_len - min_code_len) ) ;
	delete tree_root;
	tree_root=chuffman_private::load_recursive(in,min_code_len,blen + (0U == blen));

	return tree_root != 0;
}

