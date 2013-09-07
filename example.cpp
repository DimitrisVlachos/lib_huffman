
//Simple example without I/O optimizations(ie buffered IO)
//Requires my lib_bitstreams :
//https://github.com/DimitrisVlachos/lib_bitstreams

#include "chuffman.hpp"

void test_lib_huff(const char* fn) {
	bit_streams::bit_stream_writer_c out;
	bit_streams::bit_stream_reader_c in;
	huffman_c huff;
 	file_streams::file_stream_if* rd,*wr;
	std::vector<uint64_t> freq;
	

	//Open streams
	rd = new file_streams::file_stream_reader_c(fn);
	wr = new file_streams::file_stream_writer_c("decompressed.bin");
	if (!rd || !wr) {
		printf("I/O error\n");
		delete rd;
		delete wr;
		return;
	}

	//Calculate frequencies
	out.open("compressed.bin");	

	for (uint32_t i = 0;i < 256;++i)
		freq.push_back(0);
	
	for (uint32_t i = 0;i < rd->size();++i)	
		++freq[rd->read()];

	//Generate tree + store header to output stream
	huff.generate_tree_from_probabilities(freq);
	huff.store(&out);

	//Compress
	rd->seek(0);
	for (uint32_t i = 0;i < rd->size();++i)	
		huff.encode_sym(&out,rd->read());

	out.close();

	//Open compressed file and load header
	in.open("compressed.bin");
	huff.load(&in);
 
	//Validate
	bool match = true;
	rd->seek(0);
	for (uint32_t i = 0;i < rd->size();++i)	{
		uint32_t v = huff.decode_sym(&in);
		uint32_t v2 = rd->read();
		if (v != v2)  {
			printf("%c != %c\n",v,v2);
			match = false;
		}
		wr->write(v);
	}
 
	printf("Done.Files %s\n",match ? "matched!" : "did not match");

	delete rd;
	delete wr;
}

int main() {
	test_lib_huff("chuffman.cpp");
	return 0;
}


