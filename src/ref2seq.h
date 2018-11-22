//
//
//

#ifndef BWA_FASTSEQ_H
#define BWA_FASTSEQ_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "align_proc.h"

class ref2seq {
public:
    ref2seq(int blockSize_, int max_mis_, int max_readLen_, std::fstream& refFile);
    std::string getSeq(align_info& info, int readl);

private:
    int block_num;
    int blockSize, max_mis, max_readLen;
    std::fstream *reference;
    std::string ref_seg, buffer, tail; //tail for reads aligned to the end of each block
    std::map<std::string, std::string> refbase2num;
    map<char,map<char,int>> var2base;
    map<char, char> compbase;
    void getSeg();
};

ref2seq::ref2seq(int blockSize_, int max_mis_, int max_readLen_, std::fstream& refFile) {
    block_num = -1;
    blockSize = blockSize_;
    max_mis = max_mis_;
    max_readLen = max_readLen_;
    reference = &refFile;

    refbase2num["A"] = "00";
    refbase2num["C"] = "01";
    refbase2num["G"] = "10";
    refbase2num["T"] = "11";
    refbase2num["N"] = "00";

	var2base['G'] = { {0, 'C'}, {1, 'A'}, {2, 'T'} };
	var2base['C'] = { {0, 'G'}, {1, 'A'}, {2, 'T'} };
	var2base['A'] = { {0, 'G'}, {1, 'C'}, {2, 'T'} };
	var2base['T'] = { {0, 'G'}, {1, 'C'}, {2, 'A'} };

    compbase = {{'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'}};

    getSeg();
}

void ref2seq::getSeg() {
    block_num ++;
    if (ref_seg != "")
        ref_seg = ref_seg.substr(blockSize);
    while (true){
        getline(*reference, buffer);
        if (buffer[0] != '>'){
            for (int i=0;i<buffer.length();i++){
                if (buffer[i] >= 'a' &&  buffer[i] <= 'z')
                    buffer[i] -= 32;
                if (buffer[i] == 'N')
                    buffer[i] = 'A';
            }
            ref_seg += buffer;
            if (ref_seg.length() + buffer.length() >= blockSize+max_readLen || buffer == "")
                break;
        }
    }
}

std::string ref2seq::getSeq(align_info& info, int readl) {
    while (info.blockNum > block_num)
        getSeg();
    string ref_seq;
    ref_seq = ref_seg.substr(info.blockPos, readl);
    int offset = 0;
    for (int i=0; i< max_mis; i++){
        if (info.cigar_l[i] != -1){
            offset += info.cigar_l[i];
            ref_seq[offset+i] = var2base[ref_seq[offset+i]][info.cigar_v[i]];
        }
        else
            break;
    }
    if (info.isRev){
        reverse(ref_seq.begin(),ref_seq.end());
        for (int i=0; i < ref_seq.length();i++)
            ref_seq[i] = compbase[ref_seq[i]];
    }
    return ref_seq;
}


#endif //BWA_FASTSEQ_H
