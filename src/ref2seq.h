//
//
//

#ifndef BWA_FASTSEQ_H
#define BWA_FASTSEQ_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>

class ref2seq {
public:
    int blockSize;
    std::fstream *reference;
    int block_serial = 0;
    int lgst_rl = 0;
    std::string ref_seg, buffer, tail; //tail for reads aligned to the end of each block
    std::map<std::string, std::string> base2num;

    void init(int a, std::fstream& b);
    void getSeg();
    std::string getSeq(int pos, int length);
};

void ref2seq::init(int a, std::fstream& b) {
    blockSize = a;
    reference = &b;
    base2num["A"] = "00";
    base2num["a"] = "00";
    base2num["C"] = "01";
    base2num["c"] = "01";
    base2num["G"] = "10";
    base2num["g"] = "10";
    base2num["T"] = "11";
    base2num["t"] = "11";
    base2num["N"] = "00";
    base2num["n"] = "00";
}

void ref2seq::getSeg() {
    ref_seg = tail + buffer;
    while (true){
        getline(*reference, buffer);
        if (buffer[0] != '<'){
            if (ref_seg.length() + buffer.length() <= blockSize+lgst_rl){
                //全输入
                ref_seg += buffer;
            }
            else{
                //部分输入
                ref_seg += buffer.substr(0, blockSize+lgst_rl-ref_seg.length());
                buffer = buffer.substr(blockSize+lgst_rl-ref_seg.length());
                tail = buffer.substr(blockSize);
                break;
            }
        }
    }
}

std::string ref2seq::getSeq(int pos, int length) {
    if (length > lgst_rl)
        lgst_rl = length;
    return ref_seg.substr(pos, length);
}


#endif //BWA_FASTSEQ_H
