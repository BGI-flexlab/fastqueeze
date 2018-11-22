#ifndef ALIGN_PROC_H
#define ALIGN_PROC_H

#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include <fstream>

using namespace std;

typedef struct {
    int blockNum;
    int blockPos;
    bool isRev;
    int* cigar_l;
    int* cigar_v;
} align_info;

class bitProc{
public:
    int bit4int(int input); //返回表示int所需的bit数
    string int2bit(int input, int byte4output); //将a转成byte4output位的二进制赋给str，即在前填0
    string char2bit(char input); //将压缩文件中的char转回二进制
    int bit2int(string input); //将二进制转回十进制
};

int bitProc::bit4int(int input) {
    int bitNum = 1;
    while (pow(2, bitNum) < input+1) { //consider 0
        bitNum++;
    }
    return bitNum;
}

string bitProc::int2bit(int input, int byte4output){
    if (byte4output < bit4int(input)){
        cout << input << '\t' << byte4output;
        cout << "err: bytes not long enough to represent input." << endl;
        exit(-1);
    }
    if (input==0){
        string output(byte4output, '0');
        return output;
    }
    string output(byte4output - bit4int(input), '0');
    char *p = (char *) &input, c = 0, f = 0;//p指向a的首地址
    for (int o = 0; o < 4; ++o) {
        for (int i = 0; i < 8; ++i) {
            c = p[3 - o] & (1 << (7 - i));
            if (!f && !(f = c))continue;
            output += c ? "1" : "0";
        }
    }
    return output;
}

string bitProc::char2bit(char input) {
    string output;
    for (int i = 0; i < CHAR_BIT; ++i) {
        output = to_string((long long int)((input >> i) & 1)) + output;
    }
    return output;
}

int bitProc::bit2int(string input) {
    int output = 0;
    for(int i=0;i<input.length(); i++) {
        output = output * 2 + (input[i] - '0');
    }
    return output;
}

/*********** encode ***********/
class encode : public bitProc{
public:
    encode(int se_mark, int blocksize, int MaxMis_, int MaxInsr_, int MaxReadLen_);
    void parse_1(align_info &align_info1, int readLen, fstream &out);
    void parse_2(align_info &align_info2, int readLen, fstream &out);
    void end(fstream &out);

private:
    int blocksize;
    int se_mark, byte4pos, MaxMis, MaxInsr, MaxReadLen;
    int max_readLen_bit;
    int max_insr_bit;
    int max_mis_bit;

    bool init1, init2;
    int tmp_pos1, last_readLen1, last_readLen2, cigar_num;
    string buffer, cigar_l, cigar_v;
    map<int, string> cigarv2bit;

    void bufferOut(string &buffer, fstream &output); //将buffer中的内容以8bit为单位输出到output
};

encode::encode(int se_mark_, int blocksize_, int MaxMis_, int MaxInsr_, int MaxReadLen_){
    init1 = true;
    init2 = true;
    se_mark = se_mark_;
    blocksize = blocksize_;
    byte4pos = bit4int(blocksize);
    MaxMis = MaxMis_;
    MaxInsr= MaxInsr_;
    MaxReadLen = MaxReadLen_;
    max_readLen_bit = bit4int(MaxReadLen);
    max_insr_bit = bit4int(MaxInsr);
    max_mis_bit = bit4int(MaxMis);
    cigarv2bit[0]="0"; cigarv2bit[1]="10"; cigarv2bit[2]="11";
}

void encode::parse_1(align_info &align_info1, int readLen, fstream &out) {
    if (init1){ // 第一条read
        init1 = false;
        buffer += "1"; //Block有内容
        last_readLen1 = readLen;
        buffer += int2bit(readLen, max_readLen_bit);
    }
    else{
        if (readLen == last_readLen1)
            buffer += "0";
        else{
            buffer += readLen<=last_readLen1 ? "10" : "11";
            buffer += int2bit(abs(last_readLen1-readLen), max_readLen_bit);
            last_readLen1 = readLen;
        }
    }
    tmp_pos1 = align_info1.blockPos;
    buffer += int2bit(align_info1.blockPos, byte4pos);
    buffer += align_info1.isRev ? "1" : "0";

    cigar_num = 0;
    cigar_l = cigar_v = "";
    int lastl = -1, remainLen = readLen;
    for (int i=0;i<MaxMis;i++){
        if (align_info1.cigar_l[i] != -1){
            cigar_num ++;
            cigar_l += int2bit(align_info1.cigar_l[i]-lastl-1, bit4int(remainLen));
            cigar_v += cigarv2bit[align_info1.cigar_v[i]];
            lastl = align_info1.cigar_l[i];
            remainLen = readLen - align_info1.cigar_l[i]-1;
        }
        else
            break;
    }
    buffer += int2bit(cigar_num, max_mis_bit);
    buffer += cigar_l;
    buffer += cigar_v;
    bufferOut(buffer, out);
}

void encode::parse_2(align_info &align_info2, int readLen, fstream &out) {
    if (init2){ // 第一条read
        init2 = false;
        last_readLen2 = readLen;
        if (readLen == last_readLen1)
            buffer += "0";
        else{
            buffer += readLen<=last_readLen1 ? "10" : "11";
            buffer += int2bit(abs(last_readLen1-readLen), max_readLen_bit);
        }
    }
    else{
        if (readLen == last_readLen2)
            buffer += "0";
        else{
            buffer += readLen<=last_readLen2 ? "10" : "11";
            buffer += int2bit(abs(last_readLen2-readLen), max_readLen_bit);
            last_readLen2 = readLen;
        }
    }
    if (align_info2.blockPos <= tmp_pos1)
        buffer += "0" + int2bit(tmp_pos1-align_info2.blockPos, max_insr_bit); //read2在前
    else
        buffer += "1" + int2bit(align_info2.blockPos-tmp_pos1, max_insr_bit); //read2在后
    buffer += align_info2.isRev ? "1" : "0";

    cigar_num = 0;
    cigar_l = cigar_v = "";
    int lastl = -1, remainLen = readLen;
    for (int i=0;i<MaxMis;i++){
        if (align_info2.cigar_l[i] != -1){
            cigar_num ++;
            cigar_l += int2bit(align_info2.cigar_l[i]-lastl-1, bit4int(remainLen)); //debug
            cigar_v += cigarv2bit[align_info2.cigar_v[i]];
            lastl = align_info2.cigar_l[i];
            remainLen = readLen - align_info2.cigar_l[i]-1;
        }
        else
            break;
    }
    buffer += int2bit(cigar_num, max_mis_bit);
    buffer += cigar_l;
    buffer += cigar_v;
    bufferOut(buffer, out);
}

void encode::end(fstream &out) {
    if (init1){ //block内无read
        buffer += "0";
    }
    else{
        buffer += "10";
        buffer += int2bit(last_readLen1, max_readLen_bit); //Len=0 for breakpoint
        bufferOut(buffer, out);
    }
    if (buffer.length() > 0)
        buffer += string(8-buffer.length(), '0');
    bufferOut(buffer, out);
}

void encode::bufferOut(string &buffer, fstream &output)
{
    unsigned int count = 0;
    while (buffer.length() >= (count+1)*8){
        char c='\0';
        for(int i=0;i<8;i++)
        {
            if(buffer[8*count+i]=='1') c=(c<<1)|1;
            else c=c<<1;
        }
        output.write(&c, sizeof(char));
        count++;
    }
    buffer.erase(0,8*count);
}

/***********decode***********/
class decode : public bitProc{
public:
    decode(int block_size_, int max_mis_, int max_insr_, int max_readLen_);
    int parse_se(align_info &align_info1, int& readLen_, fstream &in);
    int parse_pe(align_info &align_info1, int& readLen_, fstream &in);

private:
    bool init;
    int peMark;
    int MaxMis, MaxInsr, MaxReadLen, byte4pos, readLen1, readLen2, *readLen;
    int max_readLen_bit, max_mis_bit, max_insr_bit;
    int tmp_pos1, cigar_num, offset;
    int blockNum;

    string buffer, cigar_l, cigar_v;
    string bufferIn(int length, fstream &in); //length指的是byte数，返回的是二进制字符串
    void bufferClear(); //清空buffer
};

decode::decode(int block_size_, int max_mis_, int max_insr_, int max_readLen_)
{
    init = true;
    byte4pos = bit4int(block_size_);
    MaxMis = max_mis_;
    MaxInsr = max_insr_;
    MaxReadLen = max_readLen_;
    max_mis_bit = bit4int(MaxMis);
    max_insr_bit = bit4int(MaxInsr);
    max_readLen_bit = bit4int(MaxReadLen);
    blockNum = -1;
    peMark = readLen1 = readLen2 = 0;
}

int decode::parse_se(align_info &align_info1, int& readLen_, fstream &in){
    if (init){
        init = false;
        blockNum ++;
        if (bufferIn(1, in) == "0"){ //该block为空, 即无read比对至此
            bufferClear(); //把剩余7个0提出来
            init = true;
            return 0;
        }
        else //block不为空
            readLen_ = bit2int(bufferIn(max_readLen_bit, in));
    }
    else{
        if (bufferIn(1,in) == "1"){ //readLen发生变化
            if (bufferIn(1,in) == "0")
                readLen_ -= bit2int(bufferIn(max_readLen_bit,in));
            else
                readLen_ += bit2int(bufferIn(max_readLen_bit,in));
            if (readLen_ == 0){ //block收尾
                bufferClear();
                init = true;
                return 0;
            }
        }
    }
    align_info1.blockNum = blockNum;
    align_info1.blockPos = bit2int(bufferIn(byte4pos,in));
    align_info1.isRev = (bool)atoi(bufferIn(1,in).c_str());
    cigar_num = bit2int(bufferIn(max_mis_bit, in));
    offset = 0;
    for (int i=0; i< cigar_num; i++){
        align_info1.cigar_l[i] = bit2int(bufferIn(bit4int(readLen_-offset-i),in));
        offset += align_info1.cigar_l[i];
    }
    for (int i=0; i< cigar_num; i++){
        if (bufferIn(1,in) == "0")
            align_info1.cigar_v[i] = 0;
        else{
            if (bufferIn(1,in) == "0")
                align_info1.cigar_v[i] = 1;
            else
                align_info1.cigar_v[i] = 2;
        }
    }
    if (cigar_num < 3)
        align_info1.cigar_l[cigar_num] = -1;
    return 1;
}

int decode::parse_pe(align_info &align_info1, int& readLen_, fstream &in){
    if (init){
        init = false;
        peMark = 1; //应该可以省略
        blockNum ++;
        if (bufferIn(1, in) == "0"){ //该block为空, 即无read比对至此
            bufferClear(); //把剩余7个0提出来
            init = true;
            return 0;
        }
        else{ //block不为空
            readLen1 = bit2int(bufferIn(max_readLen_bit, in));
            tmp_pos1 = align_info1.blockPos = bit2int(bufferIn(byte4pos,in));
            readLen = &readLen1;
        }
    }
    else{
        if (peMark == 0) //read1
            readLen = &readLen1;
        else
            readLen = &readLen2;
        if (*readLen == 0){ // init read2
            if (bufferIn(1,in) == "1"){
                if (bufferIn(1,in) == "0")
                    *readLen = readLen1 - bit2int(bufferIn(max_readLen_bit,in));
                else
                    *readLen = readLen1 + bit2int(bufferIn(max_readLen_bit,in));
            }
            else
                *readLen = readLen1;
        }
        else{
            if (bufferIn(1,in) == "1"){ //readLen发生变化
                if (bufferIn(1,in) == "0"){
                    *readLen -= bit2int(bufferIn(max_readLen_bit,in));
                }
                else{
                    *readLen += bit2int(bufferIn(max_readLen_bit,in));
                }
                if (*readLen == 0){ //block收尾
                    bufferClear();
                    init = true;
                    readLen1 = readLen2 = 0;
                    return 0;
                }
            }
        }
        if (peMark == 0) { //read1
            peMark ++;
            tmp_pos1 = align_info1.blockPos = bit2int(bufferIn(byte4pos, in));
        }
        else {
            peMark --;
            if (bufferIn(1, in) == "0")
                align_info1.blockPos = tmp_pos1 - bit2int(bufferIn(max_insr_bit, in));
            else
                align_info1.blockPos = tmp_pos1 + bit2int(bufferIn(max_insr_bit, in));
        }
    }
    align_info1.blockNum = blockNum;
    align_info1.isRev = (bool)atoi(bufferIn(1,in).c_str());
    cigar_num = bit2int(bufferIn(max_mis_bit, in));
    offset = 0;
    for (int i=0; i< cigar_num; i++){
        align_info1.cigar_l[i] = bit2int(bufferIn(bit4int(*readLen-offset-i),in));
        offset += align_info1.cigar_l[i];
    }
    for (int i=0; i< cigar_num; i++){
        if (bufferIn(1,in) == "0")
            align_info1.cigar_v[i] = 0;
        else{
            if (bufferIn(1,in) == "0")
                align_info1.cigar_v[i] = 1;
            else
                align_info1.cigar_v[i] = 2;
        }
    }
    if (cigar_num < 3)
        align_info1.cigar_l[cigar_num] = -1;
    readLen_ = *readLen;
    return 1;
}

string decode::bufferIn(int length, fstream &in) { //已确认, 当length大于文件剩余bit数时不会引发error
    if (length > buffer.length()){ //buffer的内容不够，需要向in读取
        int bitNum = (int)ceil((length - buffer.length())/8.0);
        char tmpIn;
        for (int i=0; i< bitNum; i++){
            in.read(&tmpIn, 1);
            buffer += char2bit(tmpIn);
        }
    }
    string output;
    output = buffer.substr(0, length);
    buffer.erase(0, length);
    return output;
}

void decode::bufferClear(){
    buffer = "";
}

#endif //ALIGN_PROC_H
