#ifndef ALIGN_PROC_H
#define ALIGN_PROC_H

#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include <fstream>

using namespace std;

#define MaxMis 3
#define MaxInsr 512
#define MaxReadLen 256 //2的倍数，用于编码readLen

typedef struct {
    int blockNum;
    int blockPos;
    bool isRev;
    int cigar_l[MaxMis];
    int cigar_v[MaxMis];
} align_info;

class bitProc{
public:
    int bit4int(int input); //返回表示int所需的bit数
    string int2bit(int input, int byte4output); //将a转成byte4output位的二进制赋给str，即在前填0
};

int bitProc::bit4int(int input) {
    int bitNum = 1;
    while (pow(2, bitNum) < input) {
        bitNum++;
    }
    return bitNum;
}

string bitProc::int2bit(int input, int byte4output){
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

/*********** encode ***********/
class encode : public bitProc{
public:
    void parse_1(align_info align_info1, int byte4pos, int readLen, fstream &out);
    void parse_2(align_info align_info2, int readLen, fstream &out);
    void end(fstream &out);

private:
    const int max_readLen_bit = bit4int(MaxReadLen); //等于8
    const int max_insr_bit = bit4int(MaxInsr); //等于9

    bool init = true;
    int tmp_pos1, last_readLen, cigar_num;

    string buffer, cigar_l, cigar_v;

    int bit4int(int input); //返回表示int所需的bit数
    string int2bit(int input); //将a转成二进制赋给str
    string int2bit(int input, int byte4output); //将a转成byte4output位的二进制赋给str，即在前填0

    void bufferOut(string &buffer, fstream &output); //将buffer中的内容以8bit为单位输出到output
};

void encode::parse_1(align_info align_info1, int byte4pos, int readLen, fstream &out) {
    if (init){ // 第一条read
        init = false;
        last_readLen = readLen;
        buffer += int2bit(readLen, max_readLen_bit);
    }
    else{
        if (readLen == last_readLen)
            buffer += "0";
        else{
            buffer += last_readLen<readLen ? "10" : "11";
            buffer += int2bit(abs(last_readLen-readLen), MaxReadLen);
            last_readLen = readLen;
        }
    }
    tmp_pos1 = align_info1.blockNum;
    buffer += int2bit(align_info1.blockNum, byte4pos);
    buffer += align_info1.isRev;

    cigar_num = 1;
    cigar_l = cigar_v = "";
    cigar_l += int2bit(align_info1.cigar_l[0], bit4int(readLen));
    cigar_v += int2bit(align_info1.cigar_v[0], 2);

    for (int i=1;i<MaxMis;i++){
        if (align_info1.cigar_l[i] != -1){
            cigar_num ++;
            cigar_l += int2bit(align_info1.cigar_l[i]-align_info1.cigar_l[i-1], bit4int(readLen-align_info1.cigar_l[i-1]));
            cigar_v += int2bit(align_info1.cigar_v[0], 2);
        }
        else
            break;
    }
    buffer += int2bit(cigar_num, MaxMis);
    buffer += cigar_l;
    buffer += cigar_v;

    bufferOut(buffer, out);
}

void encode::parse_2(align_info align_info2, int readLen, fstream &out) {
    if (readLen == last_readLen)
        buffer += "0";
    else{
        buffer += last_readLen<readLen ? "10" : "11";
        buffer += int2bit(abs(last_readLen-readLen), MaxReadLen);
        last_readLen = readLen;
    }
    if (tmp_pos1 <=align_info2.blockNum)
        buffer += "0" + int2bit(align_info2.blockNum-tmp_pos1, max_insr_bit);
    else
        buffer += "1" + int2bit(tmp_pos1-align_info2.blockNum, max_insr_bit);
    buffer += align_info2.isRev;

    cigar_num = 1;
    cigar_l = cigar_v = "";
    cigar_l += int2bit(align_info2.cigar_l[0], bit4int(readLen));
    cigar_v += int2bit(align_info2.cigar_v[0], 2);

    for (int i=1;i<MaxMis;i++){
        if (align_info2.cigar_l[i] != -1){
            cigar_num ++;
            cigar_l += int2bit(align_info2.cigar_l[i]-align_info2.cigar_l[i-1], bit4int(readLen-align_info2.cigar_l[i-1]));
            cigar_v += int2bit(align_info2.cigar_v[0], 2);
        }
        else
            break;
    }
    buffer += int2bit(cigar_num, MaxMis);
    buffer += cigar_l;
    buffer += cigar_v;

    bufferOut(buffer, out);
}

void encode::end(fstream &out) {
    buffer += string(0, 8-buffer.length());
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
    void decode(string &buffer, fstream &in); //将in的二进制数据转成align_info赋给buffer

private:
    bool init = true;
    int tmp_pos1, last_readLen, cigar_num;

    string buffer, cigar_l, cigar_v;

    void parse_se(align_info align_info1, int byte4pos, int readLen, fstream &out);
    void parse_pe(align_info align_info2, int readLen, fstream &out);
};

#endif //ALIGN_PROC_H
