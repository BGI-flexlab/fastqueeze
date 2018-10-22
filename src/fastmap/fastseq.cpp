//
//
//
#include <iostream>
#include <fstream>
#include <string>

#include "fastseq.h"
#include "encode.h"

int main(int argc, char *argv[])
{
    ifstream fin("/somewhere/example.arc");
    //先获取头部的block_size等信息

    int block_length = 1; //这个是Ref一个区块的长度，单位为碱基
    int block_num = 1000000; // 这个是单次存取的Block大小
    int pass_base = 0;
    string s;
    unsigned short tmp[block_num];
    unsigned short eightbases;

    map<string, unsigned short> base2num;
    base2num["A"] = 0;
    base2num["C"] = 1;
    base2num["G"] = 10;
    base2num["T"] = 11;
    base2num["N"] = 0;

    while(pass_base < block_length)
    {
        fin.getline(s,1);
        if(s[0] != '<'){
            for(i=0; i < s.length(); i++){
                // base2num[s[i]] 存到eightbases里 满了就转到tmp中
            }
        }

    }



}