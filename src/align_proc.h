//
// Created by cfy on 2017/6/7.
//
#ifndef DECODE2_ENCODE_H
#define DECODE2_ENCODE_H

#include <iostream>
#include <vector>
#include <map>
using namespace std;


struct result_reads
{
    string Chr_name;
    int Pos;
    string Cigar;
    int Strand_Mark;
};

struct aligent   //拆分后的数据流
{
    string Chr_name_pos[2];//染色体+状态（开始F，中间M，结尾E）
    int Pos_relative;  //相对位置，同一条染色体相对上一个reads的尾部的相对位置
    int Pieces;   //记录End的数目
    int Cigar_Mark; //是否有变异，0和1表示
    int Cigar_Num;  //记录同一条染色体中出现的第几个cigar（不带M的那种情况）
    int Position;  //记录位点
    //int reads_length;
    //int reads_Run;  //变异开始的方向为正向，记录变异的开始是否为0位置开始
    vector<int>Cigra_relative_Pos;    //保存cigar的变异相对位置，数字,0代表一个变异
    string Cigar_Variation;  //变异情况的编码
    int Strand_Mark_length;   //这个是比对的正反向，0和1表示
};

struct Chr
{
    string Mark;   //标识前后；
    int Pieces;    //个数
    int Position;   //位点
};


/*********** encode ***********/
//函数声明：
void readsIn(struct result_reads &);   //读进reads比对结果
void align_encode(int &,vector<aligent>& ,map<string,map<int,Chr> >&); //处理函数
void vec_erase(vector<aligent>&); //内存释放
void maplist(map<char,map<char,string> >&); //变异存储结构


/***********decode***********/
//函数声明
void align_decode(vector<aligent>& ,map<string,map<int,Chr> >&);   //还原reads的比对结果模式
void Reference1(map<string,string>&);

#endif //DECODE2_ENCODE_H
