#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <map>
#include <pthread.h>
#include <sys/time.h>
#include <queue>
#include <semaphore.h>
#include "bwa/bwa.h"
#include "bwa/bwamem.h"
#include "fqzcomp.h"
#include "SeqRead.hpp"

#ifdef __APPLE__
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

using namespace std;

#define MAJOR_VERS 0
#define MINOR_VERS 2
#define FORMAT_VERS 3
#define PREALIGN_NUM 2000

typedef struct {
    //int blockNum;
    uint64_t blockPos;
    bool isRev;
    int* cigar_l;
    int* cigar_v;
} align_info;

typedef struct _tagMagicParam{
    uint64_t both_strands:1;      // True if -b used
    uint64_t extreme_seq:1;       // True if -e used; 16-bit seq counters
    uint64_t multi_seq_model:1;   // True if -s<level>+; use 2 model sizes
    uint64_t fqzall:1;          //fqz compress all read
    uint64_t one_ch:1;           //fastq third line only +
    uint64_t isSE:1;            //se or pe
    uint64_t isGzip:1;          // True if gzip
    uint64_t slevel:4;         // -s level
    uint64_t qlevel:4;         // -q level
    uint64_t nlevel:4;         // -n level
    uint64_t qual_approx:2;      // 0 for lossless, >0 for lossy qual encoding
    uint64_t major_vers:4;
    uint64_t format_vers:4;
    uint64_t max_mis:4;
    uint64_t max_insr:10;
    uint64_t max_readLen:10;
    uint64_t thread_num:7;
} MagicParam;


typedef struct  _tagEncodeParam
{
    int num;
    uint64_t offset[2];
    uint64_t length[2];
    char filename[2][256];
    smem_i *pitr;
    bwaidx_t *pidx;
}EncodeParam;

typedef struct _tagDecodeParam
{
    int num;
    uint64_t offset;
    uint64_t length;
    char filename[256];
    bwaidx_t *pidx;
}DecodeParam;


unsigned char seq_nst_table[256] = {
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,  0, 12,  1, 14,255,255,  2, 11,255,255,  8,255,  5,  4,255,
        255,255,  6,  9,  3,255, 13, 10,255,  7,255,255,255,255,255,255,
        255,  0, 12,  1, 14,255,255,  2, 11,255,255,  8,255,  5,  4,255,
        255,255,  6,  9,  3,255, 13, 10,255,  7,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255
};

const int MAX_THREAD_NUM = 100;
int min_len = 17, max_iwidth = 50, max_mis = 3, lgst_num = 2, max_smem_num = 2, exp_mismatch = 1,  max_insr = 511, max_readLen = 255, thread_num = 1;
float min_alignratio = 0.5;
uint64_t g_lentharry[MAX_THREAD_NUM]={0};
bool g_show_warning = false; //是否开启异常碱基提示
fqz_params g_fqz_params;
MagicParam g_magicparam;
map<int, map<int, int>> nucleBitMap;
int g_offset_bit = 30; //偏移的位数
uint64_t g_offset_size  = 0; //偏移的数值
bool g_isgzip = false;

typedef struct _tagTask
{
    int name_len[2];
    int seq_len[2];
    char name[2][1024];
    char seq[2][1024];
    char qual[2][1024];
}Task;

bool g_bFinish = false;   //是否结束线程


int getbitnum(uint64_t data);

class BitDecode
{
public:
    BitDecode(){char_index = 8;}
    ~BitDecode(){};
    void setbuf(char *buf){m_pbuf = buf;}
    uint32_t bufferIn(int length);
    uint32_t getpos(bool isread2);
    uint32_t offset_pos;
    char *m_pbuf;

    int char_index;
    unsigned char rem_char;
};

uint32_t BitDecode::getpos(bool isread2)
{
	if(isread2)
	{
		int tmp = bufferIn(1); //判断read2是否大于read1的pos
		if(tmp == 1)
		{
			offset_pos += bufferIn(getbitnum(max_insr));
		}
		else
		{
			offset_pos -= bufferIn(getbitnum(max_insr));
		}
	}
	else
	{
		offset_pos = bufferIn(g_offset_bit);
	}
	return offset_pos;
}


uint32_t BitDecode::bufferIn(int length)
{
    uint32_t pass_bit = 0;
    uint32_t out_int = 0;
    
    while (pass_bit < length)
    {
        if (char_index == 8) 
        {
            rem_char = *m_pbuf++;
            char_index = 0;
        }
        if ((rem_char >> 7-char_index) & 1)
        {
            out_int += pow(2,(length-pass_bit-1));
        }
        ++pass_bit;
        ++char_index;
    }
    return out_int;
}

class ThreadTask
{
public:
    ThreadTask(int num);
    ~ThreadTask();
    void setTask(Task &task, int index);
    int getTask(Task &task, int index);
private:
    int m_num;
    pthread_mutex_t *m_pmutex;
    //std::vector<std::queue<Task>>  m_vec_queue;
    std::queue<Task> *m_pqueue;
};

ThreadTask::ThreadTask(int num):m_num(num)
{
    m_pmutex = new pthread_mutex_t[num];
    m_pqueue = new std::queue<Task>[num];
    //m_vec_queue.reserve(num);
    for(int i=0;i<num;i++)
    {
        pthread_mutex_init(&m_pmutex[i], 0);

        //std::queue<Task> queue;
        //m_vec_queue.emplace_back(queue);
    }
}

ThreadTask::~ThreadTask()
{
    for(int i=0;i<m_num;i++)
    {
        pthread_mutex_destroy(&m_pmutex[i]);
    } 

    delete[] m_pmutex;
    delete[] m_pqueue;
}

void ThreadTask::setTask(Task &task, int index)
{
    pthread_mutex_lock(&m_pmutex[index]);
    m_pqueue[index].push(task);
    pthread_mutex_unlock(&m_pmutex[index]);
}

int ThreadTask::getTask(Task &task, int index)
{
    pthread_mutex_lock(&m_pmutex[index]);
    if(m_pqueue[index].empty())
    {
        pthread_mutex_unlock(&m_pmutex[index]);
        if(g_bFinish)
        {
            //printf("thread id= %0x quit\n", pthread_self());
            return 1; //退出线程
        }
        else
        {
            return 2; //等待新任务
        }
    }
    task = m_pqueue[index].front();
    m_pqueue[index].pop();
    pthread_mutex_unlock(&m_pmutex[index]);
    return 0;
}

typedef struct _tagThreadParam
{
    int num;
    smem_i *pitr;
    bwaidx_t *pidx;
    ThreadTask *pthreadtask;
}ThreadParam;


bool bwtintv_lencmp(const bwtintv_t &arg1, const bwtintv_t &arg2) {     //长的SMEM排前面
    return (uint32_t) arg1.info - (uint32_t) (arg1.info >> 32) > (uint32_t) arg2.info - (uint32_t) (arg2.info >> 32);
}

int idx_compar( const void* a, const void* b ) {
    return ((const int*)a)[0] - ((const int*)b)[0];
}

void CreatBitmap(map<int, map<int, int>> &bitmap) {
    map<int, int> map_G;
    map<int, int> map_C;
    map<int, int> map_A;
    map<int, int> map_T;
    map_G[1] = 0; map_G[0] = 1; map_G[3] = 2; 
    map_C[2] = 0; map_C[0] = 1; map_C[3] = 2; 
    map_A[2] = 0; map_A[1] = 1; map_A[3] = 2; 
    map_T[2] = 0; map_T[1] = 1; map_T[0] = 2;
    bitmap[2] = map_G;
    bitmap[1] = map_C;
    bitmap[0] = map_A;
    bitmap[3] = map_T;
}

int getAlignInfoSE(char *seq, bwtintv_v *a, bwtintv_v *tmpvec[2], smem_i* func_itr, bwaidx_t *func_idx, align_info &align_info1, char *degenerate){
    int64_t rlen;
    int i, seql, base;
    int degenerate_num = 0;

    seql = strlen(seq);
    int pass_num = 0;
    int cigar_l[max_mis], cigar_v[max_mis];

    if(g_show_warning)
    {
        for (i = 0; i < seql; ++i) 
        {
            if(seq_nst_table[(int)seq[i]] == 255)
            {
                printf("WARNING: the %c is abnormal base\n", seq[i]);
            } 
        }
    }

    unsigned char seqarry[seql];
    for (i = 0; i < seql; ++i) {
        seqarry[i] = nst_nt4_table[(int)seq[i]];
    }

    int start = 0;
    vector <bwtintv_t> pvec;
    while ((start = smem_next_t(func_itr, start, seql, seqarry, a, tmpvec)) != 0) {
        for (i = 0; i < a->n; ++i) {
            bwtintv_t *p = &a->a[i];
            if ((uint32_t) p->info - (uint32_t)(p->info >> 32) >= min_len)
                pvec.emplace_back(a->a[i]);
        }
    }
    if (pvec.size() == 0)
        return 0;
    if (pvec.size() > 1)
        sort(pvec.begin(), pvec.end(), bwtintv_lencmp);

    int count = max_smem_num < pvec.size() ? max_smem_num: pvec.size();
    for (i = 0; i < count; ++i) {
        if (pvec[i].x[2] <= max_iwidth) {
            for (int k = 0; k < pvec[i].x[2]; ++k) {
                bwtint_t pos;
                int is_rev;
                pos = (bwtint_t)bns_depos(func_idx->bns, bwt_sa(func_idx->bwt, pvec[i].x[0] + k), &is_rev);
                uint8_t *rseq_l, *rseq_r;
                uint16_t missum = 0, last_missum = max_mis+1;
                int rbase, sbase;
                if (is_rev) {
                    pos -= seql - (uint32_t)(pvec[i].info >>32)-1;
                    rseq_l = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos,
                                         pos + seql - (uint32_t) (pvec[i].info), &rlen);
                    //把rseq_l和read的右截的反向互补进行比较
                    for (base = 0; base < rlen; base++) {
                        if (rseq_l[base] + seqarry[seql-base-1] != 3) {
                            if (missum == max_mis){
                                ++missum;
                                break;
                            }
                            rbase = int(rseq_l[base]);
                            sbase = int(seqarry[seql-base-1]);
                            cigar_l[missum] = base;
                            if(sbase>3)
                            {
                                cigar_v[missum] = 3;
                                degenerate[degenerate_num++] = seq[seql-base-1];
                            }
                            else
                            {
                                cigar_v[missum] = nucleBitMap[rbase][3-sbase];
                            }
                            
                            ++missum;
                        }
                    }
                    free(rseq_l);
                    if (missum <= max_mis) {
                        rseq_r = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos + seql - (uint32_t)(pvec[i].info >>32),
                                             pos + seql, &rlen);
                        //把rseq_r和read的左截的反向互补进行比较
                        for (base = 0; base < rlen; base++) {
                            if (rseq_r[base] + seqarry[(uint32_t)(pvec[i].info >>32)-base-1] != 3) {
                                if (missum == max_mis){
                                    ++missum;
                                    break;
                                }
                                rbase = int(rseq_r[base]);
                                sbase = int(seqarry[(uint32_t)(pvec[i].info >>32)-base-1]);
                                cigar_l[missum] = seql-(uint32_t)(pvec[i].info >>32)+base;
                                if(sbase > 3)
                                {
                                    cigar_v[missum] = 3;
                                    degenerate[degenerate_num++] = seq[(uint32_t)(pvec[i].info >>32)-base-1];
                                }
                                else
                                {
                                    cigar_v[missum] = nucleBitMap[rbase][3-sbase];
                                }
                                
                                ++missum;
                            }
                        }
                        free(rseq_r);
                    }
                    if (missum <= max_mis){
                        if (missum < last_missum){
                            if(pos+seql > func_idx->bns->l_pac)//超出了ref范围
                            {
                                return 0;
                            }
                            align_info1.blockPos = pos;
                            align_info1.isRev    = (bool) is_rev;
                            for (int x=0; x<missum; x++){
                                align_info1.cigar_l[x] = cigar_l[x];
                                align_info1.cigar_v[x] = cigar_v[x];
                            }
                            if (missum < max_mis)
                                align_info1.cigar_l[missum] = -1;
                            last_missum = missum;
                        }
                        ++pass_num;
                        if (missum <= exp_mismatch || pass_num >= lgst_num)
                            return pass_num;
                    }
                } else {
                    pos -= (uint32_t) (pvec[i].info >> 32);
                    rseq_l = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos,
                                         pos+(uint32_t) (pvec[i].info >> 32), &rlen);
                    //把rseq_l和read的左截进行比较
                    for (base = 0; base < rlen; base++) {
                        if (rseq_l[base] != seqarry[base]) {
                            if (missum == max_mis){
                                missum += 1;
                                break;
                            }
                            rbase = int(rseq_l[base]);
                            sbase = int(seqarry[base]);
                            cigar_l[missum] = base;
                            if(sbase>3)
                            {
                                cigar_v[missum] = 3;
                                degenerate[degenerate_num++] = seq[base];
                            }
                            else
                            {
                                cigar_v[missum] = nucleBitMap[rbase][sbase];
                            }
                            missum += 1;
                        }
                    }
                    free(rseq_l);
                    if (missum <= max_mis) {
                        rseq_r = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos+ (uint32_t)(pvec[i].info),
                                             pos + seql, &rlen);
                        //把rseq_r和read的右截进行比较
                        for (base = 0; base < rlen; base++) {
                            if (rseq_r[base] != seqarry[(uint32_t) pvec[i].info + base]) {
                                if (missum == max_mis){
                                    missum += 1;
                                    break;
                                }
                                rbase = int(rseq_r[base]);
                                sbase = int(seqarry[(uint32_t) pvec[i].info + base]);
                                cigar_l[missum] = base+(uint32_t)(pvec[i].info);
                                if(sbase>3)
                                {
                                    cigar_v[missum] = 3;
                                    degenerate[degenerate_num++] = seq[(uint32_t) pvec[i].info + base];
                                }
                                else
                                {
                                    cigar_v[missum] = nucleBitMap[rbase][sbase];
                                }
                                
                                missum += 1;
                            }
                        }
                        free(rseq_r);
                    }
                    if (missum <= max_mis){
                        if (missum < last_missum){
                            if(pos+seql > func_idx->bns->l_pac)//超出了ref范围
                            {
                                return 0;
                            }
                            align_info1.blockPos = pos;
                            align_info1.isRev    = (bool) is_rev;
                            for (int x=0; x<missum; x++){
                                align_info1.cigar_l[x] = cigar_l[x];
                                align_info1.cigar_v[x] = cigar_v[x];
                            }
                            if (missum < max_mis)
                                align_info1.cigar_l[missum] = -1;
                            last_missum = missum;
                        }
                        ++pass_num;
                        if (missum <= exp_mismatch || pass_num >= lgst_num)
                            return pass_num;
                    }
                }
            }
        }
    }

    if (pass_num)
        return pass_num;
    else
        return 0;
}

int getAlignInfoPE(char *seq1, char * seq2, bwtintv_v *a1, bwtintv_v *a2, bwtintv_v *tmpvec[2], smem_i* func_itr, bwaidx_t *func_idx, 
    align_info &align_info1, align_info &align_info2, int64_t ** idx1, int64_t ** idx2, char (*degenerate)[512]){
    int64_t rlen;
    int i, j, seql1, seql2, base, count;
    int* seql_p;

    seql1 = strlen(seq1);
    seql2 = strlen(seq2);

    if(g_show_warning)
    {
        for (i = 0; i < seql1; ++i) 
        {
            if(seq_nst_table[(int)seq1[i]] == 255)
            {
                printf("WARNING: seq1 the %c is abnormal base\n", seq1[i]);
            } 
        }

        for (i = 0; i < seql2; ++i) 
        {
            if(seq_nst_table[(int)seq2[i]] == 255)
            {
                printf("WARNING: seq2 the %c is abnormal base\n", seq2[i]);
            } 
        }
    }


    unsigned char seqarry1[seql1], seqarry2[seql2];
    unsigned char* arry_p;
    for (i = 0; i < seql1; ++i) {
        seqarry1[i] = nst_nt4_table[(int) seq1[i]];
    }
    for (i = 0; i < seql2; ++i) {
        seqarry2[i] = nst_nt4_table[(int) seq2[i]];
    }

    vector <bwtintv_t> pvec1;
    vector <bwtintv_t> pvec2;
    int start = 0, is_rev, k1 = 0, k2 = 0;
    bwtint_t pos;

    while ((start = smem_next_t(func_itr, start, seql1, seqarry1, a1, tmpvec)) != 0) {
        for (i = 0; i < a1->n; ++i) {
            bwtintv_t *p = &a1->a[i];
            if ((uint32_t) p->info - (uint32_t)(p->info >> 32) >= min_len){
                if (p->x[2] > max_iwidth)
                    break;
                pvec1.emplace_back(a1->a[i]);
            }
        }
    }
    if (pvec1.size() == 0)
        return 0;
    if (pvec1.size() > 1)
        sort(pvec1.begin(), pvec1.end(), bwtintv_lencmp);

    count = max_smem_num < pvec1.size() ? max_smem_num: pvec1.size();
    for (i = 0; i < count; ++i){
        for (j = 0; j < pvec1[i].x[2]; j++){
            pos = (bwtint_t)bns_depos(func_idx->bns, bwt_sa(func_idx->bwt, pvec1[i].x[0] + j), &is_rev);
            if (is_rev)
                pos -= seql1 - (uint32_t)(pvec1[i].info >>32)-1;
            else
                pos -= (uint32_t) (pvec1[i].info >> 32);
            idx1[k1][0] = pos;
            idx1[k1][1] = i;
            idx1[k1][2] = is_rev;
            ++k1;
        }
    }
    if (k1 > 1)
        qsort(idx1, k1, sizeof(idx1[0]), &idx_compar);

    start = 0;
    while ((start = smem_next_t(func_itr, start, seql2, seqarry2, a2, tmpvec)) != 0) {
        for (i = 0; i < a2->n; ++i) {
            bwtintv_t *p = &a2->a[i];
            if ((uint32_t) p->info - (uint32_t)(p->info >> 32) >= min_len){
                if (p->x[2] > max_iwidth)
                    break;
                pvec2.emplace_back(a2->a[i]);
            }
        }
    }
    if (pvec2.size() == 0)
        return 0;
    if (pvec2.size() > 1)
        sort(pvec2.begin(), pvec2.end(), bwtintv_lencmp);
    count = max_smem_num < pvec2.size() ? max_smem_num: pvec2.size();
    for (i = 0; i < count; ++i) {
        for (j = 0; j < pvec2[i].x[2]; j++){
            pos = (bwtint_t)bns_depos(func_idx->bns, bwt_sa(func_idx->bwt, pvec2[i].x[0] + j), &is_rev);
            if (is_rev)
                pos -= seql2 - (uint32_t)(pvec2[i].info >>32)-1;
            else
                pos -= (uint32_t) (pvec2[i].info >> 32);
            idx2[k2][0] = pos;
            idx2[k2][1] = i;
            idx2[k2][2] = is_rev;
            ++k2;
        }
    }
    if (k2 > 1)
        qsort(idx2, k2, sizeof(idx2[0]), &idx_compar);

    i = j = 0;
    while (1){
        if (idx1[i][0] < idx2[j][0]){
            if (idx2[j][0]-idx1[i][0] <= max_insr){
                    break;
            }
            else{
                ++i;
                if (i >= k1)
                    break;
            }
        }
        else{
            if (idx1[i][0]-idx2[j][0] <= max_insr){
                    break;
            }
            else{
                ++j;
                if (j >= k2)
                    break;
            }
        }
    }
    if (i>=k1 || j >= k2)
        return 0;

    uint8_t *rseq_l, *rseq_r;
    bwtint_t tmp_info;
    int rbase, sbase;
    pos = idx1[i][0];
    tmp_info = pvec1[idx1[i][1]].info;
    is_rev = idx1[i][2];
    arry_p = seqarry1;
    seql_p = &seql1;
    align_info* aligninfo_p;
    aligninfo_p = &align_info1;
    char *pseq = seq1;

    for (int k=0;k<2;++k) {
        uint16_t missum = 0;
        int degenerate_num = 0;
        if (is_rev) {
            rseq_l = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos,
                                 pos + *seql_p - (uint32_t) (tmp_info), &rlen);
            //把rseq_l和read的右截的反向互补进行比较
            for (base = 0; base < rlen; base++) {
                if (rseq_l[base] + arry_p[*seql_p - base - 1] != 3) {
                    if (missum == max_mis) {
                        ++missum;
                        break;
                    }
                    rbase = int(rseq_l[base]);
                    sbase = int(arry_p[*seql_p - base - 1]);
                    aligninfo_p->cigar_l[missum] = base;
                    if(sbase>3)
                    {
                        aligninfo_p->cigar_v[missum] = 3;
                        degenerate[k][degenerate_num++] = pseq[*seql_p - base - 1];
                    }
                    else
                    {
                        aligninfo_p->cigar_v[missum] = nucleBitMap[rbase][3 - sbase];
                    }
                    ++missum;
                }
            }
            free(rseq_l);
            if (missum <= max_mis) {
                rseq_r = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos + *seql_p - (uint32_t) (tmp_info >> 32),
                                     pos + *seql_p, &rlen);
                //把rseq_r和read的左截的反向互补进行比较
                for (base = 0; base < rlen; base++) {
                    if (rseq_r[base] + arry_p[(uint32_t) (tmp_info >> 32) - base - 1] != 3) {
                        if (missum == max_mis) {
                            ++missum;
                            break;
                        }
                        rbase = int(rseq_r[base]);
                        sbase = int(arry_p[(uint32_t) (tmp_info >> 32) - base - 1]);
                        aligninfo_p->cigar_l[missum] = *seql_p - (uint32_t) (tmp_info >> 32) + base;
                        if(sbase>3)
                        {
                            aligninfo_p->cigar_v[missum] = 3;
                            degenerate[k][degenerate_num++] = pseq[(uint32_t) (tmp_info >> 32) - base - 1];
                        }
                        else
                        {
                            aligninfo_p->cigar_v[missum] = nucleBitMap[rbase][3 - sbase];
                        }
                        
                        ++missum;
                    }
                }
                free(rseq_r);
            }
            if (missum <= max_mis) {
                if(pos+*seql_p > func_idx->bns->l_pac)//超出了ref范围
                {
                    return 0;
                }
                aligninfo_p->blockPos = pos;
                aligninfo_p->isRev = (bool) is_rev;
                if (missum < max_mis)
                    aligninfo_p->cigar_l[missum] = -1;
            } else {
                return 0;
            }
        } else {
            rseq_l = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos,
                                 pos + (uint32_t) (tmp_info >> 32), &rlen);
            //把rseq_l和read的左截进行比较
            for (base = 0; base < rlen; base++) {
                if (rseq_l[base] != arry_p[base]) {
                    if (missum == max_mis) {
                        missum += 1;
                        break;
                    }
                    rbase = int(rseq_l[base]);
                    sbase = int(arry_p[base]);
                    aligninfo_p->cigar_l[missum] = base;
                    if(sbase>3)
                    {
                        aligninfo_p->cigar_v[missum] = 3;
                        degenerate[k][degenerate_num++] = pseq[base];
                    }
                    else
                    {
                        aligninfo_p->cigar_v[missum] = nucleBitMap[rbase][sbase];
                    }
                    
                    missum += 1;
                }
            }
            free(rseq_l);
            if (missum <= max_mis) {
                rseq_r = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos + (uint32_t) (tmp_info),
                                     pos + *seql_p, &rlen);
                //把rseq_r和read的右截进行比较
                for (base = 0; base < rlen; base++) {
                    if (rseq_r[base] != arry_p[(uint32_t) tmp_info + base]) {
                        if (missum == max_mis) {
                            missum += 1;
                            break;
                        }
                        rbase = int(rseq_r[base]);
                        sbase = int(arry_p[(uint32_t) tmp_info + base]);
                        aligninfo_p->cigar_l[missum] = base + (uint32_t) (tmp_info);
                        if(sbase>3)
                        {
                            aligninfo_p->cigar_v[missum] = 3;
                            degenerate[k][degenerate_num++] = pseq[(uint32_t) tmp_info + base];
                        }
                        else
                        {
                            aligninfo_p->cigar_v[missum] = nucleBitMap[rbase][sbase];
                        }
                        
                        missum += 1;
                    }
                }
                free(rseq_r);
            }
            if (missum <= max_mis) {
                if(pos+*seql_p > func_idx->bns->l_pac)//超出了ref范围
                {
                    return 0;
                }
                aligninfo_p->blockPos = pos;
                aligninfo_p->isRev = (bool) is_rev;
                if (missum < max_mis)
                    aligninfo_p->cigar_l[missum] = -1;
            } else {
                return 0;
            }
        }
        if (k==0){
            pos = idx2[j][0];
            tmp_info = pvec2[idx2[j][1]].info;
            is_rev = idx2[j][2];
            arry_p = seqarry2;
            seql_p = &seql2;
            aligninfo_p = &align_info2;
            pseq = seq2;
        }
    }
    return 1;
}

uint64_t GetFileSize(char *path)
{
    uint64_t flength = 0;
    struct stat statbuf; 
    int ret = stat(path, &statbuf); 
    if(ret == -1) {
        return 0;
    }

    flength = statbuf.st_size;
    return flength;
}

bool DoPreAlign(smem_i* itr, bwaidx_t *idx, bool isSE, char *file1, uint64_t flength1, char *file2, uint64_t flength2)
{
    int i;
    SeqRead ssread1(file1, 0, flength1);

    align_info align_info1;
    align_info1.cigar_l = (int*)malloc(max_mis * sizeof(int));
    align_info1.cigar_v = (int*)malloc(max_mis * sizeof(int));
    bwtintv_v *matcher1 = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));

    int alignedNum = 0, prealign_num = PREALIGN_NUM; //set prealign_num for files with reads less than PREALIGN_NUM
    int seqm = 0;
    if(isSE)
    {
        for (i = 0; i < PREALIGN_NUM; i++) 
        {
            if(ssread1.getRead())
            {
                if(g_magicparam.one_ch && strlen(ssread1.name2) > 2)
                {
                    g_magicparam.one_ch = false;
                }
                char degenerate[512]={0};
                seqm = getAlignInfoSE(ssread1.seq, matcher1, NULL, itr, idx, align_info1,degenerate);
                if (seqm > 0) ++alignedNum;
            }
            else
            {
                prealign_num = i;
                break;
            }
        }
    }
    else
    {
        SeqRead ssread2(file2, 0, flength2);
        align_info align_info2;
        align_info2.cigar_l = (int*)malloc(max_mis * sizeof(int));
        align_info2.cigar_v = (int*)malloc(max_mis * sizeof(int));
        bwtintv_v *matcher2 = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
        int64_t ** idx1 = new int64_t*[max_smem_num*max_iwidth];
        for(i=0; i<max_smem_num*max_iwidth; i++)
            idx1[i] = new int64_t[3];
        int64_t ** idx2 = new int64_t*[max_smem_num*max_iwidth];
        for(i=0; i<max_smem_num*max_iwidth; i++)
            idx2[i] = new int64_t[3];


        for (i = 0; i < PREALIGN_NUM; i++) 
        {
            if(ssread1.getRead() && ssread2.getRead())
            {
                if(g_magicparam.one_ch && strlen(ssread1.name2) > 2)
                {
                    g_magicparam.one_ch = false;
                }

                char degenerate[2][512]={0};
                seqm = getAlignInfoPE(ssread1.seq, ssread2.seq, matcher1, matcher2, NULL, itr, idx, align_info1, align_info2, idx1, idx2, degenerate);
                if (seqm > 0) ++alignedNum;
            }
            else
            {
                prealign_num = i;
                break;
            }
        }

        free(matcher2->a);
        free(matcher2);
        free(align_info2.cigar_l);
        free(align_info2.cigar_v);
        for(int i=0; i<max_smem_num*max_iwidth; i++)
            delete[] idx1[i];
        delete[] idx1;
        for(int i=0; i<max_smem_num*max_iwidth; i++)
            delete[] idx2[i];
        delete[] idx2;
    }

    free(matcher1->a);
    free(matcher1);
    free(align_info1.cigar_l);
    free(align_info1.cigar_v);


    if (alignedNum < (prealign_num * min_alignratio)){
        printf("Compression begins. Only fqzcomp adopted.\n");
        return true;
    }
    return false;
}

int getbitnum(uint64_t data)
{
    int count = 0;
    while (data > 0) {
        data >>= 1;
        count++;
    }
    return count;
}

int myint2bit(uint64_t data, int limit, std::string &str)
{
    int count = 0;
    char buf[100]={0};
    while (data > 0) {
        if (data & 1) {
            buf[count] = '1';
        }
        else
        {
            buf[count] = '0';
        }

        data >>= 1;
        count++;
    }
    
    assert(limit>=count);

    str.clear();
    str.append(limit-count,'0');

    for(int i=count-1;i>=0;i--)
    {
        str.push_back(buf[i]);
    }
    return count;
}

void Cigarv2bit(std::string &buf, int ch)
{
    switch(ch)
    {
        case 0: buf.append("00"); break;
        case 1: buf.append("01"); break;
        case 2: buf.append("10"); break;
        case 3: buf.append("11"); break;
    }
}

void ToBitArry(align_info &info, int readLen, string &bitbuf)
{
	std::string stremp;
	if(info.isRev)
    {
        bitbuf.push_back('1');
    }
    else
    {
        bitbuf.push_back('0');
    }

    int cigar_num = 0;
    int lastl = 0, remainLen = readLen;
    string str_l,str_v;
    for (int i=0;i<max_mis;i++)
    {
        if (info.cigar_l[i] != -1){
            cigar_num ++;
            myint2bit(info.cigar_l[i]-lastl, getbitnum(remainLen), stremp);
            str_l.append(stremp);
            Cigarv2bit(str_v, info.cigar_v[i]);
            lastl = info.cigar_l[i];
            remainLen = readLen - info.cigar_l[i];
        }
        else
            break;
    }

    // printf("%ld %d %d\n", info.blockPos, info.isRev, cigar_num);
    // printf("%d %d %d\n",info.cigar_l[0], info.cigar_l[1],info.cigar_l[2]);
    // printf("%d %d %d\n",info.cigar_v[0], info.cigar_v[1],info.cigar_v[2]);

    myint2bit(cigar_num, 2, stremp);
    bitbuf.append(stremp);
    bitbuf.append(str_l);
    bitbuf.append(str_v);
    //printf("---%s\n", bitbuf.c_str());
}

int AlignInfoToBitArry_SE(align_info &info, int readLen, string &bitbuf)
{
    bitbuf.clear();
    std::string stremp;
    int index = info.blockPos>>g_offset_bit;
    uint64_t offset = info.blockPos & g_offset_size;
    myint2bit(offset, g_offset_bit, stremp);
    bitbuf.append(stremp);

    ToBitArry(info, readLen, bitbuf);
    return index;
}


int AlignInfoToBitArry_PE(align_info &info1, align_info &info2,int readLen, string &bitbuf1, string &bitbuf2)
{
    bitbuf1.clear();
    bitbuf2.clear();
    std::string stremp;
    int index = info1.blockPos>>g_offset_bit;
    uint64_t offset = info1.blockPos & g_offset_size;
    myint2bit(offset, g_offset_bit, stremp);
    bitbuf1.append(stremp);

    ToBitArry(info1, readLen, bitbuf1);
    
    int tmp = 0;
    if (info2.blockPos <= info1.blockPos)
    {
    	bitbuf2.push_back('0');
        tmp = info1.blockPos-info2.blockPos;
    }
    else
    {
    	bitbuf2.push_back('1');
        tmp = info2.blockPos-info1.blockPos;
    }
    
    myint2bit(tmp, getbitnum(max_insr), stremp);
    bitbuf2.append(stremp);

    ToBitArry(info2, readLen, bitbuf2);
    return index;    
}

void BitArryToBuf(const char *pbit, int len, unsigned char *pbuf, int *buflen)
{
    //printf("%s\n", pbit);
    int count = 0;
    while(len >= 8)
    {
        char c='\0';
        for(int i=0;i<8;i++)
        {
            if(*pbit++=='1') c=(c<<1)|1;
            else c=c<<1;
        }

        *pbuf++ = c;
        len -= 8;
        count ++;
    }

    if(len > 0) //剩下的补0
    {
        char c='\0';
        for(int i=0;i<len;i++)
        {
            if(*pbit++=='1') c=(c<<1)|1;
            else c=c<<1;
        }

        for(int i=len;i<8;i++)
        {
            c=c<<1;
        }

        *pbuf++ = c;
        count ++;
    }
    *buflen = count;
}

void *fqzall_encode_process(void *data)
{
    if(data == NULL) return data;

    EncodeParam *param = (EncodeParam*)data;

    fqz *pfqz = new fqz(&g_fqz_params);

    fstream out_isq; 
    char str_tmp[64]={0};
    sprintf(str_tmp, "./out_isq_%d.tmp",param->num);
    out_isq.open(str_tmp, std::ios::binary|std::ios::out);
    
    SeqRead ssread(param->filename[0], param->offset[0], param->length[0]);
    if(g_magicparam.isSE)
    {
        while(ssread.getRead())
        {
        	uint32_t name_len = strlen(ssread.name);
            uint32_t qual_len = strlen(ssread.qual);

            uint32_t tmplen = pfqz->getInLen() + name_len + qual_len + qual_len;
            if(tmplen > BLK_SIZE)
            {
            	pfqz->isq_doencode(out_isq);
            }

            pfqz->isq_addbuf(ssread.name, name_len, ssread.seq, qual_len, ssread.qual, qual_len);
        }
        pfqz->isq_doencode(out_isq);
    }
    else
    {
        SeqRead ssread2(param->filename[1], param->offset[1], param->length[1]);
        while(ssread.getRead() && ssread2.getRead())
        {
        	uint32_t name_len1 = strlen(ssread.name);
            uint32_t qual_len1 = strlen(ssread.qual);
            uint32_t name_len2 = strlen(ssread2.name);
            uint32_t qual_len2 = strlen(ssread2.qual);

        	uint32_t tmplen = pfqz->getInLen() + name_len1 + qual_len1 + name_len2 + qual_len2 + qual_len1 + qual_len2;//保证两条在一个block中
            if(tmplen > BLK_SIZE)
            {
            	pfqz->isq_doencode(out_isq);
            }
            pfqz->isq_addbuf(ssread.name, name_len1, ssread.seq, qual_len1, ssread.qual, qual_len1);
            pfqz->isq_addbuf(ssread2.name, name_len2, ssread2.seq, qual_len2, ssread2.qual, qual_len2);
        }
        pfqz->isq_doencode(out_isq);

        if(strlen(ssread.name) == 0) //f1文件已经读完
        {
            pfqz->isq_addmark(2); //添加f2分隔标识
            while(ssread2.getRead())
            {
            	uint32_t name_len2 = strlen(ssread2.name);
            	uint32_t qual_len2 = strlen(ssread2.qual);
            	uint32_t tmplen = pfqz->getInLen() + name_len2 + qual_len2 + qual_len2;
                if(tmplen > BLK_SIZE)
                {
                	pfqz->isq_doencode(out_isq);
                }
                pfqz->isq_addbuf(ssread2.name, name_len2, ssread2.seq, qual_len2, ssread2.qual, qual_len2);
            }
            pfqz->isq_doencode(out_isq);
        }
        else if(strlen(ssread2.name) == 0) //f2文件已经读完
        {
            pfqz->isq_addmark(1); //添加f1分隔标识
            do
            {
            	uint32_t name_len1 = strlen(ssread.name);
            	uint32_t qual_len1 = strlen(ssread.qual);
            	uint32_t tmplen = pfqz->getInLen() + name_len1 + qual_len1 + qual_len1;
                if(tmplen > BLK_SIZE)
                {
                	pfqz->isq_doencode(out_isq);
                }

                pfqz->isq_addbuf(ssread.name, name_len1, ssread.seq, qual_len1, ssread.qual, qual_len1);
            }while(ssread.getRead());
            pfqz->isq_doencode(out_isq);
        }
    }

    g_lentharry[param->num] = pfqz->getCompressTotalLen();
    delete param;
    delete pfqz;
    out_isq.close();
}

uint32_t parse_decode_buf(char *pbuf, char *namebuf, char *seqbuf, char *qualbuf, uint16_t *seqlen, int count)
{
    char *pstart = pbuf;
    for (int k = 0; k < count; k++) 
    {
        while ((*pbuf++ = *namebuf++) != '\n'); //id
        memcpy(pbuf, seqbuf, seqlen[k]); seqbuf += seqlen[k]; pbuf+=seqlen[k];//seq
        *pbuf++ = '\n';
        *pbuf++ = '+';
        *pbuf++ = '\n';
        memcpy(pbuf, qualbuf, seqlen[k]);qualbuf += seqlen[k]; pbuf+=seqlen[k];//qual
        *pbuf++ = '\n';
    }
    return pbuf - pstart;
}

uint32_t parse_decode_buf_more(char *pbuf, char *namebuf, char *seqbuf, char *qualbuf, uint16_t *seqlen, int count)
{
    char *pstart = pbuf;
    for (int k = 0; k < count; k++) 
    {
    	char *pname = namebuf;
        while ((*pbuf++ = *namebuf++) != '\n'); //id
        int id_len = namebuf - pname -1;
        
        memcpy(pbuf, seqbuf, seqlen[k]); seqbuf += seqlen[k]; pbuf+=seqlen[k];//seq

        *pbuf++ = '\n';
        *pbuf++ = '+';
        memcpy(pbuf, pname+1, id_len);pbuf += id_len;
        memcpy(pbuf, qualbuf, seqlen[k]);qualbuf += seqlen[k]; pbuf+=seqlen[k];//qual
        *pbuf++ = '\n';
    }
    return pbuf - pstart;
}

void fqzall_decode_SE(DecodeParam *param)
{
    fqz *pfqz = new fqz(&g_fqz_params);
    fstream in_isq;
    in_isq.open(param->filename, std::ios::binary|std::ios::in);
    in_isq.seekg(param->offset);

    fstream f1; //解压输出到文件
    char tmpath[125]={0};
    sprintf(tmpath,"./decode1_%d.tmp", param->num);
    f1.open(tmpath, std::ios::binary|std::ios::out);


    char *write_buf = NULL;
    uint32_t write_buf_len = 0;
    int read_len = 0;
    uint64_t total_len = param->length;
    while(true)
    {
        char *namebuf = NULL;
        char *seqbuf = NULL;
        char *qualbuf = NULL;
        uint16_t *seqlen = NULL;
        int count = 0;
        int mark = 0;
        if(total_len>0)
        {
        	read_len = pfqz->isq_decode(in_isq, &namebuf, &seqbuf, &qualbuf, &seqlen, &count, &mark);
        	total_len -= read_len;
        }
        else
        {
        	break;
        }

        uint32_t tmplen = count*5*seqlen[0];
        if(write_buf_len < tmplen)
        {
            write_buf_len = tmplen;
            write_buf = (char*)realloc(write_buf, write_buf_len);
        }

        uint32_t out_ind = 0;
        if(g_magicparam.one_ch)
        {
            out_ind = parse_decode_buf(write_buf, namebuf, seqbuf, qualbuf, seqlen, count);
        }
        else
        {
            out_ind = parse_decode_buf_more(write_buf, namebuf, seqbuf, qualbuf, seqlen, count);
        }

        f1.write(write_buf, out_ind);
    }
    
    f1.close();
    in_isq.close();
    delete pfqz;
    free(write_buf);
}

void fqzall_decode_PE(DecodeParam *param)
{
    fqz *pfqz = new fqz(&g_fqz_params);
    fstream in_isq;
    in_isq.open(param->filename, std::ios::binary|std::ios::in);
    in_isq.seekg(param->offset);

    fstream f1; //解压输出到文件
    char tmpath[125]={0};
    sprintf(tmpath,"./decode1_%d.tmp", param->num);
    f1.open(tmpath, std::ios::binary|std::ios::out);

    fstream f2;
    memset(tmpath,0,125);
    sprintf(tmpath,"./decode2_%d.tmp", param->num);
    f2.open(tmpath, std::ios::binary|std::ios::out);


    char *write_buf1 = NULL;
    char *write_buf2 = NULL;
    uint32_t write_buf_len = 0;
    int read_len = 0;
    uint64_t total_len = param->length;
    while(true)
    {
        char *namebuf = NULL;
        char *seqbuf = NULL;
        char *qualbuf = NULL;
        uint16_t *seqlen = NULL;
        int count = 0;
        int mark = 0;
        if(total_len>0)
        {
        	read_len = pfqz->isq_decode(in_isq, &namebuf, &seqbuf, &qualbuf, &seqlen, &count, &mark);
        	total_len -= read_len;
        }
        else
        {
        	break;
        }
//printf("fqzall_decode_PE  %d %d \n", count, read_len);
        int num = (mark != 0)? count:count/2;
        uint32_t tmplen = num*5*seqlen[0];
        if(write_buf_len < tmplen)
        {
            write_buf_len = tmplen;
            write_buf1 = (char*)realloc(write_buf1, write_buf_len);
            write_buf2 = (char*)realloc(write_buf2, write_buf_len);
        }

        uint32_t out_ind1 = 0;
        uint32_t out_ind2 = 0;

        if(mark == 1)//写入到f1
        {
        	if(g_magicparam.one_ch)
        	{
        		out_ind1 = parse_decode_buf(write_buf1, namebuf, seqbuf, qualbuf, seqlen, count);
        	}
        	else
        	{
        		out_ind1 = parse_decode_buf_more(write_buf1, namebuf, seqbuf, qualbuf, seqlen, count);
        	}
            f1.write(write_buf1, out_ind1);
        }
        else if(mark == 2)//写入到f2
        {
        	if(g_magicparam.one_ch)
        	{
				out_ind2 = parse_decode_buf(write_buf2, namebuf, seqbuf, qualbuf, seqlen, count);
        	}
        	else
        	{
        		out_ind2 = parse_decode_buf_more(write_buf2, namebuf, seqbuf, qualbuf, seqlen, count);
        	}
            f2.write(write_buf2, out_ind2);
        }
        else //分别写入f1和f2
        {
            for (int k = 0; k < count; k++) 
            {
                char *pbuf = NULL;
                char isread2 = false;
                if(k&1) //第二条read
                {
                    pbuf = &write_buf2[out_ind2];
                    isread2 = true;
                }
                else //第一条read
                {
                    pbuf = &write_buf1[out_ind1];
                }

                char *pame = namebuf;
                while ((*pbuf++ = *namebuf++) != '\n'); //id
                int id_len = namebuf - pame -1;
                memcpy(pbuf, seqbuf, seqlen[k]); seqbuf += seqlen[k]; pbuf+=seqlen[k];//seq
                *pbuf++ = '\n';
                *pbuf++ = '+';
                if(!g_magicparam.one_ch)
                {
                	memcpy(pbuf, pame+1, id_len); pbuf += id_len;
                }
                else
                {
                	*pbuf++ = '\n';
                }
                memcpy(pbuf, qualbuf, seqlen[k]);qualbuf += seqlen[k]; pbuf+=seqlen[k];//qual
                *pbuf++ = '\n';

                if(!isread2) //第一条read
                {
                    out_ind1 += pbuf - &write_buf1[out_ind1];
                }
                else
                {
                    out_ind2 += pbuf - &write_buf2[out_ind2];
                }
            }

            f1.write(write_buf1, out_ind1);
            f2.write(write_buf2, out_ind2);
        }
    }
    
    delete pfqz;
    free(write_buf1);
    free(write_buf2);
    f1.close();
    f2.close();
    in_isq.close();
}

void *fqzall_decode_process(void *data)
{
    if(data == NULL) return data;

    DecodeParam *param = (DecodeParam*)data;

    if(g_magicparam.isSE)
    {
        fqzall_decode_SE(param);
    }
    else
    {
        fqzall_decode_PE(param);
    }

    delete param;
}

void *encode_process(void *data)
{
    if(data == NULL) return data;
    EncodeParam *param = (EncodeParam*)data;
    SeqRead ssread(param->filename[0], param->offset[0], param->length[0]);
    //printf("---%0x %ld %ld %d\n", pthread_self(), param->offset[0], param->length[0], param->num);

    bwtintv_v *matcher1 = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
    fqz *pfqz = new fqz(&g_fqz_params);

    fstream out_isq; 
    char str_tmp[64]={0};
    sprintf(str_tmp, "./out_isq_%d.tmp",param->num);
    out_isq.open(str_tmp, std::ios::binary|std::ios::out);

    align_info align_info1;
    align_info1.cigar_l = (int*)malloc(max_mis * sizeof(int));
    align_info1.cigar_v = (int*)malloc(max_mis * sizeof(int));

    if (g_magicparam.isSE)
    {
        string bitbuf;
        while(ssread.getRead())
        {
            char degenerate[512]={0};
            int seqm = getAlignInfoSE(ssread.seq, matcher1, NULL, param->pitr, param->pidx, align_info1, degenerate);
            uint32_t name_len = strlen(ssread.name);
            uint32_t qual_len = strlen(ssread.qual);
            if(seqm) //比对成功
            {
                int index = AlignInfoToBitArry_SE(align_info1, strlen(ssread.seq), bitbuf);

                uint32_t tmplen = pfqz->getInLen() + name_len + qual_len + bitbuf.length();
                if(tmplen > BLK_SIZE)
                {
                	pfqz->isq_doencode_s(out_isq);
                }
                pfqz->isq_addbuf_match(ssread.name, name_len, (char*)bitbuf.c_str(), bitbuf.length(), ssread.qual, qual_len, index+1, degenerate);
            }
            else
            {
            	uint32_t tmplen = pfqz->getInLen() + name_len + qual_len + qual_len;
                if(tmplen > BLK_SIZE)
                {
                	pfqz->isq_doencode_s(out_isq);
                }
                pfqz->isq_addbuf_unmatch(ssread.name, name_len, ssread.seq, qual_len, ssread.qual, qual_len, 0);
            }       
        }
        pfqz->isq_doencode_s(out_isq);
        g_lentharry[param->num] = pfqz->getCompressTotalLen();
        //printf("%ld\n", g_lentharry[param->num]);
    }
    else
    {
        SeqRead ssread2(param->filename[1], param->offset[1], param->length[1]);
        align_info align_info2;
        align_info2.cigar_l = (int*)malloc(max_mis * sizeof(int));
        align_info2.cigar_v = (int*)malloc(max_mis * sizeof(int));

        int64_t ** idx1 = new int64_t*[max_smem_num*max_iwidth];
        for(int i=0; i<max_smem_num*max_iwidth; i++)
            idx1[i] = new int64_t[3];
        int64_t ** idx2 = new int64_t*[max_smem_num*max_iwidth];
        for(int i=0; i<max_smem_num*max_iwidth; i++)
            idx2[i] = new int64_t[3];
        bwtintv_v *matcher2 = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));

        string bitbuf1, bitbuf2;

        int count_match = 0;
        int count_unmatch = 0;
        while(ssread.getRead() && ssread2.getRead())
        {
            char degenerate[2][512]={0};
            
            int seqm = getAlignInfoPE(ssread.seq, ssread2.seq, matcher1, matcher2, NULL, param->pitr, param->pidx, align_info1, align_info2, idx1, idx2, degenerate);
            uint32_t name_len1 = strlen(ssread.name);
            uint32_t qual_len1 = strlen(ssread.qual);
            uint32_t name_len2 = strlen(ssread2.name);
            uint32_t qual_len2 = strlen(ssread2.qual);

            if (seqm > 0)
            {
            	int index = AlignInfoToBitArry_PE(align_info1, align_info2, strlen(ssread.seq), bitbuf1, bitbuf2);

            	uint32_t tmplen = pfqz->getInLen() + name_len1 + qual_len1 + name_len2 + qual_len2 + bitbuf1.length() + bitbuf2.length();//保证两条在一个block中
                if(tmplen > BLK_SIZE)
                {
                	pfqz->isq_doencode_s(out_isq);
                }

            	pfqz->isq_addbuf_match(ssread.name, name_len1, (char*)bitbuf1.c_str(), bitbuf1.length(), ssread.qual, qual_len1, index+1, degenerate[0]);
            	pfqz->isq_addbuf_match(ssread2.name, name_len2, (char*)bitbuf2.c_str(), bitbuf2.length(), ssread2.qual, qual_len2, index+1, degenerate[1]);
                count_match++;
            }
            else
            {
            	uint32_t tmplen = pfqz->getInLen() + name_len1 + qual_len1 + name_len2 + qual_len2 + qual_len1 + qual_len2;
                if(tmplen > BLK_SIZE)
                {
                	pfqz->isq_doencode_s(out_isq);
                }

            	pfqz->isq_addbuf_unmatch(ssread.name, name_len1, ssread.seq, qual_len1, ssread.qual, qual_len1, 0);
            	pfqz->isq_addbuf_unmatch(ssread2.name, name_len2, ssread2.seq, qual_len2, ssread2.qual, qual_len2, 0);
                count_unmatch++;
            }
        }
        pfqz->isq_doencode_s(out_isq);
        if(count_match < count_unmatch)
        {
            printf("%d align matching ratio too low. match=%ld unmatch=%ld\n", param->num, count_match, count_unmatch);
        }

        if(strlen(ssread.name) == 0) //f1文件已经读完
        {
            pfqz->isq_addmark(2); //添加f2分隔标识
            while(ssread2.getRead())
            {
            	uint32_t name_len2 = strlen(ssread2.name);
            	uint32_t qual_len2 = strlen(ssread2.qual);
            	uint32_t tmplen = pfqz->getInLen() + name_len2 + qual_len2 + qual_len2;
                if(tmplen > BLK_SIZE)
                {
                	pfqz->isq_doencode_s(out_isq);
                }

                pfqz->isq_addbuf_unmatch(ssread2.name, name_len2, ssread2.seq, qual_len2, ssread2.qual, qual_len2, 0);
            }
            pfqz->isq_doencode_s(out_isq);
        }
        else if(strlen(ssread2.name) == 0) //f2文件已经读完
        {
            pfqz->isq_addmark(1); //添加f1分隔标识
            do //f1已经读取了内容
            {
            	uint32_t name_len1 = strlen(ssread.name);
            	uint32_t qual_len1 = strlen(ssread.qual);
            	uint32_t tmplen = pfqz->getInLen() + name_len1 + qual_len1 + qual_len1;
                if(tmplen > BLK_SIZE)
                {
                	pfqz->isq_doencode_s(out_isq);
                }

                pfqz->isq_addbuf_unmatch(ssread.name, name_len1, ssread.seq, qual_len1, ssread.qual, qual_len1, 0);
            }while(ssread.getRead());
            pfqz->isq_doencode_s(out_isq);
        }
        

        g_lentharry[param->num] = pfqz->getCompressTotalLen();

        free(matcher2->a);
        free(matcher2);
        free(align_info2.cigar_l);
        free(align_info2.cigar_v);

        for(int i=0; i<max_smem_num*max_iwidth; i++)
            delete[] idx1[i];
        delete[] idx1;
        for(int i=0; i<max_smem_num*max_iwidth; i++)
            delete[] idx2[i];
        delete[] idx2;
    }

    free(matcher1->a);
    free(matcher1);
    free(align_info1.cigar_l);
    free(align_info1.cigar_v);

    
    delete param;
    delete pfqz;
    out_isq.close();
}

char getch(char ch, bool isrev)
{
    if(isrev)
    {
        ch = 3-ch;
    }

    switch (ch)
    {
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
    default: return 'N';
    }
}

int mapvar2base(int ch, int id)
{
    // static map<int,map<int,int>> map_var2base;
    // if(map_var2base.empty())
    // {
    //     map_var2base[2] = { {0, 1}, {1, 0}, {2, 3} };
    //     map_var2base[1] = { {0, 2}, {1, 0}, {2, 3} };
    //     map_var2base[0] = { {0, 2}, {1, 1}, {2, 3} };
    //     map_var2base[3] = { {0, 2}, {1, 1}, {2, 0} };
    // }

    switch(ch)
    {
        case 0:
        {
            switch(id)
            {
                case 0: return 2;
                case 1: return 1;
                case 2: return 3;
                case 3: return 4;
            }
        }
        case 1:
        {
            switch(id)
            {
                case 0: return 2;
                case 1: return 0;
                case 2: return 3;
                case 3: return 4;
            }
        }
        case 2:
        {
            switch(id)
            {
                case 0: return 1;
                case 1: return 0;
                case 2: return 3;
                case 3: return 4;
            }
        }
        case 3:
        {
            switch(id)
            {
                case 0: return 2;
                case 1: return 1;
                case 2: return 0;
                case 3: return 4;
            }
        }
    }
}

void decode_from_bitarry(BitDecode &bitdecode, int quallen, uint8_t order, bwaidx_t *pidx, char *outbuf, std::vector<char>::iterator& itor, bool isread2 = false)
{
    order -= 1;
	uint64_t pos = bitdecode.getpos(isread2);
    bool isrev = bitdecode.bufferIn(1);
    int cigar_num = bitdecode.bufferIn(2);
    //printf("----%ld %d %d \n", pos, isrev, cigar_num);
    int cigar_l[3] = {-1,-1,-1};
    int cigar_v[3] = {-1,-1,-1};
    int offset = 0;
    for(int i=0;i<cigar_num;i++) 
    {
        cigar_l[i] = bitdecode.bufferIn(getbitnum(quallen-offset))+offset;
        offset = cigar_l[i];
    }
    //printf("%d %d %d \n", cigar_l[0],cigar_l[1],cigar_l[2]);
    for(int i=0;i<cigar_num;i++)
    {
        cigar_v[i] = bitdecode.bufferIn(2);
    }
    //printf("%d %d %d \n", cigar_v[0],cigar_v[1],cigar_v[2]);
    uint64_t realpos = 0;
    if(g_offset_bit<32)
    {
        uint32_t tmppos = (order<<g_offset_bit) + pos;
        realpos = tmppos;
    }
    else
    {
        realpos = (order<<g_offset_bit) + pos;
    }
    
    int64_t rlen = 0;
    uint8_t *rseq = bns_get_seq(pidx->bns->l_pac, pidx->pac, realpos, realpos + quallen, &rlen);
    
    for (int i=0; i< cigar_num; i++){
        if (cigar_l[i] != -1){
            rseq[cigar_l[i]] = mapvar2base(rseq[cigar_l[i]], cigar_v[i]);
        }
        else
            break;
    }

    int j = 0;
    for(int i=0;i<quallen;i++) 
    {
        outbuf[i] = isrev ? getch(rseq[rlen-1-j],true): getch(rseq[j],false);
        if(outbuf[i] == 'N')
        {
            outbuf[i] = *itor;
            itor++;
        }
        j++;
    }
    free(rseq);
}


void decode_process_SE(DecodeParam *param)
{
	fqz *pfqz = new fqz(&g_fqz_params);
    
    fstream in_isq;
    in_isq.open(param->filename, std::ios::binary|std::ios::in);
    in_isq.seekg(param->offset);

    fstream f1; //解压输出到文件
    char tmpath[125]={0};
    sprintf(tmpath,"./decode1_%d.tmp", param->num);
    f1.open(tmpath, std::ios::binary|std::ios::out);


    char *write_buf = NULL;
    uint32_t write_buf_len = 0;
    int read_len = 0;
    uint64_t total_len = param->length;
    while(true)
    {
        char *namebuf = NULL;
        char *seqbuf = NULL;
        char *qualbuf = NULL;
		char *bitbuf = NULL;
		uint8_t *orderbuf = NULL;        
        //uint16_t *seqlen = NULL;
        uint16_t *quallen = NULL;
        std::vector<char> *pvec = NULL;
        int count = 0;
        int mark = 0;
        int vec_index = 0;
        
        if(total_len>0)
        {
        	read_len = pfqz->isq_decode_s(in_isq, &namebuf, &seqbuf, &qualbuf, &bitbuf, &orderbuf, &quallen, &count, &mark, &pvec);
        	total_len -= read_len;
        }
        else
        {
        	break;
        }

        uint32_t tmplen = count*5*quallen[0];
        if(write_buf_len < tmplen)
        {
            write_buf_len = tmplen;
            write_buf = (char*)realloc(write_buf, write_buf_len);
        }

        char buf[256]={0};
        uint32_t out_ind = 0;
        BitDecode bitdecode;
        bitdecode.setbuf(bitbuf);
        std::vector<char>::iterator itor = pvec->begin();
        for (int k = 0; k < count; k++) 
        {
        	char *pname = namebuf;
            while ((write_buf[out_ind++] = *namebuf++) != '\n');
            int id_len = namebuf - pname -1;

        	if(orderbuf[k] > 0)
        	{
				decode_from_bitarry(bitdecode, quallen[k], orderbuf[k], param->pidx, buf, itor);
                memcpy(&write_buf[out_ind], buf, quallen[k]);
                out_ind += quallen[k];
        	}
        	else
        	{
                memcpy(&write_buf[out_ind], seqbuf, quallen[k]); seqbuf += quallen[k];
                out_ind += quallen[k];
        	}


            write_buf[out_ind++] = '\n';
            write_buf[out_ind++] = '+';
            if(!g_magicparam.one_ch)//第三行有多个字符
            {
            	memcpy(&write_buf[out_ind], pname+1, id_len);out_ind += id_len;
            }
            else
            {
            	write_buf[out_ind++] = '\n';
            }
            memcpy(&write_buf[out_ind], qualbuf, quallen[k]);qualbuf += quallen[k];out_ind+=quallen[k];
            write_buf[out_ind++] = '\n';
        }

        f1.write(write_buf, out_ind);
    }
    
    f1.close();
    in_isq.close();
    delete pfqz;
    free(write_buf);
}

void decode_process_PE(DecodeParam *param)
{
	fqz *pfqz = new fqz(&g_fqz_params);
    
    fstream in_isq;
    in_isq.open(param->filename, std::ios::binary|std::ios::in);
    in_isq.seekg(param->offset);

    fstream f1; //解压输出到文件
    char tmpath[125]={0};
    sprintf(tmpath,"./decode1_%d.tmp", param->num);
    f1.open(tmpath, std::ios::binary|std::ios::out);

    fstream f2; //解压输出到文件
    char tmpath2[125]={0};
    sprintf(tmpath2,"./decode2_%d.tmp", param->num);
    f2.open(tmpath2, std::ios::binary|std::ios::out);

    char *write_buf1 = NULL;
    char *write_buf2 = NULL;
    uint32_t write_buf_len = 0;
    int read_len = 0;
    uint64_t total_len = param->length;
    while(true)
    {
        char *namebuf = NULL;
        char *seqbuf = NULL;
        char *qualbuf = NULL;
        char *bitbuf = NULL;
		uint8_t *orderbuf = NULL; 
        //uint16_t *seqlen = NULL;
        uint16_t *quallen = NULL;
        std::vector<char> *pvec = NULL;
        int count = 0;
        int mark = 0;
        if(total_len>0)
        {
        	read_len = pfqz->isq_decode_s(in_isq, &namebuf, &seqbuf, &qualbuf, &bitbuf, &orderbuf, &quallen, &count, &mark, &pvec);
        	total_len -= read_len;
        }
        else
        {
        	break;
        }

        uint32_t tmplen = count*5*quallen[0];
        if(write_buf_len < tmplen)
        {
            write_buf_len = tmplen;
            write_buf1 = (char*)realloc(write_buf1, write_buf_len);
            write_buf2 = (char*)realloc(write_buf2, write_buf_len);
        }

        char buf[256]={0};
        uint32_t out_ind1 = 0;
        uint32_t out_ind2 = 0;

        if(mark == 1)//写入到f1
        {
        	if(g_magicparam.one_ch)
        	{
        		out_ind1 = parse_decode_buf(write_buf1, namebuf, seqbuf, qualbuf, quallen, count);
        	}
        	else
        	{
        		out_ind1 = parse_decode_buf_more(write_buf1, namebuf, seqbuf, qualbuf, quallen, count);
        	}
            f1.write(write_buf1, out_ind1);
        }
        else if(mark == 2)//写入到f2
        {
        	if(g_magicparam.one_ch)
        	{
				out_ind2 = parse_decode_buf(write_buf2, namebuf, seqbuf, qualbuf, quallen, count);
        	}
        	else
        	{
        		out_ind2 = parse_decode_buf_more(write_buf2, namebuf, seqbuf, qualbuf, quallen, count);
        	}
            
            f2.write(write_buf2, out_ind2);
        }
        else //分别写入f1和f2
        {
        	BitDecode bitdecode;
        	bitdecode.setbuf(bitbuf);
            std::vector<char>::iterator itor = pvec->begin();
            for (int k = 0; k < count; k++) 
            {
                char *pbuf = NULL;
                bool isread2 = false;
                if(k&1) //第二条read
                {
                    pbuf = &write_buf2[out_ind2];
                    isread2 = true;
                }
                else //第一条read
                {
                    pbuf = &write_buf1[out_ind1];
                }

                char *pame = namebuf;
                while ((*pbuf++ = *namebuf++) != '\n'); //id
                int id_len = namebuf - pame -1;

                if(orderbuf[k] > 0)
                {
                	//memset(buf, 0, 256);
                	decode_from_bitarry(bitdecode, quallen[k], orderbuf[k], param->pidx, buf, itor, isread2);
                	memcpy(pbuf, buf, quallen[k]); pbuf+=quallen[k];//seq
                }
                else
                {
                	memcpy(pbuf, seqbuf, quallen[k]); seqbuf += quallen[k]; pbuf+=quallen[k];//seq
                }
                
                *pbuf++ = '\n';
                *pbuf++ = '+';
                if(!g_magicparam.one_ch)
                {
                	memcpy(pbuf, pame+1, id_len); pbuf += id_len;
                }
                else
                {
                	*pbuf++ = '\n';
                }
                memcpy(pbuf, qualbuf, quallen[k]);qualbuf += quallen[k]; pbuf+=quallen[k];//qual
                *pbuf++ = '\n';

                if(!isread2) //第一条read
                {
                    out_ind1 += pbuf - &write_buf1[out_ind1];
                }
                else
                {
                    out_ind2 += pbuf - &write_buf2[out_ind2];
                }
            }

            f1.write(write_buf1, out_ind1);
            f2.write(write_buf2, out_ind2);
        }
    }
    
    f1.close();
    f2.close();
    in_isq.close();
    delete pfqz;
    free(write_buf1);
    free(write_buf2);
}


void *decode_process(void *data)
{
    if(data == NULL) return data;

    DecodeParam *param = (DecodeParam*)data;

    if(g_magicparam.isSE)
    {
    	decode_process_SE(param);
    }
    else
    {
    	decode_process_PE(param);
    }

    delete param;
}

uint64_t Getoffset(uint64_t *arry, int num)
{
    int i;
    uint64_t tmp = 0;
    for(i=0;i<num ;i++)
    {
        tmp += arry[i];
    }
    return tmp;
}

bool isSimilarName(char *name1, char *name2)
{
    int len1 = strlen(name1);
    int len2 = strlen(name2);
    if(len1 != len2) return false;
    int count = 0;
    for(int i=0;i<len1;i++)
    {
        if(name1[i] != name2[i])
        {
            count++;
        }
    }

    if(count>1) return false;

    return true;
}

void AdjustPESlice(char *path1, uint64_t *slicearry1, char *path2, uint64_t *slicearry2)
{
    FILE *f1 = fopen(path1, "r");
    FILE *f2 = fopen(path2, "r");
    if(f1 == NULL || f2 == NULL) return;

    for (int i = 0; i < thread_num-1; ++i) 
    {
        if(slicearry1[i] != slicearry2[i])
        {
            char name1[1024] = {0};
            char name2[1024] = {0};

            fseek(f1, slicearry1[i], SEEK_SET);
            fgets(name1, 1024, f1);

            fseek(f2, slicearry2[i], SEEK_SET);
            fgets(name2, 1024, f2);

            if(!isSimilarName(name1,name2)) //二者切片不一致，需要调整
            {
                fseek(f2, slicearry2[i]-1024, SEEK_SET);//向前偏移一定长度，查找是否有相似的name
                int len = 0;
                for(int j=0;j<40;j++) //设定比较的上限
                {
                    fgets(name2, 1024, f2);
                    if(name2[0] == '@' && isSimilarName(name1,name2)) //找到相似的read,调整slicearry2的切片位置
                    {
                        slicearry2[i] = slicearry2[i]-1024 + len;
                        break;
                    }
                    len += strlen(name2)+1;
                }
            }
        }
    }
    fclose(f1);
    fclose(f2);
}

bool IsfirstLine(char *pbuf)
{
    int count = 0;
    if(pbuf[0] == '@')
    {
        for(int i=1;i<strlen(pbuf);i++)
        {
            switch(pbuf[i])
            {
                case '!':count++;break;
                case '"':count++;break;
//                case '#':count++;break;
                case '$':count++;break;
                case '%':count++;break;
                case '&':count++;break;
                case '\'':count++;break;
                case '(':count++;break;
                case ')':count++;break;
                case '*':count++;break;
                case '+':count++;break;
                case ',':count++;break;
                case '.':count++;break;
                case '-':count++;break;
                case '^':count++;break;
                case ';':count++;break;
                case '<':count++;break;
                case '>':count++;break;
                case '=':count++;break;
                case '?':count++;break;
                case '@':count++;break;
            }

        }
        if(count < 2)
        {
            return true;
        }
    }

    return false;
}

bool GetFileSlice(char *path, uint64_t flength, uint64_t *slicearry)
{
    if(g_isgzip) return false; //gzip格式文件不能切片
    FILE *fdna = fopen(path, "r");
    if(fdna == NULL)
        abort();
    uint64_t len = flength/thread_num;

	char buf[1024]={0}; //TODO:这个固定值的做法有出错的隐患

    int exlen = 0;
    for(int i = 0; i < thread_num-1; ++i)
    {
    	bool bfind = false;
        fseek(fdna, (i+1)*len, SEEK_SET);

        slicearry[i] = len - exlen;//去除上一个切片超过的长度
        exlen = 0;
        while(!bfind)
        {
        	do
	        {
	        	fgets(buf, 1024 , fdna);
	        	exlen += strlen(buf);
	        }while(!IsfirstLine(buf));  //认为找到第一行

	        fgets(buf, 1024 , fdna); exlen += strlen(buf);//第二行
	        fgets(buf, 1024 , fdna); exlen += strlen(buf);//第三行
	        if(buf[0] == '+')
	        {
	        	bfind = true;
	        }
	        fgets(buf, 1024 , fdna); exlen += strlen(buf);//第四行
	    }

	    slicearry[i] += exlen; //加上本次切片超过的长度
        flength -= slicearry[i];
    }

    slicearry[thread_num-1] = flength;
    fclose(fdna);
    return true;
}


bool GetFileType(char *path)
{
    FILE *f = fopen(path, "rb");
    unsigned char buf[10]={0};
    fread(buf, 1, 2, f);
    fclose(f);
    if(buf[0] == 0x1f && buf[1] == 0x8b)
    {
        return true; //是gzip格式的文件
    }
    return false;
}

void SetTaskData(SeqRead &sread, int index, Task &task)
{
    task.name_len[index] = strlen(sread.name);
    task.seq_len[index] = strlen(sread.seq);
    memcpy(task.name[index], sread.name, task.name_len[index]);
    task.name[index][task.name_len[index]] = '\0';
    memcpy(task.seq[index], sread.seq, task.seq_len[index]);
    task.seq[index][task.seq_len[index]] = '\0';
    memcpy(task.qual[index], sread.qual, task.seq_len[index]);
    task.qual[index][task.seq_len[index]] = '\0';
}


void gzip_process_fqzall_SE(ThreadParam *param, fqz *pfqz, fstream &out_isq)
{
    while (true)
    {
        Task task;
        int ret = param->pthreadtask->getTask(task, param->num);
        if(ret == 1) //退出线程
        {
            break;
        }
        else if(ret == 2)//等待任务
        {
            sleep(1);
            continue;
        }

        uint32_t name_len = task.name_len[0];
        uint32_t qual_len = task.seq_len[0];

        uint32_t tmplen = pfqz->getInLen() + name_len + qual_len + qual_len;
        if(tmplen > BLK_SIZE)
        {
            pfqz->isq_doencode(out_isq);
        }

        pfqz->isq_addbuf(task.name[0], name_len, task.seq[0], qual_len, task.qual[0], qual_len);
    }

    pfqz->isq_doencode(out_isq);
    g_lentharry[param->num] = pfqz->getCompressTotalLen();
}

void gzip_process_fqzall_PE(ThreadParam *param, fqz *pfqz, fstream &out_isq)
{
    while (true)
    {
        Task task;
        int ret = param->pthreadtask->getTask(task, param->num);
        if(ret == 1) //退出线程
        {
            break;
        }
        else if(ret == 2)//等待任务
        {
            sleep(1);
            continue;
        }

        uint32_t name_len1 = task.name_len[0];
        uint32_t qual_len1 = task.seq_len[0];
        uint32_t name_len2 = task.name_len[1];
        uint32_t qual_len2 = task.seq_len[1];

        uint32_t tmplen = pfqz->getInLen() + name_len1 + qual_len1 + name_len2 + qual_len2 + qual_len1 + qual_len2;//保证两条在一个block中
        if(tmplen > BLK_SIZE)
        {
            pfqz->isq_doencode(out_isq);
        }
        pfqz->isq_addbuf(task.name[0], name_len1, task.seq[0], qual_len1, task.qual[0], qual_len1);
        pfqz->isq_addbuf(task.name[1], name_len2, task.seq[1], qual_len2, task.qual[1], qual_len2);
    }
    
    pfqz->isq_doencode(out_isq);
    g_lentharry[param->num] = pfqz->getCompressTotalLen();
}

void gzip_process_encode_SE(ThreadParam *param, fqz *pfqz, fstream &out_isq)
{
    bwtintv_v *matcher1 = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
    align_info align_info1;
    align_info1.cigar_l = (int*)malloc(max_mis * sizeof(int));
    align_info1.cigar_v = (int*)malloc(max_mis * sizeof(int));
    string bitbuf;

    while (true)
    {
        Task task;
        int ret = param->pthreadtask->getTask(task, param->num);
        if(ret == 1) //退出线程
        {
            break;
        }
        else if(ret == 2)//等待任务
        {
            sleep(1);
            continue;
        }

        char degenerate[512]={0};
        int seqm = getAlignInfoSE(task.seq[0], matcher1, NULL, param->pitr, param->pidx, align_info1, degenerate);
        uint32_t name_len = task.name_len[0];
        uint32_t qual_len = task.seq_len[0];
        if(seqm) //比对成功
        {
            int index = AlignInfoToBitArry_SE(align_info1, qual_len, bitbuf);

            uint32_t tmplen = pfqz->getInLen() + name_len + qual_len + bitbuf.length();
            if(tmplen > BLK_SIZE-1024)
            {
                pfqz->isq_doencode_s(out_isq);
            }
            pfqz->isq_addbuf_match(task.name[0], name_len, (char*)bitbuf.c_str(), bitbuf.length(), task.qual[0], qual_len, index+1, degenerate);
        }
        else
        {
            uint32_t tmplen = pfqz->getInLen() + name_len + qual_len + qual_len;
            if(tmplen > BLK_SIZE-1024)
            {
                pfqz->isq_doencode_s(out_isq);
            }
            pfqz->isq_addbuf_unmatch(task.name[0], name_len, task.seq[0], qual_len, task.qual[0], qual_len, 0);
        }
    }

    pfqz->isq_doencode_s(out_isq);
    g_lentharry[param->num] = pfqz->getCompressTotalLen();

    free(matcher1->a);
    free(matcher1);
    free(align_info1.cigar_l);
    free(align_info1.cigar_v);
}

void gzip_process_encode_PE(ThreadParam *param, fqz *pfqz, fstream &out_isq)
{
    bwtintv_v *matcher1 = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
    align_info align_info1;
    align_info1.cigar_l = (int*)malloc(max_mis * sizeof(int));
    align_info1.cigar_v = (int*)malloc(max_mis * sizeof(int));
    align_info align_info2;
    align_info2.cigar_l = (int*)malloc(max_mis * sizeof(int));
    align_info2.cigar_v = (int*)malloc(max_mis * sizeof(int));

    int64_t ** idx1 = new int64_t*[max_smem_num*max_iwidth];
    for(int i=0; i<max_smem_num*max_iwidth; i++)
        idx1[i] = new int64_t[3];
    int64_t ** idx2 = new int64_t*[max_smem_num*max_iwidth];
    for(int i=0; i<max_smem_num*max_iwidth; i++)
        idx2[i] = new int64_t[3];
    bwtintv_v *matcher2 = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));

    string bitbuf1, bitbuf2;

    int count_match = 0;
    int count_unmatch = 0;

    while (true)
    {
        Task task;
        int ret = param->pthreadtask->getTask(task, param->num);
        if(ret == 1) //退出线程
        {
            break;
        }
        else if(ret == 2)//等待任务
        {
            sleep(1);
            continue;
        }

        char degenerate[2][512]={0};
        int seqm = getAlignInfoPE(task.seq[0], task.seq[1], matcher1, matcher2, NULL, param->pitr, param->pidx, align_info1, align_info2, idx1, idx2, degenerate);
        uint32_t name_len1 = task.name_len[0];
        uint32_t qual_len1 = task.seq_len[0];
        uint32_t name_len2 = task.name_len[1];
        uint32_t qual_len2 = task.seq_len[1];

        if (seqm > 0)
        {
            int index = AlignInfoToBitArry_PE(align_info1, align_info2, qual_len1, bitbuf1, bitbuf2);

            uint32_t tmplen = pfqz->getInLen() + name_len1 + qual_len1 + name_len2 + qual_len2 + bitbuf1.length() + bitbuf2.length();//保证两条在一个block中
            if(tmplen > BLK_SIZE)
            {
                pfqz->isq_doencode_s(out_isq);
            }

            pfqz->isq_addbuf_match(task.name[0], name_len1, (char*)bitbuf1.c_str(), bitbuf1.length(), task.qual[0], qual_len1, index+1, degenerate[0]);
            pfqz->isq_addbuf_match(task.name[1], name_len2, (char*)bitbuf2.c_str(), bitbuf2.length(), task.qual[1], qual_len2, index+1, degenerate[1]);
            count_match++;
        }
        else
        {
            uint32_t tmplen = pfqz->getInLen() + name_len1 + qual_len1 + name_len2 + qual_len2 + qual_len1 + qual_len2;
            if(tmplen > BLK_SIZE)
            {
                pfqz->isq_doencode_s(out_isq);
            }

            pfqz->isq_addbuf_unmatch(task.name[0], name_len1, task.seq[0], qual_len1, task.qual[0], qual_len1, 0);
            pfqz->isq_addbuf_unmatch(task.name[1], name_len2, task.seq[1], qual_len2, task.qual[1], qual_len2, 0);
            count_unmatch++;
        }
    }

    pfqz->isq_doencode_s(out_isq);
    g_lentharry[param->num] = pfqz->getCompressTotalLen();

    free(matcher1->a);
    free(matcher1);
    free(align_info1.cigar_l);
    free(align_info1.cigar_v);

    free(matcher2->a);
    free(matcher2);
    free(align_info2.cigar_l);
    free(align_info2.cigar_v);

    for(int i=0; i<max_smem_num*max_iwidth; i++)
        delete[] idx1[i];
    delete[] idx1;
    for(int i=0; i<max_smem_num*max_iwidth; i++)
        delete[] idx2[i];
    delete[] idx2;
}


void *gzip_process(void *data)
{
    if (!data) data;
    ThreadParam *pParam = (ThreadParam*)data;

    fqz *pfqz = new fqz(&g_fqz_params);

    fstream out_isq; 
    char str_tmp[64]={0};
    sprintf(str_tmp, "./out_isq_%d.tmp",pParam->num);
    out_isq.open(str_tmp, std::ios::binary|std::ios::out);


    if(g_magicparam.fqzall)
    {
        if(g_magicparam.isSE)
        {
            gzip_process_fqzall_SE(pParam, pfqz, out_isq);
        }
        else
        {
            gzip_process_fqzall_PE(pParam, pfqz, out_isq);
        }
    }
    else
    {
        if(g_magicparam.isSE)
        {
            gzip_process_encode_SE(pParam, pfqz, out_isq);
        }
        else
        {
            gzip_process_encode_PE(pParam, pfqz, out_isq);
        }
    }

    delete pfqz;
    delete pParam;
    out_isq.close();
}

const int WRITE_LEN = 2*1024;
bool ReadWriteData(fstream &f_s, char *writbuf, int &count)
{
    if(f_s.peek() == EOF) return false;

    f_s.getline(writbuf+count, 1024);
    count += f_s.gcount();
    writbuf[count-1]='\n';
    f_s.getline(writbuf+count, 1024);
    count += f_s.gcount();
    writbuf[count-1]='\n';
    f_s.getline(writbuf+count, 1024);
    count += f_s.gcount();
    writbuf[count-1]='\n';
    f_s.getline(writbuf+count, 1024);
    count += f_s.gcount();
    writbuf[count-1]='\n';

    return true;
}

bool isFinish(int *ret, int num)
{
    for(int i=0;i<num;i++)
    {
        if(ret[i]) return true;
    }
    return false;
}

void MergeFileForzip(string &fastq_prefix, int num, int index)
{
    char fastq_path[256]={0};
    sprintf(fastq_path,"./%s%d.fastq", fastq_prefix.c_str(), index);
    fstream out_s;
    out_s.open(fastq_path, ios::out | ios::binary);

    fstream *pf_s = new fstream[num];
    std::vector<string> vec_path;
    char *writbuf = new char[WRITE_LEN*num];

    for(int i=0;i<num;i++)
    {
        char str_tmp[64]={0};
        sprintf(str_tmp, "./decode%d_%d.tmp", index, i);
        vec_path.emplace_back(str_tmp);
        pf_s[i].open(str_tmp, ios::in | ios::binary);
    }

    int *pret = new int[num];
    do
    {
        int count = 0;
        for(int i=0;i<num;i++)
        {
            pret[i] = ReadWriteData(pf_s[i], writbuf, count);
        }
        out_s.write(writbuf, count);
    }while(isFinish(pret, num));

    for(int i=0;i<num;i++)
    {
        pf_s[i].close();
        remove(vec_path[i].c_str()); 
    }
    delete[] pf_s;
    delete[] writbuf;
    out_s.close();
}

void MergeFileForFastq(string &fastq_prefix, int num, int index)
{
    char fastq_path[256]={0};
    sprintf(fastq_path,"./%s%d.fastq", fastq_prefix.c_str(), index);
    fstream out_s;
    out_s.open(fastq_path, ios::out | ios::binary);

    for (int i = 0; i < num; i++)
    {
        char str_tmp[64]={0};
        sprintf(str_tmp, "./decode%d_%d.tmp", index, i);
        fstream f_s;
        f_s.open(str_tmp, ios::in | ios::binary);
        out_s << f_s.rdbuf();
        f_s.close();
        remove(str_tmp); 
    }

    out_s.close();
}
