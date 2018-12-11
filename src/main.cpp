#include <stdio.h>
#include <zlib.h>
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
#include <queue>
#include <pthread.h>
#include <semaphore.h>
#include <sys/time.h>
#include "bwa/bwa.h"
#include "bwa/bwamem.h"
#include "bwa/utils.h"
#include "bwa/bntseq.h"
#include "bwa/bwtaln.h"
#include "bwa/bwt.h"
#include "kseq.hpp" // for the FASTA/Q parser
#include "fqzcomp.h"
#include "ref2seq.h"
#include "align_proc.h"

using namespace std;

#define MAJOR_VERS 0
#define MINOR_VERS 1
#define FORMAT_VERS 1
#define PREALIGN_NUM 2000

/* -------------------------------------------------------------------------
 * BWA
 */
#define SEQARCTHREAD
typedef struct _tagTask
{
    int num;
    kseq seq1;
    kseq seq2;
}Task;

std::queue<Task> g_task_queue;

bool g_bFinish = false;   //是否结束线程
pthread_mutex_t g_mutex;  //读文件锁
pthread_mutex_t g_write_mutex; //写文件锁
sem_t g_sem;    //队列信号量

//Align Params
int min_len = 17, max_iwidth = 50, max_mis = 3, lgst_num = 2, max_smem_num = 2, exp_mismatch = 1, fqzall = 0, max_insr=511;
int64_t **idx1;
int64_t **idx2;

int count_blocknum(uint64_t refsize, int block_num){
    uint64_t idx = 1;
    float tmp, res=0.01;
    while (1){
        idx *= 2;
        if (idx>refsize)
            break;
    }
    while (1){
        idx /= 2;
        tmp = (float)refsize / idx;
        if (tmp > block_num)
            return (int)ceil(res);
        if (ceil(tmp)-tmp < 0.3)
            return (int)ceil(tmp);
        else{
            if (tmp-int(tmp) > res-int(res))
                res = tmp;
        }
    }
}

bool bwtintv_lencmp(const bwtintv_t &arg1, const bwtintv_t &arg2) {     //长的SMEM排前面
    return (uint32_t) arg1.info - (uint32_t) (arg1.info >> 32) > (uint32_t) arg2.info - (uint32_t) (arg2.info >> 32);
}

int idx_compar( const void* a, const void* b ) {
    return ((const int*)a)[0] - ((const int*)b)[0];
}

static map<int, map<int, int>> nucleBitMap;
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

int getAlignInfoSE(kseq &seq, bwtintv_v *a, bwtintv_v *tmpvec[2], smem_i* func_itr, bwaidx_t *func_idx, align_info &align_info1, int func_block_size){
    int64_t rlen;
    int i, seql, base;

    seql = (int) seq.seq.length();
    int pass_num = 0;
    int cigar_l[max_mis], cigar_v[max_mis];

    unsigned char seqarry[seql];
    for (i = 0; i < seql; ++i) {
        seqarry[i] = nst_nt4_table[(int) seq.seq[i]];
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
                        if ((rseq_l[base] + seqarry[seql-base-1] != 3) && (seqarry[seql-base-1] <= 3)) {
                            if (missum == max_mis){
                                ++missum;
                                break;
                            }
                            rbase = int(rseq_l[base]);
                            sbase = int(seqarry[seql-base-1]);
                            cigar_l[missum] = base;
                            cigar_v[missum] = nucleBitMap[rbase][3-sbase];
                            ++missum;
                        }
                    }
                    free(rseq_l);
                    if (missum <= max_mis) {
                        rseq_r = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos + seql - (uint32_t)(pvec[i].info >>32),
                                             pos + seql, &rlen);
                        //把rseq_r和read的左截的反向互补进行比较
                        for (base = 0; base < rlen; base++) {
                            if ((rseq_r[base] + seqarry[(uint32_t)(pvec[i].info >>32)-base-1] != 3) && (seqarry[(uint32_t)(pvec[i].info >>32)-base-1] <= 3)) {
                                if (missum == max_mis){
                                    ++missum;
                                    break;
                                }
                                rbase = int(rseq_r[base]);
                                sbase = int(seqarry[(uint32_t)(pvec[i].info >>32)-base-1]);
                                cigar_l[missum] = seql-(uint32_t)(pvec[i].info >>32)+base;
                                cigar_v[missum] = nucleBitMap[rbase][3-sbase];
                                ++missum;
                            }
                        }
                        free(rseq_r);
                    }
                    if (missum <= max_mis){
                        if (missum < last_missum){
                            align_info1.blockNum = (int) (pos / func_block_size);
                            align_info1.blockPos = (int) (pos % func_block_size);
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
                        if ((rseq_l[base] != seqarry[base]) && (seqarry[base] <= 3)) {
                            if (missum == max_mis){
                                missum += 1;
                                break;
                            }
                            rbase = int(rseq_l[base]);
                            sbase = int(seqarry[base]);
                            cigar_l[missum] = base;
                            cigar_v[missum] = nucleBitMap[rbase][sbase];
                            missum += 1;
                        }
                    }
                    free(rseq_l);
                    if (missum <= max_mis) {
                        rseq_r = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos+ (uint32_t)(pvec[i].info),
                                             pos + seql, &rlen);
                        //把rseq_r和read的右截进行比较
                        for (base = 0; base < rlen; base++) {
                            if ((rseq_r[base] != seqarry[(uint32_t) pvec[i].info + base]) && (seqarry[(uint32_t) pvec[i].info + base] <= 3)) {
                                if (missum == max_mis){
                                    missum += 1;
                                    break;
                                }
                                rbase = int(rseq_r[base]);
                                sbase = int(seqarry[(uint32_t) pvec[i].info + base]);
                                cigar_l[missum] = base+(uint32_t)(pvec[i].info);
                                cigar_v[missum] = nucleBitMap[rbase][sbase];
                                missum += 1;
                            }
                        }
                        free(rseq_r);
                    }
                    if (missum <= max_mis){
                        if (missum < last_missum){
                            align_info1.blockNum = (int) (pos / func_block_size);
                            align_info1.blockPos = (int) (pos % func_block_size);
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

int getAlignInfoPE(kseq &seq1, kseq &seq2, bwtintv_v *a1, bwtintv_v *a2, bwtintv_v *tmpvec[2], smem_i* func_itr, bwaidx_t *func_idx, align_info &align_info1, align_info &align_info2, int func_block_size){
    int64_t rlen;
    int i, j, seql1, seql2, base, count;
    int* seql_p;

    seql1 = (int) seq1.seq.length();
    seql2 = (int) seq2.seq.length();

    unsigned char seqarry1[seql1], seqarry2[seql2];
    unsigned char* arry_p;
    for (i = 0; i < seql1; ++i) {
        seqarry1[i] = nst_nt4_table[(int) seq1.seq[i]];
    }
    for (i = 0; i < seql2; ++i) {
        seqarry2[i] = nst_nt4_table[(int) seq2.seq[i]];
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
                if (idx1[i][0]/func_block_size != idx2[j][0]/func_block_size)
                    return 0;
                else
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
                if (idx1[i][0]/func_block_size != idx2[j][0]/func_block_size)
                    return 0;
                else
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

    for (int k=0;k<2;++k) {
        uint16_t missum = 0;
        if (is_rev) {
            rseq_l = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos,
                                 pos + *seql_p - (uint32_t) (tmp_info), &rlen);
            //把rseq_l和read的右截的反向互补进行比较
            for (base = 0; base < rlen; base++) {
                if ((rseq_l[base] + arry_p[*seql_p - base - 1] != 3) && (arry_p[*seql_p - base - 1] <= 3)) {
                    if (missum == max_mis) {
                        ++missum;
                        break;
                    }
                    rbase = int(rseq_l[base]);
                    sbase = int(arry_p[*seql_p - base - 1]);
                    aligninfo_p->cigar_l[missum] = base;
                    aligninfo_p->cigar_v[missum] = nucleBitMap[rbase][3 - sbase];
                    ++missum;
                }
            }
            free(rseq_l);
            if (missum <= max_mis) {
                rseq_r = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos + *seql_p - (uint32_t) (tmp_info >> 32),
                                     pos + *seql_p, &rlen);
                //把rseq_r和read的左截的反向互补进行比较
                for (base = 0; base < rlen; base++) {
                    if ((rseq_r[base] + arry_p[(uint32_t) (tmp_info >> 32) - base - 1] != 3) &&
                        (arry_p[(uint32_t) (tmp_info >> 32) - base - 1] <= 3)) {
                        if (missum == max_mis) {
                            ++missum;
                            break;
                        }
                        rbase = int(rseq_r[base]);
                        sbase = int(arry_p[(uint32_t) (tmp_info >> 32) - base - 1]);
                        aligninfo_p->cigar_l[missum] = *seql_p - (uint32_t) (tmp_info >> 32) + base;
                        aligninfo_p->cigar_v[missum] = nucleBitMap[rbase][3 - sbase];
                        ++missum;
                    }
                }
                free(rseq_r);
            }
            if (missum <= max_mis) {
                aligninfo_p->blockNum = (int) (pos / func_block_size);
                aligninfo_p->blockPos = (int) (pos % func_block_size);
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
                if ((rseq_l[base] != arry_p[base]) && (arry_p[base] <= 3)) {
                    if (missum == max_mis) {
                        missum += 1;
                        break;
                    }
                    rbase = int(rseq_l[base]);
                    sbase = int(arry_p[base]);
                    aligninfo_p->cigar_l[missum] = base;
                    aligninfo_p->cigar_v[missum] = nucleBitMap[rbase][sbase];
                    missum += 1;
                }
            }
            free(rseq_l);
            if (missum <= max_mis) {
                rseq_r = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos + (uint32_t) (tmp_info),
                                     pos + *seql_p, &rlen);
                //把rseq_r和read的右截进行比较
                for (base = 0; base < rlen; base++) {
                    if ((rseq_r[base] != arry_p[(uint32_t) tmp_info + base]) &&
                        (arry_p[(uint32_t) tmp_info + base] <= 3)) {
                        if (missum == max_mis) {
                            missum += 1;
                            break;
                        }
                        rbase = int(rseq_r[base]);
                        sbase = int(arry_p[(uint32_t) tmp_info + base]);
                        aligninfo_p->cigar_l[missum] = base + (uint32_t) (tmp_info);
                        aligninfo_p->cigar_v[missum] = nucleBitMap[rbase][sbase];
                        missum += 1;
                    }
                }
                free(rseq_r);
            }
            if (missum <= max_mis) {
                aligninfo_p->blockNum = (int) (pos / func_block_size);
                aligninfo_p->blockPos = (int) (pos % func_block_size);
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
        }
    }
    return 1;
}


/* -------------------------------------------------------------------------
 * Read modification
 */
void readModify1(string& seq, string& qual, int qualSys, int max_readLen){
    if (seq.length() > max_readLen){
        cout << "Found a read(" << seq.length() << "bp) longer than max_readLen, please increase the param" << endl;
        exit(1);
    }
    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    if (qualSys-1){ //sanger
        for (int i=0;i<seq.length();i++){
            if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T') {
                seq[i] = 'N';
                qual[i] = '!';
            }
        }
    }
    else{ //illumina
        for (int i=0;i<seq.length();i++){
            if (int(qual[i]) < 64 || (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T')) { //暂定,待修改
                seq[i] = 'N';
                qual[i] = '@';
            }
        }
    }
}

void readModify2(string& seq, string& qual, int qualSys){
    if (qualSys-1){ //sanger
        for (int i=0;i<seq.length();i++){
            if (int(qual[i]) == 33)
                seq[i] = 'N';
        }
    }
    else{ //illumina
        for (int i=0;i<seq.length();i++){
            if (int(qual[i]) == 64)
                seq[i] = 'N';
        }
    }
}

typedef struct _tagThreadParam
{
    int block_size;
    int block_num;
    smem_i *pitr;
    bwaidx_t *pidx;
    //fstream *pout_s;
    fstream *pfpOutput_s;
    //fstream *pout_iq;
    fstream *pfpOutput_iq;
    fstream *pout_isq;
    encode **pencoders;
    fqz **pfqz;
}ThreadParam;


void *task_process(void *data)
{
    if (!data)
    {
        return NULL;
    }

    ThreadParam *pParam = (ThreadParam*)data;

    bwtintv_v *matcher1 = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
    bwtintv_v *matcher2 = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
    bwtintv_v *tmpvec[2];
    tmpvec[0] = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
    tmpvec[1] = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));

    align_info align_info1;
    align_info align_info2;
    align_info1.cigar_l = (int*)malloc(max_mis * sizeof(int));
    align_info1.cigar_v = (int*)malloc(max_mis * sizeof(int));
    align_info2.cigar_l = (int*)malloc(max_mis * sizeof(int));
    align_info2.cigar_v = (int*)malloc(max_mis * sizeof(int));

    while (true)
    {
        sem_wait(&g_sem);
        pthread_mutex_lock(&g_mutex);
        if(g_task_queue.empty())
        {
            pthread_mutex_unlock(&g_mutex);
            if(g_bFinish)
            {
                printf("thread id= %0x quit\n", pthread_self());
                break; //退出线程
            }
            else
            {
                continue; //等待新任务
            }
        }
        Task task = g_task_queue.front();
        g_task_queue.pop();
        pthread_mutex_unlock(&g_mutex);

        if (task.num != 0) //有任务，开始执行
        {
            int num = task.num;
            if(num == 1) //SE
            {
                kseq &seq = task.seq1;
                if(fqzall) //匹配率太低，直接fqz压缩
                {
                    pthread_mutex_lock(&g_write_mutex);
                    pParam->pfqz[pParam->block_num]->isq_encode(seq.name, seq.seq, seq.qual, *(pParam->pout_isq));
                    pthread_mutex_unlock(&g_write_mutex);
                }
                else
                {
                    int seqm = getAlignInfoSE(seq, matcher1, tmpvec, pParam->pitr, pParam->pidx, align_info1, pParam->block_size);
                    pthread_mutex_lock(&g_write_mutex);
                    if (seqm)
                    {
                        pParam->pencoders[align_info1.blockNum]->parse_1(align_info1, seq.seq.length(), pParam->pfpOutput_s[align_info1.blockNum]);
                        pParam->pfqz[align_info1.blockNum]->iq_encode(seq.name, seq.qual, pParam->pfpOutput_iq[align_info1.blockNum]);
                    }
                    else
                    {
                        pParam->pfqz[pParam->block_num]->isq_encode(seq.name, seq.seq, seq.qual, *(pParam->pout_isq));
                    }
                    pthread_mutex_unlock(&g_write_mutex);
                }
            }
            else if(num == 2) //PE
            {
                kseq &seq1 = task.seq1;
                kseq &seq2 = task.seq2;

                if(fqzall)//匹配率太低，直接fqz压缩
                {
                    pthread_mutex_lock(&g_write_mutex);
                    pParam->pfqz[pParam->block_num]->isq_encode(seq1.name, seq1.seq, seq1.qual, *(pParam->pout_isq));
                    pParam->pfqz[pParam->block_num]->isq_encode(seq2.name, seq2.seq, seq2.qual, *(pParam->pout_isq));
                    pthread_mutex_unlock(&g_write_mutex);
                }
                else
                {
                    int seqm = getAlignInfoPE(seq1, seq2, matcher1, matcher2, tmpvec, pParam->pitr, pParam->pidx, align_info1, align_info2, pParam->block_size);

                    pthread_mutex_lock(&g_write_mutex);
                    bool bfind = false;

                    if (seqm > 0){
                        pParam->pencoders[align_info1.blockNum]->parse_1(align_info1, seq1.seq.length(), pParam->pfpOutput_s[align_info1.blockNum]);
                        pParam->pencoders[align_info2.blockNum]->parse_2(align_info2, seq2.seq.length(), pParam->pfpOutput_s[align_info2.blockNum]);
                        pParam->pfqz[align_info1.blockNum]->iq_encode(seq1.name, seq1.qual, pParam->pfpOutput_iq[align_info1.blockNum]);
                        pParam->pfqz[align_info2.blockNum]->iq_encode(seq2.name, seq2.qual, pParam->pfpOutput_iq[align_info2.blockNum]);
                    }
                    else{
                        pParam->pfqz[pParam->block_num]->isq_encode(seq1.name, seq1.seq, seq1.qual, *(pParam->pout_isq));
                        pParam->pfqz[pParam->block_num]->isq_encode(seq2.name, seq2.seq, seq2.qual, *(pParam->pout_isq));
                    }
                    pthread_mutex_unlock(&g_write_mutex);
                }
            }
        }
    }

    free(matcher1->a);
    free(matcher1);
    free(matcher2->a);
    free(matcher2);
    free(tmpvec[0]->a);
    free(tmpvec[0]);
    free(tmpvec[1]->a);
    free(tmpvec[1]);

    free(align_info1.cigar_l);
    free(align_info1.cigar_v);
    free(align_info2.cigar_l);
    free(align_info2.cigar_v);

    //printf("thread id= %0x close\n", pthread_self());
    return NULL;
}


/* -------------------------------------------------------------------------
 * Main program
 */

static void usage(int err) {
    FILE *fp = err ? stderr : stdout;

    fprintf(fp, "SeqArc v%d.%d. Yuxin Chen, Zijian Zhao, 2018\n", MAJOR_VERS, MINOR_VERS);
    fprintf(fp, "The entropy coder is derived from Fqzcomp. The aligner is derived from BWA.\n");
    fprintf(fp, "FASTA index is optional. Only Fqzcomp will be called if no index given or align ratio too low.\n");
    fprintf(fp, "For PE data, please compress each pair of files at once.\n\n");

    fprintf(fp, "To build index:\n  SeqArc -i <ref.fa>\n\n");

    fprintf(fp, "To compress:\n  SeqArc [options] [ref.fa] <input_file> [input_file2] <compress_prefix>\n");
    fprintf(fp, "    -l INT         min SMEM length to output [17]\n");
    fprintf(fp, "    -w INT         max interval size to find coordiantes [50]\n");
    fprintf(fp, "    -I INT         skip MEM mapped to over [-] places\n");
    fprintf(fp, "    -c INT         consider only the longest [2] sMEM for mapping\n");
    fprintf(fp, "    -E INT         drop out when getting a mapping result with [1] mismatch at most\n");
    fprintf(fp, "    -f INT         consider only the first [2] mapping result\n");
    fprintf(fp, "    -m INT         max mismatch to tolerate [3]\n");
    fprintf(fp, "    -B INT         maximum number of block to split reference [50]\n");
    fprintf(fp, "    -q INT         quality system, 1:illumina, 2:sanger, default as [2]\n");
    fprintf(fp, "    -s INT         max insert size between read1 and read2 [511]\n");
    fprintf(fp, "    -r INT         files with AlignRatio lower than [0.5] are processed with Fqzcomp only.\n\n");

    fprintf(fp, "    -S <level>     Sequence de novo compression level. 1-9 [3]\n");
    fprintf(fp, "                   Specifying '+' on the end (eg -s5+) will use\n");
    fprintf(fp, "                   models of multiple sizes for improved compression.\n");
    fprintf(fp, "    -N <level>     Quality compression level.  1-3 [2]\n");
    fprintf(fp, "    -n <level>     Name compression level.  1-2 [1]\n");
    fprintf(fp, "    -b             Use both strands in sequence hash table.\n");
    fprintf(fp, "    -e             Extra seq compression: 16-bit vs 8-bit counters.\n");
    fprintf(fp, "    -t INT         Thread num for multi-threading, default as [1]\n");
    fprintf(fp, "    -P             Disable multi-threadin\n\n");

    fprintf(fp, "    -X             Disable generation/verification of check sums\n\n");

    fprintf(fp, "To decompress:\n   SeqArc -d [ref.fa] <compress_prefix> <fastq_prefix>\n");

    exit(err);
}

int main(int argc, char **argv) {
    int opt, i, max_len = INT_MAX, qual_sys = 2;
    int block_num = 50, block_size;
    int max_readLen = 255;
    float min_alignratio = 0.5;
    int se_mark = 1;
    uint64_t max_intv = 0;
    int seq1l, seq2l, seqm;
    kseq seq1, seq2;
    bwtint_t k;
    gzFile fp1, fp2;
    smem_i *itr = NULL;
    const bwtintv_v *a;
    bwaidx_t *idx = NULL;

    int decompress = 0, indexing = 0;
    char *ref;
    int thread_num = 1;
    fqz_params p;

    /* Initialise and parse command line arguments */
    p.slevel = 3;
    p.qlevel = 2;
    p.nlevel = 1;
    p.both_strands = 0;
    p.extreme_seq = 0;
    p.multi_seq_model = 0;
    p.qual_approx = 0;
    p.do_threads = 1;
    p.do_hash = 1;

    while ((opt = getopt(argc, argv, "l:w:I:f:m:q:s:hdQ:S:N:bePXiB:t:n:c:E:r:")) != -1) {
        switch (opt) {
            case 'h':
                usage(0);

            case 'i':
                indexing = 1;
                break;

            case 'd':
                decompress = 1;
                break;

            case 'l':
                min_len = atoi(optarg);
                break;

            case 'w':
                max_iwidth = atoi(optarg);
                break;

            case 'I':
                max_intv = atol(optarg);
                break;

            case 'f':
                lgst_num = atoi(optarg);
                break;

            case 'm':
                max_mis = atoi(optarg);
                break;

            case 'c':
                max_smem_num = atoi(optarg);
                break;

            case 'E':
                exp_mismatch = atoi(optarg);
                break;

            case 'B':
                block_num = atoi(optarg);
                if (block_num < 2) block_num = 2;
                break;

            case 'q':
                qual_sys = atoi(optarg);
                break;

            case 's':
                max_insr = atoi(optarg);
                break;

            case 'r':
                min_alignratio = (float) atof(optarg);
                break;

            case 'S':
                char *end;
                p.slevel = strtol(optarg, &end, 10);
                if (p.slevel < 1 || p.slevel > 9)
                    usage(1);
                if (*end == '+')
                    p.multi_seq_model = 1;
                break;

            case 'N':
                p.qlevel = atoi(optarg);
                if (p.qlevel < 1 || p.qlevel > 3)
                    usage(1);
                break;

            case 'n':
                p.nlevel = atoi(optarg);
                if (p.nlevel < 1 || p.nlevel > 2)
                    usage(1);
                break;

            case 'b':
                p.both_strands = 1;
                break;

            case 'e':
                p.extreme_seq = 1;
                break;

            case 'P':
                p.do_threads = 0;
                break;

            case 'X':
                p.do_hash = 0;
                break;

            case 't':
                thread_num = atoi(optarg);
                break;

            default:
                usage(1);
        }
    }

    if (argc == optind)
        usage(1);

    if (indexing){
        ref = argv[optind];
        bwa_idx_build(ref, ref, BWTALGO_AUTO, 10000000);
        bwaidx_t *idx = bwa_idx_load_from_disk(ref, BWA_IDX_ALL);
        bwa_shm_stage(idx, ref, NULL);
        return 1;
    }
    else if (decompress) {
        bool hasFA = false;
        fstream fa, in_s, in_iq, in_isq, out1, out2;
        stringstream strInputPath;

        strInputPath.str("");
        strInputPath << argv[optind];

        if (strInputPath.str().substr(strInputPath.str().size()-2) == "fa" || (strInputPath.str().size() >= 5  && strInputPath.str().substr(strInputPath.str().size()-5) == "fasta")) { //with fasta
            hasFA = true;
            fa.open(argv[optind], std::ios::in);
            ++optind;
        }

        strInputPath.str("");
        strInputPath << "./" << argv[optind] << "_isq.arc";
        in_isq.open(strInputPath.str(), std::ios::binary|std::ios::in);

        unsigned char magic_fqz[9];
        if (9 != xget(in_isq, magic_fqz, 9)) {
            fprintf(stderr, "Abort: truncated read.\n");
            return 1;
        }
        if (memcmp(".arc", magic_fqz, 4) != 0) {
            fprintf(stderr, "Unrecognised file format.\n");
            return 1;
        }
        if (magic_fqz[4] != MAJOR_VERS || magic_fqz[5] != FORMAT_VERS) {
            fprintf(stderr, "Unsupported file format version %d.%d\n", magic_fqz[4], magic_fqz[5]);
            return 1;
        }
        fqzall = magic_fqz[6] & 1;
        se_mark = (magic_fqz[6] >> 1) & 1;
        p.slevel = magic_fqz[7] & 0x0f;
        p.qlevel = ((magic_fqz[7] >> 4) & 3);
        p.nlevel = (magic_fqz[7] >> 6);
        p.both_strands    = magic_fqz[8] & 1;
        p.extreme_seq     = magic_fqz[8] & 2;
        p.multi_seq_model = magic_fqz[8] & 4;

        if (!hasFA && !fqzall){
            fprintf(stderr, "FASTA needed\n");
            return 1;
        }
        if (p.slevel > 9 || p.slevel < 1) {
            fprintf(stderr, "Unexpected quality compression level %d\n", p.slevel);
            return 1;
        }
        if (p.qlevel > 3 || p.qlevel < 1) {
            fprintf(stderr, "Unexpected sequence compression level %d\n", p.qlevel);
            return 1;
        }
        if (p.nlevel > 2 || p.nlevel < 1) {
            fprintf(stderr, "Unexpected sequence compression level %d\n", p.qlevel);
            return 1;
        }

        strInputPath.str("");
        strInputPath << "./" << argv[optind+1] << "_1.fastq";
        out1.open(strInputPath.str(), std::ios::out);
        if (!se_mark){
            strInputPath.str("");
            strInputPath << "./" << argv[optind+1] << "_2.fastq";
            out2.open(strInputPath.str(), std::ios::out);
        }

        std::vector<string> vec_name;
        std::vector<string> vec_seq;
        std::vector<string> vec_qual;
        string str_seq;
        auto name = vec_name.begin();
        auto seq  = vec_seq.begin();
        auto qual = vec_qual.begin();
        i = 0;

        fqz *f = new fqz(&p);
        char outbuffer[BLK_SIZE*2];

        if(fqzall)
            cout << "Only fqz detected." << endl;
        else
            cout << "Fqz & alignInfo detected." << endl;

        while (true){
            if (se_mark){
                if (i == vec_qual.size()) {
                    i = 0;
                    vec_name.clear();
                    vec_seq.clear();
                    vec_qual.clear();
                    if (-1 == f->isq_decode(in_isq, vec_name, vec_seq, vec_qual))
                        break;
                    name = vec_name.begin();
                    seq  = vec_seq.begin();
                    qual = vec_qual.begin();
                }
                readModify2(*seq, *qual, qual_sys);
                out1 << *name << endl << *seq << endl << "+" << endl << *qual << endl;
                ++i;
                ++name;
                ++seq;
                ++qual;
            }
            else{
                if (i == vec_qual.size()){
                    i = 0;
                    vec_name.clear();
                    vec_seq.clear();
                    vec_qual.clear();
                    if (-1 == f->isq_decode(in_isq, vec_name, vec_seq, vec_qual))
                        break;
                    name = vec_name.begin();
                    seq  = vec_seq.begin();
                    qual = vec_qual.begin();
                }
                readModify2(*seq, *qual, qual_sys);
                out1 << *name << endl << *seq << endl << "+" << endl << *qual << endl;
                ++i;
                ++name;
                ++seq;
                ++qual;

                if (i == vec_qual.size()){
                    i = 0;
                    vec_name.clear();
                    vec_seq.clear();
                    vec_qual.clear();
                    if (-1 == f->isq_decode(in_isq, vec_name, vec_seq, vec_qual))
                        break;
                    name = vec_name.begin();
                    seq  = vec_seq.begin();
                    qual = vec_qual.begin();
                }
                readModify2(*seq, *qual, qual_sys);
                out2 << *name << endl << *seq << endl << "+" << endl << *qual << endl;
                ++i;
                ++name;
                ++seq;
                ++qual;
            }
        }
        in_isq.close();
        vec_seq.clear();
        vec_name.clear();
        vec_qual.clear();

        if (fqzall){
            delete f;
            out1.close();
            out2.close();
            return 0;
        }

        strInputPath.str("");
        strInputPath << "./" << argv[optind] << "_aligninfo.arc";
        in_s.open(strInputPath.str(), std::ios::binary|std::ios::in);

        strInputPath.str("");
        strInputPath << "./" << argv[optind] << "_iq.arc";
        in_iq.open(strInputPath.str(), std::ios::binary|std::ios::in);

        unsigned char magic_s[11];
        if (11 != xget(in_s, magic_s, 11)){
            fprintf(stderr, "Abort: truncated read.\n");
            return 1;
        }
        qual_sys = magic_s[0];
        memcpy(&block_size, magic_s+1, 4);
        block_num = magic_s[5];
        max_mis = magic_s[6];
        memcpy(&max_insr, magic_s+7, 2);
        memcpy(&max_readLen, magic_s+9, 2);

        decode seqdecoder(block_size, max_mis, max_insr, max_readLen);
        ref2seq ref2seqer(block_size, max_mis, max_readLen, fa);

        align_info align_info1;
        align_info1.cigar_l = (int*) malloc(max_mis*sizeof(int));
        align_info1.cigar_v = (int*) malloc(max_mis*sizeof(int));
        int readLen, thisblock=0;
        bool blockjump=false;

        i = 0;
        f = new fqz(&p);
        while (thisblock < block_num){
            if (se_mark){
                if (seqdecoder.parse_se(align_info1, readLen, in_s) == 1){ //request a read from vec
                    if (i == vec_qual.size()) {
                        //FQbuffer out
                        i = 0;
                        vec_name.clear();
                        vec_qual.clear();
                        f->iq_decode(in_iq, vec_name, vec_qual);
                        name = vec_name.begin();
                        qual = vec_qual.begin();
                    }
                    str_seq = ref2seqer.getSeq(align_info1, readLen);
                    readModify2(str_seq, *qual, qual_sys);
                    out1 << *name << endl << str_seq << endl << "+" << endl << *qual << endl; //换成输出到buffer
                    ++i;
                    ++name;
                    ++qual;
                    blockjump = false;
                }
                else{
                    ++thisblock;
                    if (blockjump) //the second continuous 0, means the block contains no read
                        in_iq.seekg(52, std::ios::cur);
                    else{
                        //FQbuffer out
                        f = new fqz(&p);
                        if (i < vec_qual.size())  //if align_info and iq not paired
                            return -1;
                        blockjump = true;
                    }
                }
            }
            else{
                if (seqdecoder.parse_pe(align_info1, readLen, in_s) == 0){
                    if (blockjump){
                        ++thisblock;
                        in_iq.seekg(52, std::ios::cur);
                        continue;
                    }
                    else{
                        ++thisblock;
                        f = new fqz(&p);
                        if (i < vec_qual.size()) //if align_info and iq not paired
                            return -1;
                        blockjump = true;
                        continue;
                    }
                }
                else{
                    if (i == vec_qual.size()) {
                        i = 0;
                        vec_name.clear();
                        vec_qual.clear();
                        f->iq_decode(in_iq, vec_name, vec_qual);
                        name = vec_name.begin();
                        qual = vec_qual.begin();
                    }
                    str_seq = ref2seqer.getSeq(align_info1, readLen);
                    readModify2(str_seq, *qual, qual_sys);
                    out1 << *name << endl << str_seq << endl << "+" << endl << *qual << endl;
                    ++i;
                    ++name;
                    ++qual;
                    blockjump = false;
                }

                if (seqdecoder.parse_pe(align_info1, readLen, in_s) == 1){ //request a read from vec
                    if (i == vec_qual.size()) {
                        i = 0;
                        vec_name.clear();
                        vec_qual.clear();
                        f->iq_decode(in_iq, vec_name, vec_qual);
                        name = vec_name.begin();
                        qual = vec_qual.begin();
                    }
                    str_seq = ref2seqer.getSeq(align_info1, readLen);
                    readModify2(str_seq, *qual, qual_sys);
                    out2 << *name << endl << str_seq << endl << "+" << endl << *qual << endl;
                    ++i;
                    ++name;
                    ++qual;
                }
                else{
                    cout << "jump-block shows up between read1 and read2, which is impossible." << endl;
                    return -1;
                }
            }
        }
        delete f;
        fa.close();
        in_s.close();
        in_iq.close();
        out1.close();
        out2.close();
        return 0;

    } else {
        struct timeval timeStart,timeEnd;
        gettimeofday(&timeStart, NULL);
        string tmpArg, outIndex; // 可以选择放弃comment内容

        tmpArg = argv[optind];
        if (tmpArg.substr(tmpArg.size()-2) != "fa" && tmpArg.substr(tmpArg.size()-5) != "fasta"){ //no fasta
            cout << "FASTA index missing." << endl;
            fqzall = 1;
        }
        else{
            idx = bwa_idx_load_from_shm(argv[optind]);
            if (idx == 0) {
                if ((idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1;
                bwa_shm_stage(idx, argv[optind], NULL);
            }
            itr = smem_itr_init(idx->bwt);
            smem_config(itr, 1, max_len, max_intv); //min_intv = 1
            block_num = count_blocknum(idx->bns->l_pac, block_num);
            block_size = (int) ceil(idx->bns->l_pac/block_num); //单个block的长度
            ++optind;
        }

        fp1 = xzopen(argv[optind], "r");
        FunctorZlib gzr;
        kstream<gzFile, FunctorZlib> ks1(fp1, gzr);
        kstream<gzFile, FunctorZlib> *ks2;

        if (argc - optind >= 3){
            se_mark = 0;
            fp2 = xzopen(argv[optind + 1], "r");
            ks2 = new kstream<gzFile, FunctorZlib> (fp2,gzr);
            outIndex = argv[optind + 2];
        } else
            outIndex = argv[optind + 1];

        CreatBitmap(nucleBitMap);

        gettimeofday(&timeEnd, NULL);
        double total = (timeEnd.tv_sec - timeStart.tv_sec)*1000 + (timeEnd.tv_usec -timeStart.tv_usec)*1.0/1000;
        printf("%s%d init time %f ms\n", __FUNCTION__, __LINE__, total);

        //Pre-Align
        align_info align_info1;
        align_info align_info2;
        align_info1.cigar_l = (int*)malloc(max_mis * sizeof(int));
        align_info1.cigar_v = (int*)malloc(max_mis * sizeof(int));
        align_info2.cigar_l = (int*)malloc(max_mis * sizeof(int));
        align_info2.cigar_v = (int*)malloc(max_mis * sizeof(int));

        bwtintv_v *matcher1 = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
        bwtintv_v *matcher2 = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
        bwtintv_v *tmpvec[2];
        tmpvec[0] = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
        tmpvec[1] = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));

        stringstream strOutputPath;
        fstream out_isq, out_s, out_iq; //声明id+seq+qual的输出文件指针
        strOutputPath.str("");
        strOutputPath << "./" << outIndex << "_isq.arc";
        out_isq.open(strOutputPath.str(), std::ios::binary|std::ios::out);

        idx1 = new int64_t*[max_smem_num*max_iwidth];
        for(i=0; i<max_smem_num*max_iwidth; i++)
            idx1[i] = new int64_t[3];
        idx2 = new int64_t*[max_smem_num*max_iwidth];
        for(i=0; i<max_smem_num*max_iwidth; i++)
            idx2[i] = new int64_t[3];

        int alignedNum = 0, prealign_num = PREALIGN_NUM; //set prealign_num for files with reads less than PREALIGN_NUM
        if (!fqzall) {
            for (i = 0; i < PREALIGN_NUM; i++) {
                if (se_mark) {
                    seq1l = ks1.read(seq1);
                    if (seq1l == -1){
                        prealign_num = i;
                        break;
                    }
                    seqm = getAlignInfoSE(seq1, matcher1, tmpvec, itr, idx, align_info1, block_size);
                    if (seqm > 0)
                        ++alignedNum;
                } else {
                    seq1l = ks1.read(seq1);
                    seq2l = (*ks2).read(seq2);
                    if (seq2l == -1){
                        prealign_num = i;
                        break;
                    }
                    seqm = getAlignInfoPE(seq1, seq2, matcher1, matcher2, tmpvec, itr, idx, align_info1, align_info2, block_size);
                    if (seqm > 0)
                        ++alignedNum;
                }
            }
            ks1.rewind();
            seq1.last_char = 0;
            if (!se_mark){
                (*ks2).rewind();
                seq2.last_char = 0;
            }
        }
        if (!fqzall && alignedNum < (prealign_num * min_alignratio)){
            cout << "Compression begins. Only fqzcomp adopted." << endl;
            fqzall = 1;
        }

        int level1 = fqzall | (se_mark << 1);
        int level2 = p.slevel | (p.qlevel << 4) | (p.nlevel << 6);
        int flags = p.both_strands
                    + p.extreme_seq*2
                    + p.multi_seq_model*4;
        unsigned char magic_fqz[9] = {'.', 'a', 'r', 'c',  //生成magic作为解压时参数
                                      MAJOR_VERS,
                                      FORMAT_VERS,
                                      (uint8_t)level1,
                                      (uint8_t)level2,
                                      (uint8_t)flags
        };
        if (9 != xwrite(out_isq, magic_fqz, 9)) {
            fprintf(stderr, "Abort: truncated write.2\n");
            out_isq.close();
            return 1;
        }

        gettimeofday(&timeEnd, NULL);
        total = (timeEnd.tv_sec - timeStart.tv_sec)*1000 + (timeEnd.tv_usec -timeStart.tv_usec)*1.0/1000;
        printf("%s%d load time %f ms\n", __FUNCTION__, __LINE__, total);

        unsigned char magic_s[11]{(uint8_t)qual_sys, //block_size+max_mis+max_insr+max_readLen
                                  (uint8_t)(block_size & 0xffff),
                                  (uint8_t)((block_size>>8) & 0xffff),
                                  (uint8_t)((block_size>>16) & 0xffff),
                                  (uint8_t)((block_size>>24) & 0xffff),
                                  (uint8_t)block_num,
                                  (uint8_t)max_mis,
                                  (uint8_t)(max_insr & 0xffff),
                                  (uint8_t)((max_insr>>8) & 0xffff),
                                  (uint8_t)(max_readLen & 0xffff),
                                  (uint8_t)((max_readLen>>8) & 0xffff)
        };
        if(fqzall) //预比对，匹配率低于设定值
            block_num = 0;
        else{
            strOutputPath.str("");
            strOutputPath << "./" << outIndex << "_aligninfo.arc";
            out_s.open(strOutputPath.str(), std::ios::binary|std::ios::out);
            if (11 != xwrite(out_s, magic_s, 11)) {
                fprintf(stderr, "Abort: truncated write.2\n");
                out_s.close();
                return 1;
            }

            strOutputPath.str("");
            strOutputPath << "./" << outIndex << "_iq.arc";
            out_iq.open(strOutputPath.str(), std::ios::binary|std::ios::out);
        }

        fqz* f[block_num+1]; //one more for isq
        for (i=0;i<block_num+1;i++)
            f[i] = new fqz(&p);

        fstream fpOutput_s[block_num]; //声明align_info的输出文件指针数组
        for (i=0;i<block_num;i++) {
            strOutputPath.str("");
            strOutputPath << "./s_" << i << ".tmp"; //合成文件路径
            fpOutput_s[i].open(strOutputPath.str(), std::ios::binary|std::ios::out);
        }

        fstream fpOutput_iq[block_num]; //声明id+qual的输出文件指针数组
        for (i=0;i<block_num;i++) {
            strOutputPath.str("");
            strOutputPath << "./iq_" << i << ".tmp"; //合成文件路径
            fpOutput_iq[i].open(strOutputPath.str(), std::ios::binary|std::ios::out);
        }

        encode* encoders[block_num];
        for (i = 0; i < block_num; i++)
            encoders[i] = new encode(se_mark, block_size, max_mis, max_insr, max_readLen);

#ifdef SEQARCTHREAD
        pthread_mutex_init(&g_mutex, 0);
        pthread_mutex_init(&g_write_mutex, 0);
        sem_init(&g_sem,0,0);

        ThreadParam *pParam = new ThreadParam;
        if (pParam)
        {
            pParam->block_size = block_size;
            pParam->block_num = block_num;
            pParam->pidx = idx;
            pParam->pitr = itr;
            pParam->pout_isq = &out_isq;
            pParam->pfqz = f;
            if (!fqzall)
            {
                pParam->pencoders = encoders;
                pParam->pfpOutput_iq = fpOutput_iq;
                pParam->pfpOutput_s = fpOutput_s;
            }
        }
        else
        {
            return -1;
        }

        pthread_t *tid = (pthread_t*)alloca(thread_num * sizeof(pthread_t));
        for (i = 0; i < thread_num; ++i) //创建一个线程池，等待任务
        {
            pthread_create(&tid[i], 0, task_process, pParam);
        }

        if (se_mark) //SE
        {
            while (ks1.read(seq1) >= 0)
            {
                readModify1(seq1.seq, seq1.qual, qual_sys, max_readLen);
                Task task;
                task.num = 1;
                task.seq1 = seq1;
                pthread_mutex_lock(&g_mutex);
                g_task_queue.push(task); //添加任务
                pthread_mutex_unlock(&g_mutex);
                sem_post(&g_sem);
            }
        }
        else //PE
        {
            while(ks1.read(seq1) >= 0 && (*ks2).read(seq2) >= 0)
            {
                readModify1(seq1.seq, seq1.qual, qual_sys, max_readLen);
                readModify1(seq2.seq, seq2.qual, qual_sys, max_readLen);
                Task task;
                task.num = 2;
                task.seq1 = seq1;
                task.seq2 = seq2;
                pthread_mutex_lock(&g_mutex);
                g_task_queue.push(task); //添加任务
                pthread_mutex_unlock(&g_mutex);
                sem_post(&g_sem);
            }
        }

        g_bFinish = true;

        for (i = 0; i < thread_num; ++i)
        {
            sem_post(&g_sem);
        }

        for (i = 0; i < thread_num; ++i)
        {
            pthread_join(tid[i], 0);
        }

        delete pParam;
        pthread_mutex_destroy(&g_mutex);
        pthread_mutex_destroy(&g_write_mutex);
        sem_destroy(&g_sem);
#else

        while ((seq1l = ks1.read(seq1)) >= 0) {
            readModify1(seq1.seq, seq1.qual, qual_sys, max_readLen);
            if (se_mark){ //SE
                if ((!fqzall) && (seqm = getAlignInfoSE(seq1, matcher1, tmpvec, itr, idx, align_info1, block_size))){
                    encoders[align_info1.blockNum]->parse_1(align_info1, seq1l, fpOutput_s[align_info1.blockNum]);
                    f[align_info1.blockNum]->iq_encode(seq1.name, seq1.qual, fpOutput_iq[align_info1.blockNum]); //id+qual
                }
                else{
                    f[block_num]->isq_encode(seq1.name, seq1.seq, seq1.qual, out_isq); //id+seq+qual
                }
            }
            else{ //PE
                seq2l = (*ks2).read(seq2);
                readModify1(seq2.seq, seq2.qual, qual_sys, max_readLen);
                if ((!fqzall) && (seqm = getAlignInfoPE(seq1, seq2, matcher1, matcher2, tmpvec, itr, idx, align_info1, align_info2, block_size))){
                    encoders[align_info1.blockNum]->parse_1(align_info1, seq1l, fpOutput_s[align_info1.blockNum]);
                    encoders[align_info2.blockNum]->parse_2(align_info2, seq2l, fpOutput_s[align_info2.blockNum]);
                    f[align_info1.blockNum]->iq_encode(seq1.name, seq1.qual, fpOutput_iq[align_info1.blockNum]);
                    f[align_info2.blockNum]->iq_encode(seq2.name, seq2.qual, fpOutput_iq[align_info2.blockNum]);
                }
                else{//没比上
                    f[block_num]->isq_encode(seq1.name, seq1.seq, seq1.qual, out_isq);
                    f[block_num]->isq_encode(seq2.name, seq2.seq, seq2.qual, out_isq);
                }
            }
        }
        free(matcher1->a);
        free(matcher1);
        free(matcher2->a);
        free(matcher2);
        free(tmpvec[0]->a);
        free(tmpvec[0]);
        free(tmpvec[1]->a);
        free(tmpvec[1]);
#endif
        //收尾
        string nullstr;
        for (i=0;i<block_num;i++){
            encoders[i]->end(fpOutput_s[i]);
            f[i]->iq_encode(nullstr, nullstr, fpOutput_iq[i]);
        }
        f[block_num]->isq_encode(nullstr, nullstr, nullstr, out_isq);

        if(idx)
            bwa_idx_destroy(idx);
        if(itr)
            smem_itr_destroy(itr);
        err_gzclose(fp1);
        if (!se_mark)
            err_gzclose(fp2);
        out_isq.close();

        //merge fpOutput_iq to out_iq
        for (i = 0; i < block_num; i++)
        {
            fpOutput_iq[i].close();
            stringstream str_tmp;
            str_tmp << "./iq_" << i << ".tmp";
            fstream f_iq;
            f_iq.open(str_tmp.str(), ios::in |ios::binary);
            out_iq << f_iq.rdbuf();
            f_iq.close();
            std::remove(str_tmp.str().c_str()); //delete the tmp file
        }
        if(!fqzall)
            out_iq.close();

        //merge fpOutput_s to out_s
        for (i = 0; i < block_num; i++)
        {
            fpOutput_s[i].close();
            stringstream str_tmp;
            str_tmp << "./s_" << i << ".tmp";
            fstream f_s;
            f_s.open(str_tmp.str(), ios::in | ios::binary);
            out_s << f_s.rdbuf();
            f_s.close();
            std::remove(str_tmp.str().c_str()); //delete the tmp file
        }
        if(!fqzall)
            out_s.close();

        gettimeofday(&timeEnd, NULL);
        total = (timeEnd.tv_sec - timeStart.tv_sec)*1000 + (timeEnd.tv_usec -timeStart.tv_usec)*1.0/1000;
        printf("%s%d total time %f ms\n", __FUNCTION__, __LINE__, total);
        return 0;
    }
}
