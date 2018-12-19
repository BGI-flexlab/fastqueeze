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
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h> 
#include <malloc.h>
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


pthread_mutex_t g_mutex;  //读文件锁

//Align Params
int min_len = 17, max_iwidth = 50, max_mis = 3, lgst_num = 2, max_smem_num = 2, exp_mismatch = 1,  fqzall = 0;

#define READLENGTH 8*1024

typedef struct  _tagEncodeParam
{
    uint64_t offset;
    uint64_t length;
    char filename[256];
    int num;
    smem_i *pitr;
    bwaidx_t *pidx;
}EncodeParam;

class SeqRead
{
public:
    SeqRead(EncodeParam *param);
    ~SeqRead();
    void init();
    bool getRead();

    bool m_isEof;
    char m_lastchar;
    //gzFile m_f;
    FILE *m_f;
    char *m_pbuffer;
    int index;
    int m_step;
    uint64_t m_length;
    int m_begin,m_end;
    char name[256];
    char seq[256];
    char comment[256];
    char qual[256];
};

SeqRead::SeqRead(EncodeParam *param)
{
    m_isEof = false;
    m_begin = 0;
    m_end = 0;
    //m_f = gzopen(param->filename,"r");
    m_f = fopen(param->filename,"r");
    //gzseek(m_f,param->offset,SEEK_SET);
    fseek(m_f,param->offset,SEEK_SET);
    m_length = param->length;
    m_pbuffer = (char*)malloc(READLENGTH);
}

SeqRead::~SeqRead()
{
    free(m_pbuffer);
    //err_gzclose(m_f);
    fclose(m_f);
}

void SeqRead::init()
{
    memset(name, 0, 256);
    memset(seq, 0, 256);
    memset(comment, 0, 256);
    memset(qual, 0, 256);
}

bool SeqRead::getRead()
{
    init();
    index=0;
    m_step = 1;

READ:
    if(m_begin >= m_end-1)
    {
        if(m_isEof)
        {
            return strlen(qual);
        }

        if(m_lastchar == '\n')//判断如何拼接
        {
            m_step++;
            index=0;
        }

        m_begin = 0;
        int len_read = m_length < READLENGTH ? m_length : READLENGTH;
        //m_end =  gzread(m_f, m_pbuffer,len_read);
        m_end = fread(m_pbuffer, 1, len_read, m_f);
        m_length -= m_end;
        m_lastchar = m_pbuffer[m_end-1];
        if(m_end < READLENGTH) //到达结尾
        {
            //printf("%d\n", m_length);
            m_isEof = true;
        }
    }
    
    if(m_step == 1)
    {
        while(m_pbuffer[m_begin] != '\n') //get name
        {
            if(m_begin < m_end)
            {
                name[index++] = m_pbuffer[++m_begin];
            }
            else
            {
                goto READ;
            }
        }
        if(m_begin >= m_end-1)
        {
            goto READ;
        }

        name[index-1] = '\0';
        m_step++;
        index=0;
        m_begin++;
    }



    if(m_step == 2)
    {
        while(m_pbuffer[m_begin] != '\n') //get name
        {
            if(m_begin < m_end)
            {
                seq[index++] = m_pbuffer[m_begin++];
            }
            else
            {
                goto READ;
            }
        }
        if(m_begin >= m_end-1)
        {
            goto READ;
        }
        m_step++;
        index=0;
        m_begin++;
    }

    if(m_step == 3)
    {
        while(m_pbuffer[m_begin] != '\n') //get name
        {
            if(m_begin < m_end)
            {
                comment[index++] = m_pbuffer[m_begin++];
            }
            else
            {
                goto READ;
            }
        }
        if(m_begin >= m_end-1)
        {
            goto READ;
        }
        m_step++;
        index=0;
        m_begin++;
    }

    if(m_step == 4)
    {
        while(m_pbuffer[m_begin] != '\n') //get name
        {
            if(m_begin < m_end)
            {
                qual[index++] = m_pbuffer[m_begin++];
            }
            else
            {
                goto READ;
            }
        }
        if(m_begin >= m_end-1)
        {
            goto READ;
        }
        m_begin++;
    }

    if(strlen(seq) == strlen(qual))
    {
        return true;
    }

    return false;
}

typedef struct _tagCigar
{
    unsigned short pos:13; //相对于read的pos
    unsigned short ch:3; //碱基字符
}Cigar;

typedef struct _tagAlignInfo
{
    uint64_t pos:48;//相对于参考序列的pos
    uint64_t ilenth: 12; //匹配上的长度
    uint64_t icount : 3; //差异的个数
    uint64_t Isrev : 1; //是否是反向，ture是反向
}AlignInfo;

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

int getAlignInfoSE(char *seq, bwtintv_v *a, bwtintv_v *tmpvec[2], smem_i* func_itr, bwaidx_t *func_idx, align_info &align_info1){
    int64_t rlen;
    int i, seql, base;

    seql = strlen(seq);
    int pass_num = 0;
    int cigar_l[max_mis], cigar_v[max_mis];
    for(i=0;i<max_mis;i++)
    {
        align_info1.cigar_l[i] = -1;
        align_info1.cigar_v[i] = -1;
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
                            align_info1.blockNum = 0;
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
                            align_info1.blockNum = 0;
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

int getbitnum(int data)
{
    int count = 0;
    while (data > 0) {
        data >>= 1;
        count++;
    }
    return count;
}

int myint2bit(uint32_t data, int limit, std::string &str)
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
        case 0: buf.append("0"); break;
        case 1: buf.append("10"); break;
        case 2: buf.append("11"); break;
    }
}

void AlignInfoToBitArry(align_info &info, int readLen, string &bitbuf)
{
    bitbuf.clear();
    std::string stremp;
    myint2bit(info.blockPos, 32, stremp);
    bitbuf.append(stremp);

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

    //printf("%ld %d %d\n", info.blockPos, info.isRev, cigar_num);
    //printf("%d %d %d\n",info.cigar_l[0], info.cigar_l[1],info.cigar_l[2]);
    //printf("%d %d %d\n",info.cigar_v[0], info.cigar_v[1],info.cigar_v[2]);

    myint2bit(cigar_num, 2, stremp);
    bitbuf.append(stremp);
    bitbuf.append(str_l);
    bitbuf.append(str_v);
}

void BitArryToBuf(const char *pbit, int len, char *pbuf, int *buflen)
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

void *encode_process(void *data)
{
    if(data == NULL) return data;

    //pthread_mutex_lock(&g_mutex);

    EncodeParam *param = (EncodeParam*)data;
    printf("---%0x %ld %ld %ld\n", pthread_self(), param->offset, param->length, param->num);
    SeqRead ssread(param);

    bwtintv_v *matcher1 = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
    bwtintv_v *tmpvec[2];
    tmpvec[0] = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
    tmpvec[1] = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));

    fqz_params p;
    p.slevel = 3;
    p.qlevel = 2;
    p.nlevel = 1;
    p.both_strands = 0;
    p.extreme_seq = 0;
    p.multi_seq_model = 0;
    p.qual_approx = 0;
    p.do_threads = 1;
    p.do_hash = 1;

    fqz *pfqz = new fqz(&p);

    fstream out_isq; 
    char path[64]={0};
    sprintf(path, "./out_isq_%d.arc",param->num);
    out_isq.open(path, std::ios::binary|std::ios::out);

    align_info align_info1;
    align_info1.cigar_l = (int*)malloc(max_mis * sizeof(int));
    align_info1.cigar_v = (int*)malloc(max_mis * sizeof(int));

    string bitbuf;
    while(ssread.getRead())
    {
        char outbuf[100]={0};
        int buf_len = 0;
        int seqm = getAlignInfoSE(ssread.seq, matcher1, tmpvec, param->pitr, param->pidx, align_info1);
        if(seqm) //比对成功
        {
            AlignInfoToBitArry(align_info1, strlen(ssread.seq), bitbuf);
            BitArryToBuf(bitbuf.c_str(), bitbuf.length(), outbuf, &buf_len);
            pfqz->isq_encode(ssread.name, strlen(ssread.name), outbuf, buf_len, ssread.qual,strlen(ssread.qual), out_isq);
        }
        else
        {
            pfqz->isq_encode(ssread.name, strlen(ssread.name), ssread.seq, strlen(ssread.seq), ssread.qual, strlen(ssread.qual), out_isq);
        }       
    }
    pfqz->isq_encode(NULL, 0,NULL,0, NULL,0, out_isq);
    
    delete param;
    delete pfqz;

    //pthread_mutex_unlock(&g_mutex);
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

class BitDecode
{
public:
    BitDecode(unsigned char *buf){m_pbuf = buf;char_index = 8;}
    ~BitDecode(){};
    uint32_t bufferIn(int length);
    int char_index;
    unsigned char *m_pbuf;
};


uint32_t BitDecode::bufferIn(int length)
{
    uint32_t pass_bit = 0;
    uint32_t out_int = 0;
    unsigned char rem_char;
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
            }
        }
        case 1:
        {
            switch(id)
            {
                case 0: return 2;
                case 1: return 0;
                case 2: return 3;
            }
        }
        case 2:
        {
            switch(id)
            {
                case 0: return 1;
                case 1: return 0;
                case 2: return 3;
            }
        }
        case 3:
        {
            switch(id)
            {
                case 0: return 2;
                case 1: return 1;
                case 2: return 0;
            }
        }
    }
}

typedef struct _tagDecodeParam
{
    int num;
    bwaidx_t *pidx;
}DecodeParam;

void *decode_process(void *data)//(bwaidx_t *idx)
{
    if(data == NULL) return data;

    DecodeParam *param = (DecodeParam*)data;

    fqz_params p;
    p.slevel = 3;
    p.qlevel = 2;
    p.nlevel = 1;
    p.both_strands = 0;
    p.extreme_seq = 0;
    p.multi_seq_model = 0;
    p.qual_approx = 0;
    p.do_threads = 1;
    p.do_hash = 1;

    fqz *pfqz = new fqz(&p);
    fstream in_isq;
    char path[125]={0};
    sprintf(path,"./out_isq_%d.arc", param->num);
    in_isq.open(path, std::ios::binary|std::ios::in);

    fstream original; //解压输出到文件
    char tmpath[125]={0};
    sprintf(tmpath,"./original_%d.arc", param->num);
    original.open(tmpath, std::ios::binary|std::ios::out);


    char *write_buf = NULL;
    uint32_t write_buf_len = 0;
    while(true)
    {
        char *namebuf = NULL;
        char *seqbuf = NULL;
        char *qualbuf = NULL;
        int *seqlen = NULL;
        int *quallen = NULL;
        int count = 0;
        int ret = pfqz->isq_decode(in_isq, &namebuf, &seqbuf, &qualbuf, &seqlen, &quallen, &count);
        if(ret < 0) break;

        uint32_t tmplen = count*5*quallen[0];
        if(write_buf_len < tmplen)
        {
            write_buf_len = tmplen;
            write_buf = (char*)realloc(write_buf, write_buf_len);
        }


        char buf[256]={0};
        uint32_t out_ind = 0;
        for (int k = 0; k < count; k++) 
        {
            while ((write_buf[out_ind++] = *namebuf++) != '\n');

            memset(buf, 0, 256);
            if(seqlen[k] < quallen[k])
            {
                memcpy(buf, seqbuf, seqlen[k]); seqbuf += seqlen[k];
                BitDecode bitdecode((unsigned char *)buf);
                uint64_t pos = bitdecode.bufferIn(32);
                bool isrev = bitdecode.bufferIn(1);
                int cigar_num = bitdecode.bufferIn(2);
                //printf("----%ld %d %d \n", pos, isrev, cigar_num);
                int cigar_l[3] = {-1};
                int cigar_v[3] = {-1};
                int offset = 0;
                for(int i=0;i<cigar_num;i++) 
                {
                    cigar_l[i] = bitdecode.bufferIn(getbitnum(quallen[k]-offset))+offset;
                    offset = cigar_l[i];
                }
                //printf("%d %d %d \n", cigar_l[0],cigar_l[1],cigar_l[2]);
                for(int i=0;i<cigar_num;i++)
                {
                    if (!bitdecode.bufferIn(1))
                        cigar_v[i] = 0;
                    else{
                        if (!bitdecode.bufferIn(1))
                            cigar_v[i] = 1;
                        else
                            cigar_v[i] = 2;
                    }
                }
                //printf("%d %d %d \n", cigar_v[0],cigar_v[1],cigar_v[2]);

                int64_t rlen = 0;
                uint8_t *rseq = bns_get_seq(param->pidx->bns->l_pac, param->pidx->pac, pos, pos + quallen[k], &rlen);
                
                for (int i=0; i< cigar_num; i++){
                    if (cigar_l[i] != -1){
                        rseq[cigar_l[i]] = mapvar2base(rseq[cigar_l[i]], cigar_v[i]);
                    }
                    else
                        break;
                }

                int j = 0;
                for(int i=0;i<quallen[k];i++) 
                {
                    buf[i] = isrev ? getch(rseq[rlen-1-j],true): getch(rseq[j],false);
                    j++;
                }
                free(rseq);
                memcpy(&write_buf[out_ind], buf, quallen[k]);out_ind += quallen[k];
            }
            else
            {
                memcpy(&write_buf[out_ind], seqbuf, seqlen[k]); seqbuf += seqlen[k];
                out_ind += seqlen[k];
            }

            write_buf[out_ind++] = '\n';
            write_buf[out_ind++] = '+';
            write_buf[out_ind++] = '\n';
            memcpy(&write_buf[out_ind], qualbuf, quallen[k]);qualbuf += quallen[k];out_ind+=quallen[k];
            write_buf[out_ind++] = '\n';
        }

        original.write(write_buf, out_ind);
    }
    
    original.close();
    in_isq.close();
    delete pfqz;
    free(write_buf);
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

uint64_t Getoffset(uint64_t *arry, int num)
{
    int i;
    uint64_t tmp = 0;
    for(i=0;i<num ;i++)
    {
        tmp += arry[i+1];
    }
    return tmp;
}

int main(int argc, char **argv) {
    int opt, i, max_len = INT_MAX, qual_sys = 2;
    int block_num = 50, block_size;
    int max_insr = 511, max_readLen = 255;
    float min_alignratio = 0.5;
    int se_mark = 1;
    uint64_t max_intv = 0;
    int seq1l, seq2l, seq1m, seq2m;
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

    CreatBitmap(nucleBitMap);

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

        struct timeval timeStart,timeEnd;
        gettimeofday(&timeStart, NULL);
        idx = bwa_idx_load_from_shm(argv[optind]);
        if (idx == 0) {
            if ((idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1;
            bwa_shm_stage(idx, argv[optind], NULL);
        }

        pthread_mutex_init(&g_mutex, 0);
        pthread_t *t_id = (pthread_t*)alloca(thread_num * sizeof(pthread_t));

        int i;
        for (i = 0; i < thread_num; ++i) //创建一个线程池，等待任务
        {
            DecodeParam * param = new DecodeParam;
            param->num = i;
            param->pidx = idx;
            pthread_create(&t_id[i], 0, decode_process, param);
        }

        for (i = 0; i < thread_num; ++i) 
        {
            pthread_join(t_id[i], 0);
        }

        gettimeofday(&timeEnd, NULL);
        double totald = (timeEnd.tv_sec - timeStart.tv_sec)*1000 + (timeEnd.tv_usec -timeStart.tv_usec)*1.0/1000;
        printf("%s%d decode total time %f ms\n", __FUNCTION__, __LINE__, totald);

        return 0;

    } else {
        struct timeval timeStart,timeEnd;
        gettimeofday(&timeStart, NULL);
//---------------
        // uint64_t flength=7073747447;
        // FILE *ff = fopen(argv[optind+1], "rb");
        // if(ff)
        // {
        //     fseek(ff, -4, SEEK_END);
        //     uint8_t buf[4]={0};
        //     fread(buf,1,4,ff);
        //     flength = (uint64_t)buf[0]+(uint64_t)(buf[1]<<8)+
        //                 (uint64_t)(buf[2]<<16)+(uint64_t)(buf[3]<<24);
        //     fclose(ff);
        // }

        idx = bwa_idx_load_from_shm(argv[optind]);
        if (idx == 0) {
            if ((idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1;
            bwa_shm_stage(idx, argv[optind], NULL);
        }
        itr = smem_itr_init(idx->bwt);
        smem_config(itr, 1, max_len, max_intv); //min_intv = 1

        struct stat statbuf; 
        stat(argv[optind+1],&statbuf); 
        uint64_t flength=statbuf.st_size;

        
        FILE *fdna = fopen(argv[optind+1], "r");
        uint64_t len = flength/thread_num;

        uint64_t slicearry[10] ={0};
        int exlen = 0;
        for (i = 1; i < thread_num; ++i) 
        {
            char buffer[1024] = {0};
            //gzseek(fp1,i*len,SEEK_SET);
            //gzread(fp1, buffer, 1024);
            fseek(fdna, i*len, SEEK_SET);
            fread(buffer, 1, 1024, fdna);

            std::vector<std::string> out;
            const char *delim = "\n";
            char *p = strtok(buffer, delim);
            while(p)
            {
                out.emplace_back(p);
                p = strtok(NULL, delim);
            }

            slicearry[i] =len - exlen;//去除上一个切片超过的长度
            exlen = 0;
            auto itor = out.begin();
            for(;itor != out.end();itor++)
            {
                if(itor->at(0) == '@')
                {
                    itor+=2;
                    if(itor->at(0) == '+')
                    {
                        break;
                    }
                    itor-=2;
                }
                exlen += itor->length()+1;
            }

            slicearry[i] += exlen; //加上本次切片超过的长度
            flength -= slicearry[i];
        }

        slicearry[thread_num] = flength;
        fclose(fdna);

        pthread_mutex_init(&g_mutex, 0);
        pthread_t *t_id = (pthread_t*)alloca(thread_num * sizeof(pthread_t));

        for (i = 0; i < thread_num; ++i) //创建一个线程池，等待任务
        {
            EncodeParam * test = new EncodeParam;
            test->length = slicearry[i+1];
            test->offset = Getoffset(slicearry,i); 
            strcpy(test->filename, argv[optind+1]);
            test->num = i;
            test->pidx = idx;
            test->pitr = itr;
            pthread_create(&t_id[i], 0, encode_process, test);
        }

        for (i = 0; i < thread_num; ++i) 
        {
            pthread_join(t_id[i], 0);
        }

        gettimeofday(&timeEnd, NULL);
        double totald = (timeEnd.tv_sec - timeStart.tv_sec)*1000 + (timeEnd.tv_usec -timeStart.tv_usec)*1.0/1000;
        printf("%s%d encode total time %f ms\n", __FUNCTION__, __LINE__, totald);

        return 0;
    }
}
