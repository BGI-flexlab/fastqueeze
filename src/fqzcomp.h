#ifndef SEQARC_FQZ_COMP_H
#define SEQARC_FQZ_COMP_H

#define __STDC_FORMAT_MACROS

#include "sfh.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <inttypes.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <bitset>
#include <assert.h>

#define ABS(a)   ((a)>0?(a):-(a))
#define MIN(a, b) ((a)<(b)?(a):(b))
#define MAX(a, b) ((a)>(b)?(a):(b))

#include "clr.cdr"
#include "simple_model.h"
#include "base_model.h"

#define QMAX 128
#ifdef __SSE__
    #include <xmmintrin.h>
#else
    #define _mm_prefetch(a,b)
#endif

#define BLK_SIZE 10*1024*1024
#define QBITS 12
#define QSIZE (1<<QBITS)
#define DECODE_INT(a) ((a)[0] + ((a)[1]<<8) + ((a)[2]<<16) + ((a)[3]<<24))
#define DECODE_SHORT(a) ((a)[0] + ((a)[1]<<8))

typedef struct {
    int slevel;             // -s level
    int qlevel;             // -q level
    int nlevel;             // -n level
    bool both_strands;      // True if -b used
    bool extreme_seq;       // True if -e used; 16-bit seq counters
    bool multi_seq_model;   // True if -s<level>+; use 2 model sizes
    bool do_hash;           // Generate and test check sums.
    bool fqzall;            // process data with fqzcomp only.
} fqz_params;

class FqzComp
{
public:
    FqzComp(fqz_params *p);
    ~FqzComp();
    void InitModel();
    void DestoryModel();
    int SaveModelToMem(char *modelbuf);
    void SaveModelToFile(std::fstream &s_out);
    void ReadModelFormMem(char *pbuf);
    void ReadModelFormFile(std::fstream &s_in);
    void DestoryModelbuf();
    int EncodeForModel();
    int Addbuf(char *id, int idlen, char *seq, int seqlen, char *qual, int quallen);

    int fqzall_encode(std::fstream &out);
    int align_encode(std::fstream &out);
    
    int fqzall_decode(std::fstream &in, char **namebuf, char **seqbuf, char **qualbuf, uint16_t **seqlen, int *ins,
                    int *mark);
    int align_decode(std::fstream &in, char **namebuf, char **seqbuf, char **qualbuf, char **bitbuf, uint8_t **orderbuf,
                    uint16_t **quallen, int *ins, int *mark, std::vector<char> **pvec);

    uint64_t getCompressTotalLen();
    uint32_t getInLen();
    uint32_t getReadCount();
    void isq_addmark(int mark);
    int isq_addbuf_match(char *id, int idlen, char *bit, int bitlen, char *qual, int quallen, int index,
                          char *degenerate);
    int isq_addbuf_unmatch(char *id, int idlen, char *seq, int seqlen, char *qual, int quallen, int index);
    
private:
    void encode_len(RangeCoder *rc, int len, SIMPLE_MODEL<256> &model_len1_t, 
                            SIMPLE_MODEL<256> &model_len2_t, SIMPLE_MODEL<2> &model_same_len_t);
    void encode_len_formodel(int len, SIMPLE_MODEL<256> &model_len1_t, 
                            SIMPLE_MODEL<256> &model_len2_t, SIMPLE_MODEL<2> &model_same_len_t);
    void compress_r0(bool bnewmodel, bool balign);
    void encode_name(RangeCoder *rc, char *name, int len,
                           SIMPLE_MODEL<256> *model_name_prefix_t,
                           SIMPLE_MODEL<256> *model_name_suffix_t,
                           SIMPLE_MODEL<256> *model_name_len_t,
                           SIMPLE_MODEL<128> *model_name_middle_t);
    void encode_name_formodel(char *name, int len,
                           SIMPLE_MODEL<256> *model_name_prefix_t,
                           SIMPLE_MODEL<256> *model_name_suffix_t,
                           SIMPLE_MODEL<256> *model_name_len_t,
                           SIMPLE_MODEL<128> *model_name_middle_t);
    void compress_r1(bool bnewmodel);
    void encode_seq8(RangeCoder *rc, char *seq, int len, BASE_MODEL<uint8_t> *model_seq8_t,
                          SIMPLE_MODEL<2> &seq_indicate_t, SIMPLE_MODEL<11> &seq_degenerate_t);
    void encode_seq8_formodel(char *seq, int len, BASE_MODEL<uint8_t> *model_seq8_t,
                          SIMPLE_MODEL<2> &seq_indicate_t, SIMPLE_MODEL<11> &seq_degenerate_t);
    void compress_r2(bool bnewmodel, bool balign);
    void encode_qual(RangeCoder *rc, char *qual, int len, SIMPLE_MODEL<QMAX> *model_qual_t);
    void encode_qual_formodel(char *qual, int len, SIMPLE_MODEL<QMAX> *model_qual_t);
    void compress_r3(bool bnewmodel);
    

    void decompress_r0(bool balign);
    int decode_name(RangeCoder *rc, char *name,
                           SIMPLE_MODEL<256> *model_name_prefix_t,
                           SIMPLE_MODEL<256> *model_name_suffix_t,
                           SIMPLE_MODEL<256> *model_name_len_t,
                           SIMPLE_MODEL<128> *model_name_middle_t);
    void decompress_r1(void);
    void decode_seq8(RangeCoder *rc, char *seq, int len, BASE_MODEL<uint8_t> * model_seq8_t,
                           SIMPLE_MODEL<2> &seq_indicate_t, SIMPLE_MODEL<11> &seq_degenerate_t);
    void decompress_r2(bool balign);
    void decode_qual(RangeCoder *rc, char *qual, int len, SIMPLE_MODEL<QMAX> *model_qual_t);
    void decompress_r3(void);
    int decode_len(RangeCoder *rc, SIMPLE_MODEL<256> &model_len1_t,
                    SIMPLE_MODEL<256> &model_len2_t,
                    SIMPLE_MODEL<2> &model_same_len_t);
    void isq_fqzall_decompress(char *in, int comp_len, int *out_len);
    void isq_align_decompress(char *in, int comp_len, int *out_len);

    void IntTo4Ch(int num, char *pbuf);
    void IntTo2Ch(int num, char *pbuf);
    void BitArryToBuf_s(const char *pbit, int len, char *pbuf, int *buflen);
private:
    /* --- Parameters passed into the constructor */
    int slevel, qlevel, nlevel;
    int both_strands;
    int extreme_seq;
    int multi_seq_model;
    int do_hash;
    int do_fqzall;

    int L[256];          // Sequence table lookups ACGTN->0..4

    // Input and output buffers; need to be size of BLK_SIZE
    char in_buf[BLK_SIZE];
    char out_buf[BLK_SIZE];

    int ns;
    int seq_len; // +ve if fixed per block, -ve if variable.
    char *name_buf;//[BLK_SIZE];
    char *seq_buf;//[BLK_SIZE/2];
    char *qual_buf;//[BLK_SIZE/2];
    char *name_p;
    char *seq_p;
    char *qual_p;
    uint16_t *name_len_a;//[BLK_SIZE/10];
    uint16_t *seq_len_a;//[BLK_SIZE/10];
    uint16_t *qual_len_a;//[BLK_SIZE/10];
    char *out0;//[BLK_SIZE/10]; // seq_len
    char *out1;//[BLK_SIZE]; // name
    char *out2;//[BLK_SIZE/2]; // seq
    char *out3;//[BLK_SIZE/2]; // qual
    int sz0, sz1, sz2, sz3;
    char *in_buf0, *in_buf1, *in_buf2, *in_buf3;
    uint32_t inLen, outLen;
    uint32_t pass_len;
    int uncomp_len;

    int m_seq_size;
    int m_qual_size;

/* --- Models */
    // Sequence length
    SIMPLE_MODEL<256> model_len1;
    SIMPLE_MODEL<256> model_len2;
    SIMPLE_MODEL<2> model_same_len;
    int last_len;

    bool isPEccordant;
    uint64_t m_totallen;

    char *bit_buf;
    char *bit_p;
    uint32_t bit_len;
    uint8_t *order_buf;
    int seq_count;
    char *out4, *out5, *out6;
    int sz4, sz5, sz6;
    char *in_buf4, *in_buf5, *in_buf6;
    std::vector<char> vec_degenerate;
    //SIMPLE_MODEL<5> seq_order; //block序号
    SIMPLE_MODEL<2> seq_indicate; //判断是否是简并碱基
    SIMPLE_MODEL<11> seq_degenerate;//存储比对失败的简并碱基
    //SIMPLE_MODEL<11> seq_degenerate_match;//存储比对成功的简并碱基

    std::bitset<8> m_bitset;

    // Names
    SIMPLE_MODEL<256> *model_name_prefix;
    SIMPLE_MODEL<256> *model_name_suffix;
    SIMPLE_MODEL<256> *model_name_len;
    SIMPLE_MODEL<128> *model_name_middle;

    char last_name[1024]; // Last name
    int last_name_len;    // Length of last name
    int last_p_len;       // Length of last common prefix
    int last_s_len;       // Length of last common suffix

    // Sequence
    int NS; // Number of bases of sequence context.
    BASE_MODEL<uint8_t> *model_seq8;
    BASE_MODEL<uint16_t> *model_seq16;

    // Quality
    SIMPLE_MODEL<QMAX> *model_qual;

    char *p_buf_same_len;
    char *p_buf_model_len1;
    char *p_buf_model_len2;

    char *p_buf_name_prefix;
    char *p_buf_name_suffix;
    char *p_buf_name_len;
    char *p_buf_name_middle;

    char *p_buf_seq_indicate;
    char *p_buf_seq_degenerate;
    char *p_buf_model_seq8;
    char *p_buf_model_qual;
};

#endif //SEQARC_FQZ_COMP_H