#ifndef SEQARC_FQZCOMP_H
#define SEQARC_FQZCOMP_H

#define __STDC_FORMAT_MACROS

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <inttypes.h>
#include <fstream>
#include <iostream>
#include <vector>

#include "sfh.h"
//#include "simple_model.h"

/* Keep as a power of 2 */
//#define QMAX 128
#define QMAX 64

/* Maximum length of a SOLiD sequencer entry */
#define MAX_SOLID 1024

#ifdef PTHREADS
#  include <pthread.h>
#  include <utmpx.h>
#else
#  define sched_getcpu() -1
#endif

/* Debug timing. Only works in non-threaded mode */
//#define TIMING

/*
 * SSE support to allow use of memory prefetching. It's only minor, but
 * it all helps.
 *
 * With    on -s6 -q3 -b (40million lines)
 *    encode: 2m55.691s+0m8.813s
 *    decode: 3m19.492s+0m2.720s
 *
 * Without on -s6 -q3 -b
 *    encode: 2m57.115s+0m3.260s
 *    decode: 3m46.338s+0m2.856s
 */
#ifdef __SSE__
#   include <xmmintrin.h>
#else
#   define _mm_prefetch(a,b)
#endif

#define ABS(a)   ((a)>0?(a):-(a))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

/* Range Coder:
 *
 * This is still using Eugene Shelwien's code from coders6c2.zip.
 * TODO: implement my own from scratch, although it's doubtful I'll
 * get something as efficient.
 */
//#include "../trunk/clrf.cdr"
#include "clr.cdr"
//#include "clrf256.cdr"
//#include "rc.h"

/*
 * Order 0 models, optimized for various sizes of alphabet.
 * order0_coder is the original Dmitry Shkarin code, tweaked a bit to
 * make RangeCoder a parameter and to work as a template for adjusting
 * the alphabet size.
 *
 * Simple_model is my own implementation adhering to the same API as
 * Dmitry Shkarin's original ORDER_0_CODER model. It's here for
 * purposes of adhering to the SequenceSqueeze entry rules, but
 * unfortunately it's ~5% slower.  It differs in having no escape
 * symbol (which is not a hindrance for stationary probabilities) and
 * using a different frequency updating method.
 *
 * Base coder is my own model specialising in symbols 0, 1, 2, 3, with
 * dedicated code for handling N (it encodes whichever of 0, 1, 2, 3
 * is most common and returns that value).
 *
 * All of these models are used in large arrays to implement a context,
 * ie to turn them into a high order model.
 */
#ifdef ORIG_MODEL
#    define ORDER_0_CODER SIMPLE_MODEL
#    include "order0_coder.h"
#else

#include "simple_model.h" // SIMPLE_MODEL

#endif

#include "base_model.h"       // BASE_MODEL
//#include "base_model2.h"       // BASE_MODEL

#define PREHASHED

/*
 * Defining this will use 2 sequence models, a smaller 7-mer one and a
 * larger one defined by the -s parameter. For novel kmers we encode using
 * the smaller 7-mer (and then update the full kmer version).
 *
 * This improves compression early on in a stream while the larger model
 * learns, but has no impact after a while as the larger model fills up.
 *
 * It's about 15% slower overall, but saves ~1%, depending on size of input.
 * Only a worthwhile trade when going for maximum compression.
 *
 * NOTE: now specified by -s<num>+; eg -s6+ vs -s6
 */


#define BLK_SIZE 10*1024*1024
//#define BLK_SIZE 1000000

/* QBITS is the size of the quality context, 2x quals */
#define QBITS 12
#define QSIZE (1<<QBITS)

/*
 * fqz parameter block.
 */
typedef struct {
    int slevel;         // -s level
    int qlevel;         // -q level
    int nlevel;         // -n level
    bool both_strands;      // True if -b used
    bool extreme_seq;       // True if -e used; 16-bit seq counters
    bool multi_seq_model;   // True if -s<level>+; use 2 model sizes
    int qual_approx;        // 0 for lossless, >0 for lossy qual encoding
    bool do_threads;        // Simple multi-threading enabled.
    bool do_hash;       // Generate and test check sums.
    bool fqzall;
} fqz_params;

/*
 * The fqz class itself
 */
struct fqz {
public:
    fqz();

    /* Replace with an argument struct */
    fqz(fqz_params *p);
    ~fqz();

    int encode(std::fstream &in, std::fstream &out);
    int decode(std::fstream &in, std::fstream &out);

    int iq_encode(std::string &id, std::string &qual, std::fstream &out);
    int iq_decode(std::fstream &in, std::vector<std::string> &out1, std::vector<std::string> &out2);

    int isq_encode(std::string &id, std::string &seq, std::string &qual, std::fstream &out);
    //int isq_encode(char *id, int idlen, char *seq, int seqlen, char *qual, int quallen, std::fstream &out);

    int isq_addbuf(char *id, int idlen, char *seq, int seqlen, char *qual, int quallen);
    void isq_addmark(int mark);
    int isq_doencode(std::fstream &out);

int isq_addbuf_match(char *id, int idlen, char *seq, int seqlen, char *qual, int quallen,int index);
int isq_addbuf_unmatch(char *id, int idlen, char *seq, int seqlen, char *qual, int quallen,int index);
int isq_doencode_s(std::fstream &out);
void isq_decompress_s(char *in, int comp_len, int *out_len);
int isq_decode_s(std::fstream &in, char **namebuf, char **seqbuf, char **qualbuf, char **bitbuf, uint8_t **orderbuf, uint16_t **quallen, int *ins, int *mark);

    int isq_decode(std::fstream &in, char **namebuf, char **seqbuf, char **qualbuf, uint16_t **seqlen, int *ins, int *mark);
    int isq_decode(std::fstream &in, std::vector<std::string> &out1, std::vector<std::string> &out2, std::vector<std::string> &out3);

    /* Compression metrics */
    uint64_t base_in, base_out;
    uint64_t qual_in, qual_out;
    uint64_t name_in, name_out;

    /* Entry points for pthread calls; consider as internal */
    void compress_r0();
    void compress_r1();
    void compress_r2();
    void compress_r2_s();
    void compress_r3();
    void compress_r3_s();

    void decompress_r1();
    void decompress_r2();
    void decompress_r2_s();
    void decompress_r3();
    void decompress_r3_s();

    uint64_t getCompressTotalLen();
    uint32_t getInLen();
    void test();
protected:
    /* --- Parameters passed into the constructor */
    int slevel, qlevel, nlevel;
    int both_strands;
    int extreme_seq;
    int multi_seq_model;
    int qual_approx;
    int do_threads;
    int do_hash;
    int do_fqzall;

    int L[256];          // Sequence table lookups ACGTN->0..4

    /* --- Buffers */
    // Input and output buffers; need to be size of BLK_SIZE
    char in_buf[BLK_SIZE];
    char out_buf[BLK_SIZE];
    int out_ind, in_ind; // indices into in_buf & out_buf.

    /* Block variables for the encoder to work on */
    /* FIXME: shrink these, or new[] them */
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
    uint32_t inLen, outLen, readBufMark;
    uint32_t pass_len;
    int uncomp_len;
    
    // Hashes for error detection.
    unsigned char *chk_in;
    uint32_t chk_len;
    uint32_t chksum;

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
    char *out4;
    int sz4;
    char *out5;
    int sz5;
    char *in_buf4,*in_buf5;
    SIMPLE_MODEL<4> zk_test;

    void encode_len(RangeCoder *rc, int len);
    int  decode_len(RangeCoder *rc);


    // Names
    SIMPLE_MODEL<256> model_name_prefix[256];
    SIMPLE_MODEL<256> model_name_suffix[256];
    SIMPLE_MODEL<256> model_name_len[256];
    SIMPLE_MODEL<128> model_name_middle[8192];

#define MAX_TOK 1000
    // Name lvl 2
    SIMPLE_MODEL<10>  model_name_type[MAX_TOK];
    SIMPLE_MODEL<256> model_name_alpha_len[MAX_TOK];
    SIMPLE_MODEL<256> model_name_alpha[MAX_TOK];
    SIMPLE_MODEL<256> model_name_zero[MAX_TOK];
    SIMPLE_MODEL<256> model_name_digit0[MAX_TOK];
    SIMPLE_MODEL<256> model_name_digit1[MAX_TOK];
    SIMPLE_MODEL<256> model_name_digit2[MAX_TOK];
    SIMPLE_MODEL<256> model_name_digit3[MAX_TOK];
    SIMPLE_MODEL<256> model_name_ddelta[MAX_TOK];
    SIMPLE_MODEL<256> model_name_char[MAX_TOK];

    char last_name[1024]; // Last name
    int last_name_len;    // Length of last name
    int last_p_len;       // Length of last common prefix
    int last_s_len;       // Length of last common suffix

    void encode_name(RangeCoder *rc, char *name, int len);
    int decode_name(RangeCoder *rc, char *name);

    void encode_name2(RangeCoder *rc, char *name, int len);
    int decode_name2(RangeCoder *rc, char *name);


    // Sequence
    int NS; // Number of bases of sequence context.
    BASE_MODEL<uint8_t> *model_seq8;
    BASE_MODEL<uint16_t> *model_seq16;

    //#define NS_MASK ((1<<(2*NS))-1)
#define SMALL_NS 7
#define SMALL_MASK ((1<<(2*SMALL_NS))-1)
    BASE_MODEL<uint8_t> model_seq_small[1 << (2 * SMALL_NS)];

    void encode_seq8 (RangeCoder *rc, char *seq, int len);
    void encode_seq16(RangeCoder *rc, char *seq, int len);

    void decode_seq8 (RangeCoder *rc, char *seq, int len);
    void decode_seq16(RangeCoder *rc, char *seq, int len);


    // Quality
    SIMPLE_MODEL<QMAX> *model_qual;
#define SMALL_QMASK (QSIZE-1)

    void encode_qual(RangeCoder *rc, char *qual, int len);
    void decode_qual(RangeCoder *rc, char *qual, int len);

    /* --- Main functions for compressing and decompressing blocks */
    int fq_compress(char *in,  int in_len,
                    char *out, int *out_len,
                    char **in_end, int *nseqs);

    char *fq_decompress(char *in, int comp_len, int *uncomp_len);

    char *iq_decompress(char *in, int comp_len, int *uncomp_len);
    void isq_decompress(char *in, int comp_len, int *uncomp_len);

    void load_hash_freqs(const char *fn);
};


int xget(std::fstream &in, unsigned char *in_buffer, int count);
int xwrite(std::fstream &out, unsigned char *out_buffer, int count);

#endif //SEQARC_FQZCOMP_H
