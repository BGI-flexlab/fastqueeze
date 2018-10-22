#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>

#include "fqzcomp.h"

/* -------------------------------------------------------------------------
 * Constructors and destructors   //构造函数和析构函数
 */
fqz::fqz() {
    fqz(NULL);
}

fqz::fqz(fqz_params *p) {
    if (p) {
    slevel          = p->slevel;
    qlevel          = p->qlevel;
    nlevel          = p->nlevel;
    both_strands    = p->both_strands;
    extreme_seq     = p->extreme_seq;
    multi_seq_model = p->multi_seq_model;
    qual_approx     = p->qual_approx;
    do_threads      = p->do_threads;
    do_hash         = p->do_hash; // negligible slow down
    } else {
    slevel = 3;
    qlevel = 2;
    nlevel = 1;
    both_strands = 0;
    extreme_seq = 0;
    multi_seq_model = 0;
    qual_approx = 0;
    do_threads = 1;
    do_hash = 1;
    }

    /* ACGTN* */
    for (int i = 0; i < 256; i++)
        L[i] = 0;
    L['A'] = L['a'] = 0;
        L['C'] = L['c'] = 1;
        L['G'] = L['g'] = 2;
        L['T'] = L['t'] = 3;
    
    NS = 7 + slevel;
    if (extreme_seq) {
        model_seq8  = NULL;
        model_seq16 = new BASE_MODEL<uint16_t>[1<<(2*NS)];
    } else {
        model_seq8  = new BASE_MODEL<uint8_t>[1<<(2*NS)];
        model_seq16 = NULL;
    }

    int qsize = QSIZE;
    if (qlevel > 1) qsize *= 16;
    if (qlevel > 2) qsize *= 16;
    model_qual = new SIMPLE_MODEL<QMAX>[qsize];

#ifdef PREHASHED
    /* Helps on shallow data far more than deep data */
    //load_hash_freqs("human_hash13_2strand");
    //load_hash_freqs("human_hash16");
#endif

    /* Name settings */
    memset(last_name, ' ', 1024);
    last_name_len = 0;
    last_p_len = 0;
    last_s_len = 0;

    /* Length settings */
    last_len = 0;

    name_in = name_out = 0;
    base_in = base_out = 0;
    qual_in = qual_out = 0;
}

fqz::~fqz() {
    if (model_seq8)  delete[] model_seq8;
    if (model_seq16) delete[] model_seq16;
    if (model_qual)  delete[] model_qual;
}

/*
 * Use a predefined hash table of sequence frequencies for priming the
 * model_seq contexts. If this turns out to be widely useful then
 * model_seq should be updated to allow this to work in one easy data read.
 *
 * Arguably we should also compute the stats once and then not update
 * during encoding.
 */
void fqz::load_hash_freqs(const char *fn) {
    unsigned char c4[4];
    int ctx = 0;
    FILE *fp = fopen(fn, "rb");

    fprintf(stderr, "Loading %s...", fn);
    fflush(stderr);

    if (!fp) {
        perror(fn);
        return;
    }

    while (1 == fread(c4, 4, 1, fp)) {
        int st[4];
        st[0] = c4[0] + 1;
        st[1] = c4[1] + 1;
        st[2] = c4[2] + 1;
        st[3] = c4[3] + 1;

        assert(ctx < (1<<(2*NS)));
        if (extreme_seq) {
            model_seq16[ctx++].reset(st);
        } else {
            model_seq8[ctx++].reset(st);
        }
    }

    fclose(fp);

    fprintf(stderr, "done\n");
}


/* -------------------------------------------------------------------------
 * Name model
 */
void fqz::encode_name(RangeCoder *rc, char *name, int len) {
    int p_len, s_len; // prefix and suffix length
    int i, j, k, last_char;
//    static char meta[] = {
//  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, /* 00-0f */
//  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, /* 10-1f */
//  1, 1, 1, 1, 1, 1, 1, 1,   1, 1, 1, 1, 1, 1, 1, 1, /* 20-2f */
//  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 1, 1, 1, 1, 1, 1, /* 30-3f */
//  1, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, /* 40-4f */
//  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 1, 1, /* 50-5f */
//  1, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, /* 60-6f */
//  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 1, 1, /* 70-7f */
//    };

    _mm_prefetch((const char *)&model_name_prefix[last_p_len], _MM_HINT_T0);
    _mm_prefetch((const char *)&model_name_suffix[last_s_len], _MM_HINT_T0);
    _mm_prefetch((const char *)&model_name_len[last_name_len], _MM_HINT_T0);

    // Prefix
    for (i = 0; i < len && i < last_name_len; i++) {
    if (name[i] != last_name[i])
        break;
    }
    p_len = i;

    // Suffix
    for (i = len-1, j= last_name_len-1; i >= 0 && j >= 0; i--, j--) {
    if (name[i] != last_name[j])
        break;
    }
    s_len = len-1 - i;
    if (len - s_len - p_len < 0)
    s_len = len - p_len;

    model_name_prefix[last_p_len].   encodeSymbol(rc, p_len);
    model_name_suffix[last_s_len].   encodeSymbol(rc, s_len);
    model_name_len   [last_name_len].encodeSymbol(rc, len);

    last_p_len = p_len;
    last_s_len = s_len;

    int len2 = len - s_len, lc2 = p_len ? 1 : 0;
    for (i = j = p_len, k = 0; i < len2; i++, j++, k++) {
    last_char = ((last_name[j]-32)*2 + lc2 + k*64) % 8192;
    _mm_prefetch((const char *)&model_name_middle[last_char], _MM_HINT_T0);

    //last_char = ((last_name[j]-32)*2 + lc2 + (i+j)*32) % 8192;

    model_name_middle[last_char].encodeSymbol(rc, name[i] & 0x7f);

    //if (meta[name[i]]      && name[i] != last_name[j]) j++;
    //if (meta[last_name[j]] && name[i] != last_name[j]) j--;

    //if (meta[name[i]]) k = (k+3)>>2<<2;

    if (name[i] == ' ' && last_name[j] != ' ') j++;
        if (name[i] != ' ' && last_name[j] == ' ') j--;
    if (name[i] == ':' && last_name[j] != ':') j++;
    if (name[i] != ':' && last_name[j] == ':') j--;

    if (name[i] == ':' || name[i] == ' ') k = (k+3)>>2<<2;

    lc2 = name[i] == last_name[j];
    }

    memcpy(last_name, name, len);
    last_name_len = len;
}

int fqz::decode_name(RangeCoder *rc, char *name) {
    int p_len, s_len, len; // prefix and suffix length
    int i, j, k;
    int last_char;

    _mm_prefetch((const char *)&model_name_prefix[last_p_len], _MM_HINT_T0);
    _mm_prefetch((const char *)&model_name_prefix[last_p_len], _MM_HINT_T0);
    _mm_prefetch((const char *)&model_name_prefix[last_p_len], _MM_HINT_T0);

    p_len = model_name_prefix[last_p_len].   decodeSymbol(rc);
    s_len = model_name_suffix[last_s_len].   decodeSymbol(rc);
    len   = model_name_len   [last_name_len].decodeSymbol(rc);

    last_p_len = p_len;
    last_s_len = s_len;

    for (i = 0; i < p_len; i++)
    name[i] = last_name[i];

    //fprintf(stderr, "%d: p_len = %d, s_len = %d, len = %d last='%.*s'\n",
    //column, p_len, s_len, len, last_name_len, last);

    int len2 = len - s_len, lc2 = p_len ? 1 : 0;
    for (i = j = p_len, k = 0; i < len2; i++, j++, k++) {
    unsigned char c;

    last_char = ((last_name[j]-32)*2 + lc2 + k*64) % 8192;
    //last_char = (last_name[j] + j*64) % 8192;

    c = model_name_middle[last_char].decodeSymbol(rc);
    //c = 'x';
    name[i] = c;

    if (c == ' ' && last_name[j] != ' ') j++;
        if (c != ' ' && last_name[j] == ' ') j--;
    if (c == ':' && last_name[j] != ':') j++;
    if (c != ':' && last_name[j] == ':') j--;

    if (name[i] == ':' || name[i] == ' ') k = (k+3)>>2<<2;

    lc2 = c == last_name[j];
    }

    for (j = last_name_len-s_len; i < len; i++, j++)
    name[i] = last_name[j];

    memcpy(last_name, name, len);
    last_name_len = len;

    return len;
}

/* Level 2 */
enum name_type {N_UNK = 0, N_ALPHA, N_CHAR,
        N_ZERO, N_DIGITS, N_DDELTA, N_MATCH, N_END};

void fqz::encode_name2(RangeCoder *rc, char *name, int len) {
    int i, j, k;

    static int last_token_type[1024];
    static int last_token_int[1024];
    static int last_token_str[1024];

    //fprintf(stderr, "NAME: %.*s\n", len, name);

    int ntok = 0;
    for (i = j = 0, k = 0; i < len; i++, j++, k++) {
        /* Determine data type of this segment */
        int n_type = N_ALPHA;
        if (isalpha(name[i])) {
            int s = i+1;
            while (s < len && isalpha(name[s]))
            s++;
            n_type = N_ALPHA;

            if (last_token_type[ntok] == N_ALPHA) {
                if (s-i == last_token_int[ntok] &&
                    memcmp(&name[i],
                       &last_name[last_token_str[ntok]],
                       s-i) == 0) {
                    //fprintf(stderr, "Tok %d (mat)\n", N_MATCH);
                    model_name_type[ntok].encodeSymbol(rc, N_MATCH);
                } else {
                    //fprintf(stderr, "Tok %d (alpha)\n", N_ALPHA);
                    model_name_type[ntok].encodeSymbol(rc, N_ALPHA);
                    model_name_alpha_len[ntok].encodeSymbol(rc, s-i);
                    for (int x = 0; x < s-i; x++) {
                    model_name_alpha[ntok].encodeSymbol(rc, name[i+x]);
                    }
                }
            } else {
            //fprintf(stderr, "Tok %d (alpha)\n", N_ALPHA);
                model_name_type[ntok].encodeSymbol(rc, N_ALPHA);
                model_name_alpha_len[ntok].encodeSymbol(rc, s-i);
                for (int x = 0; x < s-i; x++) {
                    model_name_alpha[ntok].encodeSymbol(rc, name[i+x]);
                }
            }

            last_token_int[ntok] = s-i;
            last_token_str[ntok] = i;
            last_token_type[ntok] = N_ALPHA;

            i = s-1;
        } else if (name[i] == '0') {
            int s = i, v;
            while (s < len && name[s] == '0')
            s++;
            v = s-i;
            n_type = N_ZERO;

            if (last_token_type[ntok] == N_ZERO) {
                if (last_token_int[ntok] == v) {
                    //fprintf(stderr, "Tok %d (mat)\n", N_MATCH);
                    model_name_type[ntok].encodeSymbol(rc, N_MATCH);
                } else {
                    //fprintf(stderr, "Tok %d (0)\n", N_ZERO);
                    model_name_type[ntok].encodeSymbol(rc, N_ZERO);
                    model_name_zero[ntok].encodeSymbol(rc, v);
                }
            } else {
            //fprintf(stderr, "Tok %d (0)\n", N_ZERO);
            model_name_type[ntok].encodeSymbol(rc, N_ZERO);
            model_name_zero[ntok].encodeSymbol(rc, v);
            }

            last_token_int[ntok] = v;
            last_token_type[ntok] = N_ZERO;

            i = s-1;
        } else if (isdigit(name[i])) {
            int s = i;
            int v = 0;
            int d = 0;
            while (s < len && isdigit(name[s]) && v < (1<<27)) {
                v = v*10 + name[s] - '0';
                s++;
            }
            n_type = N_DIGITS;

            if (last_token_type[ntok] == N_DIGITS) {
                if ((d = v - last_token_int[ntok]) == 0) {
                    //fprintf(stderr, "Tok %d (mat)\n", N_MATCH);
                    model_name_type[ntok].encodeSymbol(rc, N_MATCH);
                } else if (d < 256 && d > 0) {
                    //fprintf(stderr, "Tok %d (delta)\n", N_DDELTA);
                    model_name_type[ntok].encodeSymbol(rc, N_DDELTA);
                    model_name_ddelta[ntok].encodeSymbol(rc, d);
                } else {
                    //fprintf(stderr, "Tok %d (dig)\n", N_DIGITS);
                    model_name_type[ntok].encodeSymbol(rc, N_DIGITS);
                    model_name_digit0[ntok].encodeSymbol(rc, (v>> 0) & 0xff);
                    model_name_digit1[ntok].encodeSymbol(rc, (v>> 8) & 0xff);
                    model_name_digit2[ntok].encodeSymbol(rc, (v>>16) & 0xff);
                    model_name_digit3[ntok].encodeSymbol(rc, (v>>24) & 0xff);
                }
            } else {
            //fprintf(stderr, "Tok %d (dig)\n", N_DIGITS);
                model_name_type[ntok].encodeSymbol(rc, N_DIGITS);
                model_name_digit0[ntok].encodeSymbol(rc, (v>> 0) & 0xff);
                model_name_digit1[ntok].encodeSymbol(rc, (v>> 8) & 0xff);
                model_name_digit2[ntok].encodeSymbol(rc, (v>>16) & 0xff);
                model_name_digit3[ntok].encodeSymbol(rc, (v>>24) & 0xff);
            }

            last_token_int[ntok] = v;
            last_token_type[ntok] = N_DIGITS;

            i = s-1;
        } else {
            n_type = N_CHAR;

            if (last_token_type[ntok] == N_CHAR) {
                if (name[i] == last_token_int[ntok]) {
                    //fprintf(stderr, "Tok %d (mat)\n", N_MATCH);
                    model_name_type[ntok].encodeSymbol(rc, N_MATCH);
                } else {
                    //fprintf(stderr, "Tok %d (chr)\n", N_CHAR);
                    model_name_type[ntok].encodeSymbol(rc, N_CHAR);
                    model_name_char[ntok].encodeSymbol(rc, name[i]);
                }
            } else {
            //fprintf(stderr, "Tok %d (chr)\n", N_CHAR);
                model_name_type[ntok].encodeSymbol(rc, N_CHAR);
                model_name_char[ntok].encodeSymbol(rc, name[i]);
            }

            last_token_int[ntok] = name[i];
            last_token_type[ntok] = N_CHAR;
        }

        ntok++;
    }
    //fprintf(stderr, "Tok %d (end)\n", N_END);
    model_name_type[ntok].encodeSymbol(rc, N_END);
    
    memcpy(last_name, name, len);
    last_name_len = len;
}

int fqz::decode_name2(RangeCoder *rc, char *name) {
    enum name_type tok;
    int ntok = 0, i = 0, v;

    static int last_token_type[1024];
    static int last_token_int[1024];
    static int last_token_str[1024];

    for (;;) {
    tok = (enum name_type)model_name_type[ntok].decodeSymbol(rc);
    //fprintf(stderr, "tok=%d, last type=%d int=%d str=%d\n",
    //  tok, last_token_type[ntok], last_token_int[ntok],
    //  last_token_str[ntok]);
    if (tok == N_END)
        break;

    switch (tok) {
        /* Str delta too? */
    case N_ALPHA:
        v = model_name_alpha_len[ntok].decodeSymbol(rc);
        last_token_int[ntok] = v; // len
        last_token_str[ntok] = i;
        for (int x = 0; x < v; x++)
        // also per 'x'; per pos in name? */
            name[i++] = model_name_alpha[ntok].decodeSymbol(rc);
        last_token_type[ntok] = N_ALPHA;
        break;

    case N_CHAR:
        v = model_name_char[ntok].decodeSymbol(rc);
        name[i++] = v;
        last_token_int[ntok] = v;
        last_token_type[ntok] = N_CHAR;
        break;

    case N_ZERO:
        v = model_name_zero[ntok].decodeSymbol(rc);
        last_token_int[ntok] = v;
        for (int x = 0; x < v; x++)
            name[i++] = '0';
        last_token_type[ntok] = N_ZERO;
        break;

    case N_DIGITS: {
        char rev[100];
        int ri = 0, v0, v1, v2, v3;

        v0 = model_name_digit0[ntok].decodeSymbol(rc);
        v1 = model_name_digit1[ntok].decodeSymbol(rc);
        v2 = model_name_digit2[ntok].decodeSymbol(rc);
        v3 = model_name_digit3[ntok].decodeSymbol(rc);
        v = v0 + (v1<<8) + (v2<<16) + (v3<<24);
        last_token_int[ntok] = v;
        while (v > 0) {
            rev[ri++] = '0' + (v%10);
            v /= 10;
        }
        while (ri > 0)
            name[i++] = rev[--ri];
        last_token_type[ntok] = N_DIGITS;
        break;
    }

    case N_DDELTA: {
        char rev[100];
        int ri = 0;

        v = model_name_ddelta[ntok].decodeSymbol(rc);
        v += last_token_int[ntok];
        last_token_int[ntok] = v;
        while (v > 0) {
            rev[ri++] = '0' + (v%10);
            v /= 10;
        }
        while (ri > 0)
            name[i++] = rev[--ri];
        last_token_type[ntok] = N_DIGITS;
        break;
    }

    case N_MATCH:
        switch (last_token_type[ntok]) {
            case N_CHAR:
                name[i++] = last_token_int[ntok];
                break;

            case N_ALPHA:
                v = last_token_int[ntok];
                for (int x = 0; x < v; x++)
                    name[i++] = last_name[last_token_str[ntok]+x];
                last_token_str[ntok] = i-v;
                break;

            case N_ZERO:
                v = last_token_int[ntok];
                for (int x = 0; x < v; x++)
                    name[i++] = '0';
                break;

            case N_DIGITS: {
                char rev[100];
                int ri = 0;
                v = last_token_int[ntok];

                while (v > 0) {
                    rev[ri++] = '0' + (v%10);
                    v /= 10;
                }
                while (ri > 0)
                    name[i++] = rev[--ri];
                break;
            }
        }
        break;

    default:
        fprintf(stderr, "Unexpected name token %d\n", tok);
        return -1;
    }

    ntok++;
    }

    name[i] = '\0';
    memcpy(last_name, name, i);
    last_name_len = i;

    //fprintf(stderr, "%s\n", name);

    return i;
}

/* -------------------------------------------------------------------------
 * Sequence length model
 * last_len: alphbet's number of each line
 */
void fqz::encode_len(RangeCoder *rc, int len) {
    if (len != last_len) {
        model_same_len.encodeSymbol(rc, 0);
        model_len1.encodeSymbol(rc, len & 0xff);
        model_len2.encodeSymbol(rc, (len >> 8) & 0xff);
    } else {
        model_same_len.encodeSymbol(rc, 1);
    }
}

//decode_len return the size of each line.
int fqz::decode_len(RangeCoder *rc) {
    if (model_same_len.decodeSymbol(rc)) {
        return last_len;
    } else {
        int l1 = model_len1.decodeSymbol(rc);
        int l2 = model_len2.decodeSymbol(rc);
        last_len = l1 + (l2 << 8);
        return last_len;
    }
}

/* -------------------------------------------------------------------------
 * Sequence model
 *
 * We have 8-bit and 16-bit accumulator versions.
 * The 8-bit one is lower memory and somtimes slightly faster, but
 * marginally less optimal in compression ratios (within 1%).
 */
void fqz::encode_seq8(RangeCoder *rc, char * seq, int len) {
    int last, last2;
    int bc[4] = {(3 - 0) << (2 * NS - 2),
                 (3 - 1) << (2 * NS - 2),
                 (3 - 2) << (2 * NS - 2),
                 (3 - 3) << (2 * NS - 2)};
    const int NS_MASK = ((1 << (2 * NS)) - 1);

    /* Corresponds to a 12-mer word that doesn't occur in human genome. */
    last  = 0x007616c7 & NS_MASK;
    last2 = (0x2c6b62ff >> (32 - 2 * NS)) & NS_MASK;
    
    _mm_prefetch((const char *)&model_seq8[last], _MM_HINT_T0);

    if (multi_seq_model) {
        for (int i = 0; i < len; i++) {
            unsigned int l2 = (last << 2) & NS_MASK;
            _mm_prefetch((const char *)&model_seq8[l2+0], _MM_HINT_T0);
            //_mm_prefetch((const char *)&model_seq8[l2+3], _MM_HINT_T0);

            unsigned char  b = L[(unsigned char)seq[i]];

            /* Works OK on small files to rapidly train, but no gain on
             * large ones
             */
            if ((model_seq8 [last].getTopSym() *
                 model_seq_small[last & SMALL_MASK].getSummFreq()) >
                (model_seq_small[last & SMALL_MASK].getTopSym() *
                 model_seq8 [last].getSummFreq())) {
                 model_seq8 [last].encodeSymbol(rc, b);
            //model_seq_small[last & SMALL_MASK].updateSymbol(b);
            } else {
                model_seq_small[last & SMALL_MASK].encodeSymbol(rc, b);
                model_seq8  [last].updateSymbol(b);
            }

            last = (last*4 + b) & NS_MASK;

            /*
             * On deep data hashing both strands works well. On shallow data
             * it adds a significant CPU hit for very minimal gains (at best
             * 0.5% smaller seq portion).
             * Eg: -6=>513049628, -6b=>510382143, -8b=501978520 (Seq only)
             *
             * Pre-seeding the hash table with double-stranded human genome
             * hashes seems like a faster starting point and will help more
             * for shallow data too. However even this isn't as significant
             * as it sounds.
             * -5=>516624591, -5h=>514002730
             */
            if (both_strands) {
                int b2 = last2 & 3;
                last2 = last2/4 + ((3-b) << (2*NS-2));
                _mm_prefetch((const char *)&model_seq8[last2], _MM_HINT_T0);
                model_seq8[last2].updateSymbol(b2);
            }
        }
    } else {
        if (both_strands) {
            unsigned l2 = last2;
            for (int i = 0; i < len && i < 128; i++) {
                unsigned char  b = L[(unsigned char)seq[i]];
                l2 = l2/4 + bc[b];
                _mm_prefetch((const char *)&model_seq8[l2], _MM_HINT_T0);
            }

            for (int i = 0; i < len; i++) {
                unsigned int l2 = (last << 2) & NS_MASK;
                _mm_prefetch((const char *)&model_seq8[l2+0], _MM_HINT_T0);

                unsigned char  b = L[(unsigned char)seq[i]];
                model_seq8[last].encodeSymbol(rc, b);

                last = (last*4 + b) & NS_MASK;

                int b2 = last2 & 3;
                last2 = last2/4 + bc[b];
                model_seq8[last2].updateSymbol(b2);
            }
        } else {
            for (int i = 0; i < len; i++) {
                unsigned int l2 = (last << 2) & NS_MASK;
                _mm_prefetch((const char *)&model_seq8[l2+0], _MM_HINT_T0);

                unsigned char  b = L[(unsigned char)seq[i]];
                model_seq8[last].encodeSymbol(rc, b);

                last = ((last<<2) + b) & NS_MASK;
            }
        }
    }
}

void fqz::encode_seq16(RangeCoder *rc, char *seq, int len) {
    int last, last2;
    const int NS_MASK = ((1 << (2 * NS)) - 1);

    /* Corresponds to a 12-mer word that doesn't occur in human genome. */
    last  = 0x7616c7 & NS_MASK;
    last2 = (0x2c6b62ff >> (32 - 2 * NS)) & NS_MASK;

    for (int i = 0; i < len; i++) {
        unsigned char  b = L[(unsigned char)seq[i]];
        model_seq16[last].encodeSymbol(rc, b);
        last = (last * 4 + b) & NS_MASK;
        //volatile int p = *(int *)&model_seq16[last];
        _mm_prefetch((const char *)&model_seq16[last], _MM_HINT_T0);

        if (both_strands) {
            int b2 = last2 & 3;
            last2 = last2/4 + ((3-b) << (2 * NS - 2));
            _mm_prefetch((const char *)&model_seq16[last2], _MM_HINT_T0);
            model_seq16[last2].updateSymbol(b2);
        }
    }
}

void fqz::decode_seq8(RangeCoder *rc, char *seq, int len) {
    int last, last2;
    const char *dec = "ACGTN";
    const int NS_MASK = ((1 << (2 * NS)) - 1);

    /*
     * We can't do the same prefetch loop here as we don't know what the
     * data is yet, so we can't predict the memory addresses in model[]
     * that we're going to access.
     *
     * However we can guess it'll be one of 4 base calls, so populate the
     * cache with all 4 choices so by the time we get there it'll have been
     * loaded.
     */

    last  = 0x7616c7 & NS_MASK;
    last2 = (0x2c6b62ff >> (32 - 2 * NS)) & NS_MASK;

    if (multi_seq_model) {
        for (int i = 0; i < len; i++) {
            unsigned char b;
            unsigned int m = (last<<2) & NS_MASK;
            _mm_prefetch((const char *)&model_seq8[m+0], _MM_HINT_T0);
            //_mm_prefetch((const char *)&model_seq8[m+3], _MM_HINT_T0);

            m &= SMALL_MASK;
            _mm_prefetch((const char *)&model_seq_small[m+0], _MM_HINT_T0);
            //_mm_prefetch((const char *)&model_seq_small[m+3], _MM_HINT_T0);

            if ((model_seq8     [last             ].getTopSym() *
                 model_seq_small[last & SMALL_MASK].getSummFreq()) >
                (model_seq_small[last & SMALL_MASK].getTopSym() *
                 model_seq8     [last             ].getSummFreq())) {
                b = model_seq8[last].decodeSymbol(rc);
            //model_seq_small[last & SMALL_MASK].updateSymbol(b);
            } else {
                b = model_seq_small[last & SMALL_MASK].decodeSymbol(rc);
                model_seq8[last].updateSymbol(b);
            }
            *seq++ = dec[b];
            last = (last*4 + b) & NS_MASK;

            _mm_prefetch((const char *)&model_seq8[last], _MM_HINT_T0);

            if (both_strands) {
                int b2 = last2 & 3;
                last2 = last2/4 + ((3-b) << (2 * NS-2));
                _mm_prefetch((const char *)&model_seq8[last2], _MM_HINT_T0);
                model_seq8[last2].updateSymbol(b2);
            }
        }
    } else {
        if (both_strands) {
            for (int i = 0; i < len; i++) {
                unsigned char b;
                unsigned int m = (last<<2) & NS_MASK;
                int b2;

                /* Get next set loaded */
                _mm_prefetch((const char *)&model_seq8[m+0], _MM_HINT_T0);
                //_mm_prefetch((const char *)&model_seq8[m+3], _MM_HINT_T0);

                b = model_seq8[last].decodeSymbol(rc);

                *seq++ = dec[b];
                last = (last*4 + b) & NS_MASK;

                b2 = last2 & 3;
                last2 = last2/4 + ((3-b) << (2*NS-2));
                _mm_prefetch((const char *)&model_seq8[last2], _MM_HINT_T0);
                model_seq8[last2].updateSymbol(b2);
            }
        } else {
            for (int i = 0; i < len; i++) {
                unsigned char b;
                unsigned int m = (last<<2) & NS_MASK;

                /* Get next set loaded */
                _mm_prefetch((const char *)&model_seq8[m+0], _MM_HINT_T0);
                //_mm_prefetch((const char *)&model_seq8[m+3], _MM_HINT_T0);

                b = model_seq8[last].decodeSymbol(rc);

                *seq++ = dec[b];
                last = (last*4 + b) & NS_MASK;
            }
        }
    }
}

void fqz::decode_seq16(RangeCoder *rc, char *seq, int len) {
    int last, last2;
    const char *dec = "ACGTN";
    const int NS_MASK = ((1 << (2 * NS)) - 1);

    last  = 0x7616c7 & NS_MASK;
    last2 = (0x2c6b62ff >> (32 - 2*NS)) & NS_MASK;
    for (int i = 0; i < len; i++) {
        unsigned char b = model_seq16[last].decodeSymbol(rc);
        *seq++ = dec[b];
        last = (last*4 + b) & NS_MASK;
        _mm_prefetch((const char *)&model_seq16[last], _MM_HINT_T0);

        if (both_strands) {
            int b2 = last2 & 3;
            last2 = last2/4 + ((3-b) << (2*NS-2));
            _mm_prefetch((const char *)&model_seq16[last2], _MM_HINT_T0);
            model_seq16[last2].updateSymbol(b2);
        }
    }
}


/* -------------------------------------------------------------------------
 * Quality model
 */

#if 0
void fqz::encode_qual(RangeCoder *rc, char *seq, char *qual, int len) {
    unsigned int last = 0;
    int delta = 5;
    int i, len2 = len;
    int q1 = 0, q2 = 0;
    unsigned int X[1024];

    /* Removing "Killer Bees" */
    while (len2 > 0 && qual[len2-1] == '#')
    len2--;

    /* Prefetching & context caching. Only minor speed improvement. */
    for (i = 0; i < len2; i++) {
    unsigned char q = (qual[i] - '!') & (QMAX-1);

    X[i] = last;
    _mm_prefetch((const char *)&model_qual[last], _MM_HINT_T0);
    _mm_prefetch(64+(const char *)&model_qual[last], _MM_HINT_T0);

    // previous 2-3 bytes
    if (QBITS == 12) {
        last = ((MAX(q1, q2)<<6) + q) & ((1<<QBITS)-1);
    } else {
        last = ((last << 6) + q) & ((1<<QBITS)-1);
    }

    if (qlevel > 1) {
        last  += (q1==q2) << QBITS;
        // delta saves 3-4%, but adds 14% cpu
        delta += (q1>q)*(q1-q);
        last  += (MIN(7*8, delta)&0xf8) << (QBITS-2);
    }

    if (qlevel > 2)
        last += (MIN(i+15,127)&(15<<3))<<(QBITS+1);     // i>>3

    q2 = q1; q1 = q;
    }
    X[i] = last;

    /* The actual encoding */
    for (i = 0; i < len2; i++) {
    unsigned char q = (qual[i] - '!') & (QMAX-1);
    model_qual[X[i]].encodeSymbol(rc, q);
    }

    if (len != len2) {
    model_qual[X[i]].encodeSymbol(rc, QMAX-1); /* terminator */
    }
}
#else
void fqz::encode_qual(RangeCoder *rc, char *seq, char *qual, int len) {
    unsigned int last = 0;
    int delta = 5;
    int i, len2 = len;
    int q1 = 0, q2 = 0;

    /* Removing "Killer Bees" */
    while (len2 > 0 && qual[len2-1] == '#')
        len2--;

    for (i = 0; i < len2; i++) {
        unsigned char q = (qual[i] - '!') & (QMAX-1);

        #ifdef MULTI_QUAL_MODEL
        if (model_qual[last].bias() > model_qual[last & SMALL_QMASK].bias()) {
            model_qual[last].encodeSymbol(rc, q);
        } else {
            model_qual[last & SMALL_QMASK].encodeSymbol(rc, q);
            model_qual[last].updateSymbol(q);
        }
        #else
        model_qual[last].encodeSymbol(rc, q);
        #endif

        // previous 2-3 bytes
        if (QBITS == 12) {
            last = ((MAX(q1, q2)<<6) + q) & ((1<<QBITS)-1);
        } else {
            last = ((last << 6) + q) & ((1<<QBITS)-1);
        }

        if (qlevel > 1) {
            last  += (q1==q2) << QBITS;
            // delta saves 3-4%, but adds 14% cpu
            delta += (q1>q)*(q1-q);
            last  += (MIN(7*8, delta)&0xf8) << (QBITS-2);
        }

        if (qlevel > 2)
            //last += (MIN(i+7,127)&(7<<4))<<(QBITS);   // i>>4
            last += (MIN(i+15,127)&(15<<3))<<(QBITS+1);     // i>>3
            //last += (MIN(i+31,127)&(31<<2))<<(QBITS+2); // i>>2

        _mm_prefetch((const char *)&model_qual[last], _MM_HINT_T0);
        q2 = q1; q1 = q;

        assert(last < (QSIZE*16));
    }

    if (len != len2)
        model_qual[last].encodeSymbol(rc, QMAX-1); /* terminator */
}
#endif

/*
 * This attempts to encode qualities within a fixed distance, provided it
 * does appear to genuingly improve compression.
 *
 * Raw qual = 7722157 2.00 bpb
 * +/- 1    = 4092457 1.06 bpb
 * +/- 2    = 2853429 .739 bpb
 * +/- 3    = 2043044
 * (after first 1000 reads)
 * +/- 1    = 4758744 1.23 bpb
 * +/- 2    = 3354701 .869 bpb
 * +/- 3    = 2534989
 */

void fqz::decode_qual(RangeCoder *rc, char *qual, int len) {
    unsigned int last = 0;
    int i;
    int delta = 5;
    int q1 = 0, q2 = 0;

    for (i = 0; i < len; i++) {
        unsigned char q = model_qual[last].decodeSymbol(rc);

        if (q == QMAX-1) {
            while (i < len)
                qual[i++] = '#';
        }
        else {
            qual[i] = q + '!';
            if (QBITS == 12) {
                last = ((MAX(q1, q2)<<6) + q) & ((1<<QBITS)-1);
            }
            else {
                last = ((last << 6) + q) & ((1<<QBITS)-1);
            }

            if (qlevel > 1) {
                last  += (q1==q2) << QBITS;
                delta += (q1>q)*(q1-q);
                last  += (MIN(7*8, delta)&0xf8) << (QBITS-2);
            }

            if (qlevel > 2)
            //last += (MIN(i+7,127)&(7<<4))<<(QBITS);   // i>>4
                last += (MIN(i+15,127)&(15<<3))<<(QBITS+1);     // i>>3
            //last += (MIN(i+31,127)&(31<<2))<<(QBITS+2); // i>>2

            _mm_prefetch((const char *)&model_qual[last], _MM_HINT_T0);
            q2 = q1; q1 = q;
        }
    }
}

/* --------------------------------------------------------------------------
 * Compression functions.
 */


#ifdef TIMING
static long c1 = 0, c2 = 0, c3 = 0;
#endif

/* pthread enty points */
static void *fq_compress_r0(void *v) {
    ((fqz *)v)->compress_r0();
    return NULL;
}

static void *fq_compress_r1(void *v) {
    //fprintf(stderr, "r1 start on %d\n", sched_getcpu());
    ((fqz *)v)->compress_r1();
    //fprintf(stderr, "r1 end\n");
    return NULL;
}

static void *fq_compress_r2(void *v) {
    //fprintf(stderr, "r2 start on %d\n", sched_getcpu());
    ((fqz *)v)->compress_r2();
    //fprintf(stderr, "r2 end\n");
    return NULL;
}

static void *fq_compress_r3(void *v) {
    //fprintf(stderr, "r3 start on %d\n", sched_getcpu());
    ((fqz *)v)->compress_r3();
    //fprintf(stderr, "r3 end\n");
    return NULL;
}

/* Compute the block check sum */
void fqz::compress_r0() {
    //chksum = (do_hash && !qual_approx) ? sfhash(chk_in, chk_len) : 0;
    chksum = 0;
}

/* Sequence length & name */
void fqz::compress_r1() {
    char *name_p = name_buf;
    RangeCoder rc;

#ifdef TIMING
    clock_t c = clock();
#endif

    rc.output(out1);
    rc.StartEncode();
    if (nlevel == 1) {
        for (int i = 0; i < ns; i++) {
            encode_name(&rc, name_p, name_len_a[i]);
            name_p += name_len_a[i];
        }
    } else {
        for (int i = 0; i < ns; i++) {
            encode_name2(&rc, name_p, name_len_a[i]);
            name_p += name_len_a[i];
        }
    }
    rc.FinishEncode();

    sz1 = rc.size_out();
    name_in  += name_p - name_buf;
    name_out += sz1;

#ifdef TIMING
    c1 += clock() - c;
#endif
}

/* Sequence itself */
void fqz::compress_r2() {
    char *seq_p  = seq_buf;
    RangeCoder rc;

#ifdef TIMING
    clock_t c = clock();
#endif

    rc.output(out2);
    rc.StartEncode();
    for (int i = 0; i < ns; i++) {
        if (extreme_seq)
            encode_seq16(&rc, seq_p, seq_len_a[i]);
        else
            encode_seq8(&rc, seq_p, seq_len_a[i]);
        seq_p  += seq_len_a[i];
    }
    rc.FinishEncode();

    sz2 = rc.size_out();
    base_in  += seq_p - seq_buf;
    base_out += sz2;

#ifdef TIMING
    c2 += clock() - c;
#endif
}

/* Quality values */
void fqz::compress_r3() {
    char *qual_p = qual_buf;
    char *seq_p = seq_buf;
    RangeCoder rc;

#ifdef TIMING
    clock_t c = clock();
#endif

    rc.output(out3);
    rc.StartEncode();
    for (int i = 0; i < ns; i++) {
        encode_qual(&rc, seq_p, qual_p, seq_len_a[i]);
        qual_p += seq_len_a[i];
        seq_p  += seq_len_a[i];
    }
    rc.FinishEncode();

    sz3 = rc.size_out();
    qual_in  += qual_p - qual_buf;
    qual_out += sz3;

#ifdef TIMING
    c3 += clock() - c;
#endif
}

/*
 * Reads from in[0..in_len-1] and writes a compressed copy to out, setting
 * *out_len to the returned size. The caller needs to ensure out is large
 * enough.
 *
 * We can only compress entire sequences and in[] may end in a partial
 * sequence. Hence *in_end is set to a pointer to the end of the processed
 * buffer.
 *
 * Also returns *nseqs filled out.
 *
 * Returns total compressed length on success, with data in out[]
 *        -1 on failure
 */
int fqz::fq_compress(char *in,  int in_len,
                     char *out, int *out_len,
                     char **in_end, int *nseqs) {
    int end = 0, end_hash = 0;
    int i, j, k;
    //static char not_nl[256];
    static int not_nl[256];

    char *name_p = name_buf;
    char *seq_p  = seq_buf;
    char *qual_p = qual_buf;

    ns = 0;

    for (i = 0; i < 256; i++)
        not_nl[i] = 1;
    not_nl['\r'] = not_nl['\n'] = 0;

    /* Parse and separate into name, seq, qual buffers */
    seq_len = 0;

    // Safe method;
    for (i = k = 0; i < in_len; ) {
        char *name, *seq, *qual;

        /* Name */
        if (in[i] != '@')
            return -1;

        name = &in[i+1];
        j = i;
        in[k++] = in[i];
        i++;
        //while (i < in_len && in[i] != '\n' && in[i] != '\r')
        while (i < in_len && not_nl[(uc)in[i]])
            in[k++] = *name_p++ = in[i++];
        name_len_a[ns] = i-j-1;

        if (in[i] == '\r') i++;
        in[k++] = in[i];
        if (++i >= in_len)
            break;

        /* Sequence */
        seq = seq_p;
        //for (j = i; i < in_len && in[i] != '\n' && in[i] != '\r'; i++)
        for (j = i; i < in_len && not_nl[(uc)in[i]]; i++)
            in[k++] = *seq_p++ = in[i];
        seq_len_a[ns] = i-j;

        if (in[i] == '\r') i++;
        in[k++] = in[i];
        if (++i >= in_len)
            break;

        /* +name, assume to be identical to @name */
        if (in[i] != '+')
            return -1;
        in[k++] = in[i];

        //for (; i < in_len && in[i] != '\n' && in[i] != '\r'; i++)
        for (; i < in_len && not_nl[(uc)in[i]]; i++)
            ;
        if (in[i] == '\r') i++;
        in[k++] = in[i];
        if (++i >= in_len)
            break;

        /* Quality */
        qual = &in[i];
        static int is_N[256]={
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /*  0 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* 10 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* 20 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* 30 */
        1, 0, 1, 0,  1, 1, 1, 0,  1, 1, 1, 1,  1, 1, 1, 1, /* 40 */
        1, 1, 1, 1,  0, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* 50 */
        1, 0, 1, 0,  1, 1, 1, 0,  1, 1, 1, 1,  1, 1, 1, 1, /* 60 */
        1, 1, 1, 1,  0, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* 70 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* 80 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* 90 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* A0 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* B0 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* C10 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* D0 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* E0 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* F0 */
        };
        //for (j = 0; i < in_len && in[i] != '\n' && in[i] != '\r'; i++, j++) {
        for (j = 0; i < in_len && not_nl[(uc)in[i]]; i++, j++) {
            if (is_N[(unsigned char)seq[j]]) in[i] = '!'; // Ensure N is qual 0
            if (in[i] == '!' && !is_N[(unsigned char)seq[j]]) in[i] = '"';
            in[k++] = *qual_p++ = in[i];
        }

        if (in[i] == '\r') i++;
        in[k++] = in[i];
        if (++i > in_len)
            break;

        end = i; end_hash = k;

        if (seq_len == 0)
            seq_len = seq_len_a[ns];
        else if (seq_len != seq_len_a[ns])
            seq_len = -1;

        ns++;

        if (i == k)
            break; // Well behaved code. Go to the faster non-checking mode
    }

    // Faster method with no \r or +<name> checking
    for (; i < in_len; ) {
        char *name, *seq, *qual;

        /* Name */
        if (in[i] != '@')
            return -1;

        name = &in[i+1];
        j = i;
        i++;
        while (i < in_len && in[i] != '\n')
            *name_p++ = in[i++];
        name_len_a[ns] = i-j-1;

        if (++i >= in_len)
            break;

        /* Sequence */
        seq = seq_p;
        for (j = i; i < in_len && in[i] != '\n'; i++)
            *seq_p++ = in[i];
        seq_len_a[ns] = i-j;

        if (++i >= in_len)
            break;

        /* +name, assume to be identical to @name */
        if (in[i] != '+')
            return -1;

        for (; i < in_len && in[i] != '\n'; i++)
            ;
        if (++i >= in_len)
            break;

        /* Quality */
        qual = &in[i];
        static int is_N[256]={
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /*  0 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* 10 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* 20 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* 30 */
        1, 0, 1, 0,  1, 1, 1, 0,  1, 1, 1, 1,  1, 1, 1, 1, /* 40 */
        1, 1, 1, 1,  0, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* 50 */
        1, 0, 1, 0,  1, 1, 1, 0,  1, 1, 1, 1,  1, 1, 1, 1, /* 60 */
        1, 1, 1, 1,  0, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* 70 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* 80 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* 90 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* A0 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* B0 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* C10 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* D0 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* E0 */
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, /* F0 */
        };
        for (j = 0; i < in_len && in[i] != '\n'; i++, j++) {
            if (is_N[(unsigned char)seq[j]]) in[i] = '!'; // Ensure N is qual 0
            if (in[i] == '!' && !is_N[(unsigned char)seq[j]]) in[i] = '"';
            *qual_p++ = in[i];
        }

        if (++i > in_len)
            break;

        end = i; end_hash = i;

        if (seq_len == 0)
            seq_len = seq_len_a[ns];
        else if (seq_len != seq_len_a[ns])
            seq_len = -1;

        ns++;
    }

    *in_end = in + end;
    *nseqs = ns;

    /* Note: must be after seq==N qual editing code above */
    chk_in = (uc *) in;
    chk_len = end_hash;

    /* Encode seq len, we have a dependency on this for seq/qual */
    //fprintf(stderr, "-----\n");
    RangeCoder rc;
    rc.output(out0);
    rc.StartEncode();
    for (int i = 0; i < ns; i++) {
        //fprintf(stderr, "Encode %d: %d\n", i, seq_len_a[i]);
        encode_len(&rc, seq_len_a[i]);
    }
    rc.FinishEncode();
    sz0 = rc.size_out();

#if 1
    /* Encode the 3 buffers in parallel */
#ifdef PTHREADS
    if (do_threads) {
    pthread_t t0, t1, t2, t3;
    pthread_create(&t0, NULL, fq_compress_r0, (void*)this);
    pthread_create(&t1, NULL, fq_compress_r1, (void*)this);
    pthread_create(&t2, NULL, fq_compress_r2, (void*)this); 
    pthread_create(&t3, NULL, fq_compress_r3, (void*)this);

    pthread_join(t0, NULL);
    pthread_join(t1, NULL);
    pthread_join(t2, NULL);
    pthread_join(t3, NULL);
    } else {
    compress_r0();
    compress_r1();
    compress_r2();
    compress_r3();
    }
#else
    compress_r0();
    compress_r1();
    compress_r2();
    compress_r3();
#endif
#endif

    //fprintf(stderr, "hashes %08x %08x %08x\n", name_hash, seq_hash, qual_hash);
    
    /* Concatenate compressed output into a single block */
    char *out_p = out;
    *out_p++ = (chksum >>  0) & 0xff;
    *out_p++ = (chksum >>  8) & 0xff;
    *out_p++ = (chksum >> 16) & 0xff;
    *out_p++ = (chksum >> 24) & 0xff;

    *out_p++ = (end >>  0) & 0xff;  /* Uncompressed size */
    *out_p++ = (end >>  8) & 0xff;
    *out_p++ = (end >> 16) & 0xff;
    *out_p++ = (end >> 24) & 0xff;

    *out_p++ = (ns  >>  0) & 0xff;  /* Number of sequences */
    *out_p++ = (ns  >>  8) & 0xff;
    *out_p++ = (ns  >> 16) & 0xff;
    *out_p++ = (ns  >> 24) & 0xff;

    *out_p++ = (sz0 >>  0) & 0xff;  /* Size of 4 range-coder blocks */
    *out_p++ = (sz0 >>  8) & 0xff;
    *out_p++ = (sz0 >> 16) & 0xff;
    *out_p++ = (sz0 >> 24) & 0xff;

    *out_p++ = (sz1 >>  0) & 0xff;
    *out_p++ = (sz1 >>  8) & 0xff;
    *out_p++ = (sz1 >> 16) & 0xff;
    *out_p++ = (sz1 >> 24) & 0xff;

    *out_p++ = (sz2 >>  0) & 0xff;
    *out_p++ = (sz2 >>  8) & 0xff;
    *out_p++ = (sz2 >> 16) & 0xff;
    *out_p++ = (sz2 >> 24) & 0xff;

    *out_p++ = (sz3 >>  0) & 0xff;
    *out_p++ = (sz3 >>  8) & 0xff;
    *out_p++ = (sz3 >> 16) & 0xff;
    *out_p++ = (sz3 >> 24) & 0xff;

    memcpy(out_p, out0, sz0); out_p += sz0;
    memcpy(out_p, out1, sz1); out_p += sz1;
    memcpy(out_p, out2, sz2); out_p += sz2;
    memcpy(out_p, out3, sz3); out_p += sz3;

    *out_len = out_p - out;

    return *out_len;
}

int fqz::isq_compress(char *in,  int in_len,
                     char *out, int *out_len,
                     char **in_end, int *nseqs) {
    return *out_len;
}

int fqz::iq_compress(char *in,  int in_len,
                      char *out, int *out_len,
                      char **in_end, int *nseqs) {
    return *out_len;
}


/*
 * A blocking read that refuses to return truncated reads. //拒绝返回被删节的reads
 */
int xget(std::fstream &in, unsigned char *in_buffer, int count) {
    in.read((char *) in_buffer, count);
    if (in.bad() || in.fail() && !in.eof())
        return -1;
    return (int) in.gcount();
}

int xwrite(std::fstream &out, unsigned char *out_buffer, int count) {
    out.write((const char *) out_buffer, count);
    if (out.bad() || out.fail()) {
        return -1;
    }
    return count;
}

/*
 * 把数据传进buffer，如果buffer的大小超过blockSize，或者传入空字符串就编码输出
 */
int fqz::se_iq_encode(std::string &id, std::string &qual, std::fstream &out) {
    if (inLen+id.length()+qual.length()+2<=BLK_SIZE && id != ""){
        inLen += id.length()+qual.length()+2; // 2 for '\n'
        id = id.substr(1); //delete the '@'
        name_p = (char*)id.data();
        name_p += (int)id.length();
        name_len_a[ns] = (int)id.length();
        seq_len_a[ns] = (int)qual.length();
        qual_p = (char*)qual.data();
        qual_p += (int)qual.length();
        if (seq_len == 0)
            seq_len = (int)qual.length();
        else if (seq_len != qual.length())
            seq_len = -1;
        ns ++;
    }
    else{
        //编码
        RangeCoder rc;
        rc.output(out0);
        rc.StartEncode();
        for (int i = 0; i < ns; i++) {
            encode_len(&rc, seq_len_a[i]);
        }
        rc.FinishEncode();
        sz0 = rc.size_out();

        compress_r0();
        compress_r1();
        compress_r3();

        char *out_p0 = out_buf+4;
        char *out_p = out_p0;
        *out_p++ = (chksum >>  0) & 0xff;
        *out_p++ = (chksum >>  8) & 0xff;
        *out_p++ = (chksum >> 16) & 0xff;
        *out_p++ = (chksum >> 24) & 0xff;

        *out_p++ = (inLen >>  0) & 0xff;  /* Uncompressed size */
        *out_p++ = (inLen >>  8) & 0xff;
        *out_p++ = (inLen >> 16) & 0xff;
        *out_p++ = (inLen >> 24) & 0xff;

        *out_p++ = (ns  >>  0) & 0xff;  /* Number of sequences */
        *out_p++ = (ns  >>  8) & 0xff;
        *out_p++ = (ns  >> 16) & 0xff;
        *out_p++ = (ns  >> 24) & 0xff;

        *out_p++ = (sz0 >>  0) & 0xff;  /* Size of 4 range-coder blocks */
        *out_p++ = (sz0 >>  8) & 0xff;
        *out_p++ = (sz0 >> 16) & 0xff;
        *out_p++ = (sz0 >> 24) & 0xff;

        *out_p++ = (sz1 >>  0) & 0xff;
        *out_p++ = (sz1 >>  8) & 0xff;
        *out_p++ = (sz1 >> 16) & 0xff;
        *out_p++ = (sz1 >> 24) & 0xff;

        *out_p++ = (sz3 >>  0) & 0xff;
        *out_p++ = (sz3 >>  8) & 0xff;
        *out_p++ = (sz3 >> 16) & 0xff;
        *out_p++ = (sz3 >> 24) & 0xff;

        memcpy(out_p, out0, sz0); out_p += sz0;
        memcpy(out_p, out1, sz1); out_p += sz1;
        memcpy(out_p, out3, sz3); out_p += sz3;

        int out_len = (int)(out_p - out_p0);
        out_buf[0] = (out_len >>  0) & 0xff;
        out_buf[1] = (out_len >>  8) & 0xff;
        out_buf[2] = (out_len >> 16) & 0xff;
        out_buf[3] = (out_len >> 24) & 0xff;

        if (out_len != xwrite(out, (unsigned char*)out_buf, out_len)) {
            fprintf(stderr, "Abort: truncated write.0\n");
            return -1;
        }

        //还原
        inLen = 0;
        seq_len = 0;
        ns = 0;
        name_p = name_buf;
        qual_p = qual_buf;
    }
    if (id != ""){
        inLen += id.length()+qual.length()+2; // 2 for '\n'
        id = id.substr(1); //delete the '@'
        name_p = (char*)id.data();
        name_p += (int)id.length();
        name_len_a[ns] = (int)id.length();
        seq_len_a[ns] = (int)qual.length();
        qual_p = (char*)qual.data();
        qual_p += (int)qual.length();
        if (seq_len == 0)
            seq_len = (int)qual.length();
        else if (seq_len != qual.length())
            seq_len = -1;
        ns ++;
    }
    return 0;
}

int fqz::pe_iq_encode(std::string &id, std::string &qual, std::fstream &out) {
    if (id != ""){
        if (!readBufMark){ //为read1
            readBuf[0] = id;
            readBuf[1] = qual;
            readBufMark ++;
        }
        else{ //为read2
            readBuf[2] = id;
            readBuf[3] = qual;
            readBufMark --;
        }
    }
    if (inLen+readBuf[0].length()+readBuf[1].length()+readBuf[2].length()+readBuf[3].length()+4>BLK_SIZE || id == ""){
        //编码
        RangeCoder rc;
        rc.output(out0);
        rc.StartEncode();
        for (int i = 0; i < ns; i++) {
            encode_len(&rc, seq_len_a[i]);
        }
        rc.FinishEncode();
        sz0 = rc.size_out();

        compress_r0();
        compress_r1();
        compress_r3();

        char *out_p0 = out_buf+4;
        char *out_p = out_p0;
        *out_p++ = (chksum >>  0) & 0xff;
        *out_p++ = (chksum >>  8) & 0xff;
        *out_p++ = (chksum >> 16) & 0xff;
        *out_p++ = (chksum >> 24) & 0xff;

        *out_p++ = (inLen >>  0) & 0xff;  /* Uncompressed size */
        *out_p++ = (inLen >>  8) & 0xff;
        *out_p++ = (inLen >> 16) & 0xff;
        *out_p++ = (inLen >> 24) & 0xff;

        *out_p++ = (ns  >>  0) & 0xff;  /* Number of sequences */
        *out_p++ = (ns  >>  8) & 0xff;
        *out_p++ = (ns  >> 16) & 0xff;
        *out_p++ = (ns  >> 24) & 0xff;

        *out_p++ = (sz0 >>  0) & 0xff;  /* Size of 4 range-coder blocks */
        *out_p++ = (sz0 >>  8) & 0xff;
        *out_p++ = (sz0 >> 16) & 0xff;
        *out_p++ = (sz0 >> 24) & 0xff;

        *out_p++ = (sz1 >>  0) & 0xff;
        *out_p++ = (sz1 >>  8) & 0xff;
        *out_p++ = (sz1 >> 16) & 0xff;
        *out_p++ = (sz1 >> 24) & 0xff;

        *out_p++ = (sz3 >>  0) & 0xff;
        *out_p++ = (sz3 >>  8) & 0xff;
        *out_p++ = (sz3 >> 16) & 0xff;
        *out_p++ = (sz3 >> 24) & 0xff;

        memcpy(out_p, out0, sz0); out_p += sz0;
        memcpy(out_p, out1, sz1); out_p += sz1;
        memcpy(out_p, out3, sz3); out_p += sz3;

        int out_len = (int)(out_p - out_p0);
        out_buf[0] = (out_len >>  0) & 0xff;
        out_buf[1] = (out_len >>  8) & 0xff;
        out_buf[2] = (out_len >> 16) & 0xff;
        out_buf[3] = (out_len >> 24) & 0xff;

        if (out_len != xwrite(out, (unsigned char*)out_buf, out_len)) {
            fprintf(stderr, "Abort: truncated write.0\n");
            return -1;
        }

        //还原
        inLen = 0;
        seq_len = 0;
        ns = 0;
        name_p = name_buf;
        qual_p = qual_buf;
    }
    if (inLen+readBuf[0].length()+readBuf[1].length()+readBuf[2].length()+readBuf[3].length()+4<=BLK_SIZE){
        //pass readBuf to buffer
        inLen += readBuf[0].length()+readBuf[1].length()+2;
        readBuf[0] = readBuf[0].substr(1);
        name_p = (char*)readBuf[0].data();
        name_p += (int)readBuf[0].length();
        name_len_a[ns] = (int)readBuf[0].length();
        seq_len_a[ns] = (int)readBuf[1].length();
        qual_p = (char*)readBuf[1].data();
        qual_p += (int)readBuf[1].length();
        if (seq_len == 0)
            seq_len = (int)readBuf[1].length();
        else if (seq_len != readBuf[1].length())
            seq_len = -1;
        ns ++;

        inLen += readBuf[2].length()+readBuf[3].length()+2;
        readBuf[2] = readBuf[2].substr(1);
        name_p = (char*)readBuf[2].data();
        name_p += (int)readBuf[2].length();
        name_len_a[ns] = (int)readBuf[2].length();
        seq_len_a[ns] = (int)readBuf[3].length();
        qual_p = (char*)readBuf[3].data();
        qual_p += (int)readBuf[3].length();
        if (seq_len != readBuf[3].length())
            seq_len = -1;
        ns ++;
    }
    return 0;
}

int fqz::se_isq_encode(std::string &id, std::string &seq, std::string &qual, std::fstream &out) {
    if (inLen+id.length()+seq.length()+qual.length()+3<=BLK_SIZE && id != ""){
        inLen += id.length()+seq.length()+qual.length()+3;
        id = id.substr(1);
        name_p = (char*)id.data();
        name_p += (int)id.length();
        name_len_a[ns] = (int)id.length();
        seq_len_a[ns] = (int)qual.length();
        seq_p = (char*)seq.data();
        seq_p += (int)seq.length();
        qual_p = (char*)qual.data();
        qual_p += (int)qual.length();

        if (seq_len == 0)
            seq_len = (int)qual.length();
        else if (seq_len != qual.length())
            seq_len = -1;
        ns ++;
    }
    else{
        //编码
        RangeCoder rc;
        rc.output(out0);
        rc.StartEncode();
        for (int i = 0; i < ns; i++) {
            encode_len(&rc, seq_len_a[i]);
        }
        rc.FinishEncode();
        sz0 = rc.size_out();

        compress_r0();
        compress_r1();
        compress_r2();
        compress_r3();

        char *out_p0 = out_buf+4;
        char *out_p = out_p0;
        *out_p++ = (chksum >>  0) & 0xff;
        *out_p++ = (chksum >>  8) & 0xff;
        *out_p++ = (chksum >> 16) & 0xff;
        *out_p++ = (chksum >> 24) & 0xff;

        *out_p++ = (inLen >>  0) & 0xff;  /* Uncompressed size */
        *out_p++ = (inLen >>  8) & 0xff;
        *out_p++ = (inLen >> 16) & 0xff;
        *out_p++ = (inLen >> 24) & 0xff;

        *out_p++ = (ns  >>  0) & 0xff;  /* Number of sequences */
        *out_p++ = (ns  >>  8) & 0xff;
        *out_p++ = (ns  >> 16) & 0xff;
        *out_p++ = (ns  >> 24) & 0xff;

        *out_p++ = (sz0 >>  0) & 0xff;  /* Size of 4 range-coder blocks */
        *out_p++ = (sz0 >>  8) & 0xff;
        *out_p++ = (sz0 >> 16) & 0xff;
        *out_p++ = (sz0 >> 24) & 0xff;

        *out_p++ = (sz1 >>  0) & 0xff;
        *out_p++ = (sz1 >>  8) & 0xff;
        *out_p++ = (sz1 >> 16) & 0xff;
        *out_p++ = (sz1 >> 24) & 0xff;

        *out_p++ = (sz2 >>  0) & 0xff;
        *out_p++ = (sz2 >>  8) & 0xff;
        *out_p++ = (sz2 >> 16) & 0xff;
        *out_p++ = (sz2 >> 24) & 0xff;

        *out_p++ = (sz3 >>  0) & 0xff;
        *out_p++ = (sz3 >>  8) & 0xff;
        *out_p++ = (sz3 >> 16) & 0xff;
        *out_p++ = (sz3 >> 24) & 0xff;

        memcpy(out_p, out0, sz0); out_p += sz0;
        memcpy(out_p, out1, sz1); out_p += sz1;
        memcpy(out_p, out2, sz2); out_p += sz2;
        memcpy(out_p, out3, sz3); out_p += sz3;

        int out_len = (int)(out_p - out_p0);
        out_buf[0] = (out_len >>  0) & 0xff;
        out_buf[1] = (out_len >>  8) & 0xff;
        out_buf[2] = (out_len >> 16) & 0xff;
        out_buf[3] = (out_len >> 24) & 0xff;

        if (out_len != xwrite(out, (unsigned char*)out_buf, out_len)) {
            fprintf(stderr, "Abort: truncated write.0\n");
            return -1;
        }

        //还原
        inLen = 0;
        seq_len = 0;
        ns = 0;
        name_p = name_buf;
        seq_p  = seq_buf;
        qual_p = qual_buf;
    }
    if (id != ""){
        inLen += id.length()+seq.length()+qual.length()+3;
        id = id.substr(1);
        name_p = (char*)id.data();
        name_p += (int)id.length();
        name_len_a[ns] = (int)id.length();
        seq_len_a[ns] = (int)qual.length();
        seq_p = (char*)seq.data();
        seq_p += (int)seq.length();
        qual_p = (char*)qual.data();
        qual_p += (int)qual.length();

        if (seq_len == 0)
            seq_len = (int)qual.length();
        else if (seq_len != qual.length())
            seq_len = -1;
        ns ++;
    }
    return 0;
}

int fqz::pe_isq_encode(std::string &id, std::string &seq, std::string &qual, std::fstream &out) {
    if (id != ""){
        if (!readBufMark){ //为read1
            readBuf[0] = id;
            readBuf[1] = seq;
            readBuf[2] = qual;
            readBufMark ++;
        }
        else{ //为read2
            readBuf[3] = id;
            readBuf[4] = seq;
            readBuf[5] = qual;
            readBufMark --;
        }
    }
    if (inLen+readBuf[0].length()+readBuf[1].length()+readBuf[2].length()+readBuf[3].length()+readBuf[4].length()+readBuf[5].length()+6>BLK_SIZE || id == ""){
        //编码
        RangeCoder rc;
        rc.output(out0);
        rc.StartEncode();
        for (int i = 0; i < ns; i++) {
            encode_len(&rc, seq_len_a[i]);
        }
        rc.FinishEncode();
        sz0 = rc.size_out();

        compress_r0();
        compress_r1();
        compress_r2();
        compress_r3();

        char *out_p0 = out_buf+4;
        char *out_p = out_p0;
        *out_p++ = (chksum >>  0) & 0xff;
        *out_p++ = (chksum >>  8) & 0xff;
        *out_p++ = (chksum >> 16) & 0xff;
        *out_p++ = (chksum >> 24) & 0xff;

        *out_p++ = (inLen >>  0) & 0xff;  /* Uncompressed size */
        *out_p++ = (inLen >>  8) & 0xff;
        *out_p++ = (inLen >> 16) & 0xff;
        *out_p++ = (inLen >> 24) & 0xff;

        *out_p++ = (ns  >>  0) & 0xff;  /* Number of sequences */
        *out_p++ = (ns  >>  8) & 0xff;
        *out_p++ = (ns  >> 16) & 0xff;
        *out_p++ = (ns  >> 24) & 0xff;

        *out_p++ = (sz0 >>  0) & 0xff;  /* Size of 4 range-coder blocks */
        *out_p++ = (sz0 >>  8) & 0xff;
        *out_p++ = (sz0 >> 16) & 0xff;
        *out_p++ = (sz0 >> 24) & 0xff;

        *out_p++ = (sz1 >>  0) & 0xff;
        *out_p++ = (sz1 >>  8) & 0xff;
        *out_p++ = (sz1 >> 16) & 0xff;
        *out_p++ = (sz1 >> 24) & 0xff;

        *out_p++ = (sz2 >>  0) & 0xff;
        *out_p++ = (sz2 >>  8) & 0xff;
        *out_p++ = (sz2 >> 16) & 0xff;
        *out_p++ = (sz2 >> 24) & 0xff;

        *out_p++ = (sz3 >>  0) & 0xff;
        *out_p++ = (sz3 >>  8) & 0xff;
        *out_p++ = (sz3 >> 16) & 0xff;
        *out_p++ = (sz3 >> 24) & 0xff;

        memcpy(out_p, out0, sz0); out_p += sz0;
        memcpy(out_p, out1, sz1); out_p += sz1;
        memcpy(out_p, out2, sz2); out_p += sz2;
        memcpy(out_p, out3, sz3); out_p += sz3;

        int out_len = (int)(out_p - out_p0);
        out_buf[0] = (out_len >>  0) & 0xff;
        out_buf[1] = (out_len >>  8) & 0xff;
        out_buf[2] = (out_len >> 16) & 0xff;
        out_buf[3] = (out_len >> 24) & 0xff;

        if (out_len != xwrite(out, (unsigned char*)out_buf, out_len)) {
            fprintf(stderr, "Abort: truncated write.0\n");
            return -1;
        }

        //还原
        inLen = 0;
        seq_len = 0;
        ns = 0;
        name_p = name_buf;
        seq_p  = seq_buf;
        qual_p = qual_buf;
    }
    if (inLen+readBuf[0].length()+readBuf[1].length()+readBuf[2].length()+readBuf[3].length()+readBuf[4].length()+readBuf[5].length()+6<=BLK_SIZE){
        //pass readBuf to buffer
        inLen += readBuf[0].length()+readBuf[1].length()+readBuf[2].length()+3;
        readBuf[0] = readBuf[0].substr(1);
        name_p = (char*)readBuf[0].data();
        name_p += (int)readBuf[0].length();
        name_len_a[ns] = (int)readBuf[2].length();
        seq_len_a[ns] = (int)readBuf[2].length();
        seq_p = (char*)readBuf[1].data();
        seq_p += (int)readBuf[1].length();
        qual_p = (char*)readBuf[2].data();
        qual_p += (int)readBuf[2].length();
        if (seq_len == 0)
            seq_len = (int)readBuf[1].length();
        else if (seq_len != readBuf[1].length())
            seq_len = -1;
        ns ++;

        inLen += readBuf[3].length()+readBuf[4].length()+readBuf[5].length()+3;
        readBuf[3] = readBuf[3].substr(1);
        name_p = (char*)readBuf[3].data();
        name_p += (int)readBuf[3].length();
        name_len_a[ns] = (int)readBuf[5].length();
        seq_len_a[ns] = (int)readBuf[5].length();
        seq_p = (char*)readBuf[4].data();
        seq_p += (int)readBuf[4].length();
        qual_p = (char*)readBuf[5].data();
        qual_p += (int)readBuf[5].length();
        if (seq_len != readBuf[5].length())
            seq_len = -1;
        ns ++;
    }
    return 0;
}

/*
 * Encode an entire stream
 *
 * Returns 0 on success
 *        -1 on failure
 */
int fqz::encode(std::fstream &in, std::fstream &out) {
    int sz, blk_start = 0;
    size_t total_sz = 0;

    /*
     * Parse one block at a time. Blocks may not terminate on exact fastq
     * boundaries, so we need to know where we ended processing and move
     * that partial fastq entry back to the start for the next block.
     *
     * We write out the block size too so we can decompress block at a time.
     * This may also permits a level of parallellism in the future.
     */
    while ((sz = xget(in,(unsigned char*)in_buf+blk_start, BLK_SIZE - blk_start)) > 0) {
        char *comp, *in_end = NULL;
        int comp_len = 0, nseqs = 0;

        if (-1 == fq_compress(in_buf, sz + blk_start,
                              out_buf+4, &comp_len,
                              &in_end, &nseqs)) {
            fprintf(stderr, "Failure to parse and/or compress.\n");
            return -1;
        }
        comp = out_buf;

        if (comp == NULL) {
            fprintf(stderr, "Abort: encode_data returned NULL\n");
            return -1;
        }

        out_buf[0] = (comp_len >>  0) & 0xff;
        out_buf[1] = (comp_len >>  8) & 0xff;
        out_buf[2] = (comp_len >> 16) & 0xff;
        out_buf[3] = (comp_len >> 24) & 0xff;
        comp_len += 4;

        if (comp_len != xwrite(out, (unsigned char*)out_buf, comp_len)) {
            fprintf(stderr, "Abort: truncated write.0\n");
            return -1;
        }

        total_sz += comp_len;

        /* We maybe ended on a partial fastq entry, so start from there */
        memmove(in_buf, in_end, (sz + blk_start) - (in_end - in_buf));
        blk_start = (sz + blk_start) - (in_end - in_buf);
    }

    return 0;
}


/* --------------------------------------------------------------------------
 * Decompression functions.
 */
#define DECODE_INT(a) ((a)[0] + ((a)[1]<<8) + ((a)[2]<<16) + ((a)[3]<<24))

/* pthread enty points */
static void *fq_decompress_r1(void *v) {
    //fprintf(stderr, "r1 start on %d\n", sched_getcpu());
    ((fqz *)v)->decompress_r1();
    //fprintf(stderr, "r1 end\n");
    return NULL;
}

static void *fq_decompress_r2(void *v) {
    //fprintf(stderr, "r2 start on %d\n", sched_getcpu());
    ((fqz *)v)->decompress_r2();
    //fprintf(stderr, "r2 end\n");
    return NULL;
}

static void *fq_decompress_r3(void *v) {
    //fprintf(stderr, "r3 start on %d\n", sched_getcpu());
    ((fqz *)v)->decompress_r3();
    //fprintf(stderr, "r3 end\n");
    return NULL;
}

void fqz::decompress_r1(void) {
    RangeCoder rc;
    rc.input(in_buf1);
    rc.StartDecode();

    char *name_p = name_buf;
    if (nlevel == 1) {
        for (int i = 0; i < ns; i++) {
            *name_p++ = '@';
            name_p += decode_name(&rc, name_p);
            *name_p++ = '\n';
        }
    } else {
        for (int i = 0; i < ns; i++) {
            *name_p++ = '@';
            name_p += decode_name2(&rc, name_p);
            *name_p++ = '\n';
        }
    }
    rc.FinishDecode();
}

void fqz::decompress_r2(void) {
    RangeCoder rc;
    rc.input(in_buf2);
    rc.StartDecode();

    char *seq_p = seq_buf;
    for (int i = 0; i < ns; i++) {
        if (extreme_seq)
            decode_seq16(&rc, seq_p, seq_len_a[i]);
        else
            decode_seq8(&rc, seq_p, seq_len_a[i]);
        seq_p += seq_len_a[i];
    }
    rc.FinishDecode();
}

void fqz::decompress_r3(void) {
    RangeCoder rc;
    rc.input(in_buf3);
    rc.StartDecode();

    char *qual_p = qual_buf;
    for (int i = 0; i < ns; i++) {
        decode_qual(&rc, qual_p, seq_len_a[i]);
        qual_p += seq_len_a[i];
    }
    rc.FinishDecode();
}

/* Decompress a single block */
char *fqz::fq_decompress(char *in, int comp_len, int *out_len) {
    char *name_p, *seq_p, *qual_p;

    uint32_t chk    = DECODE_INT((unsigned char *)(in));
    //uint32_t ulen   = DECODE_INT((unsigned char *)in+4);
    uint32_t nseqs  = DECODE_INT((unsigned char *)(in+8));
    uint32_t sz0    = DECODE_INT((unsigned char *)(in+12));
    uint32_t sz1    = DECODE_INT((unsigned char *)(in+16));
    uint32_t sz2    = DECODE_INT((unsigned char *)(in+20));
    uint32_t sz3    = DECODE_INT((unsigned char *)(in+24));

    in += 28;
    ns = nseqs;

    /* Use ulen for allocating decoding buffers */

//    fprintf(stderr, "%d -> %d\n", comp_len, ulen);
//    fprintf(stderr, "   ns=%d, sz={%d, %d, %d, %d}\n",
//      nseqs, sz0, sz1, sz2, sz3);

    in_buf0 = in; in += sz0;
    in_buf1 = in; in += sz1;
    in_buf2 = in; in += sz2;
    in_buf3 = in; in += sz3;

    RangeCoder rc0;
    rc0.input(in_buf0);
    rc0.StartDecode();

    for (int i = 0; i < ns; i++)
        seq_len_a[i] = decode_len(&rc0);
    rc0.FinishDecode();

#ifdef PTHREADS
    if (do_threads && qlevel <= 3) {
    /* -q4 adds dependency between seq[] and qual[] */
    pthread_t t1, t2, t3;
    pthread_create(&t1, NULL, fq_decompress_r1, (void*)this);
    pthread_create(&t2, NULL, fq_decompress_r2, (void*)this);
    pthread_create(&t3, NULL, fq_decompress_r3, (void*)this);

    pthread_join(t1, NULL);
    pthread_join(t2, NULL);
    pthread_join(t3, NULL);
    } else {
    decompress_r1();
    decompress_r2();
    decompress_r3();
    }
#else
    decompress_r1();
    decompress_r2();
    decompress_r3();
#endif

    //fprintf(stderr, "hashes %08x %08x %08x\n", name_hash, seq_hash, qual_hash);

    /* Stick together the arrays into out_buf */
    out_ind = 0;
    name_p = name_buf;
    seq_p = seq_buf;
    qual_p = qual_buf;

    for (int i = 0; i < ns; i++) {
        /* name */
        while ((out_buf[out_ind++] = *name_p++) != '\n')
            ;

        /* seq */
        for (int j = 0; j < seq_len_a[i]; j++)
            out_buf[out_ind++] = *seq_p++;
        out_buf[out_ind++] = '\n';
        out_buf[out_ind++] = '+';
        out_buf[out_ind++] = '\n';

        /* qual */
        for (int j = 0; j < seq_len_a[i]; j++) {
            if ((out_buf[out_ind++] = *qual_p++) == '!') {
                out_buf[out_ind-4-seq_len_a[i]] = 'N';
            }
        }
        out_buf[out_ind++] = '\n';
    }


    //chksum = do_hash ? sfhash((uc *)out_buf, out_ind) : 0;
    chksum = 0;
    if (do_hash && chk && chksum != chk) {
        fprintf(stderr, "Mismatching checksums. Aborting. Rerun with -X to ignore this error.\n");
        return NULL;
    }

    *out_len = out_ind;
    return out_buf;
}

char *fqz::iq_decompress(char *in, int comp_len, int *out_len) {
    char *name_p, *qual_p;

    uint32_t chk    = DECODE_INT((unsigned char *)(in));
    //uint32_t ulen   = DECODE_INT((unsigned char *)in+4);
    uint32_t nseqs  = DECODE_INT((unsigned char *)(in+8));
    uint32_t sz0    = DECODE_INT((unsigned char *)(in+12));
    uint32_t sz1    = DECODE_INT((unsigned char *)(in+16));
    uint32_t sz3    = DECODE_INT((unsigned char *)(in+20));

    in += 24;
    ns = nseqs;

    /* Use ulen for allocating decoding buffers */

//    fprintf(stderr, "%d -> %d\n", comp_len, ulen);
//    fprintf(stderr, "   ns=%d, sz={%d, %d, %d, %d}\n",
//      nseqs, sz0, sz1, sz2, sz3);

    in_buf0 = in; in += sz0;
    in_buf1 = in; in += sz1;
    in_buf3 = in; in += sz3;

    RangeCoder rc0;
    rc0.input(in_buf0);
    rc0.StartDecode();

    for (int i = 0; i < ns; i++)
        seq_len_a[i] = decode_len(&rc0);
    rc0.FinishDecode();

#ifdef PTHREADS
    if (do_threads && qlevel <= 3) {
    /* -q4 adds dependency between seq[] and qual[] */
    pthread_t t1, t2, t3;
    pthread_create(&t1, NULL, fq_decompress_r1, (void*)this);
    pthread_create(&t2, NULL, fq_decompress_r2, (void*)this);
    pthread_create(&t3, NULL, fq_decompress_r3, (void*)this);

    pthread_join(t1, NULL);
    pthread_join(t2, NULL);
    pthread_join(t3, NULL);
    } else {
    decompress_r1();
    decompress_r2();
    decompress_r3();
    }
#else
    decompress_r1();
    decompress_r3();
#endif

    //fprintf(stderr, "hashes %08x %08x %08x\n", name_hash, seq_hash, qual_hash);

    /* Stick together the arrays into out_buf */
    out_ind = 0;
    name_p = name_buf;
    qual_p = qual_buf;

    for (int i = 0; i < ns; i++) {
        /* name */
        while ((out_buf[out_ind++] = *name_p++) != '\n')
            ;

        /* qual */
        for (int j = 0; j < seq_len_a[i]; j++) {
            if ((out_buf[out_ind++] = *qual_p++) == '!') {
                out_buf[out_ind-4-seq_len_a[i]] = 'N';
            }
        }
        out_buf[out_ind++] = '\n';
    }

    //chksum = do_hash ? sfhash((uc *)out_buf, out_ind) : 0;
    chksum = 0;
    if (do_hash && chk && chksum != chk) {
        fprintf(stderr, "Mismatching checksums. Aborting. Rerun with -X to ignore this error.\n");
        return NULL;
    }

    *out_len = out_ind;
    return out_buf;
}

char *fqz::isq_decompress(char *in, int comp_len, int *out_len) {
    char *name_p, *seq_p, *qual_p;

    uint32_t chk    = DECODE_INT((unsigned char *)(in));
    //uint32_t ulen   = DECODE_INT((unsigned char *)in+4);
    uint32_t nseqs  = DECODE_INT((unsigned char *)(in+8));
    uint32_t sz0    = DECODE_INT((unsigned char *)(in+12));
    uint32_t sz1    = DECODE_INT((unsigned char *)(in+16));
    uint32_t sz2    = DECODE_INT((unsigned char *)(in+20));
    uint32_t sz3    = DECODE_INT((unsigned char *)(in+24));

    in += 28;
    ns = nseqs;

    /* Use ulen for allocating decoding buffers */

//    fprintf(stderr, "%d -> %d\n", comp_len, ulen);
//    fprintf(stderr, "   ns=%d, sz={%d, %d, %d, %d}\n",
//      nseqs, sz0, sz1, sz2, sz3);

    in_buf0 = in; in += sz0;
    in_buf1 = in; in += sz1;
    in_buf2 = in; in += sz2;
    in_buf3 = in; in += sz3;

    RangeCoder rc0;
    rc0.input(in_buf0);
    rc0.StartDecode();

    for (int i = 0; i < ns; i++)
        seq_len_a[i] = decode_len(&rc0);
    rc0.FinishDecode();

#ifdef PTHREADS
    if (do_threads && qlevel <= 3) {
    /* -q4 adds dependency between seq[] and qual[] */
    pthread_t t1, t2, t3;
    pthread_create(&t1, NULL, fq_decompress_r1, (void*)this);
    pthread_create(&t2, NULL, fq_decompress_r2, (void*)this);
    pthread_create(&t3, NULL, fq_decompress_r3, (void*)this);

    pthread_join(t1, NULL);
    pthread_join(t2, NULL);
    pthread_join(t3, NULL);
    } else {
    decompress_r1();
    decompress_r2();
    decompress_r3();
    }
#else
    decompress_r1();
    decompress_r2();
    decompress_r3();
#endif

    //fprintf(stderr, "hashes %08x %08x %08x\n", name_hash, seq_hash, qual_hash);

    /* Stick together the arrays into out_buf */
    out_ind = 0;
    name_p = name_buf;
    seq_p = seq_buf;
    qual_p = qual_buf;

    for (int i = 0; i < ns; i++) {
        /* name */
        while ((out_buf[out_ind++] = *name_p++) != '\n')
            ;

        /* seq */
        for (int j = 0; j < seq_len_a[i]; j++)
            out_buf[out_ind++] = *seq_p++;
        out_buf[out_ind++] = '\n';

        /* qual */
        for (int j = 0; j < seq_len_a[i]; j++) {
            if ((out_buf[out_ind++] = *qual_p++) == '!') {
                out_buf[out_ind-4-seq_len_a[i]] = 'N';
            }
        }
        out_buf[out_ind++] = '\n';
    }

    //chksum = do_hash ? sfhash((uc *)out_buf, out_ind) : 0;
    chksum = 0;
    if (do_hash && chk && chksum != chk) {
        fprintf(stderr, "Mismatching checksums. Aborting. Rerun with -X to ignore this error.\n");
        return NULL;
    }

    *out_len = out_ind;
    return out_buf;
}

int fqz::iq_decode(std::fstream &in, std::string &out1, std::string &out2) {
    unsigned char len_buf[4];
    char *uncomp_buf;

    if (pass_len >= uncomp_len){
        pass_len = 0;
        if (4 == xget(in, len_buf, 4)){
            int32_t comp_len =
                    (len_buf[0] <<  0) +
                    (len_buf[1] <<  8) +
                    (len_buf[2] << 16) +
                    (len_buf[3] << 24);
            int rem_len = comp_len, in_off = 0;

            do {
                errno = 0;
                int tmp_len = xget(in, (unsigned char *) in_buf + in_off, rem_len);
                if (errno == EINTR && tmp_len == -1)
                    continue;

                if (tmp_len == -1) {
                    fprintf(stderr, "Abort: read failed, %d.\n", errno);
                    perror("foo");
                    return -1;
                }
                if (tmp_len == 0) {
                    fprintf(stderr, "Abort: truncated read, %d.\n", errno);
                    return -1;
                }
                rem_len -= tmp_len;
                in_off  += tmp_len;
            } while (rem_len);

            uncomp_buf = iq_decompress(in_buf, comp_len, &uncomp_len);
        }
    }
    if (uncomp_buf) {
        out1 = out2 = "";
        int i = 0;
        while (uncomp_buf[pass_len+i] != '\n'){
            out1 += uncomp_buf[pass_len+i];
            i++;
        }
        i++;
        while (uncomp_buf[pass_len+i] != '\n'){
            out2 += uncomp_buf[pass_len+i];
            i++;
        }
        i++;
        pass_len += i;
    } else {
        fprintf(stderr, "Failed to decompress block\n");
        return -1;
    }
    return 0;
}

int fqz::isq_decode(std::fstream &in, std::string &out1, std::string &out2, std::string &out3) {
    unsigned char len_buf[4];
    char *uncomp_buf;

    if (pass_len >= uncomp_len){
        pass_len = 0;
        if (4 == xget(in, len_buf, 4)){
            int32_t comp_len =
                    (len_buf[0] <<  0) +
                    (len_buf[1] <<  8) +
                    (len_buf[2] << 16) +
                    (len_buf[3] << 24);
            int rem_len = comp_len, in_off = 0;

            do {
                errno = 0;
                int tmp_len = xget(in, (unsigned char *) in_buf + in_off, rem_len);
                if (errno == EINTR && tmp_len == -1)
                    continue;

                if (tmp_len == -1) {
                    fprintf(stderr, "Abort: read failed, %d.\n", errno);
                    perror("foo");
                    return -1;
                }
                if (tmp_len == 0) {
                    fprintf(stderr, "Abort: truncated read, %d.\n", errno);
                    return -1;
                }
                rem_len -= tmp_len;
                in_off  += tmp_len;
            } while (rem_len);

            uncomp_buf = isq_decompress(in_buf, comp_len, &uncomp_len);
        }
    }
    if (uncomp_buf) {
        out1 = out2 = out3 = "";
        int i = 0;
        while (uncomp_buf[pass_len+i] != '\n'){
            out1 += uncomp_buf[pass_len+i];
            i++;
        }
        i++;
        while (uncomp_buf[pass_len+i] != '\n'){
            out2 += uncomp_buf[pass_len+i];
            i++;
        }
        i++;
        while (uncomp_buf[pass_len+i] != '\n'){
            out3 += uncomp_buf[pass_len+i];
            i++;
        }
        i++;
        pass_len += i;
    } else {
        fprintf(stderr, "Failed to decompress block\n");
        return -1;
    }
    return 0;
}


/*
 * Decode an entire stream
 *
 * Returns 0 on success
 *        -1 on failure
 */
int fqz::decode(std::fstream &in, std::fstream &out) {
    unsigned char len_buf[4];

    /*
     * Parse one block at a time. Blocks may not terminate on exact fastq
     * boundaries, so we need to know where we ended processing and move
     * that partial fastq entry back to the start for the next block.
     *
     * We write out the block size too so we can decompress block at a time.
     * This may also permits a level of parallellism in the future.
     */
    while (4 == xget(in, len_buf, 4)) {
        int32_t comp_len =
                (len_buf[0] <<  0) +
                (len_buf[1] <<  8) +
                (len_buf[2] << 16) +
                (len_buf[3] << 24);
        char *uncomp_buf;
        int   uncomp_len, rem_len = comp_len, in_off = 0;

        //fprintf(stderr, "Block of length %d\n", comp_len);

        do {
            errno = 0;
            int tmp_len = xget(in, (unsigned char *) in_buf + in_off, rem_len);
            if (errno == EINTR && tmp_len == -1)
                continue;

            if (tmp_len == -1) {
                fprintf(stderr, "Abort: read failed, %d.\n", errno);
                perror("foo");
                return -1;
            }
            if (tmp_len == 0) {
                fprintf(stderr, "Abort: truncated read, %d.\n", errno);
                return -1;
            }
            rem_len -= tmp_len;
            in_off  += tmp_len;
        } while (rem_len);

        uncomp_buf = fq_decompress(in_buf, comp_len, &uncomp_len);

        if (uncomp_buf) {
            if (uncomp_len != xwrite(out, (unsigned char *) uncomp_buf, uncomp_len)) {
                fprintf(stderr, "Abort: truncated write.1\n");
                return -1;
            }

        } else {
            fprintf(stderr, "Failed to decompress block\n");
            return -1;
        }
    }

    return 0;
}
