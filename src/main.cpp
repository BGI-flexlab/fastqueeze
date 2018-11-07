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

/* -------------------------------------------------------------------------
 * BWA
 */

int bwtintv_cmp(const void *arg1, const void *arg2) {     //长的SMEM排前面
    return  ((uint32_t) (*(bwtintv_t *) arg2).info - (uint32_t) ((*(bwtintv_t *) arg2).info >> 32)) -
            ((uint32_t) (*(bwtintv_t *) arg1).info - (uint32_t) ((*(bwtintv_t *) arg1).info >> 32));
}

bool sam_cmp(align_info sam1, align_info sam2){
    if (sam1.blockNum != sam2.blockNum)
        return sam1.blockNum < sam2.blockNum;
    else
        return sam1.blockPos < sam2.blockPos;
}

void CreatBitmap(map<char, map<char, int>> &bitmap) {
    map<char, int> map_G;
    map<char, int> map_C;
    map<char, int> map_A;
    map<char, int> map_T;
    map_G['1'] = 0; map_G['0'] = 1; map_G['3'] = 0;
    map_C['2'] = 0; map_C['0'] = 1; map_C['3'] = 2;
    map_A['2'] = 0; map_A['1'] = 1; map_A['3'] = 2;
    map_T['2'] = 0; map_T['1'] = 1; map_T['0'] = 2;
    bitmap['2'] = map_G;
    bitmap['1'] = map_C;
    bitmap['0'] = map_A;
    bitmap['3'] = map_T;
}

int getAlignInfo(kseq seq, smem_i* func_itr, bwaidx_t *func_idx, align_info *align_p, int func_block_size, int min_len, int max_iwidth, int max_mis, int lgst_num){
    int64_t rlen;
    int i, seql, base;

    const bwtintv_v *a;
    seql = (int) seq.seq.length();
    int pass_num = 0;
    int cigar_l[max_mis], cigar_v[max_mis];

    static map<char, map<char, int>> nucleBitMap;
    if (nucleBitMap.empty())
        CreatBitmap(nucleBitMap);

    for (i = 0; i < seql; ++i) {
        seq.seq[i] = nst_nt4_table[(int) seq.seq[i]];
    }
    smem_set_query(func_itr, seql, (uint8_t *) seq.seq.c_str());
    while ((a = smem_next(func_itr)) != 0) {
        bwtintv_t *plist[a->n];
        int short_num = 0;
        for (i = 0; i < a->n; ++i) {
            bwtintv_t *p = &a->a[i];
            if ((uint32_t) p->info - (p->info >> 32) < min_len) {
                short_num += 1;
                continue;
            } else {
                plist[i - short_num] = &a->a[i];
            }
        }
        if (a->n == short_num)
            continue;
        qsort(plist, a->n - short_num, sizeof(bwtintv_t *), bwtintv_cmp);
        for (i = 0; i < a->n - short_num; ++i) {
            if (plist[i]->x[2] <= max_iwidth) {
                for (int k = 0; k < plist[i]->x[2]; ++k) {
                    bwtint_t pos;
                    int is_rev;
                    pos = (bwtint_t)bns_depos(func_idx->bns, bwt_sa(func_idx->bwt, plist[i]->x[0] + k), &is_rev);
                    uint8_t *rseq_l, *rseq_r;
                    uint16_t missum = 0;
                    int rbase, sbase;
                    if (is_rev) {
                        pos -= seql - (uint32_t)(plist[i]->info >>32)-1;
                        rseq_l = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos,
                                             pos + seql - (uint32_t) (plist[i]->info), &rlen);
                        //把rseq_l和read的右截的反向互补进行比较
                        for (base = 0; base < rlen; base++) {
                            if ((rseq_l[base] + seq.seq[seql-base-1] != 3) && (seq.seq[seql-base-1] <= 3)) {
                                if (missum >= max_mis) break;
                                rbase = int(rseq_l[base]);
                                sbase = int(seq.seq[seql-base-1]);
                                cigar_l[missum] = base;
                                cigar_v[missum] = sbase <= 3 ? nucleBitMap[rbase][3-sbase] : nucleBitMap[rbase][sbase];
                                missum += 1;
                            }
                        }
                        if (missum <= max_mis) {
                            rseq_r = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos + seql - (uint32_t)(plist[i]->info >>32),
                                                 pos + seql, &rlen);
                            //把rseq_r和read的左截的反向互补进行比较
                            for (base = 0; base < rlen; base++) {
                                if ((rseq_r[base] + seq.seq[(uint32_t)(plist[i]->info >>32)-base-1] != 3) && (seq.seq[(uint32_t)(plist[i]->info >>32)-base-1] <= 3)) {
                                    if (missum >= max_mis) break;
                                    rbase = int(rseq_r[base]);
                                    sbase = int(seq.seq[(uint32_t)(plist[i]->info >>32)-base-1]);
                                    cigar_l[missum] = seql-(uint32_t)(plist[i]->info >>32)+base;
                                    cigar_v[missum] = sbase <= 3 ? nucleBitMap[rbase][3-sbase] : nucleBitMap[rbase][sbase];
                                    missum += 1;
                                }
                            }
                        }
                        if (missum <= max_mis){
                            (align_p+pass_num)->blockNum = (int) (pos / func_block_size);
                            (align_p+pass_num)->blockPos = (int) (pos % func_block_size);
                            (align_p+pass_num)->isRev = (bool) is_rev;
                            for (int x=0; x<missum; x++){
                                (align_p+pass_num)->cigar_l[x] = cigar_l[x];
                                (align_p+pass_num)->cigar_v[x] = cigar_v[x];
                            }
                            pass_num += 1;
                            if (pass_num >= lgst_num)
                                return pass_num;
                        }
                    } else {
                        pos -= (uint32_t) (plist[i]->info >> 32);
                        rseq_l = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos,
                                             pos+(uint32_t) (plist[i]->info >> 32), &rlen);
                        //把rseq_l和read的左截进行比较
                        for (base = 0; base < rlen; base++) {
                            if ((rseq_l[base] != seq.seq[base]) && (seq.seq[base] <= 3)) {
                                if (missum >= max_mis) break;
                                rbase = int(rseq_l[base]);
                                sbase = int(seq.seq[base]);
                                cigar_l[missum] = base;
                                cigar_v[missum] = sbase <= 3 ? nucleBitMap[rbase][3-sbase] : nucleBitMap[rbase][sbase];
                                missum += 1;
                            }
                        }
                        if (missum <= max_mis) {
                            rseq_r = bns_get_seq(func_idx->bns->l_pac, func_idx->pac, pos+ (uint32_t)(plist[i]->info),
                                                 pos + seql, &rlen);
                            //把rseq_r和read的右截进行比较
                            for (base = 0; base < rlen; base++) {
                                if ((rseq_r[base] != seq.seq[(uint32_t) plist[i]->info + base]) && (seq.seq[(uint32_t) plist[i]->info + base] <= 3)) {
                                    if (missum >= max_mis) break;
                                    rbase = int(rseq_r[base]);
                                    sbase = int(seq.seq[(uint32_t) plist[i]->info + base]);
                                    cigar_l[missum] = base+(uint32_t)(plist[i]->info);
                                    cigar_v[missum] = sbase <= 3 ? nucleBitMap[rbase][3-sbase] : nucleBitMap[rbase][sbase];
                                    missum += 1;
                                }
                            }
                        }
                        if (missum <= max_mis){
                            (align_p+pass_num)->blockNum = (int) (pos / func_block_size);
                            (align_p+pass_num)->blockPos = (int) (pos % func_block_size);
                            (align_p+pass_num)->isRev = (bool) is_rev;
                            for (int x=0; x<missum; x++){
                                (align_p+pass_num)->cigar_l[x] = cigar_l[x];
                                (align_p+pass_num)->cigar_v[x] = cigar_v[x];
                            }
                            pass_num += 1;
                            if (pass_num >= lgst_num)
                                return pass_num;
                        }
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
            if (int(qual[i]) > 74){
                cout << "Qual Sysytm wrong." << endl;
                exit(1);
            }
            if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T') {
                seq[i] = 'N';
                qual[i] = '!';
            }
        }
    }
    else{ //illumina
        for (int i=0;i<seq.length();i++){
            if (int(qual[i]) < 64){
                cout << "Qual Sysytm wrong." << endl;
                exit(1);
            }
            if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T') {
                seq[i] = 'N';
                qual[i] = '@';
            }
        }
    }
}

void readModify2(string& seq, string& quality, int qualSys){
    if (qualSys-1){ //sanger
        for (int i=0;i<seq.length();i++){
            if (int(quality[i]) == 33)
                seq[i] = 'N';
        }
    }
    else{ //illumina
        for (int i=0;i<seq.length();i++){
            if (int(quality[i]) == 64)
                seq[i] = 'N';
        }
    }
}

/* -------------------------------------------------------------------------
 * Main program
 */

static void usage(int err) {
    FILE *fp = err ? stderr : stdout;

    fprintf(fp, "SeqArc v%d.%d. Yuxin Chen, 2018\n", MAJOR_VERS, MINOR_VERS);
    fprintf(fp, "The entropy coder is derived from Fqzcomp. The aligner is derived from BWA.\n\n");

    fprintf(fp, "To build index:\n  SeqArc -i <ref.fa>\n\n");

    fprintf(fp, "To compress:\n  SeqArc [options] <ref.fa> <input_file> <output_prefix>\n\n");
    fprintf(fp, "    -l INT         min SMEM length to output [17]\n");
    fprintf(fp, "    -w INT         max interval size to find coordiantes [100]\n");
    fprintf(fp, "    -I INT         skip MEM mapped to over [-] places\n");
    fprintf(fp, "    -f INT         consider only the longest [3] MEM\n");
    fprintf(fp, "    -m INT         max mismatch to tolerate [3]\n");
    fprintf(fp, "    -B INT         number of block to split reference [50]\n");
    fprintf(fp, "    -q INT         quality system, 1:illumina, 2:sanger, default as [2]\n");
    fprintf(fp, "    -s INT         max insert size between read1 and read2 [500]\n\n");

    fprintf(fp, "    -S <level>     Sequence de novo compression level. 1-9 [3]\n");
    fprintf(fp, "                   Specifying '+' on the end (eg -s5+) will use\n");
    fprintf(fp, "                   models of multiple sizes for improved compression.\n");
    fprintf(fp, "    -N <level>     Quality compression level.  1-3 [2]\n");
    fprintf(fp, "    -n <level>     Name compression level.  1-2 [2]\n");
    fprintf(fp, "    -b             Use both strands in sequence hash table.\n");
    fprintf(fp, "    -e             Extra seq compression: 16-bit vs 8-bit counters.\n");
    fprintf(fp, "    -P             Disable multi-threading\n\n");

    fprintf(fp, "    -X             Disable generation/verification of check sums\n\n");

    fprintf(fp, "To decompress:\n   SeqArc -d <ref.fa> <compressed_prefix> <fastq_prefix>\n");

    exit(err);
}

int main(int argc, char **argv) {
    int opt, i, max_iwidth = 100, min_len = 17, max_len = INT_MAX, lgst_num = 3, qual_sys = 2;
    int block_num = 50, block_size;
    int max_mis = 3, max_insr = 511, max_readLen = 255;
    int se_mark = 1;
    uint64_t max_intv = 0;
    int seq1l, seq2l, seq1m, seq2m;
    kseq seq1, seq2;
    bwtint_t k;
    gzFile fp1, fp2;
    smem_i *itr;
    const bwtintv_v *a;
    bwaidx_t *idx;

    int decompress = 0, indexing = 0;
    char *ref;

    fqz_params p;

    /* Initialise and parse command line arguments */
    p.slevel = 3;
    p.qlevel = 2;
    p.nlevel = 2;
    p.both_strands = 0;
    p.extreme_seq = 0;
    p.multi_seq_model = 0;
    p.qual_approx = 0;
    p.do_threads = 1;
    p.do_hash = 1;

    while ((opt = getopt(argc, argv, "l:w:I:f:m:q:s:hdQ:S:N:bePXiB:")) != -1) {
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

            default:
                usage(1);
        }
    }

    if (argc == optind)
        usage(1);

    if (indexing){
        ref = argv[optind];
        bwa_idx_build(ref, ref, BWTALGO_AUTO, 10000000);
        return 1;
    }
    else if (decompress) {
        fstream fa, in_s, in_iq, in_isq, out1, out2;
        fa.open(argv[optind], std::ios::in);

        stringstream strInputPath;
        strInputPath.str("");
        strInputPath << "./" << argv[optind+1] << "_s.arc";
        in_s.open(strInputPath.str(), std::ios::binary|std::ios::in);

        strInputPath.str("");
        strInputPath << "./" << argv[optind+1] << "_iq.arc";
        in_iq.open(strInputPath.str(), std::ios::binary|std::ios::in);

        strInputPath.str("");
        strInputPath << "./" << argv[optind+1] << "_isq.arc";
        in_isq.open(strInputPath.str(), std::ios::binary|std::ios::in);

        unsigned char magic_s[11];
        if (11 != xget(in_s, magic_s, 11)){
            fprintf(stderr, "Abort: truncated read.\n");
            return 1;
        }
        qual_sys = magic_s[0];
        se_mark = magic_s[1];
        memcpy(&block_size, magic_s+2, 4);
        max_mis = magic_s[6];
        memcpy(&max_insr, magic_s+7, 2);
        memcpy(&max_readLen, magic_s+9, 2);

        decode seqdecoder(block_size, max_mis, max_insr, max_readLen);
        ref2seq ref2seqer(block_size, max_mis, max_readLen, fa);

        unsigned char magic_fqz[8];
        if (8 != xget(in_iq, magic_fqz, 8)) {
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

        p.slevel = magic_fqz[6] & 0x0f;
        p.qlevel = ((magic_fqz[6] >> 4) & 3);
        p.nlevel = (magic_fqz[6] >> 6);
        p.both_strands    = magic_fqz[7] & 1;
        p.extreme_seq     = magic_fqz[7] & 2;
        p.multi_seq_model = magic_fqz[7] & 4;

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

        fqz *f;
        f = new fqz(&p);

        strInputPath.str("");
        strInputPath << "./" << argv[optind+2] << "_1.fastq";
        out1.open(strInputPath.str(), std::ios::out);
        if (!se_mark){
            strInputPath.str("");
            strInputPath << "./" << argv[optind+2] << "_2.fastq";
            out2.open(strInputPath.str(), std::ios::out);
        }

        string name, seq, qual;
        align_info align_info1;
        int readLen;

        while (true){
            if (se_mark){
                if (-1 == f->iq_decode(in_iq, name, qual)) //不确定单用一边做休止符是否会引发异常
                    break;
                seqdecoder.parse_se(align_info1, readLen, in_s);
                seq = ref2seqer.getSeq(align_info1, readLen);
                readModify2(seq, qual, qual_sys);
                out1 << name << endl << seq << endl << "+" << endl << "@" << qual << endl;
            }
            else{
                if (-1 == f->iq_decode(in_iq, name, qual))
                    break;
                seqdecoder.parse_pe(align_info1, readLen, in_s);
                seq = ref2seqer.getSeq(align_info1, readLen);
                readModify2(seq, qual, qual_sys);
                out1 << name << endl << seq << endl << "+" << endl << "@" << qual << endl;

                if (-1 == f->iq_decode(in_iq, name, qual))
                    break;
                seqdecoder.parse_pe(align_info1, readLen, in_s);
                seq = ref2seqer.getSeq(align_info1, readLen);
                readModify2(seq, qual, qual_sys);
                out2 << name << endl << seq << endl << "+" << endl << "@" << qual << endl;
            }
        }

        fa.close();
        in_s.close();
        in_iq.close();
        in_isq.close();
        out1.close();
        out2.close();
        return 0;

    } else {
        string outIndex, seq1Name, seq2Name; // 可以选择放弃comment内容
        fp1 = xzopen(argv[optind + 1], "r");
        FunctorZlib gzr;
        kstream<gzFile, FunctorZlib> ks1(fp1, gzr);
        kstream<gzFile, FunctorZlib> *ks2;
        if (argc - optind >= 4){
            se_mark = 0;
            fp2 = xzopen(argv[optind + 2], "r");
            ks2 = new kstream<gzFile, FunctorZlib> (fp2,gzr);
            outIndex = argv[optind + 3];
        } else
            outIndex = argv[optind + 2];

        if ((idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1;
        itr = smem_itr_init(idx->bwt);
        smem_config(itr, 1, max_len, max_intv); //min_intv = 1
        block_size = (int) ceil(idx->bns->l_pac/block_num); //单个block的长度

        unsigned char magic_s[11]{(uint8_t)qual_sys, //block_size+max_mis+max_insr+max_readLen
                                  (uint8_t)se_mark,
                                  (uint8_t)block_size,
                                  (uint8_t)(block_size>>8),
                                  (uint8_t)(block_size>>16),
                                  (uint8_t)(block_size>>24),
                                  (uint8_t)max_mis,
                                  (uint8_t)max_insr,
                                  (uint8_t)(max_insr>>8),
                                  (uint8_t)max_readLen,
                                  (uint8_t)(max_readLen>>8)
        };

        int level = p.slevel | (p.qlevel << 4) | (p.nlevel << 6);
        int flags = p.both_strands
                    + p.extreme_seq*2
                    + p.multi_seq_model*4;
        unsigned char magic_fqz[8] = {'.', 'a', 'r', 'c',  //生成magic作为解压时参数
                                  MAJOR_VERS,
                                  FORMAT_VERS,
                                  (uint8_t)level,
                                  (uint8_t)flags
        };
        fqz* f[block_num+1]; //one more for isq
        for (i=0;i<block_num+1;i++)
            f[i] = new fqz(&p);

        stringstream strOutputPath;

        fstream out_s;
        strOutputPath.str("");
        strOutputPath << "./" << outIndex << "_aligninfo.arc";
        out_s.open(strOutputPath.str(), std::ios::binary|std::ios::out);

        fstream fpOutput_s[block_num]; //声明align_info的输出文件指针数组
        for (i=0;i<block_num;i++) {
            strOutputPath.str("");
            strOutputPath << "./s_" << i << ".tmp"; //合成文件路径
            fpOutput_s[i].open(strOutputPath.str(), std::ios::binary|std::ios::out);
        }

        if (11 != xwrite(out_s, magic_s, 11)) {
            fprintf(stderr, "Abort: truncated write.2\n");
            out_s.close();
            return 1;
        }

        fstream out_iq;
        strOutputPath.str("");
        strOutputPath << "./" << outIndex << "_iq.arc";
        out_iq.open(strOutputPath.str(), std::ios::binary|std::ios::out);

        fstream fpOutput_iq[block_num]; //声明id+qual的输出文件指针数组
        for (i=0;i<block_num;i++) {
            strOutputPath.str("");
            strOutputPath << "./iq_" << i << ".tmp"; //合成文件路径
            fpOutput_iq[i].open(strOutputPath.str(), std::ios::binary|std::ios::out);
        }

        fstream out_isq; //声明id+seq+qual的输出文件指针
        strOutputPath.str("");
        strOutputPath << "./" << outIndex << "_isq.arc";
        out_isq.open(strOutputPath.str(), std::ios::binary|std::ios::out);

        if (8 != xwrite(out_iq, magic_fqz, 8)) {
            fprintf(stderr, "Abort: truncated write.2\n");
            out_iq.close();
            return 1;
        }

        align_info sam1[lgst_num];
        align_info sam2[lgst_num];
        for (i=0;i<lgst_num;i++){
            sam1[i].cigar_l = (int*) malloc(max_mis*sizeof(int));
            sam1[i].cigar_v = (int*) malloc(max_mis*sizeof(int));
            sam2[i].cigar_l = (int*) malloc(max_mis*sizeof(int));
            sam2[i].cigar_v = (int*) malloc(max_mis*sizeof(int));
        }
        for (i=0;i<lgst_num;i++){
            memset(sam1[i].cigar_l,-1,max_mis*sizeof(int));
            memset(sam2[i].cigar_l,-1,max_mis*sizeof(int));
        }

        encode* encoders[block_num];
        for (i = 0; i < block_num; i++)
            encoders[i] = new encode(se_mark, block_size, max_mis, max_insr, max_readLen);

        while ((seq1l = ks1.read(seq1)) >= 0) {
            seq1Name = seq1.name + " " + seq1.comment;
            if (se_mark){ //SE
                readModify1(seq1.seq, seq1.qual, qual_sys, max_readLen);
                if (seq1m = getAlignInfo(seq1, itr, idx, sam1, block_size, min_len, max_iwidth, max_mis, lgst_num)){
                    std::sort(sam1, sam1+seq1m, sam_cmp);
                    encoders[sam1[0].blockNum]->parse_1(sam1[0], seq1l, fpOutput_s[sam1[0].blockNum]); //对sam1[0]，即比对位置最前的结果进行处理，这个操作是为了尽量使比对位置集中
                    f[sam1[0].blockNum]->iq_encode(seq1Name, seq1.qual, fpOutput_iq[sam1[0].blockNum]); //id+qual
                }
                else{
                    f[block_num]->isq_encode(seq1Name, seq1.seq, seq1.qual, out_isq); //id+seq+qual
                }
                for (i=0;i<lgst_num;i++)//把sam1清零
                    memset(sam1[i].cigar_l,-1,max_mis*sizeof(int));
            }
            else{ //PE
                seq2l = (*ks2).read(seq2);
                seq2Name = seq2.name + " " + seq2.comment;
                readModify1(seq2.seq, seq2.qual, qual_sys, max_readLen);
                if ((seq1m = getAlignInfo(seq1, itr, idx, sam1, block_size, min_len, max_iwidth, max_mis, lgst_num)) && (seq2m = getAlignInfo(seq2, itr, idx, sam2, block_size, min_len, max_iwidth, max_mis, lgst_num))){
                    std::sort(sam1, sam1+seq1m, sam_cmp);
                    std::sort(sam2, sam2+seq2m, sam_cmp);
                    int x = 0; int y = 0, find = 0;
                    while (x <= seq1m && y <= seq2m) {
                        if (sam1[x].blockNum < sam2[y].blockNum)
                            x += 1;
                        else if (sam1[x].blockNum > sam2[y].blockNum)
                            y += 1;
                        else {
                            if (abs(sam1[x].blockPos - sam2[y].blockPos) <= max_insr){
                                find = 1;
                                encoders[sam1[x].blockNum]->parse_1(sam1[x], seq1l, fpOutput_s[sam1[x].blockNum]);
                                encoders[sam2[y].blockNum]->parse_2(sam2[y], seq2l, fpOutput_s[sam2[y].blockNum]);
                                f[sam1[x].blockNum]->iq_encode(seq1Name, seq1.qual, fpOutput_iq[sam1[x].blockNum]);
                                f[sam2[y].blockNum]->iq_encode(seq2Name, seq2.qual, fpOutput_iq[sam2[y].blockNum]);
                            }
                            else if (sam1[x].blockPos < sam2[y].blockPos)
                                x += 1;
                            else
                                y += 1;
                        }
                    }
                    if (!find){ //没比上
                        f[block_num]->isq_encode(seq1Name, seq1.seq, seq1.qual, out_isq);
                        f[block_num]->isq_encode(seq2Name, seq2.seq, seq2.qual, out_isq);
                    }
                }
                else{//没比上
                    f[block_num]->isq_encode(seq1Name, seq1.seq, seq1.qual, out_isq);
                    f[block_num]->isq_encode(seq2Name, seq2.seq, seq2.qual, out_isq);
                }
                for (i=0;i<lgst_num;i++){     //把sam1和sam2清零
                    memset(sam1[i].cigar_l,0,max_mis*sizeof(int));
                    memset(sam2[i].cigar_l,-1,max_mis*sizeof(int));
                }
            }
        }
        //收尾
        string nullstr;
        for (i=0;i<block_num;i++){
            encoders[i]->end(fpOutput_s[i]);
            if (se_mark)
                f[i]->iq_encode(nullstr, nullstr, fpOutput_iq[i]);
            else
                f[i]->iq_encode(nullstr, nullstr, fpOutput_iq[i]);
        }
        f[block_num]->isq_encode(nullstr, nullstr, nullstr, out_isq);

        smem_itr_destroy(itr);
        bwa_idx_destroy(idx);
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
        out_s.close();

        return 0;
    }
}