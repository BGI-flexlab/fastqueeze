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
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
#include "bntseq.h"
#include "bwtaln.h"
#include "kseq.h" // for the FASTA/Q parser
#include "bwt.h"

KSEQ_DECLARE(gzFile)

#define MIN(x,y) ((x)<(y)?(x):(y))

int bwtintv_cmp(const void *arg1, const void *arg2) {
    return (uint32_t)(*(bwtintv_t *)arg1).info - (uint32_t)(*(bwtintv_t *)arg2).info>>32;
}

//int mis_count(uint8_t* ref, ){
//    //这里需要ref_Seq, read_Se，和是否reverse三个参数，输出是
//}

int fastsam(int argc, char *argv[])
{
    int c, i, base, min_iwidth = 20, min_len = 17, print_seq = 0, min_intv = 1, max_len = INT_MAX, lgst_num = 3, max_mis = 5, qual_sys = 2;
    uint64_t max_intv = 0;
    kseq_t *seq;
    bwtint_t k;
    gzFile fp;
    smem_i *itr;
    const bwtintv_v *a;
    bwaidx_t *idx;

    while ((c = getopt(argc, argv, "w:l:pi:I:L:")) >= 0) {
        switch (c) {
            case 'w': min_iwidth = atoi(optarg); break;
            case 'l': min_len = atoi(optarg); break;
            case 'i': min_intv = atoi(optarg); break;
            case 'I': max_intv = atol(optarg); break;
            case 'L': max_len  = atoi(optarg); break;
            case 'f': lgst_num = atoi(optarg); break;
            case 'm': max_mis = atoi(optarg); break;
            case 'q': qual_sys = atoi(optarg); break;
            default: return 1;
        }
    }

    if (optind + 1 >= argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   fastsam [options] <idxbase> <in.fq>\n\n");
        fprintf(stderr, "Options: -l INT    min SMEM length to output [%d]\n", min_len);
        fprintf(stderr, "         -w INT    max interval size to find coordiantes [%d]\n", min_iwidth);
        fprintf(stderr, "         -i INT    min SMEM interval size [%d]\n", min_intv);
        fprintf(stderr, "         -L INT    max MEM length [%d]\n", max_len);
        fprintf(stderr, "         -I INT    stop if MEM is longer than -l with a size less than INT [%ld]\n", (long)max_intv);
        fprintf(stderr, "         -f INT    consider only the longest [%d] MEM\n", lgst_num);
        fprintf(stderr, "         -m INT    max mismatch to tolerate [%d]\n", max_mis);
        fprintf(stderr, "         -q INT    quality system, 1:illumina, 2:sanger, default as [%d]\n", qual_sys);
        fprintf(stderr, "\n");
        return 1;
    }

    fp = xzopen(argv[optind + 1], "r");
    seq = kseq_init(fp);
    if ((idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1;
    itr = smem_itr_init(idx->bwt);
    smem_config(itr, min_intv, max_len, max_intv);
    while (kseq_read(seq) >= 0) {
        size_t rlen = seq->seq.l; //转存
        err_printf("SQ\t%s\t%ld", seq->name.s, seq->seq.l);
        err_putchar('\n');
        for (i = 0; i < seq->seq.l; ++i){
            seq->seq.s[i] = nst_nt4_table[(int)seq->seq.s[i]];
            if (seq->seq.s[i] >= 4){ //Ambiguous Bases的质量值转为0
                if (qual_sys == 2) seq->qual.s[i] = '!';
                else seq->qual.s[i] = '@';
            }
        }
        smem_set_query(itr, seq->seq.l, (uint8_t*)seq->seq.s);
        while ((a = smem_next(itr)) != 0) {
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
            qsort(plist, a->n-short_num, sizeof(bwtintv_t *), bwtintv_cmp);
            for (i = 0; i < MIN(a->n-short_num,lgst_num); ++i){
                err_printf("EM\t%d\t%d\t%ld", (uint32_t)(plist[i]->info>>32), (uint32_t)plist[i]->info, (long)plist[i]->x[2]);
                if (plist[i]->x[2] <= min_iwidth) {
                    for (k = 0; k < plist[i]->x[2]; ++k) {
                        bwtint_t pos;
                        int len, is_rev, ref_id;
                        len = (uint32_t)plist[i]->info - (uint32_t)(plist[i]->info>>32);
                        pos = bns_depos(idx->bns, bwt_sa(idx->bwt, plist[i]->x[0] + k), &is_rev);
//                        //测试
//                        uint8_t *rseq, *rseq_l, *rseq_r;
//                        uint16_t missum = 0;
//                        if (is_rev){
//                            pos -= len - 1;
//                            rseq_l = bns_get_seq(idx->bns->l_pac, idx->pac, pos-1-seq->seq.l+(uint32_t)(plist[i]->info), pos-1, &rlen);
//                            //把rseq_l和read的右截的反向互补进行比较，注意mismatch的坐标转化
//                            for (base=0; base<rlen; base++){
//                                if (rseq_l[base] + seq->seq.s[seq->seq.l-1-base] != 3){
//                                    missum += 1;
//
//                                }
//
//                            }
//
//                            rseq_r = bns_get_seq(idx->bns->l_pac, idx->pac, pos+len, pos+(uint32_t)(plist[i]->info), &rlen);
//                            //把rseq_r和read的左截的反向互补进行比较，注意mismatch的坐标转化
//                            for (base=0; base<rlen; base++){
//
//                            }
//                        }
//                        else{
//                            rseq_l = bns_get_seq(idx->bns->l_pac, idx->pac, pos-1-(uint32_t)(plist[i]->info>>32), pos-1, &rlen);
//                            //把rseq_l和read的左截进行比较，注意mismatch的坐标转化
//                            rseq_r = bns_get_seq(idx->bns->l_pac, idx->pac, pos+len, pos+seq->seq.l-(uint32_t)(plist[i]->info>>32), &rlen);
//                            //把rseq_r和read的右截进行比较，注意mismatch的坐标转化
//                        }
//                        //err_printf("(%d\t%d\t%d\t%d\t%d\t%d\t%d)", rseq[0], rseq[1], rseq[2], rseq[3], rseq[4], rseq[5], rseq[6]);
//                        //rseq = bns_get_seq(idx->bns->l_pac, pac, p->x[0], p->x[1], &rlen);
//                        //err_printf("(%d)", rseq[0]);
//                        //测试
//                        rseq = bns_get_seq(idx->bns->l_pac, idx->pac, 0, 10, &rlen);
                        bns_cnt_ambi(idx->bns, pos, len, &ref_id);
                        err_printf("\t%s:%c%ld", idx->bns->anns[ref_id].name, "+-"[is_rev], (long)(pos - idx->bns->anns[ref_id].offset) + 1);
                    }
                } else err_puts("\t*");
                err_putchar('\n');
            }
        }
        err_puts("//");
    }

    smem_itr_destroy(itr);
    bwa_idx_destroy(idx);
    kseq_destroy(seq);
    err_gzclose(fp);
    return 0;
}
