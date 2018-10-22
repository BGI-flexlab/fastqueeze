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

int bwtintv_cmp(const void *arg1, const void *arg2) {     //长的SMEM排前面
    return (uint32_t) (*(bwtintv_t *) arg2).info - (uint32_t) (*(bwtintv_t *) arg2).info
            >> 32 - (uint32_t) (*(bwtintv_t *) arg1).info + (uint32_t) (*(bwtintv_t *) arg1).info >> 32;
}

int var2num(int refseq, int queryseq){
    if (refseq < 4)
        return refseq + queryseq -1;
    else
        return queryseq;
}

typedef struct {
    int block; //所在block
    bwtint_t pos; //与block起始位置的距离
    int cigar_l[10], cigar_v[10];
} align_info;

int main(int argc, char *argv[])
{
    int c, i, base, min_iwidth = 20, min_len = 17, print_seq = 0, min_intv = 1, max_len = INT_MAX, lgst_num = 3, max_mis = 5, qual_sys = 2;
    int max_insr = 500;
    int se_mark = 1;
    int cigar_l[max_mis], cigar_v[max_mis];
    uint64_t max_intv = 0;
    kseq_t *seq1, *seq2;
    bwtint_t k;
    gzFile fp1, fp2;
    smem_i *itr;
    const bwtintv_v *a;
    bwaidx_t *idx;

    while ((c = getopt(argc, argv, "w:l:i:I:L:f:m:s:q:")) >= 0) {
        switch (c) {
            case 'w': min_iwidth = atoi(optarg); break;
            case 'l': min_len = atoi(optarg); break;
            case 'i': min_intv = atoi(optarg); break;
            case 'I': max_intv = atol(optarg); break;
            case 'L': max_len  = atoi(optarg); break;
            case 'f': lgst_num = atoi(optarg); break;
            case 'm': max_mis = atoi(optarg); break;
            case 's': max_insr = atoi(optarg); break;
            case 'q': qual_sys = atoi(optarg); break;
            default: return 1;
        }
    }

    if (optind + 1 >= argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   fastsam [options] <idxbase> <in_1.fq> <in_2.fq>\n\n");
        fprintf(stderr, "Options: -l INT    min SMEM length to output [%d]\n", min_len);
        fprintf(stderr, "         -w INT    max interval size to find coordiantes [%d]\n", min_iwidth);
        fprintf(stderr, "         -i INT    min SMEM interval size [%d]\n", min_intv);
        fprintf(stderr, "         -L INT    max MEM length [%d]\n", max_len);
        fprintf(stderr, "         -I INT    stop if MEM is longer than -l with a size less than INT [%ld]\n", (long)max_intv);
        fprintf(stderr, "         -f INT    consider only the longest [%d] MEM\n", lgst_num);
        fprintf(stderr, "         -m INT    max mismatch to tolerate [%d]\n", max_mis);
        fprintf(stderr, "         -q INT    quality system, 1:illumina, 2:sanger, default as [%d]\n", qual_sys);
        fprintf(stderr, "         -s INT    max insert size between read1 and read2, ");
        fprintf(stderr, "\n");
        return 1;
    }

    fp1 = xzopen(argv[optind + 1], "r");
    seq1 = kseq_init(fp1);
    if (argc - optind >= 3){
        se_mark = 0;
        fp2 = xzopen(argv[optind + 2], "r");
        seq2 = kseq_init(fp2);
    }
    if ((idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1;
    itr = smem_itr_init(idx->bwt);
    smem_config(itr, min_intv, max_len, max_intv);
    align_info alter_smem[3];
    //Seq1
    while (kseq_read(seq1) >= 0) {
        int pass_num = 0;

        size_t rlen = seq1->seq.l; //转存
        for (i = 0; i < seq1->seq.l; ++i) {
            seq1->seq.s[i] = nst_nt4_table[(int) seq1->seq.s[i]];
        }
        smem_set_query(itr, seq1->seq.l, (uint8_t *) seq1->seq.s);
        while ((a = smem_next(itr)) != 0) {   //这里表示每个smem
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
            qsort(plist, a->n - short_num, sizeof(bwtintv_t *), bwtintv_cmp);
            for (i = 0; i < a->n - short_num; ++i) {
                if (plist[i]->x[2] <= min_iwidth) {
                    for (k = 0; k < plist[i]->x[2]; ++k) {
                        bwtint_t pos;
                        int len, is_rev, ref_id;
                        len = (uint32_t) plist[i]->info - (uint32_t) (plist[i]->info >> 32);
                        pos = bns_depos(idx->bns, bwt_sa(idx->bwt, plist[i]->x[0] + k), &is_rev);
                        uint8_t *rseq, *rseq_l, *rseq_r;
                        uint16_t missum = 0;
                        if (is_rev) {
                            pos -= len - 1;
                            rseq_r = bns_get_seq(idx->bns->l_pac, idx->pac, pos + len,
                                                 pos + (uint32_t) (plist[i]->info), &rlen);
                            //把rseq_r和read的左截的反向互补进行比较
                            for (base = 0; base < rlen; base++) {
                                if (rseq_r[base] + seq1->seq.s[base] != 3) {
                                    if (missum >= max_mis) break;
                                    cigar_l[missum] = base;
                                    cigar_v[missum] = var2num(rseq_r[base], seq1->seq.s[base]);
                                    missum += 1;
                                }
                            }
                            if (missum <= max_mis) {
                                rseq_l = bns_get_seq(idx->bns->l_pac, idx->pac,
                                                     pos - seq1->seq.l + (uint32_t) (plist[i]->info), pos, &rlen);
                                //把rseq_l和read的右截的反向互补进行比较
                                for (base = 0; base < rlen; base++) {
                                    if (rseq_l[base] + seq1->seq.s[seq1->seq.l - 1 - base] != 3) {
                                        if (missum >= max_mis) break;
                                        cigar_l[missum] = (uint32_t) plist[i]->info - 1 + base;
                                        cigar_v[missum] = var2num(rseq_l[base], seq1->seq.s[seq1->seq.l - 1 - base]);
                                        missum += 1;
                                    }
                                }
                            }
                            if (missum <= max_mis)
                                //赋值;
                        } else {
                            rseq_l = bns_get_seq(idx->bns->l_pac, idx->pac, pos - (uint32_t) (plist[i]->info >> 32),
                                                 pos, &rlen);
                            //把rseq_l和read的左截进行比较
                            for (base = 0; base < rlen; base++) {
                                if (rseq_l[base] != seq1->seq.s[base]) {
                                    if (missum >= max_mis) break;
                                    cigar_l[missum] = base;
                                    cigar_v[missum] = var2num(rseq_l[base], seq1->seq.s[base]);
                                    missum += 1;
                                }
                            }
                            if (missum <= max_mis) {
                                rseq_r = bns_get_seq(idx->bns->l_pac, idx->pac, pos + len,
                                                     pos + seq1->seq.l - (uint32_t) (plist[i]->info >> 32), &rlen);
                                //把rseq_r和read的右截进行比较
                                for (base = 0; base < rlen; base++) {
                                    if (rseq_r[base] != seq1->seq.s[(uint32_t) plist[i]->info + base]) {
                                        if (missum >= max_mis) break;
                                        cigar_l[missum] = base;
                                        cigar_v[missum] = var2num(rseq_r[base], seq1->seq.s[base]);
                                        missum += 1;
                                    }
                                }
                            }
                            if (missum <= max_mis)
                                bns_cnt_ambi(idx->bns, pos, len, &ref_id);
                        }
                        if (ref_id) {
                            //return idx->bns->anns[ref_id].name, pos, is_rev, cigar_l, cigar_v;
                        }
                        //seq->qual.s seq->name.s 另作处理;
                        err_printf("\t%s:%c-%d\t%ld", idx->bns->anns[ref_id].name, "+-"[is_rev], pos,
                                   (long) (pos - idx->bns->anns[ref_id].offset) + 1);
                    }
                } else err_puts("\t*");
                err_putchar('\n');
            }
        }
        err_puts("//");
    }
    smem_itr_destroy(itr);
    bwa_idx_destroy(idx);
    kseq_destroy(seq1);
    err_gzclose(fp1);

    //Seq2
    while (kseq_read(seq2) >= 0) {
        int pass_num = 0;

        size_t rlen = seq2->seq.l; //转存
        for (i = 0; i < seq2->seq.l; ++i) {
            seq2->seq.s[i] = nst_nt4_table[(int) seq2->seq.s[i]];
        }
        smem_set_query(itr, seq2->seq.l, (uint8_t *) seq2->seq.s);
        while ((a = smem_next(itr)) != 0) {   //这里表示每个smem
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
            qsort(plist, a->n - short_num, sizeof(bwtintv_t *), bwtintv_cmp);
            for (i = 0; i < a->n - short_num; ++i) {
                if (plist[i]->x[2] <= min_iwidth) {
                    for (k = 0; k < plist[i]->x[2]; ++k) {
                        bwtint_t pos;
                        int len, is_rev, ref_id;
                        len = (uint32_t) plist[i]->info - (uint32_t) (plist[i]->info >> 32);
                        pos = bns_depos(idx->bns, bwt_sa(idx->bwt, plist[i]->x[0] + k), &is_rev);
                        uint8_t *rseq, *rseq_l, *rseq_r;
                        uint16_t missum = 0;
                        if (is_rev) {
                            pos -= len - 1;
                            rseq_r = bns_get_seq(idx->bns->l_pac, idx->pac, pos + len,
                                                 pos + (uint32_t) (plist[i]->info), &rlen);
                            //把rseq_r和read的左截的反向互补进行比较
                            for (base = 0; base < rlen; base++) {
                                if (rseq_r[base] + seq2->seq.s[base] != 3) {
                                    if (missum >= max_mis) break;
                                    cigar_l[missum] = base;
                                    cigar_v[missum] = var2num(rseq_r[base], seq2->seq.s[base]);
                                    missum += 1;
                                }
                            }
                            if (missum <= max_mis) {
                                rseq_l = bns_get_seq(idx->bns->l_pac, idx->pac,
                                                     pos - seq2->seq.l + (uint32_t) (plist[i]->info), pos, &rlen);
                                //把rseq_l和read的右截的反向互补进行比较
                                for (base = 0; base < rlen; base++) {
                                    if (rseq_l[base] + seq2->seq.s[seq2->seq.l - 1 - base] != 3) {
                                        if (missum >= max_mis) break;
                                        cigar_l[missum] = (uint32_t) plist[i]->info - 1 + base;
                                        cigar_v[missum] = var2num(rseq_l[base], seq2->seq.s[seq2->seq.l - 1 - base]);
                                        missum += 1;
                                    }
                                }
                            }
                            if (missum <= max_mis)
                            //赋值;
                        } else {
                            rseq_l = bns_get_seq(idx->bns->l_pac, idx->pac, pos - (uint32_t) (plist[i]->info >> 32),
                                                 pos, &rlen);
                            //把rseq_l和read的左截进行比较
                            for (base = 0; base < rlen; base++) {
                                if (rseq_l[base] != seq2->seq.s[base]) {
                                    if (missum >= max_mis) break;
                                    cigar_l[missum] = base;
                                    cigar_v[missum] = var2num(rseq_l[base], seq2->seq.s[base]);
                                    missum += 1;
                                }
                            }
                            if (missum <= max_mis) {
                                rseq_r = bns_get_seq(idx->bns->l_pac, idx->pac, pos + len,
                                                     pos + seq2->seq.l - (uint32_t) (plist[i]->info >> 32), &rlen);
                                //把rseq_r和read的右截进行比较
                                for (base = 0; base < rlen; base++) {
                                    if (rseq_r[base] != seq2->seq.s[(uint32_t) plist[i]->info + base]) {
                                        if (missum >= max_mis) break;
                                        cigar_l[missum] = base;
                                        cigar_v[missum] = var2num(rseq_r[base], seq2->seq.s[base]);
                                        missum += 1;
                                    }
                                }
                            }
                            if (missum <= max_mis)
                                bns_cnt_ambi(idx->bns, pos, len, &ref_id);
                        }
                        if (ref_id) {
                            //return idx->bns->anns[ref_id].name, pos, is_rev, cigar_l, cigar_v;
                        }
                        //seq->qual.s seq->name.s 另作处理;
                        err_printf("\t%s:%c-%d\t%ld", idx->bns->anns[ref_id].name, "+-"[is_rev], pos,
                                   (long) (pos - idx->bns->anns[ref_id].offset) + 1);
                    }
                } else err_puts("\t*");
                err_putchar('\n');
            }
        }
        err_puts("//");
    }
    smem_itr_destroy(itr);
    bwa_idx_destroy(idx);
    kseq_destroy(seq2);
    err_gzclose(fp1);
    return 0;
}
