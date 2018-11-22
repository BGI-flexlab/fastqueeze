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
//#define SEQARCTHREAD
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
sem_t g_sem;	//队列信号量

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

void CreatBitmap(map<int, map<int, int>> &bitmap) {
    map<int, int> map_G;
    map<int, int> map_C;
    map<int, int> map_A;
    map<int, int> map_T;
    map_G[1] = 0; map_G[0] = 1; map_G[3] = 0;
    map_C[2] = 0; map_C[0] = 1; map_C[3] = 2;
    map_A[2] = 0; map_A[1] = 1; map_A[3] = 2;
    map_T[2] = 0; map_T[1] = 1; map_T[0] = 2;
    bitmap[2] = map_G;
    bitmap[1] = map_C;
    bitmap[0] = map_A;
    bitmap[3] = map_T;
}

int getAlignInfo(kseq seq, smem_i* func_itr, bwaidx_t *func_idx, align_info *align_p, int func_block_size, int min_len, int max_iwidth, int max_mis, int lgst_num){
    int64_t rlen;
    int i, seql, base;

    const bwtintv_v *a;
    seql = (int) seq.seq.length();
    int pass_num = 0;
    int cigar_l[max_mis], cigar_v[max_mis];

    static map<int, map<int, int>> nucleBitMap;
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
                                cigar_v[missum] = nucleBitMap[rbase][3-sbase];
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
                                    cigar_v[missum] = nucleBitMap[rbase][3-sbase];
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
                                cigar_v[missum] = nucleBitMap[rbase][sbase];
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
                                    cigar_v[missum] = nucleBitMap[rbase][sbase];
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
	int min_len;
	int max_mis;
	int max_iwidth;
	int lgst_num;
	int block_num;
    int max_insr;
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

	while (!g_bFinish || !g_task_queue.empty())
	{
        //printf("thread id= %0x waitting\n", pthread_self());
        sem_wait(&g_sem);
        //printf("thread id= %0x running\n", pthread_self());
		if (!g_task_queue.empty()) //有任务，开始执行
		{
            pthread_mutex_lock(&g_mutex);
			Task task = g_task_queue.front();
            g_task_queue.pop();
            pthread_mutex_unlock(&g_mutex);

            int num = task.num;
            if(num == 1) //SE
            {
                kseq &seq = task.seq1;
                align_info *sam1 = new align_info[pParam->lgst_num];
                for (int i = 0; i < pParam->lgst_num; i++) {
                    sam1[i].cigar_l = (int*)malloc(pParam->max_mis * sizeof(int));
                    memset(sam1[i].cigar_l, -1, pParam->max_mis * sizeof(int));
                    sam1[i].cigar_v = (int*)malloc(pParam->max_mis * sizeof(int));
                }

                
                pthread_mutex_lock(&g_write_mutex);
                int seq1m = getAlignInfo(seq, pParam->pitr, pParam->pidx, sam1, pParam->block_size,
                    pParam->min_len, pParam->max_iwidth, pParam->max_mis, pParam->lgst_num);
                
                string seq1Name = seq.name + " " + seq.comment;
                
                if (seq1m)
                {
                    std::sort(sam1, sam1 + seq1m, sam_cmp);
                    pParam->pencoders[sam1[0].blockNum]->parse_1(sam1[0], seq.seq.length(), pParam->pfpOutput_s[sam1[0].blockNum]);
                    pParam->pfqz[sam1[0].blockNum]->iq_encode(seq1Name, seq.qual, pParam->pfpOutput_iq[sam1[0].blockNum]);
                }
                else 
                {
                    pParam->pfqz[pParam->block_num]->isq_encode(seq1Name, seq.seq, seq.qual, *(pParam->pout_isq)); 
                }
                pthread_mutex_unlock(&g_write_mutex);

                for (int i = 0; i < pParam->lgst_num; i++)
                {
                    free(sam1[i].cigar_l);
                    free(sam1[i].cigar_v);
                }
                delete[]sam1;
            }
            else if(num == 2) //PE
            {
                kseq &seq1 = task.seq1;
                kseq &seq2 = task.seq1;
                align_info *sam1 = new align_info[pParam->lgst_num];
                align_info *sam2 = new align_info[pParam->lgst_num];
                for (int i = 0; i < pParam->lgst_num; i++) {
                    sam1[i].cigar_l = (int*)malloc(pParam->max_mis * sizeof(int));
                    memset(sam1[i].cigar_l, -1, pParam->max_mis * sizeof(int));
                    sam1[i].cigar_v = (int*)malloc(pParam->max_mis * sizeof(int));

                    sam2[i].cigar_l = (int*)malloc(pParam->max_mis * sizeof(int));
                    memset(sam2[i].cigar_l, -1, pParam->max_mis * sizeof(int));
                    sam2[i].cigar_v = (int*)malloc(pParam->max_mis * sizeof(int));
                }

                pthread_mutex_lock(&g_write_mutex);
                int seq1m = getAlignInfo(seq1, pParam->pitr, pParam->pidx, sam1, pParam->block_size,
                    pParam->min_len, pParam->max_iwidth, pParam->max_mis, pParam->lgst_num);

                int seq2m = getAlignInfo(seq2, pParam->pitr, pParam->pidx, sam2, pParam->block_size,
                    pParam->min_len, pParam->max_iwidth, pParam->max_mis, pParam->lgst_num);

                string seq1Name = seq1.name + " " + seq1.comment;
                string seq2Name = seq2.name + " " + seq2.comment;
                bool bfind = false;

                for(int i=0;i<seq1m;i++)
                {
                    for(int j=0;j<seq2m;j++)
                    {
                        if(sam1[i].blockNum == sam2[j].blockNum &&
                           abs(sam1[i].blockPos - sam2[j].blockPos) <= pParam->max_insr) //满足比对条件
                        {
                            bfind = true;
                            pParam->pencoders[sam1[i].blockNum]->parse_1(sam1[i], seq1.seq.length(), pParam->pfpOutput_s[sam1[i].blockNum]);
                            pParam->pencoders[sam2[j].blockNum]->parse_2(sam2[j], seq2.seq.length(), pParam->pfpOutput_s[sam2[j].blockNum]);
                            pParam->pfqz[sam1[i].blockNum]->iq_encode(seq1Name, seq1.qual, pParam->pfpOutput_iq[sam1[i].blockNum]);
                            pParam->pfqz[sam2[j].blockNum]->iq_encode(seq2Name, seq2.qual, pParam->pfpOutput_iq[sam2[j].blockNum]);

                            i = seq1m; //强制结束外层循环
                            break;
                        }
                    }
                }

                if(!bfind) //没有比对上
                {
                    pParam->pfqz[pParam->block_num]->isq_encode(seq1Name, seq1.seq, seq1.qual, *(pParam->pout_isq));
                    pParam->pfqz[pParam->block_num]->isq_encode(seq2Name, seq2.seq, seq2.qual, *(pParam->pout_isq));
                }
                pthread_mutex_unlock(&g_write_mutex);

                for (int i = 0; i < pParam->lgst_num; i++)
                {
                    free(sam1[i].cigar_l);
                    free(sam1[i].cigar_v);

                    free(sam2[i].cigar_l);
                    free(sam2[i].cigar_v);
                }
                delete[]sam1;
                delete[]sam2;
            }
		}
	}
    //printf("thread id= %0x close\n", pthread_self());
}


void decode_block_num(std::fstream &in, std::vector<int> &vec_len)
{
    unsigned char len_buf[4];
    while (4 == xget(in, len_buf, 4))
    {
        int comp_len =
                (len_buf[0] <<  0) +
                (len_buf[1] <<  8) +
                (len_buf[2] << 16) +
                (len_buf[3] << 24);


        vec_len.push_back(comp_len);
        in.seekg(comp_len, std::ios::cur);
    }
    in.clear();
    in.seekg(0, std::ios::beg);
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
    fprintf(fp, "    -n <level>     Name compression level.  1-2 [1]\n");
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

    while ((opt = getopt(argc, argv, "l:w:I:f:m:q:s:hdQ:S:N:bePXiB:t:n:")) != -1) {
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
        return 1;
    }
    else if (decompress) {
        fstream fa, in_s, in_iq, in_isq, out1, out2;
        fa.open(argv[optind], std::ios::in);

        stringstream strInputPath;
        strInputPath.str("");
        strInputPath << "./" << argv[optind+1] << "_aligninfo.arc";
        in_s.open(strInputPath.str(), std::ios::binary|std::ios::in);

        strInputPath.str("");
        strInputPath << "./" << argv[optind+1] << "_iq.arc";
        in_iq.open(strInputPath.str(), std::ios::binary|std::ios::in);

        strInputPath.str("");
        strInputPath << "./" << argv[optind+1] << "_isq.arc";
        in_isq.open(strInputPath.str(), std::ios::binary|std::ios::in);

        unsigned char magic_s[12];
        if (12 != xget(in_s, magic_s, 12)){
            fprintf(stderr, "Abort: truncated read.\n");
            return 1;
        }
        qual_sys = magic_s[0];
        se_mark = magic_s[1];
        memcpy(&block_size, magic_s+2, 4);
        block_num = magic_s[6];
        max_mis = magic_s[7];
        memcpy(&max_insr, magic_s+8, 2);
        memcpy(&max_readLen, magic_s+10, 2);

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

        strInputPath.str("");
        strInputPath << "./" << argv[optind+2] << "_1.fastq";
        out1.open(strInputPath.str(), std::ios::out);
        if (!se_mark){
            strInputPath.str("");
            strInputPath << "./" << argv[optind+2] << "_2.fastq";
            out2.open(strInputPath.str(), std::ios::out);
        }

        std::vector<string> vec_name;
        std::vector<string> vec_qual;
        std::vector<string> vec_seq;
        string seq;

        align_info align_info1;
        align_info1.cigar_l = (int*) malloc(max_mis*sizeof(int));
        align_info1.cigar_v = (int*) malloc(max_mis*sizeof(int));
        int readLen, thisblock=0;
        bool blockjump=false;

//        for(i=0;i<block_num;i++)
//        {
//            char buf[100]={0};
//            sprintf(buf,"./iq_%d.tmp",i);
//            fstream ftmp;
//            ftmp.open(buf, std::ios::binary|std::ios::in);
//            fqz *pf = new fqz(&p);
//            std::vector<int> vec_len;
//            decode_block_num(ftmp, vec_len);
//            auto itor = vec_len.begin();
//            for(;itor !=vec_len.end();itor++)
//            {
//                name.clear();
//                qual.clear();
//                pf->iq_decode(ftmp, *itor, name, qual);
//                //printf("%d %d %d %d \n", i, vec_len.size(), name.size(), qual.size());
//                for(int j=0;j<qual.size();j++)
//                {
//                    //printf("%s\n", name[j].c_str());
//                    while (!seqdecoder.parse_se(align_info1, readLen, in_s));
//                    //printf("%d,%d,%d,%d,%d,%d\n", align_info1.cigar_l[0], align_info1.cigar_l[1], align_info1.cigar_l[2], align_info1.cigar_v[0], align_info1.cigar_v[1], align_info1.cigar_v[2]);
//                    seq = ref2seqer.getSeq(align_info1, readLen);
//                    //printf("%s\n", seq.c_str());
//                    //printf("%s\n", qual[j].c_str());
//                    readModify2(seq, qual[j], qual_sys);
//                    out1 << name[j] << endl << seq << endl << "+" << endl  << qual[j] << endl;
//                }
//            }
//            delete pf;
//        }
//
//        std::vector<string> vec_seq;
//        fqz *pf = new fqz(&p);
//        std::vector<int> vec_len;
//        decode_block_num(in_isq, vec_len);
//        auto itor = vec_len.begin();
//        for(;itor !=vec_len.end();itor++)
//        {
//            name.clear();
//            qual.clear();
//            vec_seq.clear();
//            pf->isq_decode(in_isq, *itor, name, vec_seq, qual);
//
//            //printf("%d %d %d %d\n", vec_len.size(), name.size(), vec_seq.size(), qual.size());
//            for(int j=0;j<qual.size();j++)
//            {
//                //printf("%s\n", name[j].c_str());
//                //printf("%s\n", vec_seq[j].c_str());
//                //printf("%s\n", qual[j].c_str());
//                out1 << name[j] << endl << vec_seq[j] << endl << "+" << endl  << qual[j] << endl;
//            }
//        }
//        delete pf;
//
//        out1.close();

        fqz *f = new fqz(&p);
        while (thisblock <= block_num){
            if (se_mark){
                if (seqdecoder.parse_se(align_info1, readLen, in_s) == 1){ //request a read from vec
                    if (vec_name.size() == 0) //if no read left in vec, call the iq_decode then
                        f->iq_decode(in_iq, vec_name, vec_qual);
                    seq = ref2seqer.getSeq(align_info1, readLen);
                    readModify2(seq, vec_qual[0], qual_sys);
                    out1 << vec_name[0] << endl << seq << endl << "+" << endl << vec_qual[0] << endl;
                    vec_name.erase(vec_name.begin());
                    vec_qual.erase(vec_qual.begin());
                    blockjump = false;
                }
                else{
                    thisblock ++;
                    if (blockjump) //the second continuous 0, means the block contains no read
                        in_iq.seekg(52, std::ios::cur);
                    else{
                        f = new fqz(&p);
                        if (vec_name.size() > 0) //it should be empty if align_info and iq paired
                            return -1;
                        blockjump = true;
                    }
                }
            }
            else{
                if (seqdecoder.parse_pe(align_info1, readLen, in_s) == 0){ //empty block
                    thisblock ++;
                    in_iq.seekg(52, std::ios::cur);
                    continue;
                }
                else{
                    if (vec_name.size() == 0) //if no read left in vec, call the iq_decode then
                        f->iq_decode(in_iq, vec_name, vec_qual);
                    seq = ref2seqer.getSeq(align_info1, readLen);
                    readModify2(seq, vec_qual[0], qual_sys);
                    out1 << vec_name[0] << endl << seq << endl << "+" << endl << vec_qual[0] << endl;
                    vec_name.erase(vec_name.begin());
                    vec_qual.erase(vec_qual.begin());
                }

                if (seqdecoder.parse_pe(align_info1, readLen, in_s) == 1){ //request a read from vec
                    if (vec_name.size() == 0) //if no read left in vec, call the iq_decode then
                        f->iq_decode(in_iq, vec_name, vec_qual);
                    seq = ref2seqer.getSeq(align_info1, readLen);
                    readModify2(seq, vec_qual[0], qual_sys);
                    out2 << vec_name[0] << endl << seq << endl << "+" << endl << vec_qual[0] << endl;
                    vec_name.erase(vec_name.begin());
                    vec_qual.erase(vec_qual.begin());
                }
                else{
                    thisblock ++;
                    f = new fqz(&p);
                    if (vec_name.size() > 0) //it should be empty if align_info and iq paired
                        return -1;
                }
            }
        }
        fa.close();
        in_s.close();
        in_iq.close();
        vec_name.clear();
        vec_qual.clear();
        vec_seq.clear();

        while (true){
            if (se_mark){
                if (vec_name.size() == 0){
                    if (-1 == f->isq_decode(in_isq, vec_name, vec_seq, vec_qual))
                        break;
                }
                readModify2(vec_seq[0], vec_qual[0], qual_sys);
                out1 << vec_name[0] << endl << vec_seq[0] << endl << "+" << endl << vec_qual[0] << endl;
                vec_name.erase(vec_name.begin());
                vec_seq.erase(vec_seq.begin());
                vec_qual.erase(vec_qual.begin());
            }
            else{
                if (vec_name.size() == 0){ // it should arise only at the end of the file
                    if (-1 == f->isq_decode(in_isq, vec_name, vec_seq, vec_qual))
                        break;
                    else
                        return -1;
                }
                readModify2(vec_seq[0], vec_qual[0], qual_sys);
                out1 << vec_name[0] << endl << vec_seq[0] << endl << "+" << endl << vec_qual[0] << endl;
                vec_name.erase(vec_name.begin());
                vec_seq.erase(vec_seq.begin());
                vec_qual.erase(vec_qual.begin());

                if (vec_name.size() == 0){
                    if (-1 == f->isq_decode(in_isq, vec_name, vec_seq, vec_qual))
                        break;
                }
                readModify2(vec_seq[0], vec_qual[0], qual_sys);
                out2 << vec_name[0] << endl << vec_seq[0] << endl << "+" << endl << vec_qual[0] << endl;
                vec_name.erase(vec_name.begin());
                vec_seq.erase(vec_seq.begin());
                vec_qual.erase(vec_qual.begin());
            }
        }
        delete f;
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

        unsigned char magic_s[12]{(uint8_t)qual_sys, //block_size+max_mis+max_insr+max_readLen
                                  (uint8_t)se_mark,
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

        if (12 != xwrite(out_s, magic_s, 12)) {
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

#ifndef SEQARCTHREAD
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
#endif

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
			pParam->max_mis = max_mis;
			pParam->min_len = min_len;
			pParam->max_iwidth = max_iwidth;
			pParam->lgst_num = lgst_num;
            pParam->block_num = block_num;
            pParam->max_insr = max_insr;
			pParam->pencoders = encoders;
			pParam->pfpOutput_iq = fpOutput_iq;
			pParam->pfpOutput_s = fpOutput_s;
			pParam->pidx = idx;
			pParam->pitr = itr;
			//pParam->pout_iq = &out_iq;
			//pParam->pout_s = &out_s;
			pParam->pout_isq = &out_isq;
			pParam->pfqz = f;
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

        for (i = 0; i < thread_num -1; ++i)
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
            			        break;
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
                    memset(sam1[i].cigar_l,-1,max_mis*sizeof(int));
                    memset(sam2[i].cigar_l,-1,max_mis*sizeof(int));
                }
            }
        }
#endif
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