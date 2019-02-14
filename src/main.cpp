#include "SeqUtil.hpp"

/* -------------------------------------------------------------------------
 * Main program
 */

static void usage(int err) {
    FILE *fp = err ? stderr : stdout;

    fprintf(fp, "SeqArc v%d.%d. Yuxin Chen, Zijian Zhao, 2019\n", MAJOR_VERS, MINOR_VERS);
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
    fprintf(fp, "    -s INT         max insert size between read1 and read2 [511]\n");
    fprintf(fp, "    -r INT         files with AlignRatio lower than [0.5] are processed with Fqzcomp only.\n\n");

    fprintf(fp, "    -S <level>     Sequence de novo compression level. 1-9 [3]\n");
    fprintf(fp, "                   Specifying '+' on the end (eg -s5+) will use\n");
    fprintf(fp, "                   models of multiple sizes for improved compression.\n");
    fprintf(fp, "    -N <level>     Quality compression level.  1-3 [2]\n");
    fprintf(fp, "    -n <level>     Name compression level.  1-2 [1]\n");
    fprintf(fp, "    -b             Use both strands in sequence hash table.\n");
    fprintf(fp, "    -e             Extra seq compression: 16-bit vs 8-bit counters.\n");
    fprintf(fp, "    -t INT         Thread num for multi-threading, default as [1]\n\n");

    fprintf(fp, "    -X             Enable generation/verification of check sums\n\n");
    fprintf(fp, "    -W             Show warning msg about abnormal base\n\n");
    fprintf(fp, "To decompress:\n   SeqArc -d [ref.fa] <compress_prefix> <fastq_prefix>\n");

    exit(err);
}


int main(int argc, char **argv) {
    int opt, i, max_len = INT_MAX;
    uint64_t max_intv = 0;
    smem_i *itr = NULL;
    bwaidx_t *idx = NULL;
    int decompress = 0, indexing = 0;


    /* Initialise and parse command line arguments */
    g_fqz_params.slevel = 3;
    g_fqz_params.qlevel = 2;
    g_fqz_params.nlevel = 1;
    g_fqz_params.both_strands = 0;
    g_fqz_params.extreme_seq = 0;
    g_fqz_params.multi_seq_model = 0;
    g_fqz_params.do_hash = 0;

    while ((opt = getopt(argc, argv, "l:w:I:f:m:q:s:hdQ:S:N:beXiB:t:n:c:E:r:W")) != -1) {
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

            case 's':
                max_insr = atoi(optarg);
                break;

            case 'r':
                min_alignratio = (float) atof(optarg);
                break;

            case 'S':
                char *end;
                g_fqz_params.slevel = strtol(optarg, &end, 10);
                if (g_fqz_params.slevel < 1 || g_fqz_params.slevel > 9)
                    usage(1);
                if (*end == '+')
                    g_fqz_params.multi_seq_model = 1;
                break;

            case 'N':
                g_fqz_params.qlevel = atoi(optarg);
                if (g_fqz_params.qlevel < 1 || g_fqz_params.qlevel > 3)
                    usage(1);
                break;

            case 'n':
                g_fqz_params.nlevel = atoi(optarg);
                if (g_fqz_params.nlevel < 1 || g_fqz_params.nlevel > 2)
                    usage(1);
                break;

            case 'b':
                g_fqz_params.both_strands = 1;
                break;

            case 'e':
                g_fqz_params.extreme_seq = 1;
                break;

            // case 'X':
            //     g_fqz_params.do_hash = 1;
            //     break;

    	    case 't':
		        thread_num = atoi(optarg);
                if(thread_num>MAX_THREAD_NUM)
                {
                    printf("Thread num can not be more than %d\n", MAX_THREAD_NUM);
                    thread_num = MAX_THREAD_NUM;
                }
		        break;

            case 'W':
                g_show_warning = true;
                break;
            default:
                usage(1);
        }
    }

    CreatBitmap(nucleBitMap);

    if (argc == optind)
        usage(1);

    if (indexing){
        char *ref = argv[optind];
        bwa_idx_build(ref, ref, BWTALGO_AUTO, 10000000);
        bwaidx_t *idx = bwa_idx_load_from_disk(ref, BWA_IDX_ALL);
        bwa_shm_stage(idx, ref, NULL);
        return 1;
    }
    else if (decompress) {
        struct timeval timeStart,timeEnd;
        gettimeofday(&timeStart, NULL);
        string tmpArg;
        bool hasFA;
        char path[256]={0};
        string fastq_prefix;

        tmpArg = argv[optind];
        if (tmpArg.length() <= 5 ? tmpArg.substr(tmpArg.size()-2) != "fa" : tmpArg.substr(tmpArg.size()-2) != "fa" && tmpArg.substr(tmpArg.size()-5) != "fasta"){ //no fasta
            hasFA = false;
            sprintf(path, "./%s.arc", optind > argc-1? "out":argv[optind]);
            fastq_prefix = (optind+1 > argc-1? "decode":argv[optind+1]);
        }
        else{
            hasFA = true;
            idx = bwa_idx_load_from_shm(argv[optind]);
            if (idx == 0) {
                if ((idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1;
                bwa_shm_stage(idx, argv[optind], NULL);
            }

            uint64_t res_len = idx->bns->l_pac; //获取参考序列的长度，小于文件的真实长度
            g_offset_bit = getbitnum(res_len)-2;

            sprintf(path, "./%s.arc", optind+1 > argc-1? "out":argv[optind+1]);
            fastq_prefix = (optind+2 > argc-1? "decode":argv[optind+2]);
        }

        fstream in_s;
        in_s.open(path, std::ios::binary|std::ios::in);
        in_s.seekg(-1024, ios::end);
        char buf[1024]={0};
        in_s.read(buf,1024);
        in_s.close();
        char *p = buf+1020;
        if(memcmp(p,".arc",4) != 0)
        {
            printf("Unrecognised file format.\n");
            return 1;
        }
        p -= sizeof(g_magicparam);
        memcpy(&g_magicparam, p, sizeof(g_magicparam));

        if(g_magicparam.major_vers != MAJOR_VERS || 
            g_magicparam.format_vers != FORMAT_VERS)
        {
            printf("Unsupported file format version %d.%d\n", g_magicparam.major_vers, g_magicparam.format_vers);
            return 1;
        }
        if (!g_magicparam.fqzall && !hasFA){
            printf("FASTA index needed.\n");
            return 1;
        }

        g_fqz_params.both_strands = g_magicparam.both_strands;
        g_fqz_params.extreme_seq = g_magicparam.extreme_seq;
        g_fqz_params.multi_seq_model = g_magicparam.multi_seq_model;
        g_fqz_params.fqzall = g_magicparam.fqzall;
        bool isSE = g_magicparam.isSE;
        g_fqz_params.slevel = g_magicparam.slevel;
        g_fqz_params.qlevel = g_magicparam.qlevel;
        g_fqz_params.nlevel = g_magicparam.nlevel;
        max_mis = g_magicparam.max_mis;
        max_insr = g_magicparam.max_insr;
        max_readLen = g_magicparam.max_readLen;
        thread_num = g_magicparam.thread_num;

        p -= thread_num*5;
        uint64_t slicearry[MAX_THREAD_NUM]={0};
        for(int i=0;i<=thread_num;i++)
        {
            unsigned char len_buf[5]={0};
            memcpy(len_buf, p, 5); p+=5; 
            slicearry[i] =
                (len_buf[0] <<  0) +
                (len_buf[1] <<  8) +
                (len_buf[2] << 16) +
                (len_buf[3] << 24) +
                (len_buf[4] << 32);
        }

        pthread_t *t_id = (pthread_t*)alloca(thread_num * sizeof(pthread_t));

        int i=0;
        for (i = 0; i < thread_num; ++i) //创建一个线程池
        {
            DecodeParam * param = new DecodeParam;
            param->num = i;
            param->pidx = idx;

            param->length = slicearry[i];
            param->offset = Getoffset(slicearry,i); 
            strcpy(param->filename, path);

            if(g_magicparam.fqzall)
            {
                pthread_create(&t_id[i], 0, fqzall_decode_process, param);
            }
            else
            {
                pthread_create(&t_id[i], 0, decode_process, param);
            }
        }

        for (i = 0; i < thread_num; ++i) 
        {
            pthread_join(t_id[i], 0);
        }

        if(g_magicparam.isGzip)
        {
            MergeFileForzip(fastq_prefix, thread_num, 1);
            if(!g_magicparam.isSE)
            {
                MergeFileForzip(fastq_prefix, thread_num, 2);
            }
        }
        else
        {
            MergeFileForFastq(fastq_prefix, thread_num, 1);
            if(!g_magicparam.isSE)
            {
                MergeFileForFastq(fastq_prefix, thread_num, 2);
            }
        }

        gettimeofday(&timeEnd, NULL);
        double totald = (timeEnd.tv_sec - timeStart.tv_sec)*1000 + (timeEnd.tv_usec -timeStart.tv_usec)*1.0/1000;
        printf("%s%d decode total time %f ms\n", __FUNCTION__, __LINE__, totald);

        return 0;

    } else {
        struct timeval timeStart,timeEnd;
        gettimeofday(&timeStart, NULL);
        string tmpArg, outIndex;

        tmpArg = argv[optind];
        if (tmpArg.substr(tmpArg.size()-2) != "fa" && tmpArg.substr(tmpArg.size()-5) != "fasta"){ //no fasta
            cout << "FASTA index missing." << endl;
            g_magicparam.fqzall = true;
        }
        else{
            g_magicparam.fqzall = false;
            idx = bwa_idx_load_from_shm(argv[optind]);
            if (idx == 0) {
                if ((idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1;
                bwa_shm_stage(idx, argv[optind], NULL);
            }
            itr = smem_itr_init(idx->bwt);
            smem_config(itr, 1, max_len, max_intv); //min_intv = 1
            ++optind;
            uint64_t res_len = idx->bns->l_pac; //获取参考序列的长度，小于文件的真实长度
            g_offset_bit = getbitnum(res_len)-2;
            g_offset_size = pow(2,g_offset_bit)-1;
        }

        g_isgzip = GetFileType(argv[optind-1]);

        uint64_t flength1 = GetFileSize(argv[optind]);
        uint64_t flength2 = GetFileSize(argv[optind+1]);

        uint64_t slicearry1[MAX_THREAD_NUM] ={0};
        uint64_t slicearry2[MAX_THREAD_NUM] ={0};

        GetFileSlice(argv[optind], flength1, slicearry1);

        bool isSE = true;
        if(flength2)
        {
            isSE = false;
            GetFileSlice(argv[optind+1], flength2, slicearry2);
            //AdjustPESlice(argv[optind+1], slicearry1, argv[optind+2], slicearry2);
            optind++;
        }

        g_magicparam.one_ch = true;
        if (!g_magicparam.fqzall)
            g_magicparam.fqzall = DoPreAlign(itr, idx, isSE, argv[optind], flength1, argv[optind+1], flength2);//执行预比对
        g_fqz_params.fqzall = g_magicparam.fqzall;

        g_magicparam.both_strands = g_fqz_params.both_strands;
        g_magicparam.extreme_seq = g_fqz_params.extreme_seq;
        g_magicparam.multi_seq_model = g_fqz_params.multi_seq_model;
        g_magicparam.isGzip = g_isgzip;
        g_magicparam.isSE = isSE;
        g_magicparam.slevel = g_fqz_params.slevel;
        g_magicparam.qlevel = g_fqz_params.qlevel;
        g_magicparam.nlevel = g_fqz_params.nlevel;
        g_magicparam.major_vers = MAJOR_VERS;
        g_magicparam.format_vers = FORMAT_VERS;
        g_magicparam.max_mis = max_mis;
        g_magicparam.max_insr = max_insr;
        g_magicparam.max_readLen = max_readLen;
        g_magicparam.thread_num = thread_num;

        char path[256]={0};
        string str_out("out");
        if(isSE)
        {
            if(argv[optind+1]) str_out = argv[optind+1];
        }
        else
        {
            if(argv[optind+2]) str_out = argv[optind+2];
        }
        sprintf(path, "./%s.arc", str_out.c_str());
        fstream out_s;
        out_s.open(path, std::ios::binary|std::ios::out);

        pthread_t *t_id = (pthread_t*)alloca(thread_num * sizeof(pthread_t));

        if(g_isgzip) //gzip文件
        {
            SeqRead ssread1(argv[optind], 0, flength1);
            ThreadTask *pthreadtask = new ThreadTask(thread_num);

            for (i = 0; i < thread_num; ++i) //创建一个线程池，等待任务
            {
                ThreadParam *pParam = new ThreadParam;
                pParam->num = i;
                pParam->pidx = idx;
                pParam->pitr = itr;
                pParam->pthreadtask = pthreadtask;
                pthread_create(&t_id[i], 0, gzip_process, pParam);
            }

            uint64_t read_num = 0;
            if(isSE)
            {
                while(ssread1.getRead())
                {
                    Task task;
                    SetTaskData(ssread1, 0, task);
                    pthreadtask->setTask(task, (read_num++)%thread_num);
                }
            }
            else
            {
                SeqRead ssread2(argv[optind+1], 0, flength2);

                while(ssread1.getRead() && ssread2.getRead())
                {
                    Task task;
                    SetTaskData(ssread1, 0, task);
                    SetTaskData(ssread2, 1, task);

                    pthreadtask->setTask(task, (read_num++)%thread_num);
                }
            }

            g_bFinish = true;
        }
        else //普通文本文件
        {
            for (i = 0; i < thread_num; ++i) //创建一个线程池，并行执行压缩
            {
                EncodeParam * test = new EncodeParam;
                test->num = i;
                test->pidx = idx;
                test->pitr = itr;

                test->length[0] = slicearry1[i];
                test->offset[0] = Getoffset(slicearry1,i); 
                strcpy(test->filename[0], argv[optind]);

                if(!isSE)
                {
                    test->length[1] = slicearry2[i];
                    test->offset[1] = Getoffset(slicearry2,i); 
                    strcpy(test->filename[1], argv[optind+1]);
                }

                if(g_magicparam.fqzall) //预比对率低，使用fqz压缩
                {
                    pthread_create(&t_id[i], 0, fqzall_encode_process, test);
                }
                else
                {
                    pthread_create(&t_id[i], 0, encode_process, test);
                }
            }
        }
        
        for (i = 0; i < thread_num; ++i) 
        {
            pthread_join(t_id[i], 0);
        }

        for (i = 0; i < thread_num; i++) //合并文件
        {
            char str_tmp[64]={0};
            sprintf(str_tmp, "./out_isq_%d.tmp", i);
            fstream f_s;
            f_s.open(str_tmp, ios::in | ios::binary);
            out_s << f_s.rdbuf();
            f_s.close();
            std::remove(str_tmp); //delete the tmp file
        }  

        for (i = 0; i < thread_num; i++) //写入压缩片段的字节数，用于解压分片
        {
            char len_buf[5]={0};
            len_buf[0] = (g_lentharry[i] >>  0) & 0xff;
            len_buf[1] = (g_lentharry[i] >>  8) & 0xff;
            len_buf[2] = (g_lentharry[i] >> 16) & 0xff;
            len_buf[3] = (g_lentharry[i] >> 24) & 0xff;
            len_buf[4] = (g_lentharry[i] >> 32) & 0xff;

            out_s.write(len_buf, 5);
        }    
        

        out_s.write((char *)&g_magicparam,sizeof(g_magicparam));
        out_s.write(".arc",4);
        out_s.close();

        gettimeofday(&timeEnd, NULL);
        double totald = (timeEnd.tv_sec - timeStart.tv_sec)*1000 + (timeEnd.tv_usec -timeStart.tv_usec)*1.0/1000;
        printf("%s%d encode total time %f ms\n", __FUNCTION__, __LINE__, totald);

        return 0;
    }
}