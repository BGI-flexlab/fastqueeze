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
    fprintf(fp, "    -M INT         Set model count when create model\n\n");

    fprintf(fp, "To decompress:\n   SeqArc -d [ref.fa] <***.arc>\n");

    exit(err);
}


int main(int argc, char **argv) {
    int opt, i, max_len = INT_MAX;
    uint64_t max_intv = 0;
    smem_i *itr = NULL;
    bwaidx_t *idx = NULL;
    int decompress = 0, indexing = 0;
    unsigned char md5 [MD5_DIGEST_LENGTH];
    int model_count = 10;
    /* Initialise and parse command line arguments */
    g_fqz_params.slevel = 3;
    g_fqz_params.qlevel = 2;
    g_fqz_params.nlevel = 1;
    g_fqz_params.both_strands = 0;
    g_fqz_params.extreme_seq = 0;
    g_fqz_params.multi_seq_model = 0;
    g_fqz_params.do_hash = 0;

    while ((opt = getopt(argc, argv, "l:w:I:f:m:q:s:hdQ:S:N:beXiB:t:n:c:E:r:WM:")) != -1) {
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
                if (thread_num > MAX_THREAD_NUM) {
                    printf("Thread num can not be more than %d\n", MAX_THREAD_NUM);
                    thread_num = MAX_THREAD_NUM;
                }
                break;

            case 'W':
                g_show_warning = true;
                break;
            case 'M':
                model_count = atoi(optarg);
                break;
            default:
                usage(1);
        }
    }

    CreatBitmap(nucleBitMap);

    if (argc == optind)
        usage(1);

    if (indexing) {
        char *ref = argv[optind];
        bwa_idx_build(ref, ref, BWTALGO_AUTO, 10000000);
        bwaidx_t *idx = bwa_idx_load_from_disk(ref, BWA_IDX_ALL);
        bwa_shm_stage(idx, ref, NULL);
        
        md5count(ref, md5);
        char path[512] = {0};
        sprintf(path, "%s.md5", ref);
        fstream out_s;
        out_s.open(path, std::ios::binary | std::ios::out);
        out_s.write((char*)md5, MD5_DIGEST_LENGTH);
        out_s.close();
        return 1;
    } else if (decompress) {
        struct timeval timeStart, timeEnd;
        gettimeofday(&timeStart, NULL);
        string tmpArg;
        bool hasFA;
        char path[256] = {0};
        string fastq_prefix;

        tmpArg = argv[optind];
        if (tmpArg.length() <= 5 ? tmpArg.substr(tmpArg.size() - 2) != "fa" :
            tmpArg.substr(tmpArg.size() - 2) != "fa" && tmpArg.substr(tmpArg.size() - 5) != "fasta") { //no fasta
            hasFA = false;
            sprintf(path, argv[optind]);
            fastq_prefix = argv[optind];
            fastq_prefix = fastq_prefix.substr(0, fastq_prefix.length()-4);
        } else {
            hasFA = true;
            idx = bwa_idx_load_from_shm(argv[optind]);
            if (idx == 0) {
                if ((idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1;
                bwa_shm_stage(idx, argv[optind], NULL);
            }

            uint64_t res_len = idx->bns->l_pac; //获取参考序列的长度，小于文件的真实长度
            g_offset_bit = getbitnum(res_len) - 2;

            sprintf(path, argv[optind + 1]);
            fastq_prefix = argv[optind + 1];
            fastq_prefix = fastq_prefix.substr(0, fastq_prefix.length()-4);
        }

        fstream in_s;
        in_s.open(path, std::ios::binary | std::ios::in);
        char buf[40960]={0};
        in_s.read(buf,40960);
        
        char *p = buf;
        if(memcmp(p,".arc",4) != 0)
        {
            printf("Unrecognised file format.\n");
            return 1;
        }
        p += 4;
        memcpy(&g_magicparam, p, sizeof(g_magicparam));
        p += sizeof(g_magicparam);
        if (!g_magicparam.fqzall){
            unsigned char md5_h [MD5_DIGEST_LENGTH];
            memcpy(md5_h, p, MD5_DIGEST_LENGTH);
            fstream md5_s;
            char md5_path[512] = {0};
            sprintf(md5_path, "%s.md5", argv[optind]);
            md5_s.open(md5_path, std::ios::binary | std::ios::in);
            md5_s.read((char*)md5, MD5_DIGEST_LENGTH);
            md5_s.close();
            //md5count(argv[optind], md5);
            if (memcmp(md5, md5_h, MD5_DIGEST_LENGTH) != 0) {
                cout << "Error: FASTA unmatched. Please check it again." << endl;
                return 1;
            }
            p += MD5_DIGEST_LENGTH;
        }
        

        if (g_magicparam.major_vers != MAJOR_VERS) {
            printf("Error: Unsupported file version %d\n", g_magicparam.major_vers); //TODO:should be downward compatible
            return 1;
        }
        if (!g_magicparam.fqzall && !hasFA) {
            printf("Error: FASTA index needed.\n");
            return 1;
        }

        g_fqz_params.both_strands = g_magicparam.both_strands;
        g_fqz_params.extreme_seq = g_magicparam.extreme_seq;
        g_fqz_params.multi_seq_model = g_magicparam.multi_seq_model;
        g_fqz_params.fqzall = g_magicparam.fqzall;
        g_fqz_params.slevel = g_magicparam.slevel;
        g_fqz_params.qlevel = g_magicparam.qlevel;
        g_fqz_params.nlevel = g_magicparam.nlevel;
        max_mis = g_magicparam.max_mis;
        max_insr = g_magicparam.max_insr;
        //max_readLen = g_magicparam.max_readLen;
        if(thread_num > g_magicparam.block_num){
            thread_num = g_magicparam.block_num;
        }
        
        int thread_block_count = g_magicparam.block_num / thread_num;
        int thread_block_left = g_magicparam.block_num % thread_num;
        int *pblockarry  = new int[g_magicparam.block_num];
        for (int i = 0; i < g_magicparam.block_num; i++) {
            unsigned char len_buf[4] = {0};
            memcpy(len_buf, p, 4);
            p += 4;
            pblockarry[i] =
                    (len_buf[0] << 0) +
                    (len_buf[1] << 8) +
                    (len_buf[2] << 16) +
                    (len_buf[3] << 24);
        }

        std::vector<uint64_t> vec_slice_length;
        int index = 0;
        for(int i=0;i<thread_num;i++){
            uint64_t length = 0;
            for(int j=0;j<thread_block_count;j++){
                length += pblockarry[index++];
            }
            if(thread_block_left-->0)
            {
                length += pblockarry[index++];                        
            }
            
            vec_slice_length.emplace_back(length);
        }
        delete[] pblockarry;
        assert(index == g_magicparam.block_num);


        unsigned char len_buf[4] = {0};
        memcpy(len_buf, p, 4);
        p += 4;
        int modelsize = (len_buf[0] << 0) +
                        (len_buf[1] << 8) +
                        (len_buf[2] << 16) +
                        (len_buf[3] << 24);
        int head_len = (int)(p - buf);
        char *pmodelbuf = new char[modelsize];
        in_s.seekg(head_len, std::ios::beg);
        in_s.read(pmodelbuf, modelsize);
        in_s.close();
        head_len += modelsize;

        pthread_t *t_id = (pthread_t *) alloca(thread_num * sizeof(pthread_t));
        int i;
        for (i = 0; i < thread_num; ++i) //创建一个线程池
        {
            DecodeParam *param = new DecodeParam;
            param->num = i;
            param->pidx = idx;
            param->length = vec_slice_length[i];
            param->offset = Getoffset(vec_slice_length, i, head_len);
            strcpy(param->filename, path);
            param->pmodelbuf = pmodelbuf;
            
            pthread_create(&t_id[i], 0, decode_process, param);
        }

        for (i = 0; i < thread_num; ++i) {
            pthread_join(t_id[i], 0);
        }

        delete[] pmodelbuf;

        MergeFileForFastq(fastq_prefix, thread_num, 1);
        if (!g_magicparam.isSE) {
            MergeFileForFastq(fastq_prefix, thread_num, 2);
        }

        gettimeofday(&timeEnd, NULL);
        double totald = (timeEnd.tv_sec - timeStart.tv_sec) * 1000 + (timeEnd.tv_usec - timeStart.tv_usec) * 1.0 / 1000;
        printf("%s%d decode total time %f ms\n", __FUNCTION__, __LINE__, totald);

        return 0;

    } else {
        if (argc - optind < 2 || argc - optind > 4) {
            cout << "Error: Wrong Param number, please check it again." << endl;
            usage(1);
        }
        struct timeval timeStart, timeEnd;
        gettimeofday(&timeStart, NULL);
        string tmpArg, outIndex;

        tmpArg = argv[optind];
        if (tmpArg.substr(tmpArg.size() - 2) != "fa" && tmpArg.substr(tmpArg.size() - 5) != "fasta") { //no fasta
            cout << "Info: FASTA index missing." << endl;
            g_magicparam.fqzall = true;
        } else {
            g_magicparam.fqzall = false;
            idx = bwa_idx_load_from_shm(argv[optind]);
            if (idx == 0) {
                if ((idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0){
                    cout << "Error: fail when loading bwa index." << endl;
                    return 1;
                }
                bwa_shm_stage(idx, argv[optind], NULL);
            }
            itr = smem_itr_init(idx->bwt);
            smem_config(itr, 1, max_len, max_intv); //min_intv = 1
            uint64_t res_len = idx->bns->l_pac; //获取参考序列的长度，小于文件的真实长度
            g_offset_bit = getbitnum(res_len) - 2;
            g_offset_size = pow(2, g_offset_bit) - 1;
            
            fstream md5_s;
            char md5_path[512] = {0};
            sprintf(md5_path, "%s.md5", argv[optind]);
            md5_s.open(md5_path, std::ios::binary | std::ios::in);
            md5_s.read((char*)md5, MD5_DIGEST_LENGTH);
            md5_s.close();

            ++optind;
        }

        g_magicparam.isGzip = GetFileType(argv[optind]);

        uint64_t flength1 = GetFileSize(argv[optind]);
        uint64_t flength2 = GetFileSize(argv[optind + 1]);
        string str_file1(argv[optind]);
        string str_file2("");
        uint64_t slicearry1[MAX_THREAD_NUM] = {0};
        uint64_t slicearry2[MAX_THREAD_NUM] = {0};

        GetFileSlice(argv[optind], flength1, slicearry1);

        bool isSE = true;
        if (flength2) {
            isSE = false;
            GetFileSlice(argv[optind + 1], flength2, slicearry2);
            str_file2 = argv[optind + 1];
        }

        g_magicparam.one_ch = true;
        if (!g_magicparam.fqzall)
            g_magicparam.fqzall = DoPreAlign(itr, idx, isSE, argv[optind], flength1, argv[optind + 1], flength2);//执行预比对
        g_fqz_params.fqzall = g_magicparam.fqzall;

        g_magicparam.both_strands = g_fqz_params.both_strands;
        g_magicparam.extreme_seq = g_fqz_params.extreme_seq;
        g_magicparam.multi_seq_model = g_fqz_params.multi_seq_model;
        g_magicparam.isSE = isSE;
        g_magicparam.slevel = g_fqz_params.slevel;
        g_magicparam.qlevel = g_fqz_params.qlevel;
        g_magicparam.nlevel = g_fqz_params.nlevel;
        g_magicparam.major_vers = MAJOR_VERS;
        g_magicparam.max_mis = max_mis;
        g_magicparam.max_insr = max_insr;
        //g_magicparam.max_readLen = max_readLen;
        g_magicparam.block_num = 0;

        char path[256] = {0};
        string str_out("out");
        if (isSE) {
            if (argv[optind + 1]) str_out = argv[optind + 1];
        } else {
            if (argv[optind + 2]) str_out = argv[optind + 2];
        }
        sprintf(path, "%s.arc", str_out.c_str());
        fstream out_s;
        out_s.open(path, std::ios::binary | std::ios::out);

        char *pmodelbuf = new char[50*1024*1024];
        int bufsize = CreateModel(model_count, str_file1, str_file2, pmodelbuf);

        pthread_t *t_id = (pthread_t *) alloca(thread_num * sizeof(pthread_t));

        if (g_magicparam.isGzip) //gzip文件
        {
            mkdir("./outtmp", 0777);
            for (i = 0; i < thread_num; ++i) //创建一个线程池，并行执行压缩
            {
                std::vector<int> vec_blocksize;
                g_vec_blocksize.emplace_back(vec_blocksize);

                EncodeParam *test = new EncodeParam;
                test->num = i;
                test->pidx = idx;
                test->pitr = itr;
                test->pmodelbuf = pmodelbuf;
                test->length[0] = 0;
                test->offset[0] = 0;
                strcpy(test->filename[0], argv[optind]);

                if (!isSE) {
                    test->length[1] = 0;
                    test->offset[1] = 0;
                    strcpy(test->filename[1], argv[optind + 1]);
                }
                pthread_create(&t_id[i], 0, gzip_encode_process, test);
            }
        }
        else //普通文本文件
        {
            for (i = 0; i < thread_num; ++i) //创建一个线程池，并行执行压缩
            {
                std::vector<int> vec_blocksize;
                g_vec_blocksize.emplace_back(vec_blocksize);

                EncodeParam *test = new EncodeParam;
                test->num = i;
                test->pidx = idx;
                test->pitr = itr;
                test->pmodelbuf = pmodelbuf;
                test->length[0] = slicearry1[i];
                test->offset[0] = Getoffset(slicearry1, i);
                strcpy(test->filename[0], argv[optind]);

                if (!isSE) {
                    test->length[1] = slicearry2[i];
                    test->offset[1] = Getoffset(slicearry2, i);
                    strcpy(test->filename[1], argv[optind + 1]);
                }
                pthread_create(&t_id[i], 0, fastq_encode_process, test);
            }
        }

        for (i = 0; i < thread_num; ++i) {
            pthread_join(t_id[i], 0);
        }

        
        for (i = 0; i < thread_num; i++) {
            g_magicparam.block_num += g_vec_blocksize[i].size();
        }

        char *psizebuf = new char[4*g_magicparam.block_num];
        char *pbuf = psizebuf;      
        if(g_magicparam.isGzip) {
            int j = 0;
            int *parry = new int[thread_num];
            memset(parry, 1 , thread_num*sizeof(int));
            while (isFinish(parry, thread_num)) {
                for (i = 0; i < thread_num; i++) {
                    if (parry[i] && j < g_vec_blocksize[i].size()) {
                        *pbuf++ = (g_vec_blocksize[i][j] >>  0) & 0xff;
                        *pbuf++ = (g_vec_blocksize[i][j] >>  8) & 0xff;
                        *pbuf++ = (g_vec_blocksize[i][j] >> 16) & 0xff;
                        *pbuf++ = (g_vec_blocksize[i][j] >> 24) & 0xff;
                    } else {
                        parry[i] = 0;
                    }
                }
                j++;
            }
            delete[] parry;
        } else {
            for (i = 0; i < thread_num; i++) //写入压缩片段的字节数，用于解压分片
            {
                for (int j = 0; j < g_vec_blocksize[i].size(); j++)
                {
                    *pbuf++ = (g_vec_blocksize[i][j] >>  0) & 0xff;
                    *pbuf++ = (g_vec_blocksize[i][j] >>  8) & 0xff;
                    *pbuf++ = (g_vec_blocksize[i][j] >> 16) & 0xff;
                    *pbuf++ = (g_vec_blocksize[i][j] >> 24) & 0xff;
                }
            }
        }

        //writing header
        out_s.write(".arc", 4);
        out_s.write((char *) &g_magicparam, sizeof(g_magicparam));

        if (!g_magicparam.fqzall)
            out_s.write((char *) md5, MD5_DIGEST_LENGTH);

        int size = pbuf - psizebuf;
        out_s.write(psizebuf, size);//写入blocksize
        delete[] psizebuf;
        char len_buf[4]={0};
        len_buf[0] = (bufsize >>  0) & 0xff;
        len_buf[1] = (bufsize >>  8) & 0xff;
        len_buf[2] = (bufsize >> 16) & 0xff;
        len_buf[3] = (bufsize >> 24) & 0xff;
        out_s.write(len_buf, 4);
        out_s.write(pmodelbuf, bufsize); //写入model
        delete[] pmodelbuf;

        if(g_magicparam.isGzip) {
            int j = 0;
            int *parry = new int[thread_num];
            memset(parry, 1 , thread_num*sizeof(int));
            while (isFinish(parry, thread_num)) {
                for (i = 0; i < thread_num; i++) {
                    if (parry[i] && j < g_vec_blocksize[i].size()) {
                        char str_tmp[64] = {0};
                        sprintf(str_tmp, "./outtmp/out_isq_%d_%d.tmp", i, j);
                        //std::cout<<str_tmp<<std::endl;
                        fstream f_s;
                        f_s.open(str_tmp, ios::in | ios::binary);
                        out_s << f_s.rdbuf();
                        f_s.close();
                    } else {
                        parry[i] = 0; //该线程生成的block已经处理完
                    }
                }
                j++;
            }
            delete[] parry;
            system("rm -rf ./outtmp");
        } else {
            for (i = 0; i < thread_num; i++) {
                char str_tmp[64] = {0};
                sprintf(str_tmp, "./out_isq_%d.tmp", i);
                fstream f_s;
                f_s.open(str_tmp, ios::in | ios::binary);
                out_s << f_s.rdbuf();
                f_s.close();
                std::remove(str_tmp);  // delete the tmp file
            }
        }

        out_s.close();

        gettimeofday(&timeEnd, NULL);
        double totald = (timeEnd.tv_sec - timeStart.tv_sec) * 1000 + (timeEnd.tv_usec - timeStart.tv_usec) * 1.0 / 1000;
        printf("%s%d encode total time %f ms\n", __FUNCTION__, __LINE__, totald);

#ifdef DEBUG
        cout << time4prepare << '\t' << time4find << '\t' << time4extend;
#endif
        return 0;
    }
}
