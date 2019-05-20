#include "fqzcomp.h"

const char *dec = "ACGTNMRYKSWHBVD";
#define MODEL_SIZE(num) ((8+(num)*3))

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

FqzComp::FqzComp(fqz_params *p){
    if (p) {
        slevel = p->slevel;
        qlevel = p->qlevel;
        nlevel = p->nlevel;
        both_strands = p->both_strands;
        extreme_seq = p->extreme_seq;
        multi_seq_model = p->multi_seq_model;
        do_hash = p->do_hash; // negligible slow down
        do_fqzall = p->fqzall;
    } else {
        slevel = 3;
        qlevel = 2;
        nlevel = 1;
        both_strands = 0;
        extreme_seq = 0;
        multi_seq_model = 0;
        do_hash = 0;
        do_fqzall = 1;
    }

    /* ACGTN* */
    for (int i = 0; i < 256; i++)
        L[i] = 0;
    L['A'] = L['a'] = 0;
    L['C'] = L['c'] = 1;
    L['G'] = L['g'] = 2;
    L['T'] = L['t'] = 3;

    L['N'] = L['n'] = 4;
    L['M'] = L['m'] = 5;
    L['R'] = L['r'] = 6;
    L['Y'] = L['y'] = 7;
    L['K'] = L['k'] = 8;
    L['S'] = L['s'] = 9;
    L['W'] = L['w'] = 10;
    L['H'] = L['h'] = 11;
    L['B'] = L['b'] = 12;
    L['V'] = L['v'] = 13;
    L['D'] = L['d'] = 14;

    NS = 7 + slevel;

    memset(last_name, ' ', 1024);
    last_name_len = 0;
    last_p_len = 0;
    last_s_len = 0;

    /* Length settings */
    last_len = 0;

    name_buf = new char[BLK_SIZE];
    seq_buf = new char[BLK_SIZE / 2];
    qual_buf = new char[BLK_SIZE];


    name_p = name_buf;
    seq_p = seq_buf;
    qual_p = qual_buf;


    inLen = 0;
    outLen = 0;
    pass_len = 0;
    uncomp_len = 0;

    name_len_a = new uint16_t[BLK_SIZE / 10];
    qual_len_a = new uint16_t[BLK_SIZE / 10];

    out0 = new char[BLK_SIZE / 2]; // seq_len
    out1 = new char[BLK_SIZE / 2]; // name
    out2 = new char[BLK_SIZE / 2]; // seq
    out3 = new char[BLK_SIZE / 2]; // qual

    m_totallen = 0;
    isPEccordant = true;

    if (!do_fqzall) {
        seq_len_a = new uint16_t[BLK_SIZE / 10];

        bit_buf = new char[BLK_SIZE / 2];
        bit_p = bit_buf;
        bit_len = 0;
        order_buf = new uint8_t[BLK_SIZE / 2];
        seq_count = 0;

        out4 = new char[BLK_SIZE / 10]; // bittobuf
        sz4 = 0;

        out5 = new char[BLK_SIZE / 10]; //orderbuf
        sz5 = 0;

        out6 = new char[10 * 1024];
        sz6 = 0;
    }

    vec_degenerate.reserve(100 * 1024);

    m_seq_size = 1 << (2 * NS);
    m_qual_size = QSIZE*16;
}

FqzComp::~FqzComp(){
    delete[] name_buf; name_buf = nullptr;
    delete[] seq_buf; seq_buf = nullptr;
    delete[] qual_buf; qual_buf = nullptr;

    delete[] name_len_a; name_len_a = nullptr;
    delete[] qual_len_a; qual_len_a = nullptr;

    delete[] out0; out0 = nullptr;
    delete[] out1; out1 = nullptr; 
    delete[] out2; out2 = nullptr;
    delete[] out3; out3 = nullptr;

    if (!do_fqzall) {
        delete[] seq_len_a;
        delete[] bit_buf;
        delete[] order_buf;
        delete[] out4;
        delete[] out5;
        delete[] out6;
    }
}

void FqzComp::InitModel(){

    model_name_prefix = new SIMPLE_MODEL<256>[256];
    model_name_suffix = new SIMPLE_MODEL<256>[256];
    model_name_len = new SIMPLE_MODEL<256>[256];
    model_name_middle = new SIMPLE_MODEL<128>[8192];

    if (extreme_seq) {
        model_seq8 = NULL;
        model_seq16 = new BASE_MODEL<uint16_t>[m_seq_size];
    } else {
        model_seq8 = new BASE_MODEL<uint8_t>[m_seq_size];
        model_seq16 = NULL;
    }

    model_qual = new SIMPLE_MODEL<QMAX>[m_qual_size];
}

void FqzComp::DestoryModel(){
    if(model_name_prefix){
        delete[] model_name_prefix;
        model_name_prefix = nullptr;
    }
    if (model_name_suffix){
        delete[] model_name_suffix;
        model_name_suffix = nullptr;
    }
    if(model_name_len){
        delete[] model_name_len;
        model_name_len = nullptr;
    }
    if (model_name_middle){
        delete[] model_name_middle;
        model_name_middle = nullptr;
    }
    if (model_seq8){
        delete[] model_seq8; 
        model_seq8 = nullptr;
    } 
    if (model_seq16){
        delete[] model_seq16; 
        model_seq16 = nullptr;
    } 
    if (model_qual){
        delete[] model_qual; 
        model_qual = nullptr;
    } 
}

int FqzComp::SaveModelToMem(char *modelbuf){
    char *pbuf = modelbuf;
    int len_test = model_same_len.save_test(pbuf);
    pbuf += len_test;
    len_test = model_len1.save_test(pbuf);
    pbuf += len_test;
    len_test = model_len2.save_test(pbuf);
    pbuf += len_test;

    for(int i=0;i<256;i++)
    {
        len_test = model_name_prefix[i].save_test(pbuf);
        pbuf += len_test;
    }

    for(int i=0;i<256;i++)
    {
        len_test = model_name_suffix[i].save_test(pbuf);
        pbuf += len_test;
    }

    for(int i=0;i<256;i++)
    {
        len_test = model_name_len[i].save_test(pbuf);
        pbuf += len_test;
    }

    for(int i=0;i<8192;i++)
    {
        len_test = model_name_middle[i].save_test(pbuf);
        pbuf += len_test;
    }

    len_test = seq_indicate.save_test(pbuf);
    pbuf += len_test;

    len_test = seq_degenerate.save_test(pbuf);
    pbuf += len_test;

    for(int i=0;i<m_seq_size;i++)
    {
        len_test = model_seq8[i].b_save_test(pbuf);
        pbuf += len_test;
    }

    for(int i=0;i<m_qual_size;i++)
    {
        len_test = model_qual[i].save_test(pbuf);
        pbuf += len_test;
    }

    return pbuf-modelbuf;
}

void FqzComp::SaveModelToFile(std::fstream &s_out){
    char tmp_buf[1048]={0};
    int len_test = model_same_len.save_test(tmp_buf);
    s_out.write(tmp_buf, len_test);

    len_test = model_len1.save_test(tmp_buf);
    s_out.write(tmp_buf, len_test);

    len_test = model_len2.save_test(tmp_buf);
    s_out.write(tmp_buf, len_test);

    for(int i=0;i<256;i++)
    {
        len_test = model_name_prefix[i].save_test(tmp_buf);
        s_out.write(tmp_buf, len_test);
    }

    for(int i=0;i<256;i++)
    {
        len_test = model_name_suffix[i].save_test(tmp_buf);
        s_out.write(tmp_buf, len_test);
    }

    for(int i=0;i<256;i++)
    {
        len_test = model_name_len[i].save_test(tmp_buf);
        s_out.write(tmp_buf, len_test);
    }

    for(int i=0;i<8192;i++)
    {
        len_test = model_name_middle[i].save_test(tmp_buf);
        s_out.write(tmp_buf, len_test);
    }

    len_test = seq_indicate.save_test(tmp_buf);
    s_out.write(tmp_buf, len_test);

    len_test = seq_degenerate.save_test(tmp_buf);
    s_out.write(tmp_buf, len_test);

    for(int i=0;i<m_seq_size;i++)
    {
        len_test = model_seq8[i].b_save_test(tmp_buf);
        s_out.write(tmp_buf, len_test);
    }

    for(int i=0;i<m_qual_size;i++)
    {
        len_test = model_qual[i].save_test(tmp_buf);
        s_out.write(tmp_buf, len_test);
    }
}

void FqzComp::ReadModelFormMem(char *pbuf){
    p_buf_same_len = pbuf;
    pbuf+=MODEL_SIZE(2);
    p_buf_model_len1 = pbuf;
    pbuf+=MODEL_SIZE(256);
    p_buf_model_len2 = pbuf;
    pbuf+=MODEL_SIZE(256);

    p_buf_name_prefix = pbuf;
    pbuf += MODEL_SIZE(256)*256;
    p_buf_name_suffix = pbuf;
    pbuf += MODEL_SIZE(256)*256;
    p_buf_name_len = pbuf;
    pbuf += MODEL_SIZE(256)*256;
    p_buf_name_middle = pbuf;
    pbuf += MODEL_SIZE(128)*8192;

    p_buf_seq_indicate = pbuf;
    pbuf+=MODEL_SIZE(2);
    p_buf_seq_degenerate = pbuf;
    pbuf+=MODEL_SIZE(11);
    p_buf_model_seq8 = pbuf;
    pbuf += 4 * m_seq_size;

    p_buf_model_qual = pbuf;
}

void FqzComp::ReadModelFormFile(std::fstream &s_in){
    p_buf_same_len = new char[MODEL_SIZE(2)];
    p_buf_model_len1 = new char[MODEL_SIZE(256)];
    p_buf_model_len2 = new char[MODEL_SIZE(256)];
    s_in.read(p_buf_same_len, MODEL_SIZE(2));
    s_in.read(p_buf_model_len1, MODEL_SIZE(256));
    s_in.read(p_buf_model_len2, MODEL_SIZE(256));

    p_buf_name_prefix = new char[MODEL_SIZE(256)*256];
    p_buf_name_suffix = new char[MODEL_SIZE(256)*256];
    p_buf_name_len = new char[MODEL_SIZE(256)*256];
    p_buf_name_middle =  new char[MODEL_SIZE(128)*8192];
    s_in.read(p_buf_name_prefix, MODEL_SIZE(256)*256);
    s_in.read(p_buf_name_suffix, MODEL_SIZE(256)*256);
    s_in.read(p_buf_name_len, MODEL_SIZE(256)*256);
    s_in.read(p_buf_name_middle, MODEL_SIZE(128)*8192);    

    p_buf_seq_indicate = new char[MODEL_SIZE(2)];
    p_buf_seq_degenerate = new char[MODEL_SIZE(11)];
    int num_seq = 4 * m_seq_size;
    p_buf_model_seq8 = new char[num_seq];  
    s_in.read(p_buf_seq_indicate, MODEL_SIZE(2));
    s_in.read(p_buf_seq_degenerate, MODEL_SIZE(11));
    s_in.read(p_buf_model_seq8, num_seq);
    

    int num_qual = MODEL_SIZE(128) * m_qual_size;
    p_buf_model_qual = new char[num_qual];
    s_in.read(p_buf_model_qual, num_qual);
}

void FqzComp::DestoryModelbuf(){
    if(p_buf_same_len){
       delete[] p_buf_same_len;
       p_buf_same_len = nullptr;
    }
    if(p_buf_model_len1){
       delete[] p_buf_model_len1;
       p_buf_model_len1 = nullptr;
    }
    if(p_buf_model_len2){
       delete[] p_buf_model_len2;
       p_buf_model_len2 = nullptr;
    }
    if(p_buf_name_prefix){
        delete[] p_buf_name_prefix;
        p_buf_name_prefix = nullptr;
    }
    if(p_buf_name_suffix){
        delete[] p_buf_name_suffix;
        p_buf_name_suffix = nullptr;
    }
    if(p_buf_name_len){
        delete[] p_buf_name_len;
        p_buf_name_len = nullptr;
    }
    if(p_buf_name_middle){
        delete[] p_buf_name_middle;
        p_buf_name_middle = nullptr;
    }
    if(p_buf_seq_indicate){
        delete[] p_buf_seq_indicate;
        p_buf_seq_indicate = nullptr;
    }
    if(p_buf_seq_degenerate){
        delete[] p_buf_seq_degenerate;
        p_buf_seq_degenerate = nullptr;
    }
    if(p_buf_model_seq8){
        delete[] p_buf_model_seq8;
        p_buf_model_seq8 = nullptr;
    }
    if(p_buf_model_qual){
        delete[] p_buf_model_qual;
        p_buf_model_qual = nullptr;
    }
}

int FqzComp::Addbuf(char *id, int idlen, char *seq, int seqlen, char *qual, int quallen){
    inLen += idlen + seqlen + quallen;
    memcpy(name_p, id, idlen);
    name_p += idlen;
    name_len_a[ns] = idlen;
    memcpy(seq_p, seq, seqlen);
    seq_p += seqlen;
    qual_len_a[ns] = seqlen;
    memcpy(qual_p, qual, quallen);
    qual_p += quallen;

    ns++;

    return ns;
}

int FqzComp::EncodeForModel(){
    if (ns == 0) {
        return -1;
    }

    compress_r0(false, false);    
    compress_r1(false);
    compress_r2(false, false);
    compress_r3(false);

    inLen = 0;
    ns = 0;
    name_p = name_buf;
    seq_p = seq_buf;
    qual_p = qual_buf;

    return 0;
}

void FqzComp::compress_r0(bool bnewmodel, bool balign){
    if(bnewmodel){
        RangeCoder rc;
        rc.output(out0);
        rc.StartEncode();

        SIMPLE_MODEL<256> model_len1_t;
        SIMPLE_MODEL<256> model_len2_t;
        SIMPLE_MODEL<2> model_same_len_t;

        model_len1_t.init_test(p_buf_model_len1);
        model_len2_t.init_test(p_buf_model_len2);
        model_same_len_t.init_test(p_buf_same_len);

        for (int i = 0; i < ns; i++) {
            encode_len(&rc, qual_len_a[i], model_len1_t, model_len2_t, model_same_len_t);
        }

        if(balign) {
            for (int i = 0; i < seq_count; i++) {
                encode_len(&rc, seq_len_a[i], model_len1_t, model_len2_t, model_same_len_t);
            }
        }

        rc.FinishEncode();
        sz0 = rc.size_out();
    }else{
        for (int i = 0; i < ns; i++) {
            encode_len_formodel(qual_len_a[i], model_len1, model_len2, model_same_len);
        }
    }
}

void FqzComp::encode_len(RangeCoder *rc, int len, SIMPLE_MODEL<256> &model_len1_t, 
                            SIMPLE_MODEL<256> &model_len2_t, SIMPLE_MODEL<2> &model_same_len_t) {
    if (len != last_len) {
        model_same_len_t.encodeSymbol(rc, 0);
        model_len1_t.encodeSymbol(rc, len & 0xff);
        model_len2_t.encodeSymbol(rc, (len >> 8) & 0xff);
    } else {
        model_same_len_t.encodeSymbol(rc, 1);
    }
}

void FqzComp::encode_len_formodel(int len, SIMPLE_MODEL<256> &model_len1_t, 
                            SIMPLE_MODEL<256> &model_len2_t, SIMPLE_MODEL<2> &model_same_len_t) {
    if (len != last_len) {
        model_same_len_t.addSymbol(0);
        model_len1_t.addSymbol(len & 0xff);
        model_len2_t.addSymbol((len >> 8) & 0xff);
    } else {
        model_same_len_t.addSymbol(1);
    }
}

void FqzComp::encode_name(RangeCoder *rc, char *name, int len,
                           SIMPLE_MODEL<256> *model_name_prefix_t,
                           SIMPLE_MODEL<256> *model_name_suffix_t,
                           SIMPLE_MODEL<256> *model_name_len_t,
                           SIMPLE_MODEL<128> *model_name_middle_t) {
    int p_len, s_len; // prefix and suffix length
    int i, j, k, last_char;

    _mm_prefetch((const char *) &model_name_prefix_t[last_p_len], _MM_HINT_T0);
    _mm_prefetch((const char *) &model_name_suffix_t[last_s_len], _MM_HINT_T0);
    _mm_prefetch((const char *) &model_name_len_t[last_name_len], _MM_HINT_T0);

    // Prefix
    for (i = 0; i < len && i < last_name_len; i++) {
        if (name[i] != last_name[i])
            break;
    }
    p_len = i;

    // Suffix
    for (i = len - 1, j = last_name_len - 1; i >= 0 && j >= 0; i--, j--) {
        if (name[i] != last_name[j])
            break;
    }
    s_len = len - 1 - i;
    if (len - s_len - p_len < 0)
        s_len = len - p_len;

    model_name_prefix_t[last_p_len].encodeSymbol(rc, p_len);
    model_name_suffix_t[last_s_len].encodeSymbol(rc, s_len);
    model_name_len_t[last_name_len].encodeSymbol(rc, len);

    last_p_len = p_len;
    last_s_len = s_len;

    int len2 = len - s_len, lc2 = p_len ? 1 : 0;
    for (i = j = p_len, k = 0; i < len2; i++, j++, k++) {
        last_char = ((last_name[j] - 32) * 2 + lc2 + k * 64) % 8192;
        _mm_prefetch((const char *) &model_name_middle_t[last_char], _MM_HINT_T0);

        //last_char = ((last_name[j]-32)*2 + lc2 + (i+j)*32) % 8192;

        model_name_middle_t[last_char].encodeSymbol(rc, name[i] & 0x7f);

        //if (meta[name[i]]      && name[i] != last_name[j]) j++;
        //if (meta[last_name[j]] && name[i] != last_name[j]) j--;

        //if (meta[name[i]]) k = (k+3)>>2<<2;

        if (name[i] == ' ' && last_name[j] != ' ') j++;
        if (name[i] != ' ' && last_name[j] == ' ') j--;
        if (name[i] == ':' && last_name[j] != ':') j++;
        if (name[i] != ':' && last_name[j] == ':') j--;

        if (name[i] == ':' || name[i] == ' ') k = (k + 3) >> 2 << 2;

        lc2 = name[i] == last_name[j];
    }

    memcpy(last_name, name, len);
    last_name_len = len;
}

void FqzComp::encode_name_formodel(char *name, int len,
                           SIMPLE_MODEL<256> *model_name_prefix_t,
                           SIMPLE_MODEL<256> *model_name_suffix_t,
                           SIMPLE_MODEL<256> *model_name_len_t,
                           SIMPLE_MODEL<128> *model_name_middle_t) {
    int p_len, s_len; // prefix and suffix length
    int i, j, k, last_char;

    _mm_prefetch((const char *) &model_name_prefix_t[last_p_len], _MM_HINT_T0);
    _mm_prefetch((const char *) &model_name_suffix_t[last_s_len], _MM_HINT_T0);
    _mm_prefetch((const char *) &model_name_len_t[last_name_len], _MM_HINT_T0);

    // Prefix
    for (i = 0; i < len && i < last_name_len; i++) {
        if (name[i] != last_name[i])
            break;
    }
    p_len = i;

    // Suffix
    for (i = len - 1, j = last_name_len - 1; i >= 0 && j >= 0; i--, j--) {
        if (name[i] != last_name[j])
            break;
    }
    s_len = len - 1 - i;
    if (len - s_len - p_len < 0)
        s_len = len - p_len;

    model_name_prefix_t[last_p_len].addSymbol(p_len);
    model_name_suffix_t[last_s_len].addSymbol(s_len);
    model_name_len_t[last_name_len].addSymbol(len);

    last_p_len = p_len;
    last_s_len = s_len;

    int len2 = len - s_len, lc2 = p_len ? 1 : 0;
    for (i = j = p_len, k = 0; i < len2; i++, j++, k++) {
        last_char = ((last_name[j] - 32) * 2 + lc2 + k * 64) % 8192;
        _mm_prefetch((const char *) &model_name_middle_t[last_char], _MM_HINT_T0);

        model_name_middle_t[last_char].addSymbol(name[i] & 0x7f);

        if (name[i] == ' ' && last_name[j] != ' ') j++;
        if (name[i] != ' ' && last_name[j] == ' ') j--;
        if (name[i] == ':' && last_name[j] != ':') j++;
        if (name[i] != ':' && last_name[j] == ':') j--;

        if (name[i] == ':' || name[i] == ' ') k = (k + 3) >> 2 << 2;

        lc2 = name[i] == last_name[j];
    }

    memcpy(last_name, name, len);
    last_name_len = len;
}

void FqzComp::compress_r1(bool bnewmodel) {
    //int name_total = 0;
    char *name_p = name_buf;
    if(bnewmodel){
        RangeCoder rc;
        rc.output(out1);
        rc.StartEncode();

        SIMPLE_MODEL<256> model_name_prefix_t[256];
        char *p1 = p_buf_name_prefix;
        for(int i=0;i<256;i++){
            model_name_prefix_t[i].init_test(p1);
            p1 += MODEL_SIZE(256);
        }
        
        SIMPLE_MODEL<256> model_name_suffix_t[256];
        char *p2 = p_buf_name_suffix;
        for(int i=0;i<256;i++){
            model_name_suffix_t[i].init_test(p2);
            p2 += MODEL_SIZE(256);
        }
        
        SIMPLE_MODEL<256> model_name_len_t[256];
        char *p3 = p_buf_name_len;
        for(int i=0;i<256;i++){
            model_name_len_t[i].init_test(p3);
            p3 += MODEL_SIZE(256);
        }
        
        SIMPLE_MODEL<128> model_name_middle_t[8192];
        char *p4 = p_buf_name_middle;
        for(int i=0;i<8192;i++){
            model_name_middle_t[i].init_test(p4);
            p4 += MODEL_SIZE(128);
        }

        last_p_len = 0;
        last_s_len = 0;
        last_name_len = 0;
        memset(last_name, ' ', 1024);
        for (int i = 0; i < ns; i++){
            encode_name(&rc, name_p, name_len_a[i], model_name_prefix_t, model_name_suffix_t, model_name_len_t, model_name_middle_t);
            name_p += name_len_a[i];
            //name_total += name_len_a[i];
        }

        rc.FinishEncode();
        sz1 = rc.size_out();
    }else{
        for (int i = 0; i < ns; i++) {
            encode_name_formodel(name_p, name_len_a[i], model_name_prefix, model_name_suffix, model_name_len, model_name_middle);
            name_p += name_len_a[i];
            //name_total += name_len_a[i];
        }
    }

    //id_totallen_original += name_total;
    //id_totallen_encode += sz1;
    //std::cout<<"name--"<<name_total<<'\t'<<sz1<<'\t'<<name_total*1.0/sz1<<std::endl;
}

void FqzComp::encode_seq8(RangeCoder *rc, char *seq, int len, BASE_MODEL<uint8_t> *model_seq8_t,
                          SIMPLE_MODEL<2> &seq_indicate_t, SIMPLE_MODEL<11> &seq_degenerate_t) {
    int last, last2;
    int bc[4] = {(3 - 0) << (2 * NS - 2),
                 (3 - 1) << (2 * NS - 2),
                 (3 - 2) << (2 * NS - 2),
                 (3 - 3) << (2 * NS - 2)};
    const int NS_MASK = ((1 << (2 * NS)) - 1);

    /* Corresponds to a 12-mer word that doesn't occur in human genome. */ //TODO: 隐患
    last = 0x007616c7 & NS_MASK;
    last2 = (0x2c6b62ff >> (32 - 2 * NS)) & NS_MASK;

    _mm_prefetch((const char *) &model_seq8_t[last], _MM_HINT_T0);

    for (int i = 0; i < len; i++) {
        unsigned char b = L[(unsigned char) seq[i]];
        if (b > 3) //简并碱基
        {
            seq_indicate_t.encodeSymbol(rc, 1);
            seq_degenerate_t.encodeSymbol(rc, b - 4);
        } else {
            seq_indicate_t.encodeSymbol(rc, 0);
            unsigned int l2 = (last << 2) & NS_MASK;
            _mm_prefetch((const char *) &model_seq8_t[l2 + 0], _MM_HINT_T0);

            model_seq8_t[last].encodeSymbol(rc, b);

            last = ((last << 2) + b) & NS_MASK;
        }
    }

}

void FqzComp::encode_seq8_formodel(char *seq, int len, BASE_MODEL<uint8_t> *model_seq8_t,
                          SIMPLE_MODEL<2> &seq_indicate_t, SIMPLE_MODEL<11> &seq_degenerate_t) {
    int last, last2;
    int bc[4] = {(3 - 0) << (2 * NS - 2),
                 (3 - 1) << (2 * NS - 2),
                 (3 - 2) << (2 * NS - 2),
                 (3 - 3) << (2 * NS - 2)};
    const int NS_MASK = ((1 << (2 * NS)) - 1);

    /* Corresponds to a 12-mer word that doesn't occur in human genome. */ //TODO: 隐患
    last = 0x007616c7 & NS_MASK;
    last2 = (0x2c6b62ff >> (32 - 2 * NS)) & NS_MASK;

    _mm_prefetch((const char *) &model_seq8_t[last], _MM_HINT_T0);

    for (int i = 0; i < len; i++) {
        unsigned char b = L[(unsigned char) seq[i]];
        if (b > 3) //简并碱基
        {
            seq_indicate_t.addSymbol(1);
            seq_degenerate_t.addSymbol(b - 4);
        } else {
            seq_indicate_t.addSymbol(0);
            unsigned int l2 = (last << 2) & NS_MASK;
            _mm_prefetch((const char *) &model_seq8_t[l2 + 0], _MM_HINT_T0);

            model_seq8_t[last].addSymbol(b);

            last = ((last << 2) + b) & NS_MASK;
        }
    }
}

void FqzComp::compress_r2(bool bnewmodel, bool balign) {
    char *seq_p = seq_buf;

    //int seq_total = 0;

    if(bnewmodel){
        RangeCoder rc;
        rc.output(out2);
        rc.StartEncode();

        BASE_MODEL<uint8_t> * model_seq8_t = new BASE_MODEL<uint8_t>[m_seq_size];
        char *p = p_buf_model_seq8;
        for(int i=0;i<m_seq_size;i++)
        {
            model_seq8_t[i].b_init_test(p);
            p += 4;
        }

        SIMPLE_MODEL<2> seq_indicate_t; 
        seq_indicate_t.init_test(p_buf_seq_indicate);
        
        SIMPLE_MODEL<11> seq_degenerate_t;
        seq_degenerate_t.init_test(p_buf_seq_degenerate);

        int count = balign ? seq_count : ns;
        uint16_t *parry = balign ? seq_len_a : qual_len_a;
        for (int i = 0; i < count; i++) {
            encode_seq8(&rc, seq_p, parry[i], model_seq8_t, seq_indicate_t, seq_degenerate_t);
            seq_p += parry[i];
            //seq_total += parry[i];
        }
        delete[] model_seq8_t;
        model_seq8_t  = nullptr;
        rc.FinishEncode();
        sz2 = rc.size_out();
    }else{
        for (int i = 0; i < ns; i++) {
            encode_seq8_formodel(seq_p, qual_len_a[i], model_seq8, seq_indicate, seq_degenerate);
            seq_p += qual_len_a[i];
            //seq_total += qual_len_a[i];
        }
    }
    

    //seq_totallen_original += seq_total;
    //seq_totallen_encode += sz2;
    //std::cout<<"seq--"<<seq_total<<'\t'<<sz2<<'\t'<<seq_total*1.0/sz2<<std::endl;
}


void FqzComp::encode_qual(RangeCoder *rc, char *qual, int len, SIMPLE_MODEL<QMAX> *model_qual_t) {
    unsigned int last = 0;
    int delta = 5;
    int i, len2 = len;
    int q1 = 0, q2 = 0;

    /* Removing "Killer Bees" */
    while (len2 > 0 && qual[len2 - 1] == '#')
        len2--;

    for (i = 0; i < len2; i++) {
        unsigned char q = qual[i];

#ifdef MULTI_QUAL_MODEL
        if (model_qual_t[last].bias() > model_qual_t[last & SMALL_QMASK].bias()) {
            model_qual_t[last].encodeSymbol(rc, q);
        } else {
            model_qual_t[last & SMALL_QMASK].encodeSymbol(rc, q);
            model_qual_t[last].updateSymbol(q);
        }
#else
        model_qual_t[last].encodeSymbol(rc, q);
#endif

        // previous 2-3 bytes
        if (QBITS == 12) {
            last = ((MAX(q1, q2) << 6) + q) & ((1 << QBITS) - 1);
        } else {
            last = ((last << 6) + q) & ((1 << QBITS) - 1);
        }

        if (qlevel > 1) {
            last += (q1 == q2) << QBITS;
            // delta saves 3-4%, but adds 14% cpu
            delta += (q1 > q) * (q1 - q);
            last += (MIN(7 * 8, delta) & 0xf8) << (QBITS - 2);
        }

        if (qlevel > 2)
            //last += (MIN(i+7,127)&(7<<4))<<(QBITS);   // i>>4
            last += (MIN(i + 15, 127) & (15 << 3)) << (QBITS + 1);     // i>>3
        //last += (MIN(i+31,127)&(31<<2))<<(QBITS+2); // i>>2

        _mm_prefetch((const char *) &model_qual_t[last], _MM_HINT_T0);
        q2 = q1;
        q1 = q;

        assert(last < m_qual_size);
    }

    if (len != len2)
        model_qual_t[last].encodeSymbol(rc, QMAX - 1); /* terminator */
}

void FqzComp::encode_qual_formodel(char *qual, int len, SIMPLE_MODEL<QMAX> *model_qual_t) {
    unsigned int last = 0;
    int delta = 5;
    int i, len2 = len;
    int q1 = 0, q2 = 0;

    /* Removing "Killer Bees" */
    while (len2 > 0 && qual[len2 - 1] == '#')
        len2--;

    for (i = 0; i < len2; i++) {
        unsigned char q = qual[i];

        model_qual_t[last].addSymbol(q);

        // previous 2-3 bytes
        if (QBITS == 12) {
            last = ((MAX(q1, q2) << 6) + q) & ((1 << QBITS) - 1);
        } else {
            last = ((last << 6) + q) & ((1 << QBITS) - 1);
        }

        if (qlevel > 1) {
            last += (q1 == q2) << QBITS;
            // delta saves 3-4%, but adds 14% cpu
            delta += (q1 > q) * (q1 - q);
            last += (MIN(7 * 8, delta) & 0xf8) << (QBITS - 2);
        }

        if (qlevel > 2)
            //last += (MIN(i+7,127)&(7<<4))<<(QBITS);   // i>>4
            last += (MIN(i + 15, 127) & (15 << 3)) << (QBITS + 1);     // i>>3
        //last += (MIN(i+31,127)&(31<<2))<<(QBITS+2); // i>>2

        _mm_prefetch((const char *) &model_qual_t[last], _MM_HINT_T0);
        q2 = q1;
        q1 = q;

        assert(last < m_qual_size);
    }

    if (len != len2)
        model_qual_t[last].addSymbol(QMAX - 1); /* terminator */
}

void FqzComp::compress_r3(bool bnewmodel) {
    char *qual_p = qual_buf;

    //int qual_total = 0;
    if(bnewmodel){
        RangeCoder rc;
        rc.output(out3);
        rc.StartEncode();

        SIMPLE_MODEL<QMAX> *model_qual_t = new SIMPLE_MODEL<QMAX>[m_qual_size];
        char *p = p_buf_model_qual;
        for(int i=0;i<m_qual_size;i++)
        {
            model_qual_t[i].init_test(p);
            p += MODEL_SIZE(128);
        }

        for (int i = 0; i < ns; i++) {
            encode_qual(&rc, qual_p, qual_len_a[i], model_qual_t);
            qual_p += qual_len_a[i];
            //qual_total+=qual_len_a[i];
        }
        delete[] model_qual_t;
        model_qual_t = nullptr;
        rc.FinishEncode();
        sz3 = rc.size_out();
    }else{
        for (int i = 0; i < ns; i++) {
            encode_qual_formodel(qual_p, qual_len_a[i], model_qual);
            qual_p += qual_len_a[i];
            //qual_total+=qual_len_a[i];
        }
    }

    //qual_totallen_original += qual_total;
    //qual_totallen_encode += sz3;
    //std::cout<<"qual--"<<qual_total<<'\t'<<sz3<<'\t'<<qual_total*1.0/sz3<<std::endl;
}

void FqzComp::IntTo4Ch(int num, char *pbuf)
{
    *pbuf++ = (num >> 0) & 0xff;
    *pbuf++ = (num >> 8) & 0xff;
    *pbuf++ = (num >> 16) & 0xff;
    *pbuf++ = (num >> 24) & 0xff;
}

void FqzComp::IntTo2Ch(int num, char *pbuf)
{
    *pbuf++ = (num >> 0) & 0xff;
    *pbuf++ = (num >> 8) & 0xff;
}

int FqzComp::fqzall_encode(std::fstream &out) 
{
    if (ns == 0) {
        return -1;
    }

    compress_r0(true, false);
    compress_r1(true);
    compress_r2(true, false);
    compress_r3(true);

    int pre_len = isPEccordant ? 4 : 8;
    char *out_p0 = out_buf + pre_len;
    char *out_p = out_p0;
    
    IntTo4Ch(ns, out_p);
    out_p += 4;

    IntTo4Ch(sz0, out_p);
    out_p += 4;

    IntTo4Ch(sz1, out_p);
    out_p += 4;

    IntTo4Ch(sz2, out_p);
    out_p += 4;

    IntTo4Ch(sz3, out_p);
    out_p += 4;

    memcpy(out_p, out0, sz0);
    out_p += sz0;
    memcpy(out_p, out1, sz1);
    out_p += sz1;
    memcpy(out_p, out2, sz2);
    out_p += sz2;
    memcpy(out_p, out3, sz3);
    out_p += sz3;

    int out_len = (int) (out_p - out_p0);
    out_buf[pre_len - 4] = (out_len >> 0) & 0xff;
    out_buf[pre_len - 3] = (out_len >> 8) & 0xff;
    out_buf[pre_len - 2] = (out_len >> 16) & 0xff;
    out_buf[pre_len - 1] = (out_len >> 24) & 0xff;
    out_len += pre_len;
    
    //printf("%d %d %d %d %d %d %d\n", ns,sz0,sz1,sz2,sz3, out_len, inLen);
    //还原
    inLen = 0;
    ns = 0;
    name_p = name_buf;
    seq_p = seq_buf;
    qual_p = qual_buf;

    if (out_len != xwrite(out, (unsigned char *) out_buf, out_len)) {
        printf("isq_encode Abort: truncated write.0\n");
        return -1;
    } else {
        m_totallen += out_len;
    }
    
    return out_len;
}

void FqzComp::decompress_r0(bool balign){
    RangeCoder rc0;
    rc0.input(in_buf0);
    rc0.StartDecode();

    SIMPLE_MODEL<256> model_len1_t;
    model_len1_t.init_test(p_buf_model_len1);
    SIMPLE_MODEL<256> model_len2_t;
    model_len2_t.init_test(p_buf_model_len2);
    SIMPLE_MODEL<2> model_same_len_t;
    model_same_len_t.init_test(p_buf_same_len);

    for (int i = 0; i < ns; i++) {
        qual_len_a[i] = decode_len(&rc0, model_len1_t, model_len2_t, model_same_len_t);
    }

    if(balign){
        for (int i = 0; i < seq_count; i++) {
            seq_len_a[i] = decode_len(&rc0, model_len1_t, model_len2_t, model_same_len_t);
        }
    }
    rc0.FinishDecode();
}

int FqzComp::decode_name(RangeCoder *rc, char *name,
                           SIMPLE_MODEL<256> *model_name_prefix_t,
                           SIMPLE_MODEL<256> *model_name_suffix_t,
                           SIMPLE_MODEL<256> *model_name_len_t,
                           SIMPLE_MODEL<128> *model_name_middle_t) {
    int p_len, s_len, len; // prefix and suffix length
    int i, j, k;
    int last_char;

    _mm_prefetch((const char *) &model_name_prefix_t[last_p_len], _MM_HINT_T0);
    _mm_prefetch((const char *) &model_name_suffix_t[last_s_len], _MM_HINT_T0);
    _mm_prefetch((const char *) &model_name_len_t[last_name_len], _MM_HINT_T0);

    p_len = model_name_prefix_t[last_p_len].decodeSymbol(rc);
    s_len = model_name_suffix_t[last_s_len].decodeSymbol(rc);
    len = model_name_len_t[last_name_len].decodeSymbol(rc);

    last_p_len = p_len;
    last_s_len = s_len;

    for (i = 0; i < p_len; i++)
        name[i] = last_name[i];

    //fprintf(stderr, "%d: p_len 9= %d, s_len = %d, len = %d last='%.*s'\n",
    //column, p_len, s_len, len, last_name_len, last);

    int len2 = len - s_len, lc2 = p_len ? 1 : 0;
    for (i = j = p_len, k = 0; i < len2; i++, j++, k++) {
        unsigned char c;

        last_char = ((last_name[j] - 32) * 2 + lc2 + k * 64) % 8192;
        //last_char = (last_name[j] + j*64) % 8192;

        c = model_name_middle_t[last_char].decodeSymbol(rc);
        //c = 'x';
        name[i] = c;

        if (c == ' ' && last_name[j] != ' ') j++;
        if (c != ' ' && last_name[j] == ' ') j--;
        if (c == ':' && last_name[j] != ':') j++;
        if (c != ':' && last_name[j] == ':') j--;

        if (name[i] == ':' || name[i] == ' ') k = (k + 3) >> 2 << 2;

        lc2 = c == last_name[j];
    }

    for (j = last_name_len - s_len; i < len; i++, j++)
        name[i] = last_name[j];

    memcpy(last_name, name, len);
    last_name_len = len;

    return len;
}

void FqzComp::decompress_r1(void) {
    RangeCoder rc;
    rc.input(in_buf1);
    rc.StartDecode();

    char *name_p = name_buf;

    SIMPLE_MODEL<256> model_name_prefix_t[256];
    char *p1 = p_buf_name_prefix;
    for(int i=0;i<256;i++){
        model_name_prefix_t[i].init_test(p1);
        p1 += MODEL_SIZE(256);
    }
    SIMPLE_MODEL<256> model_name_suffix_t[256];
    char *p2 = p_buf_name_suffix;
    for(int i=0;i<256;i++){
        model_name_suffix_t[i].init_test(p2);
        p2 += MODEL_SIZE(256);
    }
    SIMPLE_MODEL<256> model_name_len_t[256];
    char *p3 = p_buf_name_len;
    for(int i=0;i<256;i++){
        model_name_len_t[i].init_test(p3);
        p3 += MODEL_SIZE(256);
    }
    SIMPLE_MODEL<128> model_name_middle_t[8192];
    char *p4 = p_buf_name_middle;
    for(int i=0;i<8192;i++){
        model_name_middle_t[i].init_test(p4);
        p4 += MODEL_SIZE(128);
    }

    last_p_len = 0;
    last_s_len = 0;
    last_name_len = 0;
    memset(last_name, ' ', 1024);
    for (int i = 0; i < ns; i++) {
        *name_p++ = '@';
        name_p += decode_name(&rc, name_p, model_name_prefix_t, model_name_suffix_t, model_name_len_t, model_name_middle_t);
        *name_p++ = '\n';
    }

    rc.FinishDecode();
}

void FqzComp::decode_seq8(RangeCoder *rc, char *seq, int len, BASE_MODEL<uint8_t> * model_seq8_t,
                           SIMPLE_MODEL<2> &seq_indicate_t, SIMPLE_MODEL<11> &seq_degenerate_t) {
    int last, last2;
    //const char *dec = "ACGTNMRYKSWHBVD";
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

    last = 0x7616c7 & NS_MASK;
    last2 = (0x2c6b62ff >> (32 - 2 * NS)) & NS_MASK;
    
    for (int i = 0; i < len; i++) {
        unsigned char b;
        if (seq_indicate_t.decodeSymbol(rc)) //简并碱基
        {
            b = seq_degenerate_t.decodeSymbol(rc) + 4;
        } else {
            unsigned int m = (last << 2) & NS_MASK;

            /* Get next set loaded */
            _mm_prefetch((const char *) &model_seq8_t[m + 0], _MM_HINT_T0);
            //_mm_prefetch((const char *)&model_seq8[m+3], _MM_HINT_T0);
            b = model_seq8_t[last].decodeSymbol(rc);
            last = (last * 4 + b) & NS_MASK;
        }

        *seq++ = dec[b];
    }

}

void FqzComp::decompress_r2(bool balign) {
    RangeCoder rc;
    rc.input(in_buf2);
    rc.StartDecode();

    BASE_MODEL<uint8_t> * model_seq8_t = new BASE_MODEL<uint8_t>[m_seq_size];
    char *p = p_buf_model_seq8;
    for(int i=0;i<m_seq_size;i++){
        model_seq8_t[i].b_init_test(p);
        p+=4;
    }
    SIMPLE_MODEL<2> seq_indicate_t; 
    seq_indicate_t.init_test(p_buf_seq_indicate);
    SIMPLE_MODEL<11> seq_degenerate_t;
    seq_degenerate_t.init_test(p_buf_seq_degenerate);

    char *seq_p = seq_buf;
    int count = balign ? seq_count : ns;
    uint16_t *parry = balign ? seq_len_a : qual_len_a;
    for (int i = 0; i < count; i++) {
        decode_seq8(&rc, seq_p, parry[i], model_seq8_t, seq_indicate_t, seq_degenerate_t);
        seq_p += parry[i];
    }
    delete[] model_seq8_t;
    rc.FinishDecode();
}

void FqzComp::decode_qual(RangeCoder *rc, char *qual, int len, SIMPLE_MODEL<QMAX> *model_qual_t) {
    unsigned int last = 0;
    int i;
    int delta = 5;
    int q1 = 0, q2 = 0;

    for (i = 0; i < len; i++) {
        unsigned char q = model_qual_t[last].decodeSymbol(rc);

        if (q == QMAX - 1) {
            while (i < len)
                qual[i++] = '#';
        } else {
            qual[i] = q;
            if (QBITS == 12) {
                last = ((MAX(q1, q2) << 6) + q) & ((1 << QBITS) - 1);
            } else {
                last = ((last << 6) + q) & ((1 << QBITS) - 1);
            }

            if (qlevel > 1) {
                last += (q1 == q2) << QBITS;
                delta += (q1 > q) * (q1 - q);
                last += (MIN(7 * 8, delta) & 0xf8) << (QBITS - 2);
            }

            if (qlevel > 2)
                //last += (MIN(i+7,127)&(7<<4))<<(QBITS);   // i>>4
                last += (MIN(i + 15, 127) & (15 << 3)) << (QBITS + 1);     // i>>3
            //last += (MIN(i+31,127)&(31<<2))<<(QBITS+2); // i>>2

            _mm_prefetch((const char *) &model_qual_t[last], _MM_HINT_T0);
            q2 = q1;
            q1 = q;
        }
    }
}

void FqzComp::decompress_r3(void) {
    RangeCoder rc;
    rc.input(in_buf3);
    rc.StartDecode();
    
    SIMPLE_MODEL<QMAX> *model_qual_t = new SIMPLE_MODEL<QMAX>[m_qual_size];
    char *p = p_buf_model_qual;
    for(int i=0;i<m_qual_size;i++){
        model_qual_t[i].init_test(p);
        p += MODEL_SIZE(128);
    }

    char *qual_p = qual_buf;
    for (int i = 0; i < ns; i++) {
        decode_qual(&rc, qual_p, qual_len_a[i], model_qual_t);
        qual_p += qual_len_a[i];
    }
    delete[] model_qual_t;
    rc.FinishDecode();
}

int FqzComp::decode_len(RangeCoder *rc, SIMPLE_MODEL<256> &model_len1_t,
                    SIMPLE_MODEL<256> &model_len2_t,
                    SIMPLE_MODEL<2> &model_same_len_t) {
    if (model_same_len_t.decodeSymbol(rc)) {
        return last_len;
    } else {
        int l1 = model_len1_t.decodeSymbol(rc);
        int l2 = model_len2_t.decodeSymbol(rc);
        last_len = l1 + (l2 << 8);
        return last_len;
    }
}

void FqzComp::isq_fqzall_decompress(char *in, int comp_len, int *out_len) {
    ns = DECODE_INT((unsigned char *) (in));
    sz0 = DECODE_INT((unsigned char *) (in + 4));
    sz1 = DECODE_INT((unsigned char *) (in + 8));
    sz2 = DECODE_INT((unsigned char *) (in + 12));
    sz3 = DECODE_INT((unsigned char *) (in + 16));

    in += 20;

    in_buf0 = in;
    in += sz0;
    in_buf1 = in;
    in += sz1;
    in_buf2 = in;
    in += sz2;
    in_buf3 = in;
    in += sz3;

    //printf("%d %d %d %d %d \n", ns,sz0,sz1,sz2,sz3);

    decompress_r0(false);
    decompress_r1();
    decompress_r2(false);
    decompress_r3();
}

int FqzComp::fqzall_decode(std::fstream &in, char **namebuf, char **seqbuf, char **qualbuf, uint16_t **seqlen, int *ins,
                    int *mark) {
    unsigned char len_buf[4];
    if (4 != xget(in, len_buf, 4))
        return -1;

    int read_len = 4;
    if (memcmp(len_buf, "f1f1", 4) == 0) {
        *mark = 1;
        xget(in, len_buf, 4);
        read_len += 4;
    } else if (memcmp(len_buf, "f2f2", 4) == 0) {
        *mark = 2;
        xget(in, len_buf, 4);
        read_len += 4;
    }

    int rem_len =
            (len_buf[0] << 0) +
            (len_buf[1] << 8) +
            (len_buf[2] << 16) +
            (len_buf[3] << 24);
    read_len += rem_len;

    do {
        errno = 0;
        int tmp_len = xget(in, (unsigned char *) in_buf, rem_len);
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
    } while (rem_len);

    isq_fqzall_decompress(in_buf, rem_len, &uncomp_len);
    *namebuf = name_buf;
    *seqbuf = seq_buf;
    *qualbuf = qual_buf;
    *seqlen = qual_len_a;
    *ins = ns;
    return read_len;
}

uint64_t FqzComp::getCompressTotalLen() {
    return m_totallen;
}

uint32_t FqzComp::getInLen() {
    return inLen;
}

uint32_t FqzComp::getReadCount() {
    return ns;
}

void FqzComp::isq_addmark(int mark) {
    if (mark == 1) {
        memcpy(out_buf, "f1f1", 4);
    } else if (mark == 2) {
        memcpy(out_buf, "f2f2", 4);
    }

    isPEccordant = false;
}

int FqzComp::isq_addbuf_match(char *id, int idlen, char *bit, int bitlen, char *qual, int quallen, int index,
                          char *degenerate) {
    inLen += idlen + quallen;
    memcpy(name_p, id, idlen);
    name_p += idlen;
    name_len_a[ns] = idlen;
    memcpy(bit_p, bit, bitlen);
    bit_p += bitlen;
    bit_len += bitlen;
    memcpy(qual_p, qual, quallen);
    qual_p += quallen;
    qual_len_a[ns] = quallen;

    for (int i = 0; i < strlen(degenerate); i++) {
        vec_degenerate.emplace_back(degenerate[i]);
    }

    order_buf[ns] = index;
    ns++;

    return ns;
}

int FqzComp::isq_addbuf_unmatch(char *id, int idlen, char *seq, int seqlen, char *qual, int quallen, int index) {
    inLen += idlen + seqlen + quallen;
    memcpy(name_p, id, idlen);
    name_p += idlen;
    name_len_a[ns] = idlen;
    memcpy(seq_p, seq, seqlen);
    seq_p += seqlen;
    seq_len_a[seq_count] = seqlen;
    memcpy(qual_p, qual, quallen);
    qual_p += quallen;
    qual_len_a[ns] = quallen;

    order_buf[ns] = index;

    ns++;
    seq_count++;
    return ns;
}

int FqzComp::align_encode(std::fstream &out) {
    if (ns == 0) return -1;

    compress_r0(true, true);
    compress_r1(true);
    compress_r2(true, true);
    compress_r3(true);

    BitArryToBuf_s(bit_buf, bit_len, out4, &sz4);


    RangeCoder rc5;
    rc5.output(out5);
    rc5.StartEncode();
    SIMPLE_MODEL<5> seq_order_t;
    for (int i = 0; i < ns; i++) {
        seq_order_t.encodeSymbol(&rc5, order_buf[i]);
    }
    rc5.FinishEncode();
    sz5 = rc5.size_out();

    int len_degenerate = vec_degenerate.size();
    if (len_degenerate > 0) {
        RangeCoder rc6;
        rc6.output(out6);
        rc6.StartEncode();
        SIMPLE_MODEL<11> seq_degenerate_match_t;
        for (int i = 0; i < len_degenerate; i++) {
            unsigned char b = L[(unsigned char) vec_degenerate[i]];
            if (b >= 4) {
                seq_degenerate_match_t.encodeSymbol(&rc6, b - 4);
            } else {
                seq_degenerate_match_t.encodeSymbol(&rc6, b); //异常碱基设成N
                //printf("%c\n", (unsigned char) vec_degenerate[i]);
            }

        }
        rc6.FinishEncode();
        sz6 = rc6.size_out();
    } else {
        sz6 = 0;
    }

    int pre_len = isPEccordant ? 4 : 8;
    char *out_p0 = out_buf + pre_len;
    char *out_p = out_p0;

    IntTo4Ch(seq_count, out_p);
    out_p += 4;

    // IntTo4Ch(inLen, out_p);
    // out_p += 4;

    IntTo4Ch(ns, out_p);
    out_p += 4;

    IntTo4Ch(sz0, out_p);
    out_p += 4;

    IntTo4Ch(sz1, out_p);
    out_p += 4;

    IntTo4Ch(sz2, out_p);
    out_p += 4;

    IntTo4Ch(sz3, out_p);
    out_p += 4;

    IntTo4Ch(sz4, out_p);
    out_p += 4;

    IntTo4Ch(sz5, out_p);
    out_p += 4;

    IntTo2Ch(len_degenerate, out_p);
    out_p += 2;

    IntTo2Ch(sz6, out_p);
    out_p += 2;

    memcpy(out_p, out0, sz0);
    out_p += sz0;
    memcpy(out_p, out1, sz1);
    out_p += sz1;
    memcpy(out_p, out2, sz2);
    out_p += sz2;
    memcpy(out_p, out3, sz3);
    out_p += sz3;
    memcpy(out_p, out4, sz4);
    out_p += sz4;
    memcpy(out_p, out5, sz5);
    out_p += sz5;
    memcpy(out_p, out6, sz6);
    out_p += sz6;

    int out_len = (int) (out_p - out_p0);
    out_buf[pre_len - 4] = (out_len >> 0) & 0xff;
    out_buf[pre_len - 3] = (out_len >> 8) & 0xff;
    out_buf[pre_len - 2] = (out_len >> 16) & 0xff;
    out_buf[pre_len - 1] = (out_len >> 24) & 0xff;
    out_len += pre_len;
    
    //还原
    inLen = 0;
    ns = 0;
    seq_count = 0;
    bit_len = 0;
    name_p = name_buf;
    seq_p = seq_buf;
    qual_p = qual_buf;
    bit_p = bit_buf;
    vec_degenerate.clear();


    if (out_len != xwrite(out, (unsigned char *) out_buf, out_len)) {
        //fprintf(stderr, "Abort: truncated write.0\n");
        printf("isq_encode Abort: truncated write.0\n");
        return -1;
    } else {
        m_totallen += out_len;
    }
    return out_len;
}

int FqzComp::align_decode(std::fstream &in, char **namebuf, char **seqbuf, char **qualbuf, char **bitbuf, uint8_t **orderbuf,
                  uint16_t **quallen, int *ins, int *mark, std::vector<char> **pvec) {
    unsigned char len_buf[4];
    if (4 != xget(in, len_buf, 4))
        return -1;
    int read_len = 4;
    if (memcmp(len_buf, "f1f1", 4) == 0) {
        *mark = 1;
        xget(in, len_buf, 4);
        read_len += 4;
    } else if (memcmp(len_buf, "f2f2", 4) == 0) {
        *mark = 2;
        xget(in, len_buf, 4);
        read_len += 4;
    }

    int rem_len =
            (len_buf[0] << 0) +
            (len_buf[1] << 8) +
            (len_buf[2] << 16) +
            (len_buf[3] << 24);
    read_len += rem_len;

    do {
        errno = 0;
        int tmp_len = xget(in, (unsigned char *) in_buf, rem_len);
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
    } while (rem_len);

    isq_align_decompress(in_buf, rem_len, &uncomp_len);
    *namebuf = name_buf;
    *seqbuf = seq_buf;
    *qualbuf = qual_buf;
    *bitbuf = in_buf4;
    *orderbuf = order_buf;
    *quallen = qual_len_a;
    *ins = ns;
    *pvec = &vec_degenerate;
    return read_len;
}

void FqzComp::isq_align_decompress(char *in, int comp_len, int *out_len) {
    seq_count = DECODE_INT((unsigned char *) (in));
    ns = DECODE_INT((unsigned char *) (in + 4));
    sz0 = DECODE_INT((unsigned char *) (in + 8));
    sz1 = DECODE_INT((unsigned char *) (in + 12));
    sz2 = DECODE_INT((unsigned char *) (in + 16));
    sz3 = DECODE_INT((unsigned char *) (in + 20));
    sz4 = DECODE_INT((unsigned char *) (in + 24));
    sz5 = DECODE_INT((unsigned char *) (in + 28));
    int num = DECODE_SHORT((unsigned char *) (in + 32));
    sz6 = DECODE_SHORT((unsigned char *) (in + 34));

    in += 36;
    in_buf0 = in;
    in += sz0;
    in_buf1 = in;
    in += sz1;
    in_buf2 = in;
    in += sz2;
    in_buf3 = in;
    in += sz3;
    in_buf4 = in;
    in += sz4;
    in_buf5 = in;
    in += sz5;
    in_buf6 = in;
    in += sz6;

    decompress_r0(true);
    decompress_r1();
    decompress_r2(true);
    decompress_r3();

    RangeCoder rc5;
    rc5.input(in_buf5);
    rc5.StartDecode();
    SIMPLE_MODEL<5> seq_order_t;
    for (int i = 0; i < ns; i++) {
        order_buf[i] = seq_order_t.decodeSymbol(&rc5);
    }
    rc5.FinishDecode();

    vec_degenerate.clear();
    if (num > 0) {
        RangeCoder rc6;
        rc6.input(in_buf6);
        rc6.StartDecode();
        SIMPLE_MODEL<11> seq_degenerate_match_t;
        for (int i = 0; i < num; i++) {
            unsigned char b = seq_degenerate_match_t.decodeSymbol(&rc6) + 4;
            vec_degenerate.emplace_back(dec[b]);
        }
        rc6.FinishDecode();
    }
}

void FqzComp::BitArryToBuf_s(const char *pbit, int len, char *pbuf, int *buflen) {
    int count = 0;
    while (len > 8) {
        for (int i = 0; i < 8; i++) {
            m_bitset.set(7 - i, *pbit++ == '0' ? 0 : 1);
        }

        len -= 8;
        *pbuf++ = m_bitset.to_ulong();
        count++;
    }
    if (len > 0) {
        m_bitset.reset();
        for (int i = 0; i < len; i++) {
            m_bitset.set(7 - i, *pbit++ == '0' ? 0 : 1);
        }

        *pbuf++ = m_bitset.to_ulong();
        count++;
    }

    *buflen = count;
}