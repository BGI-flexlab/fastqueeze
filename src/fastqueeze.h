#ifndef FASTQUEEZE_H
#define FASTQUEEZE_H

#include "libzpaq.h"
#include <vector>
#include <string>

namespace fastqueeze {

typedef pthread_t ThreadID; // job ID type

class CompressJob;

class fqcompressor {
    int threads;
    libzpaq::StringBuffer id_buffer[3], id[3], qs[4];
    int id_prefix;
    int32_t id_cnt;
    std::vector<std::string> qs_raw;
    CompressJob* job;
    std::vector<ThreadID> tid;
    ThreadID wid;

    void split_qs(const std::string& center, const float s, libzpaq::Writer& out);
public:
    fqcompressor(int _threads);
    void end();
    void qs_add(const char* s) { qs_raw.push_back(s); }
    void qs_compress();
    void get_qs(libzpaq::StringBuffer& sb, int i) { sb.swap(qs[i]); }
    void id_add(const char* s);
    void end_id_add();
    void get_id(libzpaq::StringBuffer& sb, int i) { sb.swap(id[i]); }
};

class ExtractJob;

class fqdecompressor {
    int threads;
    std::vector<std::string> qs_raw[4];
    std::vector<std::string> id_raw;
    int id_compressed_cnt;
    ExtractJob* job;
    std::vector<ThreadID> tid;
public:
    fqdecompressor(int _threads);
    void qs_add(libzpaq::StringBuffer& q, int i);
    void start();
    void end();
    void get_qs(std::vector<std::string>& ans);
    void id_add(libzpaq::StringBuffer& sb, int i);
    void get_id(std::vector<std::string>& ans) { ans.swap(id_raw); }
};

}

#endif  // FASTQUEEZE_H
