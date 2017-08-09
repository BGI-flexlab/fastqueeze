#include "fastqueeze.h"
#include "libzpaq.h"
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <pthread.h>
#include <cstdio>
#include <iostream>

using namespace libzpaq;
using namespace std;

namespace fastqueeze {

typedef void* ThreadReturn;                                // job return type
void run(ThreadID& tid, ThreadReturn(*f)(void*), void* arg)// start job
    {pthread_create(&tid, NULL, f, arg);}
void join(ThreadID tid) {pthread_join(tid, NULL);}         // wait for job
typedef pthread_mutex_t Mutex;                             // mutex type
void init_mutex(Mutex& m) {pthread_mutex_init(&m, 0);}     // init mutex
void lock(Mutex& m) {pthread_mutex_lock(&m);}              // wait for mutex
void release(Mutex& m) {pthread_mutex_unlock(&m);}         // release mutex
void destroy_mutex(Mutex& m) {pthread_mutex_destroy(&m);}  // destroy mutex

class Semaphore {
public:
    Semaphore() {sem=-1;}
    void init(int n) {
        assert(n>=0);
        assert(sem==-1);
        pthread_cond_init(&cv, 0);
        pthread_mutex_init(&mutex, 0);
        sem=n;
    }
    void destroy() {
        assert(sem>=0);
        pthread_mutex_destroy(&mutex);
        pthread_cond_destroy(&cv);
    }
    int wait() {
        assert(sem>=0);
        pthread_mutex_lock(&mutex);
        int r=0;
        if (sem==0) r=pthread_cond_wait(&cv, &mutex);
        assert(sem>0);
        --sem;
        pthread_mutex_unlock(&mutex);
        return r;
    }
    void signal() {
        assert(sem>=0);
        pthread_mutex_lock(&mutex);
        ++sem;
        pthread_cond_signal(&cv);
        pthread_mutex_unlock(&mutex);
    }
private:
    pthread_cond_t cv;  // to signal FINISHED
    pthread_mutex_t mutex; // protects cv
    int sem;  // semaphore count
};

// A CompressJob is a queue of blocks to compress and write to the archive.
// Each block cycles through states EMPTY, FILLING, FULL, COMPRESSING,
// COMPRESSED, WRITING. The main thread waits for EMPTY buffers and
// fills them. A set of compressThreads waits for FULL threads and compresses
// them. A writeThread waits for COMPRESSED buffers at the front
// of the queue and writes and removes them.

// Buffer queue element
struct CJ {
    enum {EMPTY, FULL, COMPRESSING, COMPRESSED, WRITING} state;
    StringBuffer in;       // uncompressed input
    StringBuffer out;      // compressed output
    string filename;       // to write in filename field
    string comment;        // if "" use default
    string method;         // compression level or "" to mark end of data
    Semaphore full;        // 1 if in is FULL of data ready to compress
    Semaphore compressed;  // 1 if out contains COMPRESSED data
    Writer* out2;
    CJ(): state(EMPTY) {}
};

// Instructions to a compression job
class CompressJob {
public:
    Mutex mutex;           // protects state changes
private:
    int job;               // number of jobs
    CJ* q;                 // buffer queue
    unsigned qsize;        // number of elements in q
    int front;             // next to remove from queue
    Semaphore empty;       // number of empty buffers ready to fill
    Semaphore compressors; // number of compressors available to run
public:
    friend ThreadReturn compressThread(void* arg);
    friend ThreadReturn writeThread(void* arg);
    CompressJob(int threads, int buffers): job(0), q(0), qsize(buffers), front(0) {
        q=new CJ[buffers];
        if (!q) throw std::bad_alloc();
        init_mutex(mutex);
        empty.init(buffers);
        compressors.init(threads);
        for (int i=0; i<buffers; ++i) {
            q[i].full.init(0);
            q[i].compressed.init(0);
        }
    }
    ~CompressJob() {
        for (int i=qsize-1; i>=0; --i) {
            q[i].compressed.destroy();
            q[i].full.destroy();
        }
        compressors.destroy();
        empty.destroy();
        destroy_mutex(mutex);
        delete[] q;
    }      
    void write(StringBuffer& s, const char* filename, string method, Writer& f, const char* comment=0);
    vector<int> csize;  // compressed block sizes
};

// Write s at the back of the queue. Signal end of input with method=""
void CompressJob::write(StringBuffer& s, const char* fn, string method, Writer& f, const char* comment) {
    for (unsigned k=(method=="")?qsize:1; k>0; --k) {
        empty.wait();
        lock(mutex);
        unsigned i, j;
        for (i=0; i<qsize; ++i) {
            if (q[j=(i+front)%qsize].state==CJ::EMPTY) {
                q[j].filename=fn?fn:"";
                q[j].comment=comment?comment:"";
                q[j].method=method;
                q[j].in.resize(0);
                q[j].in.swap(s);
                q[j].out2 = &f;
                q[j].state=CJ::FULL;
                q[j].full.signal();
                break;
            }
        }
        release(mutex);
        assert(i<qsize);  // queue should not be full
    }
}

// Compress data in the background, one per buffer
ThreadReturn compressThread(void* arg) {
    CompressJob& job=*(CompressJob*)arg;
    int jobNumber=0;
    try {

        // Get job number = assigned position in queue
        lock(job.mutex);
        jobNumber=job.job++;
        assert(jobNumber>=0 && jobNumber<int(job.qsize));
        CJ& cj=job.q[jobNumber];
        release(job.mutex);

        // Work until done
        while (true) {
            cj.full.wait();
            lock(job.mutex);

            // Check for end of input
            if (cj.method=="") {
                cj.compressed.signal();
                release(job.mutex);
                return 0;
            }

            // Compress
            assert(cj.state==CJ::FULL);
            cj.state=CJ::COMPRESSING;
            release(job.mutex);
            job.compressors.wait();
            // libzpaq::compressBlock(&cj.in, &cj.out, cj.method.c_str(),
            //     cj.filename.c_str(), cj.comment=="" ? 0 : cj.comment.c_str());
            libzpaq::compressBlock(&cj.in, &cj.out, cj.method.c_str(), cj.filename.c_str(), 0);
            cj.in.resize(0);
            lock(job.mutex);
            cj.state=CJ::COMPRESSED;
            cj.compressed.signal();
            job.compressors.signal();
            release(job.mutex);
        }
    }
    catch (std::exception& e) {
        lock(job.mutex);
        fflush(stdout);
        fprintf(stderr, "job %d: %s\n", jobNumber+1, e.what());
        release(job.mutex);
        exit(1);
    }
    return 0;
}

inline int32_t btoi(StringBuffer& sb) {
    int32_t x = 0;
    for(int i = 0; i < 4; ++i) x = x * 256 + sb.get();
    return x;
}

inline void write_len(Writer& out, int32_t len) {
    out.put(len>>24); out.put(len>>16&255); out.put(len>>8&255); out.put(len&255);
}

// Write compressed data in the background
ThreadReturn writeThread(void* arg) {
    CompressJob& job=*(CompressJob*)arg;
    try {

        // work until done
        while (true) {

            // wait for something to write
            CJ& cj=job.q[job.front];  // no other threads move front
            cj.compressed.wait();

            // Quit if end of input
            lock(job.mutex);
            if (cj.method=="") {
                release(job.mutex);
                return 0;
            }

            // Write
            assert(cj.state==CJ::COMPRESSED);
            cj.state=CJ::WRITING;
            job.csize.push_back(cj.out.size());
            if (cj.out2 && cj.out.size()>0) {
                release(job.mutex);
                assert(cj.out.c_str());
                const char* p=cj.out.c_str();
                int64_t n=cj.out.size();
                const int64_t N=1<<30;
                write_len(*cj.out2, n);
                while (n>N) {
                    cj.out2->write(p, N);
                    p+=N;
                    n-=N;
                }
                cj.out2->write(p, n);
                lock(job.mutex);
            }
            cj.out.resize(0);
            cj.state=CJ::EMPTY;
            job.front=(job.front+1)%job.qsize;
            job.empty.signal();
            release(job.mutex);
        }
    }
    catch (std::exception& e) {
        fflush(stdout);
        fprintf(stderr, "zpaq exiting from writeThread: %s\n", e.what());
        exit(1);
    }
    return 0;
}

int numberOfProcessors() {
    int rc=0;  // result
    // Count lines of the form "processor\t: %d\n" in /proc/cpuinfo
    // where %d is 0, 1, 2,..., rc-1
    FILE *in=fopen("/proc/cpuinfo", "r");
    if (!in) return 1;
    std::string s;
    int c;
    while ((c=getc(in))!=EOF) {
        if (c>='A' && c<='Z') c+='a'-'A';  // convert to lowercase
        if (c>' ') s+=c;  // remove white space
        if (c=='\n') {  // end of line?
            if (s.size()>10 && s.substr(0, 10)=="processor:") {
                c=atoi(s.c_str()+10);
                if (c==rc) ++rc;
            }
            s="";
        }
    }
    fclose(in);
    if (rc<1) rc=1;
    if (sizeof(char*)==4 && rc>2) rc=2;
    return rc;
}

const float shreshold[3] = {0.6f, 0.45f, 0.25f};
const size_t buffer_size = 1 << 26;
const char METHOD[] = "56,230,0";

void get_qs_center(string& center, const vector<string>& sample) {
    int32_t cnt[128];
    for(size_t j = 0; j < center.size(); ++j) {
        memset(cnt, 0, sizeof cnt);
        for(size_t i = 0; i < center.size(); ++i) {
            ++cnt[sample[i][j]];
        }
        center[j] = max_element(cnt, cnt + 128) - cnt;
    }
}

inline int same_qs(const string& center, const string& s) {
    int cnt = 0;
    for(size_t i = 0; i < center.size(); ++i) {
        if(s[i] == center[i]) ++cnt;
    }
    return cnt;
}

fqcompressor::fqcompressor(int _threads): threads(_threads), id_prefix(0), id_cnt(0) {
    if(threads < 1) threads = numberOfProcessors();
    tid.resize(threads*2-1);
    job = new CompressJob(threads, tid.size());
    for (unsigned i=0; i<tid.size(); ++i) run(tid[i], compressThread, job);
    run(wid, writeThread, job);
}

void fqcompressor::end() {
    StringBuffer _;
    job->write(_, 0, "", _);  // signal end of input
    for (unsigned i=0; i<tid.size(); ++i) join(tid[i]);
    join(wid);
    delete job;
    write_len(id[0], id_cnt);
}

void fqcompressor::split_qs(const string& center, const float s, Writer& out) {
    int limit = floor(s*center.size()), cnt = 0;
    StringBuffer sb;
    write_len(sb, 0);
    for(size_t i = 0; i < qs_raw.size(); ++i) {
        if(same_qs(center, qs_raw[i]) <= limit) {
            qs_raw[cnt++].swap(qs_raw[i]);
            sb.put('\n');
        }
        else {
            sb.write(qs_raw[i].c_str(), center.size());
            sb.put('\n');
            if(sb.size() > buffer_size) {
                job->write(sb, "", METHOD, out);
                write_len(sb, i+1);
            }
        }
    }
    if(sb.size() > 4) {
        job->write(sb, "", METHOD, out);
    }
    qs_raw.resize(cnt);
}

void fqcompressor::qs_compress() {
    StringBuffer sb;
    write_len(qs[3], qs_raw[0].size());
    for(int i = 0; i < 3; ++i) {
        vector<string> sample(qs_raw.begin(), qs_raw.size() < 10000 ? qs_raw.end() : qs_raw.begin()+10000);
        string center(sample[0].size(), 0);
        get_qs_center(center, sample);
        write_len(qs[i], qs_raw.size());
        split_qs(center, shreshold[i], qs[i]);
    }
    write_len(qs[3], qs_raw.size());
    write_len(sb, 0);
    for(size_t i = 0; i < qs_raw.size(); ++i) {
        sb.write(qs_raw[i].c_str(), qs_raw[0].size());
        sb.put('\n');
        if(sb.size() > buffer_size) {
            job->write(sb, "", METHOD, qs[3]);
            write_len(sb, i+1);
        }
    }
    if(sb.size() > 4) {
        job->write(sb, "", METHOD, qs[3]);
    }
    qs_raw.resize(0);
}

void fqcompressor::id_add(const char* s) {
    ++id_cnt;
    if(id_prefix == 0) {
        id_prefix = strlen(s) - 1;
        while(s[id_prefix] != 'C') --id_prefix;
        id[0].write(s, id_prefix);
        id[0].put(0);
        id[0].put(s[strlen(s)-1]);
        for(int i = 0; i < 3; ++i) write_len(id_buffer[i], 0);
    }
    for(int i = 1; i < 4; ++i) {
        id_buffer[0].put(s[id_prefix+i]);
        id_buffer[1].put(s[id_prefix+4+i]);
    }
    int32_t res = 0;
    for(int i = id_prefix+9; s[i] != '/'; ++i) {
        res = res * 10 + s[i] - '0';
    }
    write_len(id_buffer[2], res);
    for(int i = 0; i < 3; ++i) {
        if(id_buffer[i].size() > buffer_size) {
            job->write(id_buffer[i], "", METHOD, id[i]);
            write_len(id_buffer[i], id_cnt);
        }
    }
}

struct Block {
    StringBuffer* in;
    vector<string>& out;
    enum {READY, WORKING, GOOD, BAD} state;
    int id;
    Block(vector<string>& _out, int _id = -1): state(READY), out(_out), id(_id), in(new StringBuffer) {}

    void operator = (const Block& b) {
        in = b.in; out = b.out; state = b.state; id = b.id;
    }
};

struct ExtractJob {         // list of jobs
    Mutex mutex;              // protects state
    int job;                  // number of jobs started
    vector<Block> block;      // list of data blocks to extract
    int id_prefix;
    char pairend;
    size_t q_len;
    ExtractJob(): job(0) {
        init_mutex(mutex);
    }
    ~ExtractJob() {
        destroy_mutex(mutex);
    }
};

// Decompress blocks in a job until none are READY
ThreadReturn decompressThread(void* arg) {
    ExtractJob& job=*(ExtractJob*)arg;
    char pairend = job.pairend;
    int id_prefix = job.id_prefix;
    
    int jobNumber=0;

    // Get job number
    lock(job.mutex);
    jobNumber=++job.job;
    release(job.mutex);

    // Look for next READY job.
    int next=0;  // current job
    while (true) {
        lock(job.mutex);
        for (unsigned i=0; i<=job.block.size(); ++i) {
            unsigned k=i+next;
            if (k>=job.block.size()) k-=job.block.size();
            if (i==job.block.size()) {  // no more jobs?
                release(job.mutex);
                return 0;
            }
            Block& b=job.block[k];
            if (b.state==Block::READY) {
                b.state=Block::WORKING;
                release(job.mutex);
                next=k;
                break;
            }
        }
        Block& b=job.block[next];

        // Decompress
        StringBuffer out;
        decompress(b.in, &out);
        delete b.in;
        b.in = NULL;

        // Write
        if(b.id == -1) {
            size_t os = btoi(out), q_len = job.q_len;
            char c;
            while((c=out.get()) != EOF) {
                if(c != '\n') {
                    b.out[os].resize(q_len);
                    b.out[os][0] = c;
                    for(size_t i = 1; i < q_len; ++i) b.out[os][i] = out.get();
                    out.get();
                }
                ++os;
            }
            continue;
        }

        size_t os = btoi(out);
        char buf[11];
        if(b.id == 2) {
            while(out.remaining()) {
                int32_t x = btoi(out), tot = 0;
                while(x) {
                    buf[tot++] = x % 10 + '0';
                    x /= 10;
                }
                if(tot == 0) tot = 1, buf[0] = '0';
                b.out[os].resize(id_prefix+tot+11);
                for(int i = 1; i <= tot; ++i) {
                    b.out[os][id_prefix+8+i] = buf[tot-i];
                }
                for(int i = 0; i < id_prefix; ++i) {
                    b.out[os][i] = b.out[0][i];
                }
                b.out[os][id_prefix] = 'C';
                b.out[os][id_prefix+4] = 'R';
                b.out[os][id_prefix+8] = '_';
                b.out[os][b.out[os].size()-2] = '/';
                b.out[os][b.out[os].size()-1] = pairend;
                ++os;
            }
        }
        else {
            while(out.remaining()) {
                for(int i = 1; i < 4; ++i) {
                    b.out[os][id_prefix+i+b.id*4] = out.get();
                }
                ++os;
            }
        }
        
    } // end while true

    // Last block
    return 0;
}

void fqcompressor::end_id_add() {
    for(int i = 0; i < 3; ++i) {
        if(id_buffer[i].size() > 4) {
            job->write(id_buffer[i], "", METHOD, id[i]);
        }
    }
}

fqdecompressor::fqdecompressor(int _threads): threads(_threads), job(new ExtractJob) {
    if(threads < 1) threads = numberOfProcessors();
}

void fqdecompressor::qs_add(StringBuffer& q, int i) {
    if(i == 3) job->q_len = btoi(q);
    qs_raw[i].resize(btoi(q));
    int32_t len;
    do {
        len = btoi(q);
        job->block.push_back(Block(qs_raw[i]));
        job->block[job->block.size()-1].in->read_buffer(q, len);
    } while(q.remaining());
}

void fqdecompressor::start() {
    for(size_t i = 0; i < id_raw.size(); ++i) {
        id_raw[i].resize(job->id_prefix+16);
    }
    tid.resize(threads);
    for (unsigned i=0; i<tid.size(); ++i) run(tid[i], decompressThread, job);
}

void fqdecompressor::end() {
    for (unsigned i=0; i<tid.size(); ++i) join(tid[i]);
    delete job;
}

void fqdecompressor::get_qs(vector<string>& ans) {
    for(size_t i = 2; ; --i) {
        size_t a, b;
        for(a = 0, b = 0; a < qs_raw[i].size(); ++a) {
            if(qs_raw[i][a].size() == 0) qs_raw[i][a].swap(qs_raw[i+1][b++]);
        }
        if(i == 0) break;
    }
    ans.swap(qs_raw[0]);
}

void fqdecompressor::id_add(StringBuffer& sb, int i) {
    StringBuffer temp;
    temp.swap(sb);
    size_t _ = 0;
    if(i == 0) {
        if(id_raw.size() == 0) id_raw.resize(1);
        char c;
        while((c = temp.get())) {
            id_raw[0].push_back(c);
        }
        job->pairend = temp.get();
        job->id_prefix = id_raw[0].size();
        _ = 4;
    }
    int32_t len;
    do {
        len = btoi(temp);
        job->block.push_back(Block(id_raw, i));
        job->block[job->block.size()-1].in->read_buffer(temp, len);
    } while(temp.remaining() != _);
    if(i == 0) id_raw.resize(btoi(temp));
}

}
