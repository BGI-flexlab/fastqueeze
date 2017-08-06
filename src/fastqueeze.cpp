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
typedef pthread_t ThreadID;                                // job ID type
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
    void write(StringBuffer& s, const char* filename, string method, Writer *f, const char* comment=0);
    vector<int> csize;  // compressed block sizes
};

// Write s at the back of the queue. Signal end of input with method=""
void CompressJob::write(StringBuffer& s, const char* fn, string method, Writer *f, const char* comment) {
    for (unsigned k=(method=="")?qsize:1; k>0; --k) {
        empty.wait();
        lock(mutex);
        unsigned i, j;
        for (i=0; i<qsize; ++i) {
            if (q[j=(i+front)%qsize].state==CJ::EMPTY) {
                q[j].filename=fn?fn:"";
                q[j].comment=comment?comment:"jDC\x01";
                q[j].method=method;
                q[j].in.resize(0);
                q[j].in.swap(s);
                q[j].out2 = f;
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

void write_len(Writer& out, int32_t len) {
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

vector<string> id;
vector<string> qs;
const float shreshold[3] = {0.6f, 0.45f, 0.25f};
const size_t buffer_size = 1 << 26;
const char METHOD[] = "56,230,0";
int threads;
CompressJob* job;
vector<ThreadID> tid;
ThreadID wid;

void c_init(int _threads) {
    threads = _threads;
    id.clear();
    qs.clear();

    if(threads < 1) threads = numberOfProcessors();
    tid.resize(threads*2-1);
    job = new CompressJob(threads, tid.size());
    for (unsigned i=0; i<tid.size(); ++i) run(tid[i], compressThread, job);
    run(wid, writeThread, job);
}

void c_destroy() {
    StringBuffer _;
    job->write(_, 0, "", NULL);  // signal end of input
    for (unsigned i=0; i<tid.size(); ++i) join(tid[i]);
    join(wid);
}

void d_init(int _threads) {
    threads = _threads;
    if(threads < 1) threads = numberOfProcessors();
}

void d_destroy() {

}

void id_add(const char* s) {
    id.push_back(s);
}

void qs_add(const char* s) {
    qs.push_back(s);
}

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

void split_qs(const string& center, vector<string>& qs, const float s, Writer* out) {
    int limit = floor(s*center.size()), cnt = 0;
    StringBuffer sb;
    write_len(sb, 0);
    for(size_t i = 0; i < qs.size(); ++i) {
        if(same_qs(center, qs[i]) <= limit) {
            qs[cnt++].swap(qs[i]);
            sb.put('\n');
        }
        else {
            sb.write(qs[i].c_str(), center.size());
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
    qs.resize(cnt);
}

void qs_compress(vector<StringBuffer*>& q) {
    StringBuffer sb;
    write_len(*q[3], qs[0].size());
    for(int i = 0; i < 3; ++i) {
        vector<string> sample(qs.begin(), qs.size() < 10000 ? qs.end() : qs.begin()+10000);
        string center(sample[0].size(), 0);
        get_qs_center(center, sample);
        write_len(*q[3], qs.size());
        split_qs(center, qs, shreshold[i], q[i]);
    }
    write_len(*q[3], qs.size());
    write_len(sb, 0);
    for(size_t i = 0; i < qs.size(); ++i) {
        sb.write(qs[i].c_str(), qs[0].size());
        sb.put('\n');
        if(sb.size() > buffer_size) {
            job->write(sb, "", METHOD, q[3]);
            write_len(sb, i+1);
        }
    }
    if(sb.size() > 4) {
        job->write(sb, "", METHOD, q[3]);
    }
}

size_t q_len;

struct Block {
    StringBuffer* in;
    vector<string>& out;
    enum {READY, WORKING, GOOD, BAD} state;
    Block(vector<string>& _out): state(READY), out(_out) {
        in = new StringBuffer;
        in->resize(0);
    }

    void operator = (const Block& b) {
        in = b.in; out = b.out; state = b.state;
    }
};

struct ExtractJob {         // list of jobs
    Mutex mutex;              // protects state
    int job;                  // number of jobs started
    vector<Block> block;      // list of data blocks to extract
    ExtractJob(): job(0) {
        init_mutex(mutex);
    }
    ~ExtractJob() {
        destroy_mutex(mutex);
    }
};

inline int32_t btoi(StringBuffer& sb) {
    int32_t x = 0;
    for(int i = 0; i < 4; ++i) x = x * 256 + sb.get();
    return x;
}

// Decompress blocks in a job until none are READY
ThreadReturn decompressThread(void* arg) {
    ExtractJob& job=*(ExtractJob*)arg;
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
        size_t os = btoi(out);
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
        
    } // end while true

    // Last block
    return 0;
}

void qs_decompress(vector<StringBuffer*>& q, vector<string>& ans) {
    ExtractJob job;
    vector<string> res[4];
    q_len = btoi(*q[3]);
    for(int i = 0; i < 4; ++i) res[i].resize(btoi(*q[3]));
    for(int i = 0; i < 4; ++i) {
        StringBuffer raw_q;
        int32_t len;
        do {
            len = btoi(*q[i]);
            job.block.push_back(Block(res[i]));
            job.block[job.block.size()-1].in->read_buffer(*q[i], len);
        } while(q[i]->remaining());
    }
    tid.resize(threads);
    for (unsigned i=0; i<tid.size(); ++i) run(tid[i], decompressThread, &job);
    for (unsigned i=0; i<tid.size(); ++i) join(tid[i]);
    for(size_t i = 2; ; --i) {
        for(size_t a = 0, b = 0; a < res[i].size(); ++a) {
            if(res[i][a].size() == 0) res[i][a].swap(res[i+1][b++]);
        }
        if(i == 0) break;
    }
    ans.swap(res[0]);
}

}
