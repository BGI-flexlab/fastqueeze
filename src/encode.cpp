// rangemapper.cpp : Defines the entry point for the console application.
// (C) 2010, Andrew Polar under GPL ver. 3.
// Released  Oct, 2010.
//
//   LICENSE
//
//   This program is free software; you can redistribute it and/or
//   modify it under the terms of the GNU General Public License as
//   published by the Free Software Foundation; either version 3 of
//   the License, or (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful, but
//   WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//   General Public License for more details at
//   Visit <http://www.gnu.org/copyleft/gpl.html>.
//
//	This is experimental arithmetic or range coder. It is designed as close to
//  theoretical range narrowing concept as possible. It is inspired by
//  Matt Mahoney binary arithmetic coder used in FPAQ and ZPAQ. His concept
//  of binary coding is simply generalized on any size alphabet.
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <string.h>

///////Data generation functions///////////////////////////////
template <class T>
double calculate_entropy(T* data, int data_size) {
    int min, max, counter, buffer_size;
    int* buffer;
    double entropy;
    double log2 = log(2.0);
    double prob;

    min = data[0];
    max = data[0];
    for (counter=0; counter<data_size; ++counter) {
        if (data[counter] < min)
            min = data[counter];
        if (data[counter] > max)
            max = data[counter];
    }

    buffer_size = max - min + 1;
    buffer = (int*)malloc(buffer_size * sizeof(int));
    memset(buffer, 0x00, buffer_size * sizeof(int));
    for (counter=0; counter<data_size; ++counter) {
        ++buffer[data[counter]-min];
    }

    entropy = 0.0;
    for (counter=0; counter<buffer_size; ++counter) {
        if (buffer[counter] > 0) {
            prob = (double)(buffer[counter]);
            prob /= (double)(data_size);
            entropy += log(prob)/log2*prob;
        }
    }
    entropy *= -1.0;

    if (buffer)
        free(buffer);
    return  entropy;
}

static __inline double rounds(double x) {
    if ((x - floor(x)) >= 0.5)
        return ceil(x);
    else
        return floor(x);
}

double makeData(int* data, int data_size, int min, int max, int redundancy_factor) {
    int counter, cnt, high, low;
    double buf;

    if (redundancy_factor <= 1)
        redundancy_factor = 1;

    if (max <= min)
        max = min + 1;

    srand((unsigned)(time(0)));
    for (counter=0; counter<data_size; counter++) {
        buf = 0;
        for (cnt=0; cnt<redundancy_factor; cnt++) {
            buf += (double)(rand());
        }
        data[counter] = ((int)(buf))/redundancy_factor;
    }
    low  = data[0];
    high = data[0];
    for (counter=0; counter<data_size; counter++) {
        if (data[counter] > high)
            high = data[counter];
        if (data[counter] < low)
            low = data[counter];
    }
    for (counter=0; counter<data_size; counter++) {
        buf = (double)(data[counter] - low);
        buf /= (double)(high - low);
        buf *= (double)(max - min);
        buf = rounds(buf);
        data[counter] = (int)(buf) + min;
    }
    return calculate_entropy(data, data_size);
}
//end of data generation functions

//global data buffer for round trip
unsigned char* g_buffer = 0;
int current_byte = 0;

//input/output functions
static __inline void writeByte(unsigned char byte) {
    g_buffer[current_byte++] = byte;
}

static __inline unsigned char readByte() {
    return g_buffer[current_byte++];
}

static __inline int findInterval(int* cm, int size, int point) {
    int index = -1;
    int left  = 0;
    int right = size-2;
    int cnt = 0;
    while (true) {
        int mid = (right + left)>>1;
        if (point >= cm[mid] && point < cm[mid+1]) {
            index = mid;
            break;
        }
        if (point >= cm[mid+1]) left = mid + 1;
        else right = mid;
        if (cnt++ >= size) break;
    }
    return index;
}

//This is the main class RangeMapper. It can be used in adaptive data
//processing. It processes ranges that passed to encoder and decoder.
//In adaptive coder the ranges can be computed dynamically and depend
//on the context.
//Class is written for generic case: alphabets from 2 up to 10000 can be
//processed without changes in code. In case of special data coder can
//be optimized and run at least twice faster.
class RangeMapper {
public:
    RangeMapper(int range_size) {
        LOW = 0;
        MID = 0;
        RANGE_SIZE = range_size;
        RANGE_LIMIT = ((unsigned long long)(1) << RANGE_SIZE);
        BYTES_IN_BUFFER = (64 - RANGE_SIZE) / 8;
        SHIFT = (BYTES_IN_BUFFER - 1) * 8;
        MASK  = ((long long)(1)<<(BYTES_IN_BUFFER * 8)) - 1;
        HIGH = MASK;
    }
    ~RangeMapper() {}
    void encodeRange(int cmin, int cmax);
    void decodeRange(int cmin, int cmax);
    int getMidPoint();
    void flush();
    void init();
private:
    void updateModel(int cmin, int cmax);
    unsigned long long LOW, HIGH, MID, RANGE_LIMIT, MASK;
    unsigned char RANGE_SIZE, SHIFT, BYTES_IN_BUFFER;
};

void RangeMapper::updateModel(int cmin, int cmax) {
    unsigned long long range = HIGH - LOW;
    HIGH = LOW + ((range * cmax) >> RANGE_SIZE);
    LOW += ((range * cmin) >> RANGE_SIZE) + 1;
}

int RangeMapper::getMidPoint() {
    return (int)(((MID - LOW) << RANGE_SIZE) / (HIGH - LOW));
}

void RangeMapper::encodeRange(int cmin, int cmax) {
    updateModel(cmin, cmax);
    if ((HIGH - LOW) < RANGE_LIMIT) HIGH = LOW; //preventing narrowing range
    while (((LOW ^ HIGH) >> SHIFT) == 0) {
        writeByte((unsigned char)(LOW >> SHIFT));
        LOW <<= 8;
        HIGH = (HIGH << 8) | 0xff;
    }
    HIGH &= MASK;
    LOW  &= MASK;
}

void RangeMapper::decodeRange(int cmin, int cmax) {
    updateModel(cmin, cmax);
    if ((HIGH - LOW) < RANGE_LIMIT)  HIGH = LOW;
    while (((LOW ^ HIGH) >> SHIFT) == 0) {
        LOW <<= 8;
        HIGH = (HIGH << 8) | 0xff;
        MID =  (MID << 8)  | readByte();
    }
    HIGH &= MASK;
    LOW  &= MASK;
    MID  &= MASK;
}

void RangeMapper::flush() {
    LOW += 1;
    for (int i=0; i<BYTES_IN_BUFFER; i++) {
        writeByte((unsigned char)(LOW >> SHIFT));
        LOW <<= 8;
    }
}

void RangeMapper::init() {
    for (int i=0; i<BYTES_IN_BUFFER; ++i) {
        MID = (MID << 8) + readByte();
    }
}

void makeRanges(int* data, int data_size, int* cm, int alphabet_size, int PRECISION) {
    //we make ranges for data
    int* freq = (int*)malloc(alphabet_size * sizeof(int));
    memset(freq, 0x00, alphabet_size * sizeof(int));
    for (int i=0; i<data_size; ++i) {
        ++freq[data[i]];
    }

    cm[0] = 0;
    for (int i=0; i<alphabet_size; ++i) {
        cm[i+1] = cm[i] + freq[i];
    }

    int total = cm[alphabet_size];
    int upper_limit = (1<<PRECISION) - 2;
    for (int i=0; i<alphabet_size + 1; ++i) {
        cm[i] = (int)((long long)(cm[i]) * (long long)(upper_limit) / (long long)(total));
    }
    cm[alphabet_size+1] = (1<<PRECISION) - 1;
    //ranges are ready

    //correction of ranges
    for (int i=0; i<alphabet_size; ++i) {
        if (cm[i+1] <= cm[i]) cm[i+1] = cm[i] + 1;
    }
    for (int i=alphabet_size; i>=0; --i) {
        if (cm[i] >= cm[i+1]) cm[i] = cm[i+1] - 1;
    }
    //end of correction
    if (freq) free(freq);
}

void makeLookupTable(int* cm, int size, short* lookup) {
    for (int i=0; i<size-1; ++i) {
        for (int j=cm[i]; j<cm[i+1]; ++j) {
            lookup[j] = i;
        }
    }
}

int main() {
    //parameters that can be changed for testing
    //the RANGE_SIZE_IN_BITS depends on your alphabet size
    const int alphabet_size = 256; //Please notice the size of the alphabet.
    int data_size = 4000000;
    int RANGE_SIZE_IN_BITS = 14;
    ////////////////////////////////

    //we make original data
    printf("Wait, the data is being made ...\n");
    int* data = (int*)malloc(data_size * sizeof(int));
    double entropy = makeData(data, data_size, 0, alphabet_size - 1, 10);
    int expected_size = (int)(entropy * (double)(data_size) / 8.0);
    //data is made

    //make ranges for data
    const int cm_size = alphabet_size + 2;
    int cm[cm_size];
    makeRanges(data, data_size, cm, alphabet_size, RANGE_SIZE_IN_BITS);
    //ranges are completed


    //encode
    g_buffer = (unsigned char*)malloc(expected_size * 2);
    clock_t start_encoding = clock();
    RangeMapper* rm_encode = new RangeMapper(RANGE_SIZE_IN_BITS);
    for (int i=0; i<data_size; ++i) {
        rm_encode->encodeRange(cm[data[i]], cm[data[i]+1]);
    }
    std::cout << data[30];
    rm_encode->encodeRange(cm[alphabet_size], cm[alphabet_size+1]); //end of data marker
    rm_encode->flush();
    delete rm_encode;

    clock_t end_encoding = clock();
    printf("Time for encoding %2.3f sec.\n", (double)(end_encoding - start_encoding)/CLOCKS_PER_SEC);
    int actual_size = current_byte;
    //end encoding

    //decode
    clock_t start_decoding = clock();
    short* lookup = (short*)malloc(((1<<RANGE_SIZE_IN_BITS) + 1) * sizeof(short));
    makeLookupTable(cm, alphabet_size + 2, lookup);
    bool isOK = true;
    RangeMapper* rm_decode = new RangeMapper(RANGE_SIZE_IN_BITS);
    current_byte = 0;
    rm_decode->init();
    int i=0;
    while (true) {
        int midpoint = rm_decode->getMidPoint();
        //next is binary search algorithm that does not need having lookup array
        //int index = findInterval(cm, alphabet_size + 2, midpoint);
        //this is lookup table that expedites execution, either of these functions works
        int index = lookup[midpoint]; //midpoint presumed being within correct boundaries
        if (index == alphabet_size) break; //end of data marker
        if (index != data[i++]) {
            printf("Data mismatch %d element %d\n", i, index);
            isOK = false;
            break;
        }
        rm_decode->decodeRange(cm[index], cm[index+1]);
    }
    if (i != data_size) isOK = false;
    delete rm_decode;
    if (lookup) free(lookup);
    clock_t end_decoding = clock();
    printf("Time for decoding %2.3f sec.\n\n", (double)(end_decoding - start_decoding)/CLOCKS_PER_SEC);
    //end decoding

    printf("Expected size %d\n", expected_size);
    printf("Actual size   %d, \n\n", actual_size);

    if (isOK) printf("Round trip is OK\n");
    else printf("Data mismatch\n");

    if (g_buffer) free(g_buffer);
    if (data) free(data);

    return 0;
}
