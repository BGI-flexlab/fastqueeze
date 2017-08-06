#ifndef FASTQUEEZE_H
#define FASTQUEEZE_H

#include "libzpaq.h"
#include <vector>

namespace fastqueeze {

void c_init(int _threads = 1);
void c_destroy();

void d_init(int _threads = 1);
void d_destroy();

void id_add(const char* s);
void id_compress();

void qs_add(const char* s);
void qs_compress(std::vector<libzpaq::StringBuffer*>& q);

void qs_decompress(std::vector<libzpaq::StringBuffer*>& q, std::vector<std::string>& ans);

}

#endif  // FASTQUEEZE_H
