# coding=utf-8
# !/usr/bin/env python

#
# Sort large FASTQ files in a minimum amount of memory
# Forked from https://github.com/spiralout/external-sort
#

import os
import argparse


class FileSplitter(object):
    BLOCK_FILENAME_FORMAT = 'block_{0}.dat'

    def __init__(self, filename, line_unit):
        self.filename = filename
        self.line_unit = line_unit
        self.block_filenames = []

    def write_block(self, data, block_number):
        filename = self.BLOCK_FILENAME_FORMAT.format(block_number)
        block = open(filename, 'w')
        block.write(data)
        block.close()
        self.block_filenames.append(filename)

    def get_block_filenames(self):
        return self.block_filenames

    def split(self, block_size, sort_key=None):
        block = open(self.filename, 'r')
        i = 0

        while True:
            lines = self.read_src(block, block_size)

            if not lines:
                break

            if sort_key is None:
                lines.sort()
            else:
                lines.sort(key=sort_key)

            self.write_block(''.join(lines), i)
            i += 1

    def read_src(self, srcfile, block_size):
        i = 0
        lines = []

        while i < block_size:
            read = ''

            for j in range(self.line_unit):
                read += srcfile.readline()
            if not read:
                break

            lines.append(read)
            i += len(read)

        return lines

    def cleanup(self):
        map(lambda f: os.remove(f), self.block_filenames)


class NWayMerge(object):
    def __init__(self):
        self.sort_index = []
        self.init = True

    def select(self, choice):
        if self.init:
            self.init = False
            self.sort_index = sorted(choice.items(), key=lambda x: x[1], reverse=True)
            return self.sort_index.pop()[0]

        else:
            if choice[1]:
                for i in range(len(self.sort_index)):
                    if choice[1] > self.sort_index[i][1]:
                        self.sort_index.insert(i, choice)
                        return self.sort_index.pop()[0]
                return choice[0]

            else:
                return self.sort_index.pop()[0]


class FilesArray(object):
    def __init__(self, files, line_unit):
        self.files = files
        self.line_unit = line_unit
        self.empty = 0
        self.vacancy = 0
        self.num_buffers = len(files)
        self.buffers = {i: '' for i in range(self.num_buffers)}
        self.init = True

    def get_alter(self):
        if self.init:
            self.init = False
            return self.buffers
        else:
            return self.vacancy, self.buffers[self.vacancy]

    def refresh(self):
        if self.init:
            for i in range(self.num_buffers):
                for j in range(self.line_unit):
                    self.buffers[i] += self.files[i].readline()
        else:
            for j in range(self.line_unit):
                self.buffers[self.vacancy] += self.files[self.vacancy].readline()

            if not self.buffers[self.vacancy]:
                self.empty += 1

        if self.empty == self.num_buffers:
            return False

        return True

    def unshift(self, index):
        value = self.buffers[index]
        self.buffers[index] = ''
        self.vacancy = index

        return value


class FileMerger(object):
    def __init__(self, merge_strategy):
        self.merge_strategy = merge_strategy

    def merge(self, filenames, outfilename, buffer_size, line_unit):
        outfile = open(outfilename, 'w', buffer_size)
        buffers = FilesArray(self.get_file_handles(filenames, buffer_size), line_unit)

        while buffers.refresh():
            min_index = self.merge_strategy.select(buffers.get_alter())
            outfile.write(buffers.unshift(min_index))

    @staticmethod
    def get_file_handles(filenames, buffer_size):
        files = {}

        for i in range(len(filenames)):
            files[i] = open(filenames[i], 'r', buffer_size)

        return files


class ExternalSort(object):
    def __init__(self, block_size, line_unit):
        self.block_size = block_size
        self.line_unit = line_unit

    def sort(self, filename, sort_key=None):
        num_blocks = self.get_number_blocks(filename)
        splitter = FileSplitter(filename, self.line_unit)
        splitter.split(self.block_size, sort_key)

        merger = FileMerger(NWayMerge())
        buffer_size = int(self.block_size / (num_blocks + 1) / 0.4)
        merger.merge(splitter.get_block_filenames(), filename + '.out', buffer_size, self.line_unit)

        splitter.cleanup()

    def get_number_blocks(self, filename):
        return os.stat(filename).st_size / self.block_size


def parse_memory(string):
    if string[-1].lower() == 'k':
        return int(string[:-1]) * 1024 * 0.4
    elif string[-1].lower() == 'm':
        return int(string[:-1]) * 1024 * 1024 * 0.4
    elif string[-1].lower() == 'g':
        return int(string[:-1]) * 1024 * 1024 * 1024 * 0.4
    else:
        return int(string) * 0.4


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m',
                        '--mem',
                        help='amount of memory to use for sorting [100M]',
                        default='500M')
    parser.add_argument('-l',
                        '--lineunit',
                        help='number of lines processed as a unit, [4] for FQ file, 1 for regular file',
                        default='4')
    parser.add_argument('filename',
                        metavar='<filename>',
                        nargs='+',
                        help='name of file to sort')
    args = parser.parse_args()

    sorter = ExternalSort(parse_memory(args.mem), int(args.lineunit))
    
    for i in args.filename:
        sorter.sort(i)


if __name__ == '__main__':
    main()
