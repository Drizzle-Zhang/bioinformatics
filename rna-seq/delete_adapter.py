#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: delete_adapter.py
# @time: 2018/10/23 20:52

import sys


def delete_adapter(file1, num_del):
    with open(file1, 'r') as r_f:
        with open(file1[:-6] + '.fq', 'w') as w_f:
            for line in r_f:
                if line[0] != '@' and line[0] != '+':
                    post_line = line[int(num_del):]
                    w_f.write(post_line)
                else:
                    w_f.write(line)

    return


if __name__ == '__main__':
    delete_adapter(sys.argv[1], sys.argv[2])
