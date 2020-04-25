#!/usr/bin/env python3

from collections import Counter

with open("words", "r") as f:
    words = [x.strip() for x in f.readlines()]

words = [x for x in words if len(x) == 6]

s_words = [''.join(sorted(x)) for x in words]

alpha_words = sorted(words)

count = Counter(s_words)

maximum = count.most_common(1)[0][0]

print(maximum)

print(next(w for w in alpha_words if ''.join(sorted(w)) == maximum))
