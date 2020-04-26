#
# Yusuf Bham
# Period 3
#

from pprint import pprint as e

def lev_same_len(s1, s2):
    if len(s1) != len(s2):
        raise ValueError(f"Called on two strings of non-equal distance. {s1} {s2}")

    return sum(i != j for (i, j) in zip(s1, s2))

def lev_one_off(s1, s2):
    if abs(len(s1) - len(s2)) != 1:
        raise ValueError("Called on two strings with length not differing by 1.")

    small, big = sorted((s1, s2), key=len)

    minimum = float("inf")

    for ind in range(len(big)):
        gapped = small[:ind] + "-" + small[ind:]

        dist = lev_same_len(gapped, big)

        minimum = min(dist, minimum)

    return minimum


def lev_recur(orig_s1, orig_s2):

    small, big = sorted((orig_s1, orig_s2), key=len)

    def recur(s1, s2):
        if len(s1) == len(s2):
            return lev_same_len(s1, s2)

        minimum = float("inf")

        for ind in range(len(s2)):
            gapped = s1[:ind] + "-" + s1[ind:]
            dist = recur(gapped, s2)
            minimum = min(dist, minimum)

        return minimum

    return recur(small, big)


def main():
    # print("cat, bat:", lev_same_len("cat", "bat"))
    # print("dog, cat:", lev_same_len("dog", "cat"))

    # print("kitchen, kitten:", lev_one_off("kitchen", "kitten"))
    # print("meat, eat:", lev_one_off("meat", "eat"))
    # print("eat, eats:", lev_one_off("eat", "eats"))

    print("kitchen, kitten:", lev_recur("kitchen", "kitten"))
    print("meat, eat:", lev_recur("meat", "eat"))
    print("eat, eats:", lev_recur("eat", "eats"))
    print("anagrams, alarms:", lev_recur("anagrams", "alarms"))

    lev_dp("Saturday", "Sunday")

if __name__ == "__main__":
    main()
