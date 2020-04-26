#
# Yusuf Bham
# Period 3
#

def lev_same_len(s1, s2):
    if len(s1) != len(s2):
        raise ValueError(f"Called on two strings of non-equal distance. {s1} {s2}")

    return sum(i != j for (i, j) in zip(s1, s2))


def lev_recur(orig_s1, orig_s2):
    def recur(s1, s2):
        if len(s1) == len(s2):
            return lev_same_len(s1, s2)

        return min(recur(s1[:ind] + "-" + s1[ind:], s2) for ind in range(len(s2)))

    return recur(*sorted((orig_s1, orig_s2), key=len))


def lev_dp(s1, s2):
    # s1 x s2
    mat = [[0] * (len(s2) + 1) for _ in range(len(s1) + 1)]

    for ind in range(len(s1) + 1):
        mat[ind][0] = ind

    for ind in range(len(s2) + 1):
        mat[0][ind] = ind

    for ind_i, i in enumerate(mat):
        for ind_j, j in enumerate(i):
            if ind_i == 0 or ind_j == 0:
                continue

            cost = 0 if s1[ind_i - 1] == s2[ind_j - 1] else 1

            mat[ind_i][ind_j] = min(
                mat[ind_i - 1][ind_j    ] + 1,
                mat[ind_i    ][ind_j - 1] + 1,
                mat[ind_i - 1][ind_j - 1] + cost,
            )

    return mat[len(s1) - 1][len(s2) - 1]


def main():
    w1, w2 = input("w1: "), input("w2: ")
    print(lev_dp(w1, w2))


if __name__ == "__main__":
    main()
