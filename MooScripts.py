#!/Anaconda/python
__author__ = 'Keagan Moo'
def euclidd(list1, list2):
    iter = len(list1)
    sum = 0
    for i in range(iter):
        gap = (list1[i] - list2[i]) ** 2
        sum = sum + gap
    return sum ** (0.5)

def weightt(t, N, wgtp, pThet, pThetp, pTurb, piTurb):
    if t == 0:
        return 1
    elif t > 0:
        sum = 0
        for i in range(N):
            bit = wgtp * pTurb(pThetp, pThet)
            sum = sum + bit
        return piTurb(pThet)/sum
    else:
        return -666