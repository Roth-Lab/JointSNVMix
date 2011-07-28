'''
Created on 2011-07-27

@author: andrew
'''
import unittest


class Test(unittest.TestCase):


    def testName(self):
        pass
def test():
    # found these tests on this ticket. (thanks).
    # http://projects.scipy.org/scipy/ticket/956
    # these values were taken from R as a means to test the code in that ticket.
    tablist = [
            ([[100, 2], [1000, 5]], (2.505583993422285e-001, 1.300759363430016e-001)),
            ([[2, 100], [5, 1000]], (2.505583993422285e-001, 1.300759363430016e-001)),
            ([[2, 7], [8, 2]], (8.586235135736206e-002, 2.301413756522114e-002)),
            ([[5, 1], [10, 10]], (4.725646047336584e+000, 1.973244147157190e-001)),
            ([[5, 15], [20, 20]], (3.394396617440852e-001, 9.580440012477637e-002)),
            ([[5, 16], [20, 25]], (3.960558326183334e-001, 1.725864953812994e-001)),
            ([[10, 5], [10, 1]], (2.116112781158483e-001, 1.973244147157190e-001)),
            ([[10, 5], [10, 0]], (0.000000000000000e+000, 6.126482213438734e-002)),
            ([[5, 0], [1, 4]], ('inf', 4.761904761904762e-002)),
            ([[0, 5], [1, 4]], (0.000000000000000e+000, 1.000000000000000e+000)),
            ([[5, 1], [0, 4]], ('inf', 4.761904761904758e-002)),
            ([[0, 1], [3, 2]], (0.000000000000000e+000, 1.000000000000000e+000))
            ]

    
    for table, ab in tablist:
        p = pvalue(table[0][0], table[0][1], table[1][0], table[1][1])
        print table, p
        assert abs(p.two_tail - ab[1]) < 0.1, (table, ab, p)


def test_speed():
    cdef int i
    import time
    t = time.time()
    N = 5000
    for i in range(N):
        p = pvalue(160, 40, 60, 404)
    t = time.time() - t
    print "iterations/sec:", float(N) / t

    t = time.time()
    N = 5
    a = np.zeros(N, np.uint)
    print pvalue_npy(a + 160, a + 40, a + 60, a + 404)
    t = time.time() - t
    print "npy iterations/sec:", float(N) / t

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()