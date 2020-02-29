import numpy
import time


n = 25600
a = numpy.random.random((n,n)) + numpy.random.random((n,n)) * 1j

b1 = numpy.random.random((n,n*2))

b2 = numpy.random.random((n*2,n))

ntimes = 1

t0 = time.time()
for i in range(ntimes): aa = a + a
print("complex: ",(time.time()-t0)/ntimes)


t0 = time.time()
for i in range(ntimes): bb1 = b1 + b1
print("float1: ",(time.time()-t0)/ntimes)

t0 = time.time()
for i in range(ntimes): bb2 = b2 + b2
print("float2: ",(time.time()-t0)/ntimes)