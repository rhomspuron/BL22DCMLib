cpdef long calc_e(long a):
    if a <=1:
       return 1
    return a* calc_e(a-1)

