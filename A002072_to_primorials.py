# A002072_to_primorials.py
"""
Use the A002072 data to calculate the iota value for each
L * (L+1) value in its primorial interval.

By Ken Clements, June 25, 2026
"""

from math import isqrt, log10
from sympy import prime, primorial, primefactors, primepi

A002072 = [ 0, 1, 8, 80, 4374, 9800,
      123200, 336140, 11859210, 11859210, 177182720,
      1611308699, 3463199999, 63927525375, 421138799639, 1109496723125,
      1453579866024, 20628591204480, 31887350832896, 31887350832896, 119089041053696,
      2286831727304144, 9591468737351909375, 9591468737351909375, 9591468737351909375, 9591468737351909375,
      9591468737351909375, 19316158377073923834000, 19316158377073923834000, 19316158377073923834000, 19316158377073923834000, 
      19316158377073923834000, 19316158377073923834000, 19316158377073923834000, 124225935845233319439173, 124225935845233319439173,
      124225935845233319439173, 27129807647978258459761875
]
# LxL1 = L_r * (L_r + 1)
LxL1 = [ L * (L+1) for L in A002072]

primorials = [0]

for r in range(1, 100):
    P_r = primorial(r)
    primorials.append(P_r)
    if P_r > LxL1[-1]:
      break
minIota = [0]

print(f"\nFor the values of L in OEIS Sequence A002072, here is the list of where")
print(f"    the product L times (L+1) fits into the primorial intervals.")
print(f"Iota is the prime index of the GPF of the product, and when greater")
print(f"    than the index of the primorial interval, precludes the existance of")
print(f"    any prime-complete product of two consecutive integers in that interval.")

for r in range(1, len(primorials)-1):
    P_r = primorials[r]
    slot = 1
    while slot < len(LxL1) and LxL1[slot] < P_r: slot += 1
    if slot >= len(LxL1): break
    minIota.append(slot-1)
# print(minIota)    
print("\n") 
pidx, lidx, last_n, r = 1, 1, 0, 1
while pidx < len(primorials) and lidx < len(LxL1):
    if LxL1[lidx] >= last_n:
      r = lidx
      pf = primefactors(LxL1[lidx])
      iota = primepi(pf[-1])
      print(f"LxL1[{lidx:2}] ={LxL1[lidx]:55}, log10={log10(LxL1[lidx]):6.2f} iota = {iota}")
      last_n = LxL1[lidx]
      lidx += 1
      while lidx < len(LxL1) and LxL1[lidx] == last_n:
        pf = primefactors(LxL1[lidx])
        iota = primepi(pf[-1])
        print(f"LxL1[{lidx:2}] ={LxL1[lidx]:55}, log10={log10(LxL1[lidx]):6.2f} iota = {iota}")
        last_n = LxL1[lidx]
        lidx += 1      
    while pidx < len(primorials) and lidx < len(LxL1) and primorials[pidx] < LxL1[lidx]:
      print(f" P_r[{pidx:2}] ={primorials[pidx]:55}, log10={log10(primorials[pidx]):6.2f}")
      pidx += 1
print(f" P_r[{pidx:2}] ={primorials[pidx]:55}, log10={log10(primorials[pidx]):6.2f}")
print(f"\nEnd of Program")      
        


"""
Program output:

For the values of L in OEIS Sequence A002072, here is the list of where
    the product L times (L+1) fits into the primorial intervals.
Iota is the prime index of the GPF of the product, and when greater
    than the index of the primorial interval, precludes the existance of
    any prime-complete product of two consecutive integers in that interval.


LxL1[ 1] =                                                      2, log10=  0.30 iota = 1
P_r[ 1]  =                                                      2, log10=  0.30
P_r[ 2]  =                                                      6, log10=  0.78
P_r[ 3]  =                                                     30, log10=  1.48
LxL1[ 2] =                                                     72, log10=  1.86 iota = 2
P_r[ 4]  =                                                    210, log10=  2.32
P_r[ 5]  =                                                   2310, log10=  3.36
LxL1[ 3] =                                                   6480, log10=  3.81 iota = 3
P_r[ 6]  =                                                  30030, log10=  4.48
P_r[ 7]  =                                                 510510, log10=  5.71
P_r[ 8]  =                                                9699690, log10=  6.99
LxL1[ 4] =                                               19136250, log10=  7.28 iota = 4
LxL1[ 5] =                                               96049800, log10=  7.98 iota = 5
P_r[ 9]  =                                              223092870, log10=  8.35
P_r[10]  =                                             6469693230, log10=  9.81
LxL1[ 6] =                                            15178363200, log10= 10.18 iota = 6
LxL1[ 7] =                                           112990435740, log10= 11.05 iota = 7
P_r[11]  =                                           200560490130, log10= 11.30
P_r[12]  =                                          7420738134810, log10= 12.87
LxL1[ 8] =                                        140640873683310, log10= 14.15 iota = 8
LxL1[ 9] =                                        140640873683310, log10= 14.15 iota = 8
P_r[13]  =                                        304250263527210, log10= 14.48
P_r[14]  =                                      13082761331670030, log10= 16.12
LxL1[10] =                                      31393716443781120, log10= 16.50 iota = 10
P_r[15]  =                                     614889782588491410, log10= 17.79
LxL1[11] =                                    2596315725084381300, log10= 18.41 iota = 11
LxL1[12] =                                   11993754236536800000, log10= 19.08 iota = 12
P_r[16]  =                                   32589158477190044730, log10= 19.51
P_r[17]  =                                 1922760350154212639070, log10= 21.28
LxL1[13] =                                 4086728500635196416000, log10= 21.61 iota = 13
P_r[18]  =                               117288381359406970983270, log10= 23.07
LxL1[14] =                               177357888561798925329960, log10= 23.25 iota = 14
LxL1[15] =                              1230982978626222406488750, log10= 24.09 iota = 15
LxL1[16] =                              2112894426911803369434600, log10= 24.32 iota = 16
P_r[19]  =                              7858321551080267055879090, log10= 24.90
LxL1[17] =                            425538775081570245763274880, log10= 26.63 iota = 17
P_r[20]  =                            557940830126698960967415390, log10= 26.75
LxL1[18] =                           1016803143140225112266579712, log10= 27.01 iota = 18
LxL1[19] =                           1016803143140225112266579712, log10= 27.01 iota = 18
LxL1[20] =                          14182199699089010382996314112, log10= 28.15 iota = 20
P_r[21]  =                          40729680599249024150621323470, log10= 28.61
P_r[22]  =                        3217644767340672907899084554130, log10= 30.51
LxL1[21] =                        5229599349004857113477606876880, log10= 30.72 iota = 21
P_r[23]  =                      267064515689275851355624017992790, log10= 32.43
P_r[24]  =                    23768741896345550770650537601358310, log10= 34.38
P_r[25]  =                  2305567963945518424753102147331756070, log10= 36.36
LxL1[22] =                 91996272539599030715854727695564800000, log10= 37.96 iota = 22
LxL1[23] =                 91996272539599030715854727695564800000, log10= 37.96 iota = 22
LxL1[24] =                 91996272539599030715854727695564800000, log10= 37.96 iota = 22
LxL1[25] =                 91996272539599030715854727695564800000, log10= 37.96 iota = 22
LxL1[26] =                 91996272539599030715854727695564800000, log10= 37.96 iota = 22
P_r[26]  =                232862364358497360900063316880507363070, log10= 38.37
P_r[27]  =              23984823528925228172706521638692258396210, log10= 40.38
P_r[28]  =            2566376117594999414479597815340071648394470, log10= 42.41
P_r[29]  =          279734996817854936178276161872067809674997230, log10= 44.45
LxL1[27] =          373113974448203123099782895727610333479834000, log10= 44.57 iota = 27
LxL1[28] =          373113974448203123099782895727610333479834000, log10= 44.57 iota = 27
LxL1[29] =          373113974448203123099782895727610333479834000, log10= 44.57 iota = 27
LxL1[30] =          373113974448203123099782895727610333479834000, log10= 44.57 iota = 27
LxL1[31] =          373113974448203123099782895727610333479834000, log10= 44.57 iota = 27
LxL1[32] =          373113974448203123099782895727610333479834000, log10= 44.57 iota = 27
LxL1[33] =          373113974448203123099782895727610333479834000, log10= 44.57 iota = 27
LxL1[34] =        15432083136624024515389371619608911236566363102, log10= 46.19 iota = 33
LxL1[35] =        15432083136624024515389371619608911236566363102, log10= 46.19 iota = 33
LxL1[36] =        15432083136624024515389371619608911236566363102, log10= 46.19 iota = 33
P_r[30]  =        31610054640417607788145206291543662493274686990, log10= 46.50
P_r[31]  =      4014476939333036189094441199026045136645885247730, log10= 48.60
P_r[32]  =    525896479052627740771371797072411912900610967452630, log10= 50.72
LxL1[37] =    736026463016299604294737333092019573986740163277500, log10= 50.87 iota = 37
P_r[33]  =  72047817630210000485677936198920432067383702541010310, log10= 52.86

"""