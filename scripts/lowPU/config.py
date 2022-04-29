
## common configuration for lowPU analysis


import decimal
def drange(x, y, jump):
    while x < y:
        yield float(x)
        x += decimal.Decimal(jump)

bins_recoil_reco = [0, 5, 10, 15, 20, 30, 40, 50, 60, 75, 90, 10000]
bins_recoil_gen = [0.0, 10.0, 20.0, 40.0, 60.0, 90.0, 10000]

bins_recoil_qT = list(drange(0, 30, 0.5)) + list(range(30, 60, 2)) + list(range(60, 100, 5)) + list(range(100, 210, 10)) + [10000]