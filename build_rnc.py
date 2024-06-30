## simple numerical random number converter from inv(CDF)
## Barry Y. Li and Tim Duong (2024)

import sys
import numpy as np

def generate_random_numbers(filename):
    # 1. parameters
    lower = 0                                       # sampling range lower bound
    upper = 16                                      # sampling range upper bound
    n = int(1e7)                                    # number of random numbers to draw
    randy = np.random.rand(n)
    x = np.arange(lower, upper + 1e-4, 1e-4)

    # 2. define PDF
    def f(x):
        return np.exp(-((x-5) / 2)**2)

    # 3. trapezoidal integral function
    def numint(x, f, i):
        return np.trapz(f[:i+1], x[:i+1])

    # 4. compute cdf numerically
    intf = np.zeros(len(x) - 1)
    for i in range(len(x) - 1):
        intf[i] = 0.5 * (f(x[i+1]) + f(x[i])) * (x[i+1] - x[i])

    cdf_f = np.cumsum(intf)
    # 5. normalize the cdf
    f = f(x) / cdf_f[-1]
    intf = np.zeros(len(x) - 1)
    for i in range(len(x) - 1):
        intf[i] = 0.5 * (f[i+1] + f[i]) * (x[i+1] - x[i])
    cdf_f = np.cumsum(intf)

    # 6. convert random numbers based on the cdf
    def closest_value(y, x):
        findind = 0
        endind = len(y) - 1
        while endind - findind > 1:
            midind = (endind + findind) // 2
            if y[midind] >= x:
                endind = midind
            else:
                findind = midind
        if endind - findind == 1 and abs(y[endind] - x) < abs(y[findind] - x):
            findind = endind
        return findind

    new_vec = np.zeros(n)
    for i in range(n):
        c = closest_value(cdf_f, randy[i])
        new_vec[i] = x[c]

    # 7. Write new_vec to a .log file
    np.savetxt(filename, new_vec)

    print(f"Random numbers saved to {filename}.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 rng.py <filename>")
    else:
        filename = sys.argv[1]
        generate_random_numbers(filename)
