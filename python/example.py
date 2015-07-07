import ECCSAMPLES
import random

# Number of samples to draw
n = int(1e4)
transit = 1
verbose = 0

# Define alpha and beta
alpha = 0.867
beta = 3.030

# Generate samples
samples=[[0 for j in range(2)] for i in range(n)]
for i in range(n):
  xe = random.uniform(0, 1)
  xw = random.uniform(0, 1)
  ew = ECCSAMPLES.ecc_sample(transit,verbose,alpha,beta,xe,xw)
  samples[i][0] = ew[0]
  samples[i][1] = ew[1]
  
# Export samples
with open("ECCSAMPLES.dat", 'w') as f:
    f.writelines(' '.join(str(j) for j in i) + '\n' for i in samples)
