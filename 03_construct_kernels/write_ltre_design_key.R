# File for designating treatment combinations for LTRE
# Permanent version that can be read in for future use
# Will add other read-ins as needed

ltre.trt.key = data.frame(
  trt.phen.idx = 1:7,
  trt.phen = c('drought', 'control', 'irrigated', 'control', 'drought', 'irrigated', 'control'),
  trt.rate = c('drought', 'drought', 'irrigated', 'irrigated', 'control', 'control', 'control')
)

write.csv(ltre.trt.key, na = '', row.names = FALSE, file = '03_construct_kernels/ltre_treatment_key.csv')
