library(dplyr)
library(tidyr)
library(purrr)

# Growth + survival subkernel
gs.control = read.csv('03_construct_kernels/out/deterministic_growsurv_kernel_phen_ltre.csv') %>%
  filter(trt.phen.idx %in% 7) %>%
  select(-trt.phen.idx)
# Flowering + reproduction subkernel
fr.control = read.csv('03_construct_kernels/out/deterministic_reprod_kernel_phen_ltre.csv') %>%
  filter(trt.phen.idx %in% 7) %>%
  select(-trt.phen.idx)

kernel = merge(
  # Survival + growth subkernel
  gs.control,
  # Reproductive subkernels
  fr.control,
  by.x = c('size.prev', 'size.cur'), by.y = c('size.prev', 'size.nex'),
  suffixes = c('.g', '.r')
)

p.germ.cur = 0.000
cur.digit = 1000
tol = 1e-7

p.germ.df = data.frame(p.germ = p.germ.cur + (0:9)/cur.digit)

lambda.cur = kernel %>%
  merge(p.germ.df) %>%
  # Combine growth and survival entries into single kernel entry
  mutate(
    p.size.cur = pred.surv * (pv.grow.size * (1 - prob.flower) + pf.grow.size * prob.flower) + (p.size.cur * p.germ)
  ) %>%
  split(~ p.germ) %>%
  map_vec(
    \(m) m %>%
      select(size.prev, size.cur, p.size.cur) %>%
      arrange(size.prev, size.cur) %>%
      pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
      select(-c(size.cur)) %>%
      as.matrix() %>%
      (\(m) Re(eigen(m)$values[1]))
  )

while (all(abs(lambda.cur - 1) > tol)) {
  
  print(p.germ.cur)
  
  p.germ.cur = p.germ.cur + ((max(which(lambda.cur < 1)) - 1) /cur.digit)
  cur.digit = cur.digit * 10
  
  p.germ.df = data.frame(p.germ = p.germ.cur + (0:9)/cur.digit)
  
  lambda.cur = kernel %>%
    merge(p.germ.df) %>%
    # Combine growth and survival entries into single kernel entry
    mutate(
      p.size.cur = pred.surv * (pv.grow.size * (1 - prob.flower) + pf.grow.size * prob.flower) + (p.size.cur * p.germ)
    ) %>%
    split(~ p.germ) %>%
    map_vec(
      \(m) m %>%
        select(size.prev, size.cur, p.size.cur) %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(size.cur)) %>%
        as.matrix() %>%
        (\(m) Re(eigen(m)$values[1]))
    )
  
}

p.germ = lambda.cur %>%
  (\(v) names(v)[c(max(which(v < 1)), min(which(v > 1)))])() %>%
  as.numeric() %>%
  mean()

p.germ
# 0.004632985

kernel %>%
  mutate(
    p.size.cur = pred.surv * (pv.grow.size * (1 - prob.flower) + pf.grow.size * prob.flower) + (p.size.cur * p.germ)
  ) %>%
  (
    \(m) m %>%
      select(size.prev, size.cur, p.size.cur) %>%
      arrange(size.prev, size.cur) %>%
      pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
      select(-c(size.cur)) %>%
      as.matrix() %>%
      (\(m) Re(eigen(m)$values[1]))
  )()
