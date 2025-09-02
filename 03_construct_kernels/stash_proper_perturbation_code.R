## =================================== ##
## =================================== ##
## == Doing reproductive phenology  == ##
## =================================== ##
## =================================== ##

perturb.list[[5]] = ltre.backbone %>%
  rename(year = Year) %>%
  mutate(
    # Probability of flowering
    # (not used in this script, but used for phen-growth trade-off)
    prob.flower = 1 - predict(
      u_s_s.ty, newdata = .,
      allow.new.levels = TRUE, re.form = ~ 0, type = 'zprob'
    ),
    # Umbel count
    phen.umbels = predict(
      u_s_s.ty, newdata = .,
      allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
    )
  ) %>%
  rename(Year = year) %>%
  mutate(
    # Model predictions for seed set on linear (link) scale for averaging
    seeds.zinf.linear =  predict(
      s_st.p_s.u.p, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      newparams = s_st.p_s.u.p$fit$par %>%
        (function(x) {
          x[11] <- x[11] + delta
          return(x)
        }),
      type = 'zlink'
    ),
    seeds.seed.linear =  predict(
      s_st.p_s.u.p, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      type = 'link'
    )
  ) %>%
  # Taking out the size-zinf terms for 2021 - very extrapolatory, affects averages too much
  mutate(seeds.zinf.linear = ifelse(Year %in% 2021, NA, seeds.zinf.linear)) %>%
  group_by(size, size.nex, trt, trt.phen, trt.phen.idx, phen.c, phen.umbels, prob.flower) %>%
  summarise(across(c(seeds.zinf.linear, seeds.seed.linear), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(
    # Transform from linear scale to response scale
    seeds.per.umbel = (1 / (1 + exp(seeds.zinf.linear))) * exp(seeds.seed.linear),
    # Total umbels per plant
    seeds.total = seeds.per.umbel * phen.umbels
  ) %>%
  select(-c(seeds.zinf.linear, seeds.seed.linear)) %>%
  mutate(
    # Mean recruit size
    recr.mean = predict(
      r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .
    ),
    # Get the number of seeds produced for each size grouping
    p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sigma.recr)
  ) %>%
  # Rename column
  rename(size.prev = size) %>%
  mutate(
    param = 'phen.succ',
    # here: the *phenology slope* for the zero inflation model * phen
    orig.parval = s_st.p_s.u.p$fit$par[17] * phen.c
  )


# 6: phen effects on seed set (linear term)

perturb.list[[6]] = ltre.backbone %>%
  rename(year = Year) %>%
  mutate(
    # Probability of flowering
    # (not used in this script, but used for phen-growth trade-off)
    prob.flower = 1 - predict(
      u_s_s.ty, newdata = .,
      allow.new.levels = TRUE, re.form = ~ 0, type = 'zprob'
    ),
    # Umbel count
    phen.umbels = predict(
      u_s_s.ty, newdata = .,
      allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
    )
  ) %>%
  rename(Year = year) %>%
  mutate(
    # Model predictions for seed set on linear (link) scale for averaging
    # Model predictions for seed set on linear (link) scale for averaging
    seeds.zinf.linear =  predict(
      s_st.p_s.u.p, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      type = 'zlink'
    ),
    seeds.seed.linear =  predict(
      s_st.p_s.u.p, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      newparams = s_st.p_s.u.p$fit$par %>%
        (function(x) {
          x[1] <- x[1] + delta
          return(x)
        }),
      type = 'link'
    )
  ) %>%
  # Taking out the size-zinf terms for 2021 - very extrapolatory, affects averages too much
  mutate(seeds.zinf.linear = ifelse(Year %in% 2021, NA, seeds.zinf.linear)) %>%
  group_by(size, size.nex, trt, trt.phen, trt.phen.idx, phen.c, phen.umbels, prob.flower) %>%
  summarise(across(c(seeds.zinf.linear, seeds.seed.linear), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(
    # Transform from linear scale to response scale
    seeds.per.umbel = (1 / (1 + exp(seeds.zinf.linear))) * exp(seeds.seed.linear),
    # Total umbels per plant
    seeds.total = seeds.per.umbel * phen.umbels
  ) %>%
  select(-c(seeds.zinf.linear, seeds.seed.linear)) %>%
  mutate(
    # Mean recruit size
    recr.mean = predict(
      r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .
    ),
    # Get the number of seeds produced for each size grouping
    p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sigma.recr)
  ) %>%
  # Rename column
  rename(size.prev = size) %>%
  mutate(
    param = 'phen.seed',
    orig.parval = s_st.p_s.u.p$fit$par[8] * phen.c
  )

perturb.df = do.call(rbind, perturb.list) %>%
  select(-c(phen.c, phen.umbels, seeds.per.umbel, seeds.total, recr.mean, trt, trt.phen))

write.csv(
  perturb.df,
  file = '03_construct_kernels/out/deterministic_repr_coef_perturbation_phen_try2.csv',
  row.names = FALSE, na = ''
)

## =================================== ##
## =================================== ##
## == Doing growth phenology ========= ##
## =================================== ##
## =================================== ##

outputs[[3]] = ltre.backbone %>%
  # Predicted survival
  mutate(
    pred.surv = predict(
      newdata = .,
      object = s_s, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Predicted growth
  mutate(
    # Model with no phenology
    pred.grow.mean = predict(
      newdata = .,
      object = g_st.ty, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Model with phenology
  mutate(phen.grow.mean = pred.grow.mean + delta + (phen.c * phen.effect)) %>%
  # # Take the average of the growth kernel across years
  # group_by(size.prev, size.cur, trt, trt.phen, phen.c, pred.surv, pred.grow.mean) %>%
  # summarise(phen.grow.mean = mean(phen.grow.mean)) %>%
  # ungroup() %>%
  # Predicted distribution of sizes in next time step
  mutate(
    pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = gv.sd),
    pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = gv.sd)
  ) %>%
  # Add in perturbation information
  mutate(
    perturb.param = 'phen.grow',
    orig.par.val = (phen.c) * phen.effect
  )

# Bind them all together
outputs.all = do.call(rbind, outputs)

write.csv(
  outputs.all %>% 
    select(
      size.prev, size.cur, trt.phen.idx, 
      pred.surv, pv.grow.size, pf.grow.size, perturb.param, orig.par.val
    ),
  file = '03_construct_kernels/out/deterministic_grow_coef_perturbation_phen_v2.csv',
  row.names = FALSE
)

## =================================== ##
## =================================== ##
## == Phen diffs for LTRE ============ ##
## =================================== ##
## =================================== ##


obsv.phen.diffs = rbind(
  gs.obsv.pert %>% distinct(trt.phen.idx, param, orig.parval), 
  fr.obsv.pert %>% distinct(trt.phen.idx, param, orig.parval)
) %>%
  merge(trt.phen.ltre.key %>% select(-mean.phen)) %>%
  select(-c(trt.phen.idx)) %>%
  rename(trt = trt.rate, rate = param) %>%
  # get rid of the phen 
  filter(grepl('phen', rate)) %>%
  # Get differences between treatments
  pivot_wider(names_from = trt.phen, values_from = orig.parval) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  select(-c(drought, irrigated, control)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast.phen', values_to = 'pardiff', values_drop_na = TRUE)
  
phen.ltre = merge(midp.phen.sens, obsv.phen.diffs) %>%
  mutate(contrib = pardiff * sv)  
  
  
control.ltre.all = rbind(
  # --- Observed treatment effects
  obsv.trt.ltre %>%
    # give me LTRE values for the control dates and remove column
    filter(trt.phen %in% 'control') %>%
    select(-trt.phen) %>%
    # marker for type of observation
    mutate(varb = 'beta', samp = 'obsv', type = 'trt'),
  # --- Observed phenology effects (within treatment)
  obsv.phen.ltre %>%
    # give me LTRE values where the reference date is the control
    # and remove unnecessary column
    filter(trt %in% 'control') %>%
    select(-trt) %>%
    # Rename column for column agreement
    rename(contrast = contrast.phen) %>%
    mutate(varb = 'alpha', samp = 'obsv', type = 'phen')
) %>%
  mutate(ltre.varb = paste0(varb, '[', rate, ']'))

control.ltre.all %>%
  mutate(
    contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control')
  ) %>%
  ggplot(aes(x = ltre.varb)) +
  # geom_col_pattern(
  #   aes(y = contrib, fill = contrast, pattern = varb),
  #   colour = 'gray22',
  #   pattern_colour = 'gray22', pattern_fill = 'gray22',
  #   pattern_density = 0.025
  # ) +
  geom_col(
    aes(y = contrib, fill = contrast), colour = 'gray22'
  ) +
  # geom_segment(aes(xend = ltre.varb, y = lo, yend = hi), linewidth = 1.2) +
  scale_x_discrete(
    labels = scales::label_parse(),
    limits = c(
      'alpha[grow]', 'alpha[succ]', 'alpha[seed]',
      'beta[grow]', 'beta[flow]', 'beta[seed]', 'beta[recr]'
    ),
    guide = guide_axis(n.dodge = 2)
  ) +
  # scale_pattern_manual(values = c('stripe', 'crosshatch')) +
  scale_fill_manual(values = c('goldenrod', 'dodgerblue')) +
  facet_wrap(~ contr.pretty) +
  labs(x = '', y = expression(paste('Contribution to ', Delta, lambda))) +
  guides(fill = 'none', pattern = 'none') +
  theme(
    panel.background = element_blank(),
    strip.background = element_part_rect(fill = 'white', side = 'b', colour = 'gray22'),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 7)
  )

# Looks good to me.

# What does the Delta lambda comparison look like?

merge(
  phen.ltre %>%
    group_by(trt, contrast.phen) %>%
    summarise(contrib.sum = sum(contrib)),
  obsv.lambda %>%
    select(lambda, trt, trt.phen) %>%
    pivot_wider(names_from = trt.phen, values_from = lambda) %>%
    mutate(d.c = drought - control, i.c = irrigated - control) %>%
    select(-c(drought, irrigated, control)) %>%
    pivot_longer(c(d.c, i.c), names_to = 'contrast.phen', values_to = 'lambda.diff', values_drop_na = TRUE)
)
# but with the original parameterization this is spot on lmao oly moly
# holy guacamole!!!!!!

# ayy there we go