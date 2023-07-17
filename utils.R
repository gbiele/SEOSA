source("data/get_item_name.R")


model_summary_tbl = function(fit,grep_string = "b_|sd_|bs_|bcs_|sigma") {
  footnote_text = 
    paste(
      c(paste("model: ",fit$formula[1]),
        "b_ = fix effect parameters",
        "sd_ standard deviation of hierarchical parameters",
        "ess = effective sample size"),
      collapse = "\n")
  
  sample_diags = 
    fit$fit %>% 
    as_draws() %>% 
    subset_draws(grep_string, regex = TRUE) %>% 
    summarize_draws() %>% 
    data.table() %>% 
    .[, c("mad","median") := NULL] %>% 
    setnames("variable","parameter") 
  
  if (fit$family$family == "acat")
    sample_diags[, parameter := gsub("b_Intercept","Threshold",parameter)]
  
  sample_diags %>% 
    kable(caption = "Summary of model paramters",
          digits = c(1,1,2,1,1,3,0,0)) %>% 
    kable_styling(full_width = FALSE) %>% 
    kable_classic() %>% 
    add_footnote(footnote_text, notation = "none")
}

print_effect = function(x, digits = 0, short = FALSE, BF.thresh = 0, contr = "S>C") {
  correctBF = function(d) ifelse(d > 10e5,">10^5", ifelse(d < 10e-5,"<10^-5", paste0("=",d)))
  correctP = function(d) ifelse(d > .999,">0.999", ifelse(d < .001,"<0.001", paste0("=",d)))
  
  
  x = as.numeric(x)
  vals = c(
    mean(x),
    quantile(x,probs = c(.025,.975)),
    mean(x>BF.thresh),
    round(sum(x>BF.thresh)/sum(x<BF.thresh),1)
  )  
  vals[1:3] = round(vals[1:3], digits = digits)
  vals[4] = round(vals[4], digits = 2)
  vals[5] = round(vals[5], digits = 2)
  return(
    ifelse(
      short == FALSE,
      paste0(
        vals[1],
        " (CI=[",vals[2],", ",vals[3],"]",
        "; BF", correctBF(vals[5]),
        "; P(", contr, ")", correctP(vals[4]),")"
      ),
      paste0(vals[1]," (",vals[2],", ",vals[3],")"))
  )
}



get_mci = function(value, digits=0, get.P = TRUE) {
  get_P = function(value) {
    x = mean(value>0)
    if (x > .999) {
      return(">.999")
    } else if (x < .001) {
      return("<.001")
    } else {
      return(x)
    }
  }
  frmt = gsub("D",digits,"%.Df [%.Df, %.Df]")
  Pr = paste0("; ",ifelse(get.P == TRUE, get_P(value),""))
  return(paste0(sprintf(frmt,
                        mean(value),
                        quantile(value,probs = .025),
                        quantile(value,.975)),
                Pr))
}

print_PSs = function(x) {
  paste0(round(x[1,mean],1),
         " (CI=[", x[1,CI],"], BF=",x[1,BF],")")
}

make_tbl = function(dt, caption) {
  tmp = 
    dt %>% 
    .[, .(`stats` = get_mci(est)), by = .(effect)] %>% 
    .[, effect := gsub("_vs_"," - ", effect)] %>%
    .[, effect := gsub("_in_"," in ", effect)] %>% 
    .[, effect := gsub("in_","", effect)]
  tmp %>% 
    setnames(names(tmp), c("Effect~a~", "mean (CrI)~b~")) %>% 
    kable(caption = caption) %>% 
    kable_styling(full_width = FALSE) %>% 
    kable_classic() 
}

# selector = function(dt) {
#   n.oxy2 = 
#     unique(dt[,.(participant,session,Drug)]) %>% 
#     .[, N.oxy := sum(Drug == "oxycodone"), by = .(participant)] %>% 
#     .[N.oxy == 2, participant] %>% 
#     unique()
#   dt = 
#     dt %>% 
#     .[, N.sessions := length(unique(session)), by = .(participant)] %>% 
#     .[pilot == F & !grepl("ADHD", comment) & N.sessions == 4] %>% 
#     .[pilot == F & !grepl("ADHD", comment) & participant %in% n.oxy2] 
# }

#Marie and Martin tested this "new selector" May 9th 23
selector = function(dt) {
  dt = 
    dt %>% 
    .[pilot == F & !grepl("ADHD", comment) & !grepl("PSY", comment) & !grepl("failed drug test", comment) ]
}


coarsen = function(dt, var, bins = 20, range = c(0,100), breaks = NULL) {
  var.c = paste0(var,".c")
  if (is.null(breaks)) {
    breaks = seq(range[1],range[2],length.out = bins + 1)
  }
  my_cut = function(x) {
    return(
      cut(x,
          breaks = breaks, 
          include.lowest = TRUE, 
          right = FALSE,
          ordered_result = TRUE)
    )
  }
  
  dt %>% 
    .[, (var.c) := my_cut(get(var))]
  
  c.levels = 
    my_cut(range[1]:range[2]) %>% 
    levels()
  
  var.L = paste0(var,".L")
  dt[, (var.L) := which(get(var.c) == c.levels),by = 1:nrow(dt)]
  
  return(dt)
}

get_ftt = function(effect_obtained,target, time, bias) {
  if (any(bias < 0)) {
    ftt = 
      ifelse(any(effect_obtained >= target),
             min(time[effect_obtained >= target]),
             Inf)
  } else {
    ftt = 
      ifelse(any(effect_obtained <= target),
             min(time[effect_obtained <= target]),
             Inf)
  }
  return(ftt)
}

#Here is a function to create plots of ratings across stages
plot_by_stage = function(dt, stages = NULL, item_names = NULL, color.var = NULL, lty.var = NULL) {
  if (is.null(color.var) & is.null(lty.var))
    stop("Specify at least either color or lty.")
  if (is.null(stages))
    stages = unique(dt$stage)
  if (is.null(item_names))
    stop("Specify at least one item_name")
  add.factet = FALSE
  if(is.list(item_names)) {
    item.categorie = data.table(
      state = names(unlist(item_names)) %>% gsub("[0-9]","",.),
      item_name = unlist(item_names)
    ) %>% 
      setkeyv("item_name")
    add.facet = TRUE
    item_names = item.categorie$item_name
    dt = merge(dt,item.categorie)
  } else {
    dt[,state := NA]
  }
  
  p.data = 
    dt %>% 
    .[item_name %in% item_names & stage %in% stages] %>% 
    .[, response.demeaned := response - mean(response), by = c("participant","item_name","state","Sex")] %>% 
    .[, .(m = mean(response),
          sd = sd(response.demeaned),
          N = length(unique(participant))), 
      by = c(color.var,lty.var,"stage","state")] %>%
    .[, se := sd/sqrt(N-1)] %>% 
    .[, `:=`(lower = m - qnorm(.975)*se, upper = m + qnorm(.975))]
  
  p = 
    p.data %>% 
    ggplot(aes_string(x = "stage", y = "m", color = color.var, lty = lty.var, fill = color.var)) +
    annotate("rect",xmin = 5.5, xmax = 5.75, ymin = 0, ymax = 100, color = NA, fill = "blue", alpha = .05) +
    annotate("rect",xmin = 6.5, xmax = 6.75, ymin = 0, ymax = 100, color = NA, fill = "blue", alpha = .05) +
    annotate("rect",xmin = 10.25, xmax = 10.75, ymin = 0, ymax = 100, color = NA, fill = "red", alpha = .1) +
    geom_line() + 
    geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = .25) + 
    ylab("response") +
    coord_cartesian(ylim = c(0,100))
  
  if (min(stages)<5)
    p = p + 
    annotate("rect",xmin = 3.25, xmax = 3.75, ymin = 0, ymax = 100, color = NA, fill = "red", alpha = .1) +
    annotate("rect",xmin = 4.25, xmax = 4.75, ymin = 0, ymax = 100, color = NA, fill = "blue", alpha = .1) 
  
  if (add.facet == TRUE)
    p = p + facet_wrap(~state)
  
  p = p + theme(legend.position = "right")
  
  return(p)
}


# contrast samples for effect obtained
# As a general procedure, we predict from model parameters 
# for each participant expected outcomes under each condition
# and calculate contrasts as comparisons between conditions
make_plot_samples = function(fit, dt, outcome.var = NULL) {
  if (is.null(outcome.var)) {
    stop("The function requires explicit declaration of an outcome variable 'outcome.var'")
  }
  if (is.null(dt)) {
    stop("The function requires explicit declaration of the original data set")
  }
  
  fn = paste0("brmsP1/obtained/pp_",abs(as_draws_df(fit)[,1][[1]][1]),".Rdata")
  
  if (file.exists(fn)) {
    load(fn)
  } else {
    # generate data for posterior predictions
    # to remove effects of session, we replicate data such that
    # we generated predicted values for each participant and condition four times
    # one time for each session
    new_data = do.call(rbind,lapply(1:4, function(x) data.table(fit$data)[, session := x]))
    # predict probabilities that response falls into levels of the coarsened outcome variable
    pp = posterior_epred(fit, newdata = new_data)
    
    var = outcome.var %>%
      gsub("\\.c","",.)
    
    # mean response for each level of the coarsened categorical variable in oberved data
    bin.means =
      dt[, .(m = mean(get(var))), by = c(outcome.var)] %>%
      setkeyv(outcome.var) %>%
      .[,m]
    # calculate expected response on orginal scale as weighted mean
    # of average response for each level (bin), with level probabilities as weights
    tmp =
      apply(pp, 1, function(x1) { # loop over iterations
        apply(x1,1, function(x2) { # loop over participants
          fsum(bin.means,w=x2) # expected score per participant
        })
      })
    
    ## calculate the effect obtained by state.condition and sex
    dt.samples = 
      new_data %>% 
      .[, .(state.condition, Sex, session)] %>% 
      cbind(tmp) %>% 
      melt(id.vars = c("state.condition", "Sex", "session"),
           variable.name = "iter") %>% 
      .[, .(value = mean(value)), by = .(state.condition, Sex, iter)] %>% 
      .[, iter := as.numeric(gsub("V","",iter))]
    
    dt.samples = rbind(
      dt.samples,
      dt.samples %>% 
        .[, .(value = mean(value)), by = .(state.condition,iter)] %>% 
        .[, Sex := "all"]
    )
    
    ## contrast between stress (state.conditions), split by sex
    dt.samples.d = 
      new_data %>% 
      .[, .(state.condition, Sex, session)] %>% 
      cbind(tmp) %>% 
      melt(id.vars = c("state.condition", "Sex", "session"),
           variable.name = "iter") %>% 
      .[, .(value = mean(value)), by = .(state.condition, Sex, iter)] %>% 
      dcast(Sex + iter ~ state.condition, value.var = "value") %>% 
      .[, value := stress - control] %>% 
      .[, iter := as.numeric(gsub("V","",iter))]
    
    ## add contrasts
    ## - average effect across men and women
    ## - sex difference
    dt.samples.d = rbind(
      dt.samples.d[, .(Sex,iter,value)],
      dt.samples.d %>% 
        .[, .(value = mean(value)), by = .(iter)] %>% 
        .[, Sex := "all"],
      dt.samples.d %>% 
        dcast(iter ~ Sex, value.var = "value") %>% 
        .[, .(value = men-women, Sex = "sex diff: m-f"), by = .(iter)]
    )
    
    dt.samples.d %>% 
      .[, grp := ifelse(Sex %in% c("women","men"), "by_sex","all")]
    save(dt.samples,dt.samples.d, file = fn)
  }
  
  return(list(
    dt.samples = dt.samples,
    dt.samples.d = dt.samples.d
  ))
  
}


## plot histogram fo contrast for plot_results_origscale
diff.plotter = function(dts, sx) {
  dts[, Sex := factor(Sex, levels = sx)]
  clr.values = sex_colors[sx]
  clr.values[is.na(clr.values)] = "grey"
  
  annotations <- data.frame(
    Sex = sx,
    text = c(sx),
    X = c(Inf),
    Y =  c(Inf),
    x_adjust = c(1),
    y_adjust = c(1))
  
  p = 
    dts %>% 
    .[Sex %in% sx] %>% 
    ggplot(aes(y = value)) + 
    stat_halfeye(aes(fill = Sex, alpha = after_stat(y>0)), color  ="black") + 
    scale_fill_manual(values = clr.values) + 
    scale_alpha_manual(values = c(.2,1)) +
    theme(legend.position = "none") +
    ylab("Stress effect") + 
    coord_cartesian(ylim = quantile(dts[!is.na(Sex),value], c(.0005,.9995))) +
    geom_hline(yintercept = 0, lty = 2) + 
    facet_wrap(~Sex, ncol = 1) + 
    theme(panel.border = element_blank(),
          text = element_text(size = 17),
          axis.line.y = element_line(color="black"),
          axis.title.x = element_blank(),
          axis.ticks.length.x = unit(0,"cm"),
          axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text = element_blank())
    p = p + 
      geom_text(aes(x = X, y=Y,hjust=x_adjust,vjust=y_adjust,label=text),
                data = annotations, inherit.aes = FALSE)
    return(p)
}

# Function to plot results (and more)
plot_results_origscale = function(fit, dt, outcome.var = NULL, plot.sex_avg = "yes") {
  
  if (is.null(outcome.var)) {
    stop("The function requires explicit declaration of an outcome variable 'outcome.var'")
  }
  if (is.null(dt)) {
    stop("The function requires explicit declaration of the original data set")
  }
  
  var = outcome.var %>%
    gsub("\\.c","",.)
  
  plotting_samples = make_plot_samples(fit, dt, outcome.var)
  
  dt.samples = plotting_samples[["dt.samples"]] %>% 
    .[, scl := ifelse(state.condition == "stress",.25, -.25)] # control direction (left/right) of posterior density plot
  dt.samples.d = plotting_samples[["dt.samples.d"]] 
  
  dt.samples %>% 
    .[, Sex := gsub("man","men",Sex)] %>% 
    .[Sex == "all", Sex := "total"]
  
  dt.samples.d %>% 
    .[, Sex := gsub("man","men",Sex)] %>% 
    .[Sex == "all", Sex := "total"]
  
  my.clrs = palette()[c(8,2,4)]
  names(my.clrs) = c("total","women","men")
  if(plot.sex_avg == "no") {
    dt.samples = dt.samples[Sex != "total"]
    dt.samples.d = dt.samples.d[Sex != "total"]
    my.clrs = my.clrs[-1]
  } else if (plot.sex_avg == "only") {
    dt.samples = dt.samples[Sex == "total"]
    dt.samples.d = dt.samples.d[Sex == "total"]
  } else {
    dt.samples[, grp := ifelse(Sex == "total","a","b")]
    dt.samples.d[, grp := ifelse(Sex == "total","a","b")]
  }
  
  ## wanting by stress
  
  stats.scale = dt.samples %>% 
    na.omit() %>%
    .[, .(mean = mean(value),
          lower = quantile(value,.025),
          upper = quantile(value, .975)),
      by = .(state.condition, Sex)] %>% 
    .[, CI := paste(round(c(lower,upper),2), collapse = ", "), by = .(Sex,state.condition)]
  
  stats.scale[, grp := ifelse(Sex == "all", "all", "by sex")]
  
  if (plot.sex_avg != "only") {
  ## main effect of sex
  tbl_me_sex = 
    dt.samples %>% 
    dcast(state.condition + iter ~ Sex, value.var = "value") %>% 
    .[, main_woman := men-women] %>% 
    .[, .(state.condition, iter, main_woman)] %>% 
    dcast(iter ~ state.condition, value.var = "main_woman")  %>% 
    .[, overall := (control+stress)/2] %>% 
    melt(id.vars = "iter", variable.name = "state condition") %>% 
    .[, .(effect_woman = get_mci(value,get.P = TRUE)), by = c("state condition")] %>% 
    kable(caption = "Effect of sex, calculated as effect obtained in men - women.") %>% 
    kable_styling(full_width = FALSE) %>% 
    kable_classic() %>% 
    add_footnote("Numbers are means and CIs: Upper and lower bound of 95% credible intervals",
                 notation = "none")
  } else  {
    tbl_me_sex = NULL
  }
  p1 =
    dt.samples %>%
    na.omit() %>%
    .[, sc := as.numeric(factor(state.condition))] %>% 
    ggplot(aes(x = sc, y = value, color = Sex, group = factor(iter):factor(Sex))) + 
    geom_line(data = dt.samples[iter %in% 1:150], alpha = .15)  +
    stat_halfeye(aes(group = Sex, fill = Sex, scale = scl, slab_alpha = after_stat(f)),
                 height = .5, show.legend = c(slab_alpha = FALSE, size = FALSE)) +
    stat_pointinterval(aes(group = Sex)) + 
    theme(legend.position = "none", axis.line = element_line(),text = element_text(size = 20)) +
    ylab("% effect obtained (0-125)") +
    xlab("") +
    scale_col_sex 
  yrange = layer_scales(p1)$y$range$range
  p1 = p1 + 
    guides(y = guide_axis_truncated(
      trunc_lower = yrange[1]+diff(yrange)/2000)) +
    scale_x_continuous(breaks = c(1,2), labels = c("control","stress")) + 
    coord_cartesian(xlim = c(.8,2.2))
  
  if (plot.sex_avg == "yes")
    p1 = p1 + facet_wrap(~grp) + theme(strip.text = element_blank())
  
  ## wanting stress contrast
  contr.stats = dt.samples.d %>% 
    .[, .(mean = mean(value),
          lower = quantile(value,.025),
          upper = quantile(value, .975),
          `P > 0` = mean(value>0),
          BF = round(mean(value > 0)/mean(value < 0),1)),
      by = .(Sex, grp)] %>% 
    .[, CI := paste(round(c(lower,upper),2), collapse = ", "), by = .(Sex)] 
  
  # estimated effect obtained
  css = function(value) {sprintf("%.0f (%.0f, %.0f)",
                                 mean(value),
                                 quantile(value,probs = .025),
                                 quantile(value,.975))}
  simple.stats = 
    dt.samples %>% 
    .[, .(`value` = css(value)), by = .(state.condition,Sex)] %>% 
    dcast(Sex ~state.condition, value.var = "value")
  

  if (plot.sex_avg == "yes") {
    p2 = diff.plotter(dt.samples.d,c("men","total","women"))
  } else if (plot.sex_avg == "no") {
    p2 = diff.plotter(dt.samples.d,c("men","women"))
  } else {
    p2 = diff.plotter(dt.samples.d,c("total"))
  }
  
  p = (p1 | p2) + plot_layout(widths = c(4,1))
  
  return(list(stats = contr.stats, simple.stats = simple.stats, stats.scale = stats.scale, plot = p, p1 = p1, p2 = p2,
              dt.samples = dt.samples, dt.samples.d = dt.samples.d, tbl_me_sex = tbl_me_sex))
}

plot_results_origscale2 = function(fit) {
  outcome.var = "effect_obtained.c"
  var = outcome.var %>% gsub("\\.c","",.)
  
  fn = paste0("brmsP1/pp_",outcome.var,"_",abs(as_draws_df(fit)[,1][[1]][1]),".Rdata")
  if (file.exists(fn)) {
    load(fn)
  } else {
    new_data = do.call(rbind,lapply(1:4, function(x) data.table(fit$data)[, session := x]))
    pp = posterior_epred(fit, newdata = new_data)
    
    bin.means =
      my_data_beh %>% 
      .[, .(m = mean(get(var))), by = c(outcome.var)] %>%
      setkeyv(outcome.var) %>%
      .[,m]
    tmp =
      apply(pp, 1, function(x1) { # loop over iterations
        apply(x1,1, function(x2) { # loop over participants
          sum(x2*bin.means) # expected score per participant
        })
      })
    
    ## wanting by stress, split by sex
    dt.samples = 
      new_data %>% 
      .[, .(state.condition, session)] %>% 
      cbind(tmp) %>% 
      melt(id.vars = c("state.condition", "session"),
           variable.name = "iter") %>% 
      .[, .(value = fmean(value)), by = .(state.condition, iter)] %>% 
      .[, iter := as.numeric(gsub("V","",iter))]
    
    ## wanting stress contrast, split by sex
    dt.samples.d = 
      new_data %>% 
      .[, .(state.condition, session)] %>% 
      cbind(tmp) %>% 
      melt(id.vars = c("state.condition", "session"),
           variable.name = "iter") %>% 
      .[, .(value = fmean(value)), by = .(state.condition, iter)] %>% 
      dcast(iter ~ state.condition, value.var = "value") %>% 
      .[, value := stress - control] %>% 
      .[, iter := as.numeric(gsub("V","",iter))]
    
    save(dt.samples,dt.samples.d, file = fn)
  }
  # needed for prior predictive
  dt.samples = dt.samples %>% na.omit()
  dt.samples.d = dt.samples.d %>% na.omit()
  
  # control direction (left/right) of posterior density plot
  dt.samples[, scl := ifelse(state.condition == "stress",.25, -.25)]
  
  p1 =
    dt.samples %>% 
    na.omit() %>%
    ggplot(aes(x = state.condition, y = value)) + 
    geom_line(data = dt.samples[iter %in% 1:150], aes(group = factor(iter)), alpha = .15)  +
    stat_halfeye(aes(scale = scl, slab_alpha = after_stat(f)),
                 height = .5, show.legend = c(slab_alpha = FALSE, size = FALSE)) +
    stat_pointinterval() + 
    theme(legend.position = "none", axis.line = element_line(),text = element_text(size = 20)) +
    guides(y = guide_axis_truncated(
      trunc_lower = floor(fmin(dt.samples[iter %in% 1:150,value])/10)*10,
      trunc_upper = c(120))) + 
    coord_cartesian(xlim = c(1.25,1.75)) +
    ylab("Drug self-administration (% effect obtained)") +
    xlab("State condition") 
  
  contr.stats = 
    dt.samples.d %>% 
    .[, .(mean = fmean(value),
          lower = round(fquantile(value,.025)),
          upper = round(fquantile(value, .975)),
          `P > 0` = fmean(value>0),
          BF = round(fmean(value > 0)/mean(value < 0),1)),
      by = .()] %>% 
    .[, CI := paste(round(c(lower,upper),2), collapse = ", ")] 
  
  simple.stats = 
    dt.samples %>% 
    .[, .(`m (CrI)` = sprintf("%.0f (%.0f, %.0f)",
                              fmean(value),
                              fquantile(value,probs = .025),
                              fquantile(value,.975))),
      by = .(state.condition)]
  
  p2 = 
    diff.plotter(dt.samples.d[, Sex := "all"],"all") + 
    ylab("Effect of stress on drug wanting (difference % obtained)")
  
  p = p1 | p2
  
  return(list(stats = contr.stats, simple.stats = simple.stats, plot = p, p1 = p1, p2 = p2,
              dt.samples = dt.samples, dt.samples.d = dt.samples.d))
  
}


plot_buffer_effect = function(fit,the_data,var) {
  tmp = fit$fit %>% as_draws()
  fn = paste0("brmsP1/buffer/contr_",var,tmp[1],".Rdata")
  rm(tmp)
  gc()
  if (!file.exists(fn)) {
    new_data = do.call(rbind,lapply(1:4, function(x) data.table(fit$data)[, session := x]))
    epred.data = fit$data %>% data.table() %>% .[, response.c := NULL]
    # average response on 0-100 scale in response categories
    obs.ratings = the_data[,.(r = mean(response)), by = .(response.c)][,r] %>% sort()
    batches = (fit$fit@sim$iter-fit$fit@sim$warmup)*4/500
    pp_basic_l = vector(mode = "list", length = batches)
    pp_l = vector(mode = "list", length = batches)
    for(j in 1:batches) { # circumvent memory problems
      my_draw_ids = as.integer((1:500)+((j-1)*500))
      epred = posterior_epred(fit, newdata = new_data, draw_ids = my_draw_ids)
      # calculate expected score per row in data
      epred.r = 
        apply(epred, 1, function(a) apply(a, 1, function(p) sum(p*obs.ratings)))
      rm(epred)
      gc()
      
      pp_basic_l[[j]] = 
        cbind(epred.data, epred.r) %>% 
        melt(id.vars = names(epred.data), variable.name = "iter") %>% 
        .[, .(value = collapse::fmean(value)), by = .(Sex,Drug,Stress, reminder,iter)] %>% 
        .[, iter := NULL] %>% 
        .[, iter := my_draw_ids, by = .(Sex,Drug,Stress, reminder)] %>% 
        dcast(iter ~ Sex + Drug + Stress + reminder, value.var = "value", sep = ":") %>% 
        # stress-effect = post-pre difference
        .[, .(m_pp_plac_contr = `men:placebo:control:post` - `men:placebo:control:pre`,
              m_pp_plac_stress = `men:placebo:stress:post` - `men:placebo:stress:pre`,
              m_pp_oxy_contr = `men:oxycodone:control:post` - `men:oxycodone:control:pre`,
              m_pp_oxy_stress = `men:oxycodone:stress:post` - `men:oxycodone:stress:pre`,
              f_pp_plac_contr = `women:placebo:control:post` - `women:placebo:control:pre`,
              f_pp_plac_stress = `women:placebo:stress:post` - `women:placebo:stress:pre`,
              f_pp_oxy_contr = `women:oxycodone:control:post` - `women:oxycodone:control:pre`,
              f_pp_oxy_stress = `women:oxycodone:stress:post` - `women:oxycodone:stress:pre`),
          by = .(iter)]
      
      tmp = 
        pp_basic_l[[j]] %>% 
        # stress-effects by sex and drug
        .[, .(m_plac_sr = m_pp_plac_stress - m_pp_plac_contr,
              m_oxy_sr = m_pp_oxy_stress - m_pp_oxy_contr,
              f_plac_sr = f_pp_plac_stress - f_pp_plac_contr,
              f_oxy_sr = f_pp_oxy_stress - f_pp_oxy_contr), by = .(iter)] %>% 
        # state-condition X stress interaction 
        .[, .(men = m_oxy_sr - m_plac_sr,
              women = f_oxy_sr - f_plac_sr), by = .(iter)] %>% 
        .[, avg := (men+women)/2] %>% 
        .[, `women-men` := (women-men)] %>% 
        melt(id.var = "iter", value.name = "buffer_drug_effect", variable.name = "Sex")
      pp_l[[j]] = tmp
      rm(tmp)
      gc()
    }
    
    pp = do.call(rbind,pp_l)
    pp_basic = do.call(rbind,pp_basic_l)
    rm(pp_l,pp_basic_l)
    gc()
    save(pp,pp_basic, file = fn)
  } else {
    load(fn)
  }
  
  tmp =
    pp_basic %>% 
    melt(id.vars = "iter") %>% 
    .[, variable := gsub("_pp","",variable)] %>% 
    .[, c("Sex","Drug","state.condition") := tstrsplit(variable,"_")] %>% 
    .[, Sex := ifelse(Sex == "m","men","women")] %>% 
    .[, Drug := ifelse(Drug == "plac","placebo","oxycodone")] %>% 
    .[, state.condition := ifelse(state.condition == "stress","stress","control")] %>% 
    .[, variable := NULL]
  
  obs_diff = 
    the_data %>% 
    dcast(participant + Drug + state.condition + Sex ~ stage, value.var = "response") %>% 
    .[, .(diff = mean(`6`-`5`)), by = .(Drug, state.condition, Sex)] 
  
  ppc = 
    obs_diff %>% 
    ggplot(aes(x = state.condition, y=diff, color = Drug)) + 
    geom_point(position = position_dodge(0.2)) + 
    facet_wrap(~Sex) + 
    geom_line(aes(group=Drug),position = position_dodge(0.2)) + 
    stat_interval(data = tmp[, diff := value],alpha = .25,
                  position = position_dodge(0.2)) + 
    ylab("post-pre reminder difference")
  
  pp_cells =
    tmp %>% 
    dcast(iter + Drug + Sex ~ state.condition) %>% 
    .[, value := stress-control] %>% 
    dcast(iter + Sex ~ Drug) %>% 
    .[, .(stress = (oxycodone+placebo)/2,
          `stress in placebo` = placebo,
          `stress in oxycodone` = oxycodone,
          buffer = oxycodone-placebo), by = .(iter,Sex)] %>% 
    melt(id.vars = c("iter","Sex")) %>% 
    dcast(iter + variable ~ Sex) %>% 
    .[, total := (men+women)/2] %>% 
    .[, `Sex difference` := women-men] %>% 
    melt(id.vars = c("iter","variable"), variable.name = "Sex") %>% 
    setnames("variable","contrast")
  
  stress_effects = 
    pp_basic %>% 
    .[, .(m_plac_sr = m_pp_plac_stress - m_pp_plac_contr,
          m_oxy_sr = m_pp_oxy_stress - m_pp_oxy_contr,
          f_plac_sr = f_pp_plac_stress - f_pp_plac_contr,
          f_oxy_sr = f_pp_oxy_stress - f_pp_oxy_contr), by = .(iter)] %>% 
    .[, iter := 1:nrow(pp_basic)] %>% 
    melt(id.var = "iter", value.name = "stress_effect") %>% 
    .[, c("Sex","Drug","x") := tstrsplit(variable,"_")] %>% 
    .[, Sex := ifelse(Sex == "m","men","women")] %>% 
    .[, Drug := ifelse(Drug == "plac","placebo","oxycodone")] %>% 
    .[, x := NULL]
  
  stress.plot =
    stress_effects %>% 
    ggplot(aes(x = Sex, y = stress_effect, color = Drug)) +
    stat_pointinterval(position = position_dodge(.2)) + 
    geom_hline(yintercept = 0, col = "red", lty = 2) + 
    scale_col_drug
  
  buffer.plot =
    pp %>% 
    .[, Sex := gsub("man","men",Sex)] %>% 
    .[Sex == "women-men", Sex := "Sex difference"] %>% 
    .[Sex == "avg", Sex := "total"] %>% 
    .[, Sex := factor(Sex, levels = c("men","women","total","Sex difference"))] %>% 
    ggplot(aes(x = buffer_drug_effect, fill = Sex)) + 
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = c(-5,5), lty = 2) +
    stat_halfeye(alpha = .5) + 
    facet_wrap(~Sex, scales = "free") +
    ylab("") +
    gg_no_y_axis +
    scale_col_sex.b
  
  buffer.stats = 
    pp_cells %>% 
    .[, .(stats = get_mci(value)), by = .(Sex,contrast)] %>% 
    dcast(Sex ~ contrast, value.var = "stats")
  
  return(list(buffer.stats = buffer.stats,
              stress.plot = stress.plot,
              buffer.plot = buffer.plot,
              ppc = ppc))
}

get_timing_data = function() {
  #load("data/Timings_from_Martin.Rdata")
  #TS = Timing_summary %>% data.table()
  #setkeyv(TS,"Activity")
  #save(TS,file="data/Times_from_Martin.Rdata")
  #rm(Timing_summary)
  load("data/Times_from_Martin.Rdata")
  tx = TS[Activity %in% c("state_prepmainmath","drug_1_admin","drug_2_admin", "reminder_1","reminder_2"),
          .(Activity, Activity_mid, Activity_start, Activity_end)] %>% 
    setnames(c("Activity_start", "Activity_end","Activity_mid"),c("xmin","xmax","time")) %>% 
    .[!is.na(xmin), xmid := NA] %>% 
    .[, .(time = mean(time), xmin = min(xmin), xmax = max(xmax)), by = .(Activity)]
  
  return(list(tx = tx, TS = TS))
}

get_timing_data() %>% list2env(envir = .GlobalEnv)

plot_stages = function(dt, my_stages = 1:10, by.Sex = FALSE, ymin.min = 0, by.Drug = FALSE, 
                       fill = "state.condition", shape = "state.condition", lty = "state.condition") {
  bw = NULL
  if (by.Sex != FALSE) bw = "Sex"
  wi = c("state.condition","stage")
  if (by.Drug != FALSE) wi = c(wi,"Drug")
  pdata = 
    seWithin(
      data = dt[stage %in% my_stages][, state.condition := factor(state.condition)][,stage := factor(stage)],
      measurevar = "response",
      withinvars = wi,
      betweenvars = bw,
      idvar = "participant",
      conf.interval = .95
    ) %>% 
    data.table() %>% 
    .[, Activity := paste0("Q",stage)] %>% setkeyv("Activity") %>% 
    .[TS, time := Activity_start] %>% 
    .[, `:=`(lower = response-2*se, upper = response+2*se)]
  
  
  gg.aes = aes(x = time, y = response, shape = .data[[shape]], color = .data[[fill]], 
                 fill = .data[[fill]], lty = .data[[shape]])
  scale_col = scale_col_stress
  if (by.Drug == TRUE) {
    scale_col = scale_col_drug
  }
  
  
  tx[, `:=`(ymin = min(c(ymin.min,pdata$lower),na.rm = TRUE)*0.5,
            ymax = min(c(100,max(pdata$upper)),na.rm = TRUE)*1.5,
            response = mean(pdata$response), 
            state.condition = pdata$state.condition[1],
            Drug = "placebo",
            Sex = "women")]
  tx = rbind(copy(tx)[, Sex := "men"],tx[, Sex := "women"])
  p = 
    pdata %>% 
    ggplot(gg.aes) + 
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax, group = Activity),
              data = tx, color = NA, fill = "gray",
              alpha = .5) + 
    geom_line() + 
    geom_point() + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2, color = NA) + 
    scale_x_continuous(breaks = seq(-60,80,20)) + 
    scale_y_continuous(expand=c(0,0)) +
    geom_linerange(aes(x = time, ymin = ymin, ymax = ymax), color = "black",
                   data = tx[grepl("drug",Activity)], lty = 2) + 
    theme(panel.grid.major.y = element_blank(),
          axis.line = element_line(),
          legend.key.size = unit(1,"line")) +
    coord_cartesian(ylim = c(min(ymin.min,pdata$lower, na.rm = TRUE),max(pdata$upper))) + 
    theme(axis.line = element_line()) +
    scale_col
  
    yrange = c(min(pdata$lower),max(pdata$upper))
    
    if (is.na(ymin.min))
      p = p + 
      guides(y = guide_axis_truncated(
        trunc_upper = yrange[2],
        trunc_lower = yrange[1]+diff(yrange)/500))
    p
}

# marie added this to look at drug effects quickly (also changed the analysis script a bit)
plot_stages_drug = function(dt, my_stages = 1:10, by.Drug = FALSE) {
  
  if (by.Drug == FALSE) {
    dt.by = c("stage","state.condition")
    gg.aes = aes(x = stage, y = m, lty = state.condition, color = state.condition, fill = state.condition)
    scale_col = scale_col_stress
  } else {
    dt.by = c("stage","Drug","state.condition")
    gg.aes = aes(x = stage, y = m, color = Drug, lty = state.condition, fill = Drug)
    scale_col = scale_col_drug
  }
  
  dt %>% 
    .[stage %in% my_stages] %>% 
    .[,.(m = median(response, na.rm = TRUE),
         se = sd(response, na.rm = TRUE)/sqrt(length(unique(participant)))), 
      by = dt.by] %>% 
    .[, `:=`(lower = m-2*se, upper = m+2*se)] %>% 
    ggplot(gg.aes) + 
    geom_line() + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2, color = NA) + 
    scale_x_continuous(breaks = my_stages, labels = stage.names[my_stages]) + 
    scale_col
}

make_the_data = function(dt, phase = "induction", n.post.stages = 1) {
  if (phase == "induction") {
    pre.stage = 3
    stages = pre.stage:(pre.stage + n.post.stages) %>% as.character()
  } else if (phase == "reminder") {
    pre.stage = c(5,7)
    stages = c(5:8)  %>% as.character()
  }  else if (phase == "reminder1") {
    pre.stage = c(5)
    stages = c(5:6)  %>% as.character()
  } else if (phase == "drug") {
    pre.stage = c(4)
    stages = c(4:5)  %>% as.character()
  } else if (phase == "stressdrug") {
    pre.stage = 1:3
    stages = c(1:3,5)  %>% as.character()
  }
  
  
  stage.nms = paste0(c("pre.","post."),phase)
  dt = 
    dt %>% 
    .[stage %in% stages] %>% 
    .[, orig.stage := stage] %>% 
    .[, stage := as.character(stage)] %>% 
    .[, stage := ifelse(stage %in% pre.stage,stage.nms[1],stage.nms[2])] %>% 
    .[, stage := factor(stage, levels = stage.nms)]
  
  tbl = 
    dt[, .(N = length(unique(participant))), by = .(stage, Drug, state.condition)] %>% 
    dcast(Drug + state.condition ~ stage, value.var = "N") %>% 
    kable(caption = "Number of particicipants per condition by condition and stage") %>% 
    kable_styling(full_width = FALSE) %>% 
    kable_classic()
  
  attr(dt,"N") = tbl
  
  return(dt)
}

plot_pre_post = function(dt, phase = "induction", by.Sex = FALSE) {
  
  if (phase == "drug") {
    p.aes = aes(x = x, y = response, color = Drug,
                group = factor(participant):factor(session))
    p.wrap = facet_wrap(~Drug)
    if (by.Sex == TRUE)
      p.wrap = facet_grid(Sex ~ Drug)
    scale_col = scale_col_drug
  } else if (phase == "drug:induction") {
    p.aes = aes(x = x, y = response, color = Drug,
                group = factor(participant):factor(session):factor(state.condition))
    p.wrap = facet_wrap(~state.condition)
    if (by.Sex == TRUE)
      p.wrap = facet_grid(Sex ~ Drug)
    scale_col = scale_col_drug
  } else {
    p.aes = aes(x = x, y = response, color = state.condition,
                group = factor(participant):factor(session))
    p.wrap = facet_wrap(~state.condition)
    if (by.Sex == TRUE)
      p.wrap = facet_grid(Sex ~ state.condition)
    scale_col = scale_col_sex
  }
  
  p = 
    dt %>% 
    .[, .(response = mean(response)),
      by = .(stage,Sex,state.condition, participant, session, Drug)] %>% 
    .[, x := ifelse(grepl("pre",stage),1,2)] %>% 
    .[, x := x + runif(1,-.1,.1), by = .(participant, session)] %>% 
    ggplot(p.aes) + 
    geom_line(alpha = .25) + 
    geom_point(alpha = .5) +
    scale_x_continuous(breaks = c(1,2), labels = paste(c("pre","post"),phase)) +
    xlab("") + 
    p.wrap + 
    ylim(0,100) + 
    scale_col
  
  return(p)
} 

plot_drug_effect_raw = function(the_data) {
  the_data %>%  
    .[, .(response = mean(response)),
      by = .(stage,Sex,state.condition, participant, session, Drug)] %>% 
    dcast(participant + Sex + Drug + state.condition ~ stage, value.var = "response")  %>% 
    .[, stage_diff := post.drug - pre.drug] %>% 
    dcast(participant + Sex + state.condition ~ Drug, value.var = "stage_diff")  %>%
    .[, drug_effect := oxycodone - placebo] %>% 
    ggplot(aes(x = drug_effect, fill = Sex)) + 
    facet_wrap(~state.condition) +
    geom_density(alpha = .25, color = NA) + 
    geom_vline(xintercept = 0) +
    ggtitle("drug_effect = rating_in_oxycodone-rating_in_placebo") + 
    scale_col_sex
}

plot_stresseff_raw = function(the_data, by.Sex = TRUE) {
  
  if (by.Sex == FALSE) {
    dt.by = c("stage","state.condition", "participant", "session", "Drug")
    dt.cast1 = participant + Drug + state.condition ~ stage
    dt.cast2 = participant + Drug ~ state.condition
    by.tbl = c("Drug")
    gg.aes = aes(x = stress_effect, fill = "1")
    scale_col = list(
      scale_fill_manual(values = c("grey25")), 
      scale_color_manual(values = c("grey25"))  
    )
  } else {
    dt.by = c("stage","Sex","state.condition", "participant", "session", "Drug")
    dt.cast1 = participant + Sex + Drug + state.condition ~ stage
    dt.cast2 = participant + Sex + Drug ~ state.condition
    by.tbl = c("Drug","Sex")
    gg.aes = aes(x = stress_effect, fill = Sex)
    scale_col = scale_col_sex
  }
  
  postpre = unique(the_data$stage) %>% as.character() %>% sort()
  pdata = 
    the_data %>%  
    .[, .(response = mean(response)), by = dt.by] %>% 
    dcast(dt.cast1, value.var = "response")  %>% 
    .[, stage_diff := get(postpre[2]) - get(postpre[1])] %>% 
    dcast(dt.cast2, value.var = "stage_diff")  %>%
    .[, stress_effect := control - stress] 
  
  tbl = pdata %>% 
    seWithin(
      measurevar = "stress_effect",
      withinvars = by.tbl,
      idvar = "participant")
  
  p = 
    pdata %>% 
    ggplot(gg.aes) + 
    facet_wrap(~Drug) +
    geom_density(alpha = .25, color = NA) + 
    geom_vline(xintercept = 0) +
    ggtitle("stress effect = rating_in_stress-rating_in_control") + 
    scale_col + 
    xlab("stress effect")
  
  attr(p,"table") = tbl
  
  return(p)
}

plot_stresseff_raw2 = function(the_data, by.Sex = TRUE) {
  
  if (by.Sex == FALSE) {
    dt.by = c("stage","state.condition", "participant", "session", "Drug")
    dt.by2 = c("Drug")
    dt.cast1 = participant + Drug + state.condition ~ stage
    dt.cast2 = participant + Drug ~ state.condition
    gg.aes = aes(x = Drug, y = stress_effect, color = Drug, fill = Drug)
    scale_col = list(
      scale_fill_manual(values = c("grey25")), 
      scale_color_manual(values = c("grey25"))  
    )
  } else {
    dt.by = c("stage","Sex","state.condition", "participant", "session", "Drug")
    dt.by2 = c("Sex", "Drug")
    dt.cast1 = participant + Sex + Drug + state.condition ~ stage
    dt.cast2 = participant + Sex + Drug ~ state.condition
    gg.aes = aes(x = Drug, y = stress_effect, fill = Drug, color = Drug)
    scale_col = scale_col_drug
  }
  
  postpre = unique(the_data$stage) %>% as.character() %>% sort()
  pdata.i =
    the_data %>%  
    .[, .(response = mean(response, na.rm = TRUE)),
      by = .(stage,Sex,state.condition, participant, session, Drug)] %>% 
    dcast(dt.cast1, value.var = "response")  %>% 
    .[, stage_diff := get(postpre[2]) - get(postpre[1])]  %>% 
    dcast(dt.cast2, value.var = "stage_diff") %>% 
    .[, stress_effect := control - stress]
  pdata.stats = 
    pdata.i %>% 
    .[, .(stress_effect = mean(stress_effect, na.rm = TRUE),
          sd = sd(control - stress, na.rm = TRUE),
          N = length(unique(participant))), by  = dt.by2] %>% 
    .[, se := sd/sqrt(N)] %>% 
    .[, lower := stress_effect - 2*se] %>% 
    .[, upper := stress_effect + 2*se] 
  
  p = 
    pdata.stats %>% 
    ggplot(gg.aes) + 
    geom_bar(stat = "identity", alpha = .75, color = NA) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), linewidth = 1, color = "black", width = .2) + 
    geom_quasirandom(data = pdata.i) +
    scale_col_drug + 
    guides(fill="none") + 
    ylab("Effect of stress induction")
  if (by.Sex == TRUE) p = p + facet_wrap(~Sex)
  return(p)
  
}

plot_prepost_contrast = function(fit,the_data,outcome, my_contrast = NULL) {
  stopifnot(!is.null(my_contrast))
  my_effect = paste0(my_contrast,"_effect")
  tmp = fit %>% as_draws_df() %>% .[1,1] %>% as.numeric() %>% round(5)*100000
  fn = paste0("brmsP1/",outcome,"_contrasts_",
              paste(sort(unique(the_data$stage)), collapse = "_"),"_",
              "_",tmp,".Rdata")
  if (file.exists(fn)) {
    load(fn)
  } else {
    # names of pre and post manipulation stages
    postpre = unique(the_data$stage) %>% as.character() %>% sort()
    # posterior expectations, here probabilities for response categories
    new_data = do.call(rbind,lapply(1:4, function(x) data.table(fit$data)[, session := x]))
    epred = posterior_epred(fit, newdata = new_data)
    # average response on 0-100 scale in response categories
    obs.ratings = the_data[,.(r = mean(response)), by = .(response.c)][,r] %>% sort()
    
    # calculate expected score per row in data
    epred.r = 
      apply(epred, 1, function(a) apply(a, 1, function(p) sum(p*obs.ratings)))
    
    # posterior predictive (sanity) check
    # plot(the_data$response, rowMeans(epred.r))
    
    # get relevant condition info from the data
    if (my_contrast == "stress") {
      epred.data = 
        fit$data %>% data.table() %>% 
        .[, .(induced.state, stage, Sex,participant, session)] %>% 
        .[, state.condition := ifelse(any(induced.state == "stress"),"stress","control"),
          by = .(participant,session)] %>% 
        .[, induced.state := NULL] %>% 
        .[, condition := state.condition]
      c.1 = "stress"
      c.2 = "control"
    } else if (my_contrast == "drug") {
      epred.data = 
        fit$data %>% data.table() %>% 
        .[, .(Drug, stage, Sex, participant, session)] %>% 
        .[, condition := Drug]
      c.1 = "oxycodone"
      c.2 = "placebo"
    }
    
    # merge posterior predictions with condition info from the data
    # calculate stress effects and 
    pp = 
      cbind(epred.data, epred.r) %>% 
      melt(id.vars = names(epred.data), variable.name = "iter") %>% # everything in long format
      dcast(participant + Sex + condition + iter ~ stage, # pre/post manipulation stages in columns
            value.var = "value", fun.aggregate = mean)
    
    pp = 
      pp %>% 
      .[, stage_diff := get(postpre[1]) - get(postpre[2])] %>% # calculate pre/post manipulation scores
      dcast(participant + Sex + iter ~ condition, value.var = "stage_diff") %>%  # manipulation conditions in columns
      .[, effect := get(c.1) - get(c.2)] %>% # calculate effect of manipulation
      .[, .(effect = mean(effect),
            c.1 = mean(get(c.1)),
            c.2 = mean(get(c.2))),
        by = .(Sex,iter)]  %>% # effect of manipulation by sex
      setnames(c("c.1","c.2","effect"), c(c.1,c.2,my_effect)) %>% 
      melt(id.vars = c("Sex","iter"), variable.name = "contrast")
    
    # sex difference in effect
    sex.diff = 
      pp %>% 
      dcast(iter + contrast ~ Sex, value.var = "value") %>% 
      .[, value := women-men] %>% 
      .[,.(iter,contrast,value)] %>% 
      .[, Sex := "women-men"]
    
    # average effect over sexes
    sex.avg = 
      pp %>% 
      dcast(iter + contrast ~ Sex, value.var = "value") %>% 
      .[, value := (women+men)/2] %>% 
      .[,.(iter,contrast,value)] %>% 
      .[, Sex := "avg"]
    
    # put all effects in one dt
    pp = 
      rbind(pp, sex.diff, sex.avg) %>% 
      .[, Sex := factor(Sex, levels = c("men","women","women-men","avg"))]
    save(pp,file = fn)
  }
  ## rename some values
  new.nms = c("men","women","total","Sex difference")
  pp %>% 
    .[, Sex := gsub("man","men",Sex)] %>% 
    .[Sex == "avg", Sex := "total"] %>% 
    .[Sex == "women-men", Sex := "Sex difference"] %>% 
    .[, Sex := factor(Sex, levels = new.nms)]
  
  post = pp[contrast == my_effect]
  p =
    post %>% 
    ggplot(aes(x = value, color = Sex, fill = Sex)) + 
    geom_vline(xintercept = 0) +
    stat_halfeye(alpha = .5) + 
    facet_wrap(~Sex, scales = "free") +
    ylab("") + 
    xlab(paste(my_effect,"for",outcome)) +
    gg_no_y_axis +
    scale_col_sex.b
  my_caption = 
    ifelse(my_contrast == "stress",
           "Modeled post minus pre state induction on xx ratings by condition (control, stress) and the stress effect by sex.",
           "Modeled post minus pre drug administration on xx ratings by condition (placebo, oxycodone) and the drug effect by sex.")
  my_caption = gsub("xx",outcome,my_caption)
  
  tbl = 
    pp %>% 
    .[, .(stats = get_mci(value, get.P = FALSE)), by = .(Sex,contrast)] %>% 
    .[, contrast := factor(contrast,levels = c("placebo","control","oxycodone","stress",my_effect))] %>% 
    dcast(Sex~contrast, value.var = "stats") %>% 
    kable(caption = my_caption) %>% 
    kable_styling(full_width = FALSE) %>% 
    kable_classic()
  
  summary_tbl = 
    pp %>% 
    .[, .(stats = get_mci(value, get.P = ifelse(grepl("-",Sex),TRUE,FALSE))), by = .(Sex,contrast)] %>% 
    .[, contrast := factor(contrast,levels = c("placebo","control","oxycodone","stress",my_effect))] %>% 
    dcast(contrast~Sex, value.var = "stats") %>% 
    .[contrast == my_effect] %>% 
    .[,new.nms,with = FALSE]
  summary_tbl = cbind(Outcome = outcome, summary_tbl)
  
  return(list(plot = p,post_stress_effect = post, table = tbl, summary_table = summary_tbl))
}

plot_drug = function(the_data) {
  the_data %>% 
    .[, .(r = mean(response)), by = .(participant, Sex, state.condition, Drug)] %>% 
    dcast(participant + Sex + state.condition ~ Drug, value.var = "r") %>% 
    .[, drug_effect := oxycodone-placebo] %>% 
    ggplot(aes(x = drug_effect, fill = Sex)) + 
    facet_wrap(~state.condition) +
    geom_density(alpha = .25, color = NA) + 
    geom_vline(xintercept = 0) +
    ggtitle("drug_effect = rating_in_oxycodone-rating_in_placebo") + 
    scale_col_sex
}

plot_drug_lines = function(the_data) {
  if (length(unique(the_data$item_name)) > 1) {
    p.data = the_data[, response := mean(response), by = .(participant,Drug,stage,state.condition)]
  } else {
    p.data = the_data
  }
    
  p.data %>% 
    .[, Drug := factor(Drug, levels = c("placebo","oxycodone"))] %>% 
    ggplot(aes(x = Drug, y = response, color = Sex, group = factor(participant):factor(stage))) + 
    geom_line(alpha = .5) +
    facet_wrap(~state.condition) + 
    scale_col_sex
}

get_drug_contrasts = function(fit, the_data) {
  fn = paste0("brmsP1/drug/contrasts_",the_data$item_name[1],"_",fit$fit@sim$samples[[1]][1,1],".Rdata")
  if (file.exists(fn)) {
    load(fn)
  } else {
    new_data = do.call(rbind,lapply(1:4, function(x) data.table(fit$data)[, session := x]))
    epred = posterior_epred(fit, newdata = new_data)
    obs.ratings = 
      the_data[,.(r = mean(response)),
               by = .(response.c)][,r] %>% sort()
    
    epred.r = apply(epred, 1, function(a) apply(a, 1, function(p) sum(p*obs.ratings)))
    epred.data = 
      fit$data %>% data.table() %>% 
      .[, .(state.condition, Sex,participant, session, Drug)]
    ppx = cbind(epred.data, epred.r) %>% 
      melt(id.vars = names(epred.data), variable.name = "iter") %>% 
      .[, .(value = mean(value)), by = .(state.condition, Sex, Drug, iter)] %>% 
      .[, .(mCI = get_mci(value,get.P = FALSE)), by = .(state.condition, Sex,Drug)] %>% 
      dcast(Sex + state.condition ~ Drug, value.var = "mCI")
    pp = cbind(epred.data, epred.r) %>% 
      melt(id.vars = names(epred.data), variable.name = "iter") %>% 
      .[, .(value = mean(value)), by = .(state.condition, Sex, Drug, iter)] %>% 
      dcast(iter ~ Drug + state.condition + Sex, value.var = "value") %>% 
      # calculate drug effects:
      .[, .(women_in_control = oxycodone_control_women - placebo_control_women,
            men_in_control = oxycodone_control_men - placebo_control_men,
            women_in_stress = oxycodone_stress_women - placebo_stress_women,
            men_in_stress = oxycodone_stress_men - placebo_stress_men),
        by = .(iter)] %>% 
      .[, `:=`(in_women = (women_in_control + women_in_stress)/2,
               in_men = (men_in_control + men_in_stress)/2,
               in_stress = (women_in_stress+men_in_stress)/2,
               in_control = (women_in_control+men_in_control)/2)] %>% 
      .[, `:=`(women_vs_men = in_women-in_men,
               stress_vs_ctrl = in_stress-in_control,
               stress_vs_ctrl_in_men = men_in_stress-men_in_control,
               stress_vs_ctrl_in_women = women_in_stress-women_in_control,
               women_vs_men_in_ctrl = women_in_control-men_in_control,
               women_vs_men_in_stress = women_in_stress-men_in_stress)] %>% 
      .[, three_way := women_vs_men_in_ctrl-women_vs_men_in_stress] %>% 
      melt(id.vars = "iter", variable.name = "effect", value.name = "est")
    attr(pp,"cell_stats") = ppx
    save(pp,file = fn)
  }
  attr(pp,"outcome") = paste(unique(the_data$item_name),collapse = ", ")
  pp[, effect := gsub("man","men",effect)]
  attr(pp,"cell_stats")[, Sex  := gsub("man","men",Sex )]
  return(pp)
}


drug_plotter = function(post.contrasts) {
  dplotter = function(dt, title) {
    dt %>% ggplot(aes(x = est)) + 
      stat_halfeye() + 
      facet_wrap(~effect, scale = "free_x") + 
      geom_vline(xintercept = 0) + 
      ggtitle(title)
  }
  ox = attr(post.contrasts,"outcome")
  p1 = post.contrasts %>% 
    .[effect %in% c("in_women","in_men")] %>% 
    .[, effect := gsub("in_","",effect)] %>% 
    ggplot(aes(x = est, fill = effect)) +
    stat_halfeye(alpha = .5) + 
    geom_vline(xintercept = 0) + 
    ggtitle(paste0("Drug effect by sex on '",ox,"'")) + 
    scale_col_sex + 
    ylab("") +
    gg_no_y_axis
  p2 = post.contrasts %>% 
    .[effect %in% c("in_control","in_stress")] %>% 
    .[, effect := gsub("in_","",effect)] %>% 
    ggplot(aes(x = est, fill = effect)) +
    stat_halfeye(alpha = .5) + 
    geom_vline(xintercept = 0) + 
    ggtitle(paste0("Drug effect by stress on '",ox,"'")) + 
    scale_col_stress + 
    ylab("") +
    gg_no_y_axis
  
  pA = p1|p2
  
  p3 = post.contrasts[grepl("vs",effect)] %>% 
    .[, effect := gsub("_vs_"," - ", effect)] %>%
    .[, effect := gsub("_in_"," in ", effect)] %>% 
    .[, effect := gsub("_"," in ", effect)] %>% 
    dplotter(paste0("Contrasts and interactions for '",ox,"'")) + 
    theme(strip.text = element_text(size = 6), axis.text = element_text(size = 6)) + 
    ylab("") +
    gg_no_y_axis
  p4 = post.contrasts[grepl("way",effect)] %>% 
    dplotter("Three way interaction:\nwoman_vs_man_in_ctrl-\nwoman_vs_man_in_stress") +
    theme(plot.title = element_text(size = 10)) + 
    ylab("") +
    gg_no_y_axis
  
  pB = (p3|p4) + plot_layout(widths = c(2,1))
  
  return(
    (pA / pB) + plot_layout(heights = c(1.2,2))
  )
}

get_drug_effects.stats = function(post.contrasts, my_outcome) {
  old.nms = c("in_men","in_women","total","women_vs_men")
  new.nms = c("Men","Women","Total","Sex difference")
  tmp = 
    post.contrasts %>% 
    .[effect %in% c("in_men","in_women","women_vs_men")] %>% 
    dcast(iter ~ effect, value.var = "est") %>%
    .[, total := (in_women+in_men)/2] %>% 
    melt(id.vars = "iter") %>% 
    .[,.(stats = get_mci(value, get.P = ifelse(grepl("vs",variable),TRUE,FALSE))), by = .(variable)] %>% 
    dcast(.~variable) %>% 
    setnames(old.nms,new.nms) %>% 
    .[, new.nms, with = FALSE]
  
  return(cbind(Outcome = my_outcome,tmp))
}

get_drug_effects.stats.b = function(post.contrasts, my_outcome) {
  tmp = 
    post.contrasts %>% 
    .[grepl("^men_in|^women_in", effect)] %>% 
    .[,.(stats = get_mci(est, get.P = FALSE)), by = .(effect)] %>% 
    .[, c("Sex","state.condition") := tstrsplit(effect,"_in_")] %>% 
    .[, stats := gsub("\\]|;","",stats)] %>% 
    .[, c("drug_effect","lower","upper") := tstrsplit(stats,", | \\[")] %>% 
    .[, `:=`(drug_effect = as.numeric(drug_effect), lower = as.numeric(lower), upper = as.numeric(upper))] %>% 
    .[, c("effect","stats") := NULL] %>% 
    .[, Outcome := my_outcome]
  return(tmp)
}

add_cumthresh_prior = function(my_prior = c(),responses) {
  get_thresh_prior_mean = function(x) {
    t_starts = x %>% table() %>% prop.table() %>% cumsum() 
    return(t_starts[1:(length(unique(x))-1)] %>% qlogis())
  }
  t_starts = get_thresh_prior_mean(responses)
  for (k in 1:length(t_starts)) {
    txt = paste0("prior(normal(",t_starts[k],",1), class = Intercept, coef = ",k,")")
    my_prior = my_prior + eval(parse(text = txt))
  }
  return(my_prior)
}

get_acat_thresh_priors = function(y,my_prior,sd = 2) {
  cat.counts = table(y)
  cat.counts = cat.counts[cat.counts > 0]
  n.cat = length(cat.counts)
  N = sum(cat.counts)
  t_starts = log((head(cat.counts,n.cat-1)/N)/(tail(cat.counts,n.cat-1)/N)) %>% round(2)
  for (k in 1:(n.cat-1)) {
    txt = paste0("prior(normal(",t_starts[k],",",sd,"), class = Intercept, coef = ",k,")")
    my_prior = my_prior + eval(parse(text = txt))
  }
  return(my_prior)
}

fit_ordinal_model = function(the_data, the_formula, fn, family = cumulative()) {
  if (family$family == "cumulative") {
    # get good priors for thresholds
    t_starts = get_thresh_prior_mean(the_data[,response.c])
    my_prior = prior(normal(0,2), class = b)
    for (k in 1:length(t_starts)) {
      txt = paste0("prior(normal(",t_starts[k],",1), class = Intercept, coef = ",k,")")
      my_prior = my_prior + eval(parse(text = txt))
    }
    control = list(adapt_delta = .8)
    warmup = 1000
    iter = 2000
  } else if (family$family == "acat") {
    my_prior = prior(normal(0,2), class = b) + prior(normal(0,2), class = "sd")
    my_prior = get_acat_thresh_priors(the_data$response.c,my_prior)
    control = list(adapt_delta = .9)
    warmup = 1000
    iter = 2000
  }
  # run model
  fit = brm(
    the_formula,
    family = family,
    data = the_data,
    prior = my_prior,
    backend = "cmdstanr",
    control=control,
    threads = threading(2, grainsize = 168),
    iter = iter,
    warmup = warmup)
  save(fit,file = fn)
  return(fit)
}
select_data = function(items) {
  my_data.Q[item_name %in% items, 
            .(participant, session, stage, response, state.condition, 
              induced.state, Sex, Drug, item_name)] %>% 
    .[, session := factor(session)]
}