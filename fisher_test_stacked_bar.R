library(ggstatsplot)

# df = your df, cat_var_1 and 2 = categorical variables to compare

test <- fisher.test(table(df$cat_var_1, df$cat_var2))

ggbarstats(
  df, cat_var1, cat_var2,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
  )
)