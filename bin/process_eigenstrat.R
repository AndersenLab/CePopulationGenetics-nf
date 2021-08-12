
eval <- read.table("~/transfer/eigenstrat/eigenstrat_outliers_removed.eval") %>%
  dplyr::rename(Eigenvalue = V1) %>%
  dplyr::mutate(PC =  paste0("PC",c(1:n())),
                nPC = 1:n(),
                VE = round(Eigenvalue/sum(Eigenvalue)*100,digits=2)) %>%
  dplyr::mutate(C_VE = cumsum(VE))

evec <- read.table("~/transfer/eigenstrat/eigenstrat_outliers_removed.evac")

n.pcs <- ncol(evec)-2

colnames(evec) <- c("Strain",
                    paste0("PC",c(1:50)),
                    "Population")

processed_pcs <- evec %>%
  tidyr::gather(PC, Value, -Strain, -Population) %>%
  dplyr::left_join(., eval, by = "PC")

# plot eigenvalues

processed_pcs %>%
  dplyr::distinct(PC, VE, .keep_all = T) %>%
  ggplot() +
  geom_point(aes(x = nPC, y = Eigenvalue)) +
  geom_point(aes(x = nPC, y = VE), color = "red") +
  theme_bw(18) + 
  ylim(c(0,20)) +
  labs(x = "Principal Component")

pcx <- 1
pcy <- 2

x_pc <- dplyr::filter(processed_pcs, nPC == pcx) %>%
  dplyr::select(Strain, xPC = Value)
y_pc <- dplyr::filter(processed_pcs, nPC == pcy) %>%
  dplyr::select(Strain, yPC = Value) %>%
  dplyr::left_join(.,x_pc)

ggplot(y_pc) +
  aes(x = xPC, y = yPC) +
  geom_point() +
  theme_bw(18) +
  theme(panel.grid = element_blank()) +
  labs(x = glue::glue("Principal Component {pcx}"),
       y = glue::glue("Principal Component {pcy}"))


