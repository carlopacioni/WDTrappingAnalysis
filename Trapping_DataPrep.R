xls <- file.path(data_path,"TrappingSimulation.xlsx")
ID_sheet <- "Scats"
Trap_sheet <- "Trap"
Coords_sheet <- "Coords"
matches_sheet <- "Matches"

# IDs
IDs <- read_excel(xls, sheet = ID_sheet)
# subset to Mandfield
IDs <- IDs[IDs$Location == "Mansfield", ]
IDs <- IDs$ID

# Trap matrix
trap_matrix <- as.data.frame(read_excel(xls, sheet = Trap_sheet, skip = 1))

# Fixing issue with read_excel reading the dates as integer
names(trap_matrix)[-1] <- format(as.Date(as.numeric(tail(names(trap_matrix), -1)), 
                                         origin = "1899-12-30"), format="%d-%m-%Y")
names(trap_matrix)[1] <- "Sample_ID"

# Coords
coords_trap <- read_excel(xls, sheet = Coords_sheet)
coords_scats <- read_excel(xls, sheet = ID_sheet)

tmp <- coords_scats[, c("Long", "Lat")]
coordinates(tmp) <- c("Long", "Lat" )
proj4string(tmp) <- CRS("+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

coord_dec <- spTransform(tmp, CRS("+proj=longlat +datum=WGS84"))
coord_dec <-as.data.frame(coord_dec)
coords_scats$Long <- coord_dec$Long
coords_scats$Lat <- coord_dec$Lat

# Matches
Matches <- read_excel(xls, sheet = matches_sheet)
Matches$Date_chr <- as.character(format(Matches$Date, "%d-%m-%Y"))

# caphist df
nind <- length(IDs)
caphist <- as.integer(as.matrix(trap_matrix[, -1], ncol=ncol(trap_matrix)-1))
caphist_df <- data.frame(Y=caphist, 
                         Trap=rep(trap_matrix[, 1], ncol(trap_matrix)-1),
                         Day=rep(names(trap_matrix)[-1], each=nrow(trap_matrix)))
caphist_df <- caphist_df[!is.na(caphist_df$Y),]
caphist_m <- t(replicate(nind, caphist_df$Y))
dimnames(caphist_m) <- list(IDs, paste(caphist_df$Trap, caphist_df$Day, sep="_"))

# The for loop insert matches (trapped) animals, and fill in with NA after that day.
#------------------------------------------------------------------------------#
#   NEED TO REVISIT THIS BECAUSE IF AN ANIMAL IS TRAPPED IN ONE TRAP IS NOT 
#   AVAILABLE ANYMORE TO BE TRAPPED IN THE OTHERS, SO POSSIBLY ALL THE OTHERS 
#   SHOULD BE SET TO na ON THE DAY THEY ARE TRAPPED
#------------------------------------------------------------------------------#
for(i in seq_len(nrow(Matches))) {
  caphist_m[Matches$Matched_to[i], paste(Matches$Trap[i], Matches$Date_chr[i], sep="_")] <- 1
  d <- max(grep(dimnames(caphist_m)[[2]], pattern = Matches$Date_chr[i])) + 1
  caphist_m[Matches$Matched_to[i], d:ncol(caphist_m)] <- NA
}

nobs <- ncol(caphist_m) - apply(is.na(caphist_m), 1, sum)

# time
attack_date <- as.integer(
  as.Date(
    coords_scats[coords_scats$ID %in% row.names(caphist_m), ][["Date"]], 
    "%d-%m-%Y")
)
trap_date <- sapply(strsplit(colnames(caphist_m), "_"), "[", 2)
trap_date_m <- matrix(
  as.Date(rep(trap_date, length(IDs)), "%d-%m-%Y"), 
  nrow = length(IDs), byrow = TRUE)
time_diff <- ((trap_date_m - attack_date) + 1) / time_period

# Distance
attack_coords <- coords_scats[coords_scats$ID %in% row.names(caphist_m), c("Long", "Lat")]
trap_names <- sapply(strsplit(colnames(caphist_m), "_"), "[", 1)
trap_locs <- merge(data.frame(Site=trap_names, nr=seq_len(length(trap_names))), 
                   coords_trap)
trap_locs <- trap_locs[order(trap_locs$nr),]
d <- vector("list", length = nrow(attack_coords))
for(i in seq_len(nrow(attack_coords))) {
  d[[i]] <- calc.latlong.dist(attack_coords[i,], trap_locs[, c("lon", "lat")])
}
dist <- matrix(unlist(d), nrow = length(IDs), byrow = TRUE)


# Reformat Survival data as for JAGS example
surv_df <- data.frame(ID=rep(IDs, each=length(unique(trap_names))),
                      ind=rep(seq_along(IDs), each=length(unique(trap_names))),
                      Trap=rep(unique(trap_names), length(IDs)),
                      censored=0,
                      TimeTrap=NA,
                      t.cen=(rep(max(as.integer(as.Date(trap_date, "%d-%m-%Y"))) - attack_date,
                                each=length(unique(trap_names))) + 1) / time_period,
                      t.start=(rep(max(as.integer(as.Date(trap_date, "%d-%m-%Y"))) - attack_date,
                                   each=length(unique(trap_names))) + 2) / time_period)

# fixe numb of days for traps closed earlier
surv_df[surv_df$Trap=='MT09', 't.cen'] <- 
  ((as.integer(max(as.Date(trap_date, "%d-%m-%Y")[trap_names== 'MT09'])) - attack_date) + 1)/time_period
surv_df[surv_df$Trap=='MT10', 't.cen'] <- 
  ((as.integer(max(as.Date(trap_date, "%d-%m-%Y")[trap_names== 'MT10'])) - attack_date) + 1)/time_period

for(i in seq_len(nrow(Matches))) {
  surv_df[surv_df$ID==Matches$Matched_to[i] & surv_df$Trap==Matches$Trap[i], 
          c('censored', 'TimeTrap', 't.start')] <- 
    c(1, ((as.integer(as.Date(Matches$Date[i])) - 
        attack_date[which(coords_scats$ID == Matches$Matched_to[i])]) + 1)/time_period, NA)
}
surv_df$t.comb <- surv_df$t.cen
surv_df[!is.na(surv_df$TimeTrap), "t.comb"] <- surv_df[!is.na(surv_df$TimeTrap), "TimeTrap"]

# dist as vector
attack_coords_IDs <- cbind(IDs, attack_coords)
dist_vect <- vector("numeric")
for(r in seq_len(nrow(surv_df))) {
  dist_vect[r] <- calc.latlong.dist(attack_coords_IDs[attack_coords_IDs$IDs==surv_df[r, ][["ID"]], c("Long", "Lat")],
                                    coords_trap[coords_trap$Site == surv_df[r, ][["Trap"]], c("lon", "lat")])
  
}

ind <- as.numeric(unclass(as.factor(surv_df$ID)))
nb0 <- length(unique(ind))