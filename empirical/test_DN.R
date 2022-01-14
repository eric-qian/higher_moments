# Preliminaries ----
library('steadyICA')
library('readxl')
set.seed(2022)

# Settings
n = 3     
R = 5000  # Number of permutations

# Read shocks, from ICA procedure
shockPath = 'figures2019/initSetting=GlobalSearch_shocks.xls'
sheets    = excel_sheets(shockPath)

# Store results
Results = list()


for(jSheet in 1:length(sheets)){  # Loop through each sheet
  
  # Read shocks
  shocks = read_excel(shockPath, sheet = sheets[jSheet])
  S      = as.matrix( shocks[,2:4])
  
  pPerm = permTest(S, group = 1:n, R, FUN = 'gmultidcov', 
                   symmetric = FALSE, alpha=1)
  
  # Store results
  Results[[jSheet]] = list(shocks = shocks, S = S, pPerm = pPerm, type = sheets[jSheet])
}

# Store output
out        = as.data.frame(sapply(Results, function(x) x$type))
names(out) = 'type'
out$pPerm  = sapply(Results, function(x) x$pPerm)

write.csv(out, 'figures2019/initSetting=GlobalSearch_shocksPermTest.csv')

