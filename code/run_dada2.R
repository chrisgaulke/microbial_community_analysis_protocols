library(dada2)
library(argparser)

main <- function(){

p <- arg_parser("run_dada2.R")
p <- add_argument(p, "in_dir", help="the file path to the directory containing input fastq files")
p <- add_argument(p, "filt_dir", help="a directory for filtered files. If this does not exist it will be created.")
p <- add_argument(p, "out_dir", help="the file path to an output directory")
p <- add_argument(p, "--nprocs", help="number of cores", default=1)
p <- add_argument(p, "--tax_path",
					 help="path to taxonomy",
					 default= "/nfs3/Sharpton_Lab/tmp/data/public/databases/taxonomy/dada2_silva128/silva_nr_v128_train_set.fa.gz")

args <- parse_args(p)


###########################
#                         #
#    Filter and trim      #
#                         #
###########################

# Filename parsing
#path <- "/path/to/FWD" # CHANGE ME to the directory containing your demultiplexed fastq files
path <- args$in_dir
#filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
filtpath <- args$filt_dir
fns <- list.files(path, pattern="fastq.gz") # CHANGE if different file extensions

# Filtering
#might need to revisit these params
#need to add option to assign number of threads

filterAndTrim(file.path(path,fns), file.path(filtpath,fns), 
              truncLen=240, maxEE=1, truncQ=11, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=args$nprocs)


###########################
#                         #
# infer sequence variants #
#                         #
###########################

filts <- list.files(filtpath, pattern="fastq.gz", full.names=TRUE) # CHANGE if different file extensions

#probably need to modify this

sample.names <- sapply(strsplit(basename(filts), "_"), `[`, 1) # Assumes filename = sample_XXX.fastq.gz

names(filts) <- sample.names

# Learn error rates
set.seed(100)

#need to add option to assign number of threads
err <- learnErrors(filts, nreads = 1e6, multithread=args$nprocs, randomize=TRUE)
# Infer sequence variants
dds <- vector("list", length(sample.names))
names(dds) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err=err, multithread=args$nprocs)
}
# Construct sequence table 
st.all <- makeSequenceTable(dds)


###########################
#                         #
#   Chimera Filtering     #
#                         #
###########################


seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=args$nprocs)
# Assign taxonomy
tax <- assignTaxonomy(seqtab, args$tax_path, multithread=TRUE)

#do I want "species" level classification too?

# Write to disk 
#do I really want this to be just an R object or should I print it to tab? 

fn <- paste0(args$out_dir, "/", "sequence_table.rds")
saveRDS(seqtab, fn)

fn <- paste0(args$out_dir, "/", "sequence_table.txt")
write.table(seqtab, fn, quote=F, sep = "\t")  


fn <- paste0(args$out_dir, "/", "tax_final.rds")
saveRDS(tax, fn) 

fn <- paste0(args$out_dir, "/", "tax_final.txt")
write.table(tax, fn, quote=F, sep = "\t")

}
main()
