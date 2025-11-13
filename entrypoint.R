#!/usr/bin/Rscript
library(optparse)
library(tictoc)
library(processx)
library(httr2)
library(jsonlite)
library(curl)

printf = function(str, ...) {
  cat(sprintf(str, ...))
}

`%_%` = paste0

argv = commandArgs(trailingOnly=TRUE)
if(length(argv) == 0) {
	printf(
		"Usage: \n" %_%
		"bam2bed merge|bed|mnpflex [--only | --force] --name samplename [options]\n" %_%
		"By default, specifying a later step in the merge→bed→mnpflex pipeline\n" %_%
		"will also run the earlier steps, unless their output files already exist.\n" %_%
		"    --only     : Only try to run the specified step.\n" %_%
		"    --force    : Redo the specified step(s) even if the output seems to exist.\n" %_%
		"    --name     : Sample name to use in the output\n" %_%
		"    --threads  : Number of threads to use in samtools and modkit [default 8]\n" %_%
		"    --tag      : Tag to mark output files with, default is blank\n" %_%
		"\nThe different stages also take options:\n" %_%
		"merge:\n" %_%
		"    --Nbam     : Number of bamfiles to use [default 8000]\n" %_%
		"BED:\n" %_%
		"    --mbases   : Which base codes to retain in the output (e.g. 'hm' for 5mc and 5hmC) [default 'm']\n" %_%
		"mnpflex:\n" %_%
		"    --username : MNP-Flex username (required)\n" %_%
		"    --password : MNP-Flex password (required)\n"
		)
	q(save="no", status=22)
}

tic("Everything")

verb = argv[1]
option_list <- list(
	make_option("--only", action="store_true", default=FALSE, 
	  help="Only try to run the given step (useful if you have e.g. a merged bam but not original data)"
	),
	make_option("--force", action="store_true", default=FALSE, 
	  help="Rerun the given step(s) even if the output appears to exist"
	),
	make_option("--name", type="character", help="Sample name for use in output. No spaces."),
	make_option("--tag", type="character", default="", help="Tag to mark output files with"),
	make_option("--threads", type="integer", default=8, help="Threads to use in modkit and samtool"),
	make_option("--Nbam", type="integer", default=8000, help="Number of BAM files to merge"),
	make_option("--username", type="character", default="", help="MNP-Flex username"),
	make_option("--password", type="character", default="", help="MNP-Flex password")
)

parser = OptionParser(option_list=option_list, usage="%prog [--only | --force], --name samplename [options]")
args = parse_args(parser, argv[-1])

if(!"name" %in% names(args)) {
	printf("Missing --name argument\n")
	q(save="no", status=22)
}

do_merge   = verb %in% c("merge", "bed", "mnpflex") & !(verb != "merge" & args$only)
do_bed  = verb %in% c("bed", "mnpflex") & !(verb != "bed" & args$only)
do_mnpflex = (verb == "mnpflex")

printf("SampleID: \t\x1b[1m%s\x1b[0m\n", args$name)
printf("Merge:   %s\t%s\n", do_merge, ifelse(do_merge & args$force, "[Overwrite]", ""))
printf("BED   :  %s\t%s\n", do_bed,ifelse(do_bed & args$force, "[Overwrite]", ""))
printf("MNPFlex: %s\t%s\n", do_mnpflex, ifelse(do_mnpflex & args$force, "[Overwrite]", ""))

merge_dir   = file.path("/work/out/", args$name, "merged_bam")
bed_dir     = file.path("/work/out/", args$name, "BED")
mnpflex_dir = file.path("/work/out/", args$name, "mnpflex")

if(do_merge)   dir.create(merge_dir,   recursive=TRUE, showWarnings=FALSE)
if(do_bed)     dir.create(bed_dir,  recursive=TRUE, showWarnings=FALSE)
if(do_mnpflex) dir.create(mnpflex_dir, recursive=TRUE, showWarnings=FALSE)

if (args$tag != "") {
	bamfile = file.path(merge_dir,   sprintf("%s_%s_merged.bam", args$name, args$tag))
	bedfile = file.path(bed_dir,     sprintf("%s_%s_pileup.bed", args$name, args$tag))
	mnpfile = file.path(bed_dir,     sprintf("%s_%s_mnpflex.bed", args$name, args$tag))
	report  = file.path(mnpflex_dir, sprintf("%s_%s_report.pdf", args$name, args$tag))
} else {
	bamfile = file.path(merge_dir,   sprintf("%s_merged.bam", args$name))
	bedfile = file.path(bed_dir,     sprintf("%s_pileup.bed", args$name))
	mnpfile = file.path(bed_dir,     sprintf("%s_mnpflex.bed", args$name))
	report  = file.path(mnpflex_dir, sprintf("%s_report.pdf", args$name))
}


merge_bam = function(Nbam, outfile, force, tag) {
	if(!force & file.exists(outfile)) {
		printf("Merge: Skipping, %s already exists and --force not specified\n", outfile)
		return(invisible())
	}
	printf("Merging max %d BAM files to \x1b[1m%s\x1b[0m\n", Nbam, outfile)
	printf("Scanning for BAM files in in/*/bam_pass/ and in/bam_pass/ :\n")
	bamfiles = Sys.glob("/work/in/bam_pass/*.[bB][aA][mM]")
	bamfiles = c(bamfiles, Sys.glob("/work/in/*/bam_pass/*.[bB][aA][mM]"))
	printf("Found %d BAM files\n", length(bamfiles))
	bamfiles = bamfiles[1:min(Nbam, length(bamfiles))]
	
	args = c("merge", "-l", "4", sprintf("-@%d", args$threads), 
	         "--write-index", "--verbosity", "255", "-o", outfile)
	if(force) args=c(args, "-f")
	args = c(args, bamfiles)
	printf("Starting samtools merge on %d files. This could take a while.\n", length(bamfiles))
	tic("samtools")
	discard = run("/usr/bin/samtools", args=args, echo=TRUE)
	toc()
	printf("Merge: Done.\n")
}

make_bed = function(bamfile, bedfile, mnpfile, force) {
	if(file.exists(bedfile) & !force) {
		printf("BED: Skipping, %s already exists and --force not specified\n", bedfile)
		return(invisible())
	} 
	printf("Running modkit pileup on \x1b[1m%s\x1b[0m\n to create \x1b[1m%s\x1b[0m\n\n", bamfile, bedfile)
	printf("Starting modkit pileup. This may take a while.\n")
	tic("Modkit")
	run("modkit", args=c("pileup", "--cpg", "--ignore", "h", "-t", args$threads, bamfile, bedfile))
	toc()
	printf("BED: Modkit done.\n")
	
	printf("Formatting \x1b[1m%s\x1b[0m\n to \x1b[1m%s\x1b[0m\n\n", bedfile, mnpfile)
	tic("make_mnpflex.R")
	run("Rscript", c("make_mnpflex.R", "-i", bedfile, "-o", mnpfile, 
	    "-r", "MNP-flex-methylationepic-v-1-0-b5-manifest-file-MGMT-complete.sorted.bed")
	)
	toc()
	printf("BED: MNP-flex prefiltering/formatting done. \nDone.\n")
}

run_mnpflex = function(mnpfile, reportfile, name, tag, username, password, force) {
	if(username=="" | password=="") {
		printf("Both --username and --password must be set to submit to MNP Flex\n")
		return(invisible())
	}
	
	if(! file.exists(mnpfile)) {
		printf("MNP file missing - did you run with --only without having run 'bed' for this name + tag?\n")
		return(invisible())
	}
	
	if(file.exists(reportfile) & !force) {
		printf("MNP Flex: Skipping, %s already exists and --force not specified\n", reportfile)
		return(invisible())
	}
	
	upname = sprintf("%s_%s_%s", name, tag, format(Sys.time(), "%Y-%m-%dT%H%M"))
	printf("Submitting to MNP Flex as \x1b[1m%s\x1b[0m\n", upname)
	
	req = request("https://mnp-flex.org") |>
	      req_options(ssl_verifypeer = 0)

	# Log in and get a bearer token
	token = req |> 
	    req_url_path("/api/v1/auth/token") |>
	    req_body_form(username=username, password=password) |>
	    req_perform() |>
	    resp_body_json() |>
	    getElement("access_token")
	
	# Upload our samples
	job = req |> 
	    req_url_path("/api/v1/samples") |>
	    req_auth_bearer_token(token) |>
	    req_url_query(sample_name = upname, disclaimer_confirmed='true') |>
	    req_body_multipart(file = form_file(path=mnpfile)) |>
	    req_method("PUT") |>
	    req_perform() |>
	    resp_body_json()
	
	printf("Waiting for MNP Flex to process the sample\n")
	done = FALSE
	i=0
	while(!done) {
		Sys.sleep(1)
		i = i+1
		res = req |>
		    req_url_path("/api/v1/samples/" %_% job$id) |>
		    req_auth_bearer_token(token) |>
		    req_perform() |>
		    resp_body_json()

		analysis_status = res$bed_file_sample$analysis_status
		printf("Waiting: %d sec, status is %s\n", i, analysis_status)
		if(analysis_status == "Analysis error") {
		    printf("Analysis error!\n")
		    return(invisible())
		}
		done <- analysis_status == "done"
	}
	
	# Download the report
	report = req |>
	    req_url_path("/api/v1/samples/download_sample_result/" %_% job$id) |>
	    req_auth_bearer_token(token) |>
	    req_perform() |>
	    resp_body_raw()
	writeBin(report, reportfile)

	printf("MNP Flex done; report saved as \x1b[1m%s\x1b[0m.\n", reportfile)

}

if(do_merge)   merge_bam(args$Nbam, bamfile, args$force, args$tag)
if(do_bed)     make_bed(bamfile, bedfile, mnpfile, args$force)
if(do_mnpflex) run_mnpflex(mnpfile, report, args$name, args$tag, args$username, args$password, args$force)

toc()
