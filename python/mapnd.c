#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "../bseq.h"
#include "../minimap.h"
#include "../mmpriv.h"
#include "../ketopt.h"

#define MM_VERSION "2.20-r1061"

static ko_longopt_t long_options[] = {
	{ "bucket-bits",    ko_required_argument, 300 },
	{ "mb-size",        ko_required_argument, 'K' },
	{ "seed",           ko_required_argument, 302 },
	{ "no-kalloc",      ko_no_argument,       303 },
	{ "print-qname",    ko_no_argument,       304 },
	{ "no-self",        ko_no_argument,       'D' },
	{ "print-seeds",    ko_no_argument,       306 },
	{ "max-chain-skip", ko_required_argument, 307 },
	{ "min-dp-len",     ko_required_argument, 308 },
	{ "print-aln-seq",  ko_no_argument,       309 },
	{ "splice",         ko_no_argument,       310 },
	{ "cost-non-gt-ag", ko_required_argument, 'C' },
	{ "no-long-join",   ko_no_argument,       312 },
	{ "sr",             ko_no_argument,       313 },
	{ "frag",           ko_required_argument, 314 },
	{ "secondary",      ko_required_argument, 315 },
	{ "cs",             ko_optional_argument, 316 },
	{ "end-bonus",      ko_required_argument, 317 },
	{ "no-pairing",     ko_no_argument,       318 },
	{ "splice-flank",   ko_required_argument, 319 },
	{ "idx-no-seq",     ko_no_argument,       320 },
	{ "end-seed-pen",   ko_required_argument, 321 },
	{ "for-only",       ko_no_argument,       322 },
	{ "rev-only",       ko_no_argument,       323 },
	{ "heap-sort",      ko_required_argument, 324 },
	{ "all-chain",      ko_no_argument,       'P' },
	{ "dual",           ko_required_argument, 326 },
	{ "max-clip-ratio", ko_required_argument, 327 },
	{ "min-occ-floor",  ko_required_argument, 328 },
	{ "MD",             ko_no_argument,       329 },
	{ "lj-min-ratio",   ko_required_argument, 330 },
	{ "score-N",        ko_required_argument, 331 },
	{ "eqx",            ko_no_argument,       332 },
	{ "paf-no-hit",     ko_no_argument,       333 },
	{ "split-prefix",   ko_required_argument, 334 },
	{ "no-end-flt",     ko_no_argument,       335 },
	{ "hard-mask-level",ko_no_argument,       336 },
	{ "cap-sw-mem",     ko_required_argument, 337 },
	{ "max-qlen",       ko_required_argument, 338 },
	{ "max-chain-iter", ko_required_argument, 339 },
	{ "junc-bed",       ko_required_argument, 340 },
	{ "junc-bonus",     ko_required_argument, 341 },
	{ "sam-hit-only",   ko_no_argument,       342 },
	{ "chain-gap-scale",ko_required_argument, 343 },
	{ "alt",            ko_required_argument, 344 },
	{ "alt-drop",       ko_required_argument, 345 },
	{ "mask-len",       ko_required_argument, 346 },
	{ "rmq",            ko_optional_argument, 347 },
	{ "help",           ko_no_argument,       'h' },
	{ "max-intron-len", ko_required_argument, 'G' },
	{ "version",        ko_no_argument,       'V' },
	{ "min-count",      ko_required_argument, 'n' },
	{ "min-chain-score",ko_required_argument, 'm' },
	{ "mask-level",     ko_required_argument, 'M' },
	{ "min-dp-score",   ko_required_argument, 's' },
	{ "sam",            ko_no_argument,       'a' },
	{ 0, 0, 0 }
};

static inline int64_t mm_parse_num2(const char *str, char **q)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9, ++p;
	else if (*p == 'M' || *p == 'm') x *= 1e6, ++p;
	else if (*p == 'K' || *p == 'k') x *= 1e3, ++p;
	if (q) *q = p;
	return (int64_t)(x + .499);
}

static inline int64_t mm_parse_num(const char *str)
{
	return mm_parse_num2(str, 0);
}

static inline void yes_or_no(mm_mapopt_t *opt, int flag, int long_idx, const char *arg, int yes_to_set)
{
	if (yes_to_set) {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag |= flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag &= ~flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	} else {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag &= ~flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag |= flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	}
}

int parse_option(mm_mapopt_t *opt, mm_idxopt_t *ipt, char *_opt_str)
{
	_opt_str = strdup(_opt_str);
	int argc = 1;
	char *argv[1024] = {0};

	argv[argc++] = strtok(_opt_str," ");
	while (argv[argc - 1] != NULL)
	{
		// printf ("argv[%d]: %s\n",argc -1, argv[argc - 1]);
		argv[argc++] = strtok(NULL, " ");
	}
	argc --;
	
	const char *opt_str = "2aSDw:k:K:t:r:f:Vv:g:G:I:d:XT:s:x:Hcp:M:n:z:A:B:O:E:m:N:Qu:R:hF:LC:yYPo:e:U:";
	ketopt_t o = KETOPT_INIT;
	int c, old_best_n = -1;
	char *fnw = 0, *s;

	mm_verbose = 0;
	mm_set_opt(0, ipt, opt);

	ipt->batch_size = 0x7fffffffffffffffL;
	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) { // test command line options and apply option -x/preset first
		if (c == 'x') {
			if (mm_set_opt(o.arg, ipt, opt) < 0) {
				fprintf(stderr, "[ERROR] unknown preset '%s'\n", o.arg);
				return 1;
			}
		} else if (c == ':') {
			fprintf(stderr, "[ERROR] missing option argument\n");
			return 1;
		} else if (c == '?') {
			fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[o.i - 1]);
			return 1;
		}
	}
	o = KETOPT_INIT;

	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) {
		if (c == 'w') ipt->w = atoi(o.arg);
		else if (c == 'k') ipt->k = atoi(o.arg);
		else if (c == 'H') ipt->flag |= MM_I_HPC;
		// else if (c == 'd') fnw = o.arg; // the above are indexing related options, except -I
		else if (c == 't') opt->n_threads = atoi(o.arg);
		else if (c == 'v') mm_verbose = atoi(o.arg);
		else if (c == 'g') opt->max_gap = (int)mm_parse_num(o.arg);
		else if (c == 'G') mm_mapopt_max_intron_len(opt, (int)mm_parse_num(o.arg));
		else if (c == 'F') opt->max_frag_len = (int)mm_parse_num(o.arg);
		else if (c == 'N') old_best_n = opt->best_n, opt->best_n = atoi(o.arg);
		else if (c == 'p') opt->pri_ratio = atof(o.arg);
		else if (c == 'M') opt->mask_level = atof(o.arg);
		else if (c == 'c') opt->flag |= MM_F_OUT_CG | MM_F_CIGAR;
		else if (c == 'D') opt->flag |= MM_F_NO_DIAG;
		else if (c == 'P') opt->flag |= MM_F_ALL_CHAINS;
		else if (c == 'X') opt->flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN; // -D -P --no-long-join --dual=no
		else if (c == 'a') opt->flag |= MM_F_OUT_SAM | MM_F_CIGAR;
		else if (c == 'Q') opt->flag |= MM_F_NO_QUAL;
		else if (c == 'Y') opt->flag |= MM_F_SOFTCLIP;
		else if (c == 'L') opt->flag |= MM_F_LONG_CIGAR;
		else if (c == 'y') opt->flag |= MM_F_COPY_COMMENT;
		else if (c == 'T') opt->sdust_thres = atoi(o.arg);
		else if (c == 'n') opt->min_cnt = atoi(o.arg);
		else if (c == 'm') opt->min_chain_score = atoi(o.arg);
		else if (c == 'A') opt->a = atoi(o.arg);
		else if (c == 'B') opt->b = atoi(o.arg);
		else if (c == 's') opt->min_dp_max = atoi(o.arg);
		else if (c == 'C') opt->noncan = atoi(o.arg);
		// else if (c == 'I') ipt->batch_size = mm_parse_num(o.arg);
		else if (c == 'K') opt->mini_batch_size = mm_parse_num(o.arg);
		else if (c == 'e') opt->occ_dist = mm_parse_num(o.arg);
		// else if (c == 'R') rg = o.arg;
		// else if (c == 'h') fp_help = stdout;
		else if (c == '2') opt->flag |= MM_F_2_IO_THREADS;
		else if (c == 'o') {
			if (strcmp(o.arg, "-") != 0) {
				if (freopen(o.arg, "wb", stdout) == NULL) {
					fprintf(stderr, "[ERROR]\033[1;31m failed to write the output to file '%s'\033[0m: %s\n", o.arg, strerror(errno));
					exit(1);
				}
			}
		}
		else if (c == 300) ipt->bucket_bits = atoi(o.arg); // --bucket-bits
		else if (c == 302) opt->seed = atoi(o.arg); // --seed
		else if (c == 303) mm_dbg_flag |= MM_DBG_NO_KALLOC; // --no-kalloc
		else if (c == 304) mm_dbg_flag |= MM_DBG_PRINT_QNAME; // --print-qname
		else if (c == 306) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_SEED, opt->n_threads = 1; // --print-seed
		else if (c == 307) opt->max_chain_skip = atoi(o.arg); // --max-chain-skip
		else if (c == 339) opt->max_chain_iter = atoi(o.arg); // --max-chain-iter
		else if (c == 308) opt->min_ksw_len = atoi(o.arg); // --min-dp-len
		else if (c == 309) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_ALN_SEQ, opt->n_threads = 1; // --print-aln-seq
		else if (c == 310) opt->flag |= MM_F_SPLICE; // --splice
		else if (c == 312) opt->flag |= MM_F_NO_LJOIN; // --no-long-join
		else if (c == 313) opt->flag |= MM_F_SR; // --sr
		else if (c == 317) opt->end_bonus = atoi(o.arg); // --end-bonus
		else if (c == 318) opt->flag |= MM_F_INDEPEND_SEG; // --no-pairing
		else if (c == 320) ipt->flag |= MM_I_NO_SEQ; // --idx-no-seq
		else if (c == 321) opt->anchor_ext_shift = atoi(o.arg); // --end-seed-pen
		else if (c == 322) opt->flag |= MM_F_FOR_ONLY; // --for-only
		else if (c == 323) opt->flag |= MM_F_REV_ONLY; // --rev-only
		else if (c == 327) opt->max_clip_ratio = atof(o.arg); // --max-clip-ratio
		else if (c == 328) opt->min_mid_occ = atoi(o.arg); // --min-occ-floor
		else if (c == 329) opt->flag |= MM_F_OUT_MD; // --MD
		else if (c == 331) opt->sc_ambi = atoi(o.arg); // --score-N
		else if (c == 332) opt->flag |= MM_F_EQX; // --eqx
		else if (c == 333) opt->flag |= MM_F_PAF_NO_HIT; // --paf-no-hit
		// else if (c == 334) opt->split_prefix = o.arg; // --split-prefix
		else if (c == 335) opt->flag |= MM_F_NO_END_FLT; // --no-end-flt
		else if (c == 336) opt->flag |= MM_F_HARD_MLEVEL; // --hard-mask-level
		else if (c == 337) opt->max_sw_mat = mm_parse_num(o.arg); // --cap-sw-mat
		else if (c == 338) opt->max_qlen = mm_parse_num(o.arg); // --max-qlen
		else if (c == 340) ipt->junc_bed = strdup(o.arg); // --junc-bed
		else if (c == 341) opt->junc_bonus = atoi(o.arg); // --junc-bonus
		else if (c == 342) opt->flag |= MM_F_SAM_HIT_ONLY; // --sam-hit-only
		else if (c == 343) opt->chain_gap_scale = atof(o.arg); // --chain-gap-scale
		else if (c == 344) ipt->alt_list = strdup(o.arg); // --alt
		else if (c == 345) opt->alt_drop = atof(o.arg); // --alt-drop
		else if (c == 346) opt->mask_len = mm_parse_num(o.arg); // --mask-len
		else if (c == 330) {
			fprintf(stderr, "[WARNING] \033[1;31m --lj-min-ratio has been deprecated.\033[0m\n");
		} else if (c == 314) { // --frag
			yes_or_no(opt, MM_F_FRAG_MODE, o.longidx, o.arg, 1);
		} else if (c == 315) { // --secondary
			yes_or_no(opt, MM_F_NO_PRINT_2ND, o.longidx, o.arg, 0);
		} else if (c == 316) { // --cs
			opt->flag |= MM_F_OUT_CS | MM_F_CIGAR;
			if (o.arg == 0 || strcmp(o.arg, "short") == 0) {
				opt->flag &= ~MM_F_OUT_CS_LONG;
			} else if (strcmp(o.arg, "long") == 0) {
				opt->flag |= MM_F_OUT_CS_LONG;
			} else if (strcmp(o.arg, "none") == 0) {
				opt->flag &= ~MM_F_OUT_CS;
			} else if (mm_verbose >= 2) {
				fprintf(stderr, "[WARNING]\033[1;31m --cs only takes 'short' or 'long'. Invalid values are assumed to be 'short'.\033[0m\n");
			}
		} else if (c == 319) { // --splice-flank
			yes_or_no(opt, MM_F_SPLICE_FLANK, o.longidx, o.arg, 1);
		} else if (c == 324) { // --heap-sort
			yes_or_no(opt, MM_F_HEAP_SORT, o.longidx, o.arg, 1);
		} else if (c == 326) { // --dual
			yes_or_no(opt, MM_F_NO_DUAL, o.longidx, o.arg, 0);
		} else if (c == 347) { // --rmq
			yes_or_no(opt, MM_F_RMQ, o.longidx, o.arg, 1);
		} else if (c == 'S') {
			opt->flag |= MM_F_OUT_CS | MM_F_CIGAR | MM_F_OUT_CS_LONG;
			if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING]\033[1;31m option -S is deprecated and may be removed in future. Please use --cs=long instead.\033[0m\n");
		} else if (c == 'V') {
			puts(MM_VERSION);
			return 0;
		} else if (c == 'r') {
			opt->bw = (int)mm_parse_num2(o.arg, &s);
			if (*s == ',') opt->bw_long = (int)mm_parse_num2(s + 1, &s);
		} else if (c == 'U') {
			opt->min_mid_occ = strtol(o.arg, &s, 10);
			if (*s == ',') opt->max_mid_occ = strtol(s + 1, &s, 10);
		} else if (c == 'f') {
			double x;
			char *p;
			x = strtod(o.arg, &p);
			if (x < 1.0) opt->mid_occ_frac = x, opt->mid_occ = 0;
			else opt->mid_occ = (int)(x + .499);
			if (*p == ',') opt->max_occ = (int)(strtod(p+1, &p) + .499);
		} else if (c == 'u') {
			if (*o.arg == 'b') opt->flag |= MM_F_SPLICE_FOR|MM_F_SPLICE_REV; // both strands
			else if (*o.arg == 'f') opt->flag |= MM_F_SPLICE_FOR, opt->flag &= ~MM_F_SPLICE_REV; // match GT-AG
			else if (*o.arg == 'r') opt->flag |= MM_F_SPLICE_REV, opt->flag &= ~MM_F_SPLICE_FOR; // match CT-AC (reverse complement of GT-AG)
			else if (*o.arg == 'n') opt->flag &= ~(MM_F_SPLICE_FOR|MM_F_SPLICE_REV); // don't try to match the GT-AG signal
			else {
				fprintf(stderr, "[ERROR]\033[1;31m unrecognized cDNA direction\033[0m\n");
				return 1;
			}
		} else if (c == 'z') {
			opt->zdrop = opt->zdrop_inv = strtol(o.arg, &s, 10);
			if (*s == ',') opt->zdrop_inv = strtol(s + 1, &s, 10);
		} else if (c == 'O') {
			opt->q = opt->q2 = strtol(o.arg, &s, 10);
			if (*s == ',') opt->q2 = strtol(s + 1, &s, 10);
		} else if (c == 'E') {
			opt->e = opt->e2 = strtol(o.arg, &s, 10);
			if (*s == ',') opt->e2 = strtol(s + 1, &s, 10);
		}else if (c != 'x'){
			fprintf(stderr, "Unrecognized parameter: %c %s", c, o.arg);
			return 1;
		}
	}

	if ((opt->flag & MM_F_SPLICE) && (opt->flag & MM_F_FRAG_MODE)) {
		fprintf(stderr, "[ERROR]\033[1;31m --splice and --frag should not be specified at the same time.\033[0m\n");
		return 1;
	}
	if (!fnw && !(opt->flag&MM_F_CIGAR))
		ipt->flag |= MM_I_NO_SEQ;
	if (mm_check_opt(ipt, opt) < 0)
		return 1;
	if (opt->best_n == 0) {
		fprintf(stderr, "[WARNING]\033[1;31m changed '-N 0' to '-N %d --secondary=no'.\033[0m\n", old_best_n);
		opt->best_n = old_best_n, opt->flag |= MM_F_NO_PRINT_2ND;
	}

	if ((opt->flag & MM_F_SR) && argc - o.ind > 3) {
		fprintf(stderr, "[ERROR] incorrect input: in the sr mode, please specify no more than two query files.\n");
		return 1;
	}
	free(_opt_str);
	return 0;
}

mm_idx_t *build_index(char *file, char *options){
	mm_mapopt_t opt;
	mm_idxopt_t ipt;
	parse_option(&opt, &ipt, options);
	mm_idx_reader_t *idx_rdr = mm_idx_reader_open(file, &ipt, NULL);
	if (idx_rdr == 0) {
		fprintf(stderr, "[ERROR] failed to open file '%s': %s\n", file, strerror(errno));
		return NULL;
	}
	mm_idx_t *mi = mm_idx_reader_read(idx_rdr, opt.n_threads);
	mm_idx_reader_close(idx_rdr);
	return mi;
}

void destroy_index(mm_idx_t *mi){
	if (mi){
		mm_idx_destroy(mi);
		mi = NULL;
	}
}

typedef struct {
	char *hdr;
	char *aln;
} mm_out;

void destroy_out(mm_out *out){
	if (out){
		if (out->hdr) {
			free (out->hdr);
			out->hdr = NULL;
		}
		if (out->aln) {
			free (out->aln);
			out->aln = NULL;
		}
		free (out);
		out = NULL;
	}
}

mm_out *iter_map(mm_idx_t *mi, char *file, char *options){
	return NULL;
}

mm_out *multi_map(mm_idx_t *mi, char **file, int file_count, char *options){
	mm_mapopt_t opt;
	mm_idxopt_t ipt;
	parse_option(&opt, &ipt, options);
	opt.py_mode = 1;
	if (mi->k != ipt.k || mi->w != ipt.w || (mi->flag&MM_I_HPC) != (ipt.flag&MM_I_HPC)){
		fprintf(stderr, "Indexing parameters (-k, -w or -H) overridden by parameters used in the prebuilt index.\n");
		return NULL;
	}

	int i, ret;
	if ((opt.flag & MM_F_CIGAR) && (mi->flag & MM_I_NO_SEQ)) {
		fprintf(stderr, "[ERROR] the prebuilt index doesn't contain sequences.\n");
		mm_idx_destroy(mi);
		return NULL;
	}

	mm_out *out = (mm_out *) calloc(1, sizeof(mm_out));
	if (opt.flag & MM_F_OUT_SAM) {
		out->hdr = get_sam_hdr(mi, NULL, MM_VERSION, options);
		if (!out->hdr) {
			mm_idx_destroy(mi);
			destroy_out(out);
			return NULL;
		}
	}

	mm_mapopt_update(&opt, mi);
	if (ipt.junc_bed) {
		mm_idx_bed_read(mi, ipt.junc_bed, 1);
		free (ipt.junc_bed);
	}
	if (ipt.alt_list) {
		mm_idx_alt_read(mi, ipt.alt_list);
		free (ipt.alt_list);
	}

	ret = 0;
	if (!(opt.flag & MM_F_FRAG_MODE)) {
		for (i = 0; i < file_count; ++i) {
			ret = mm_map_file(mi, file[i], &opt, opt.n_threads);
			if (ret < 0) break;
		}
	} else {
		ret = mm_map_file_frag(mi, file_count, (const char**)file, &opt, opt.n_threads);
	}
	if (ret < 0) {
		fprintf(stderr, "ERROR: failed to map the query file\n");
		destroy_out(out);
		return NULL;
	}
	out->aln = opt.out;
	return out;
}

mm_out *map(mm_idx_t *mi, char *file, char *options){
	char *files[] = {file};
	return multi_map(mi, files, 1, options);
}

int32_t get_mid_occ(mm_idx_t *mi, char *options) {
	mm_mapopt_t opt;
	mm_idxopt_t ipt;
	parse_option(&opt, &ipt, options);
	return mm_idx_cal_max_occ(mi, opt.mid_occ_frac);
}

int main(int argc, char *argv[]){
	char *options = "--sam-hit-only -Y -a -x map-hifi -t 10";
	char *ref = "/data/project/huj/14_grandsv/01_data/human_g1k_v37.chr1.fasta";
	char *read = "/data/project/huj/14_grandsv/00_grandsv/chr1.hifi.bam.fasta.dup.fa";
	mm_idx_t *idx = build_index(ref, options);
	mm_out *out = map(idx, read, options);

	if (out->hdr) printf("%s\n", out->hdr);
	if (out->aln) printf("%s", out->aln);
	destroy_index(idx);
	destroy_out(out);
	return 0;
}
