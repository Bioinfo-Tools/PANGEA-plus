#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>


#define NCBITC_FILE_DMP_GI_TAXID	"gi_taxid_nucl.dmp"
#define	NCBITC_FILE_DMP_NODES		"nodes.dmp"
#define	NCBITC_FILE_DMP_NAMES		"names.dmp"

#define NCBITC_FILE_BIN_GI_TAXID	"gi_taxid_nucl.dmp.bin"
#define	NCBITC_FILE_BIN_NODES		"nodes.dmp.bin"
#define	NCBITC_FILE_BIN_NAMES		"names.dmp.bin"

#undef NCBITC_WITH_COMMENTS
#undef NCBITC_WITH_GIINDEX

enum ncbitc_name_class {
	NCBITC_CLASS_ACRONYM,
	NCBITC_CLASS_ANAMORPH,
	NCBITC_CLASS_AUTHORITY,
	NCBITC_CLASS_BLAST_NAME,
	NCBITC_CLASS_COMMON_NAME,
	NCBITC_CLASS_EQUIVALENT_NAME,
	NCBITC_CLASS_GENBANK_ACRONYM,
	NCBITC_CLASS_GENBANK_ANAMORPH,
	NCBITC_CLASS_GENBANK_COMMON_NAME,
	NCBITC_CLASS_GENBANK_SYNONYM,
	NCBITC_CLASS_INCLUDES,
	NCBITC_CLASS_IN_PART,
	NCBITC_CLASS_MISNOMER,
	NCBITC_CLASS_MISSPELLING,
	NCBITC_CLASS_SCIENTIFIC_NAME,
	NCBITC_CLASS_SYNONYM,
	NCBITC_CLASS_TELEOMORPH,
	NCBITC_CLASS_UNPUBLISHED_NAME,
};

enum ncbitc_rank {
	NCBITC_RANK_CLASS,
	NCBITC_RANK_FAMILY,
	NCBITC_RANK_FORMA,
	NCBITC_RANK_GENUS,
	NCBITC_RANK_INFRACLASS,
	NCBITC_RANK_INFRAORDER,
	NCBITC_RANK_KINGDOM,
	NCBITC_RANK_NO_RANK,
	NCBITC_RANK_ORDER,
	NCBITC_RANK_PARVORDER,
	NCBITC_RANK_PHYLUM,
	NCBITC_RANK_SPECIES,
	NCBITC_RANK_SPECIES_GROUP,
	NCBITC_RANK_SPECIES_SUBGROUP,
	NCBITC_RANK_SUBCLASS,
	NCBITC_RANK_SUBFAMILY,
	NCBITC_RANK_SUBGENUS,
	NCBITC_RANK_SUBKINGDOM,
	NCBITC_RANK_SUBORDER,
	NCBITC_RANK_SUBPHYLUM,
	NCBITC_RANK_SUBSPECIES,
	NCBITC_RANK_SUBTRIBE,
	NCBITC_RANK_SUPERCLASS,
	NCBITC_RANK_SUPERFAMILY,
	NCBITC_RANK_SUPERKINGDOM,
	NCBITC_RANK_SUPERORDER,
	NCBITC_RANK_SUPERPHYLUM,
	NCBITC_RANK_TRIBE,
	NCBITC_RANK_VARIETAS
};

#define NCBITC_MAX_SHORT_NAME	32
#define NCBITC_MAX_LONG_NAME	64
#define NCBITC_MAX_HUGE_NAME	256
#define NCBITC_LINE_SIZE		512

/**
nodes.dmp
---------

This file represents taxonomy nodes. The description for each node includes 
the following fields:

	tax_id					-- node id in GenBank taxonomy database
 	parent tax_id				-- parent node id in GenBank taxonomy database
 	rank					-- rank of this node (superkingdom, kingdom, ...) 
 	embl code				-- locus-name prefix; not unique
 	division id				-- see division.dmp file
 	inherited div flag  (1 or 0)		-- 1 if node inherits division from parent
 	genetic code id				-- see gencode.dmp file
 	inherited GC  flag  (1 or 0)		-- 1 if node inherits genetic code from parent
 	mitochondrial genetic code id		-- see gencode.dmp file
 	inherited MGC flag  (1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
 	GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
 	hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
 	comments				-- free-text comments and citations
 */

struct nodes_dmp {
	int tax_id;
	int parent_tax_id;
	char rank;
	char embl_code[3];
	short division_id;
	char inherited_div_flag;
	short genetic_code_id;
	char inherited_GC_flag;
	int mitochondrial_genetic_code_id;
	char inherited_MGC_flag;
	char GenBank_hidden_flag;
	char hidden_subtree_root_flag;
#ifdef NCBITC_WITH_COMMENTS
	char comments[255];
#endif
};

/**

names.dmp
---------
Taxonomy names file has these fields:

	tax_id					-- the id of node associated with this name
	name_txt				-- name itself
	unique name				-- the unique variant of this name if name not unique
	name class				-- (synonym, common name, ...)
*/

struct names_dmp {
	int tax_id;
	char name_txt[NCBITC_MAX_LONG_NAME];
	char unique_name[NCBITC_MAX_LONG_NAME];
	char name_class[NCBITC_MAX_LONG_NAME];
};

struct gi_taxid_nucl_dmp {
#ifdef NCBITC_WITH_GIINDEX
	int gi;
#endif
	int tax_id;
};

/* Flag set by ‘--verbose’. */
static int verbose_flag = 0;

int ncbitc_bin_to_str_class(int id, char *str)
{
	switch (id) {
	case NCBITC_CLASS_ACRONYM:
		sprintf(str, "acronym");
		break;

	case NCBITC_CLASS_ANAMORPH:
		sprintf(str, "anamorph");
		break;

	case NCBITC_CLASS_AUTHORITY:
		sprintf(str, "authority");
		break;

	case NCBITC_CLASS_BLAST_NAME:
		sprintf(str, "blast name");
		break;

	case NCBITC_CLASS_COMMON_NAME:
		sprintf(str, "common name");
		break;

	case NCBITC_CLASS_EQUIVALENT_NAME:
		sprintf(str, "equivalent name");
		break;

	case NCBITC_CLASS_GENBANK_ACRONYM:
		sprintf(str, "genbank acronym");
		break;

	case NCBITC_CLASS_GENBANK_ANAMORPH:
		sprintf(str, "genbank anamorph");
		break;

	case NCBITC_CLASS_GENBANK_COMMON_NAME:
		sprintf(str, "genbank common name");
		break;

	case NCBITC_CLASS_GENBANK_SYNONYM:
		sprintf(str, "genbank synonym");
		break;

	case NCBITC_CLASS_INCLUDES:
		sprintf(str, "includes");
		break;

	case NCBITC_CLASS_IN_PART:
		sprintf(str, "in-part");
		break;

	case NCBITC_CLASS_MISNOMER:
		sprintf(str, "misnomer");
		break;

	case NCBITC_CLASS_MISSPELLING:
		sprintf(str, "misspelling");
		break;

	case NCBITC_CLASS_SCIENTIFIC_NAME:
		sprintf(str, "scientific name");
		break;

	case NCBITC_CLASS_SYNONYM:
		sprintf(str, "synonym");
		break;

	case NCBITC_CLASS_TELEOMORPH:
		sprintf(str, "teleomorph");
		break;

	case NCBITC_CLASS_UNPUBLISHED_NAME:
		sprintf(str, "unpublished name");
		break;

	default:
		sprintf(str, "invalid id");
		return -1;
	}

	return 0;
};

int ncbitc_str_to_bin_class(char *str, int *id)
{
	if (strncmp(str, "acronym", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_ACRONYM;
	else if (strncmp(str, "anamorph", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_ANAMORPH;
	else if (strncmp(str, "authority", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_AUTHORITY;
	else if (strncmp(str, "blast name", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_BLAST_NAME;
	else if (strncmp(str, "common name", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_COMMON_NAME;
	else if (strncmp(str, "equivalent name", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_EQUIVALENT_NAME;
	else if (strncmp(str, "genbank acronym", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_GENBANK_ACRONYM;
	else if (strncmp(str, "genbank anamorph", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_GENBANK_ANAMORPH;
	else if (strncmp(str, "genbank common name", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_GENBANK_COMMON_NAME;
	else if (strncmp(str, "genbank synonym", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_GENBANK_SYNONYM;
	else if (strncmp(str, "includes", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_INCLUDES;
	else if (strncmp(str, "in-part", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_IN_PART;
	else if (strncmp(str, "misnomer", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_MISNOMER;
	else if (strncmp(str, "misspelling", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_MISSPELLING;
	else if (strncmp(str, "scientific name", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_SCIENTIFIC_NAME;
	else if (strncmp(str, "synonym", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_SYNONYM;
	else if (strncmp(str, "teleomorph", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_TELEOMORPH;
	else if (strncmp(str, "unpublished name", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_CLASS_UNPUBLISHED_NAME;
	else {
		*id = -1;
	}

	return *id;
}

int ncbitc_bin_to_str_rank(int id, char *str)
{
	switch (id) {
	case NCBITC_RANK_CLASS:
		sprintf(str, "class");
		break;

	case NCBITC_RANK_FAMILY:
		sprintf(str, "family");
		break;

	case NCBITC_RANK_FORMA:
		sprintf(str, "forma");
		break;

	case NCBITC_RANK_GENUS:
		sprintf(str, "genus");
		break;

	case NCBITC_RANK_INFRACLASS:
		sprintf(str, "infraclass");
		break;

	case NCBITC_RANK_INFRAORDER:
		sprintf(str, "infraorder");
		break;

	case NCBITC_RANK_KINGDOM:
		sprintf(str, "kingdom");
		break;

	case NCBITC_RANK_NO_RANK:
		sprintf(str, "no rank");
		break;

	case NCBITC_RANK_ORDER:
		sprintf(str, "order");
		break;

	case NCBITC_RANK_PARVORDER:
		sprintf(str, "parvorder");
		break;

	case NCBITC_RANK_PHYLUM:
		sprintf(str, "phylum");
		break;

	case NCBITC_RANK_SPECIES:
		sprintf(str, "species");
		break;

	case NCBITC_RANK_SPECIES_GROUP:
		sprintf(str, "species group");
		break;

	case NCBITC_RANK_SPECIES_SUBGROUP:
		sprintf(str, "species subgroup");
		break;

	case NCBITC_RANK_SUBCLASS:
		sprintf(str, "subclass");
		break;

	case NCBITC_RANK_SUBFAMILY:
		sprintf(str, "subfamily");
		break;

	case NCBITC_RANK_SUBGENUS:
		sprintf(str, "subgenus");
		break;

	case NCBITC_RANK_SUBKINGDOM:
		sprintf(str, "subkingdom");
		break;

	case NCBITC_RANK_SUBORDER:
		sprintf(str, "suborder");
		break;

	case NCBITC_RANK_SUBPHYLUM:
		sprintf(str, "subphylum");
		break;

	case NCBITC_RANK_SUBSPECIES:
		sprintf(str, "subspecies");
		break;

	case NCBITC_RANK_SUBTRIBE:
		sprintf(str, "subtribe");
		break;

	case NCBITC_RANK_SUPERCLASS:
		sprintf(str, "superclass");
		break;

	case NCBITC_RANK_SUPERFAMILY:
		sprintf(str, "superfamily");
		break;

	case NCBITC_RANK_SUPERKINGDOM:
		sprintf(str, "superkingdom");
		break;

	case NCBITC_RANK_SUPERORDER:
		sprintf(str, "superorder");
		break;

	case NCBITC_RANK_SUPERPHYLUM:
		sprintf(str, "superphylum");
		break;

	case NCBITC_RANK_TRIBE:
		sprintf(str, "tribe");
		break;

	case NCBITC_RANK_VARIETAS:
		sprintf(str, "varietas");
		break;

	default:
		sprintf(str, "invalid id");
		return -1;
	}

	return 0;
}

int ncbitc_str_to_bin_rank(char *str, int *id)
{
	if (strncmp(str, "class", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_CLASS;
	else if (strncmp(str, "family", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_FAMILY;
	else if (strncmp(str, "forma", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_FORMA;
	else if (strncmp(str, "genus", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_GENUS;
	else if (strncmp(str, "infraclass", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_INFRACLASS;
	else if (strncmp(str, "infraorder", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_INFRAORDER;
	else if (strncmp(str, "kingdom", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_KINGDOM;
	else if (strncmp(str, "no rank", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_NO_RANK;
	else if (strncmp(str, "order", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_ORDER;
	else if (strncmp(str, "parvorder", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_PARVORDER;
	else if (strncmp(str, "phylum", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_PHYLUM;
	else if (strncmp(str, "species", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_SPECIES;
	else if (strncmp(str, "species group", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_SPECIES_GROUP;
	else if (strncmp(str, "species subgroup", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_SPECIES_SUBGROUP;
	else if (strncmp(str, "subclass", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_SUBCLASS;
	else if (strncmp(str, "subfamily", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_SUBFAMILY;
	else if (strncmp(str, "subgenus", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_SUBGENUS;
	else if (strncmp(str, "subkingdom", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_SUBKINGDOM;
	else if (strncmp(str, "suborder", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_SUBORDER;
	else if (strncmp(str, "subphylum", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_SUBPHYLUM;
	else if (strncmp(str, "subspecies", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_SUBSPECIES;
	else if (strncmp(str, "subtribe", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_SUBTRIBE;
	else if (strncmp(str, "superclass", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_SUPERCLASS;
	else if (strncmp(str, "superfamily", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_SUPERFAMILY;
	else if (strncmp(str, "superkingdom", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_SUPERKINGDOM;
	else if (strncmp(str, "superorder", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_SUPERORDER;
	else if (strncmp(str, "superphylum", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_SUPERPHYLUM;
	else if (strncmp(str, "tribe", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_TRIBE;
	else if (strncmp(str, "varietas", NCBITC_MAX_SHORT_NAME) == 0)
		*id = NCBITC_RANK_VARIETAS;
	else {
		*id = -1;
	}

	return *id;
}

int ncbitc_print_node_entry(struct nodes_dmp *node)
{
	char rank[NCBITC_MAX_SHORT_NAME];
/*	
	if (node->tax_id == 0) {
		printf("0\n");
		return 0;
	}
*/
	ncbitc_bin_to_str_rank(node->rank, rank);

	printf
	    ("%d | %d | %s | %s | %d | %d | %d | %d | %d | %d | %d | %d | %s |\n",
	     node->tax_id, node->parent_tax_id, rank, node->embl_code,
	     node->division_id, node->inherited_div_flag, node->genetic_code_id,
	     node->inherited_GC_flag, node->mitochondrial_genetic_code_id,
	     node->inherited_MGC_flag, node->GenBank_hidden_flag,
	     node->hidden_subtree_root_flag,
#ifdef NCBITC_WITH_COMMENTS
	     node->comments
#else
	     ""
#endif
	    );

	return 0;
}

int ncbitc_cleanup_str(char *str, int max)
{
	int len;

	len = strnlen(str, max);

	if (len <= 2) {
		str[0] = '\0';
	} else {
		str[0] = ' ';
		str[len - 1] = '\0';
	}

	return 0;
}

int ncbitc_read_node_entry(char *line, struct nodes_dmp *node)
{
	char rank[NCBITC_MAX_SHORT_NAME];
	char embl_code[NCBITC_MAX_SHORT_NAME];
	char comments[255];

	sscanf(line,
	       "%d | %d |%[^|]|%[^|]| %d | %d | %d | %d | %d | %d | %d | %d |%[^|]|",
	       &node->tax_id, &node->parent_tax_id, rank, embl_code,
	       (int *)&node->division_id, (int *)&node->inherited_div_flag,
	       (int *)&node->genetic_code_id, (int *)&node->inherited_GC_flag,
	       (int *)&node->mitochondrial_genetic_code_id,
	       (int *)&node->inherited_MGC_flag,
	       (int *)&node->GenBank_hidden_flag,
	       (int *)&node->hidden_subtree_root_flag, comments);

	ncbitc_cleanup_str(rank, 32);
	ncbitc_str_to_bin_rank(rank + 1, (int *)&node->rank);

	ncbitc_cleanup_str(embl_code, 32);
	if (embl_code[0] != '\0') {
		node->embl_code[0] = embl_code[1];
		node->embl_code[1] = embl_code[2];
		node->embl_code[2] = '\0';
	} else {
		node->embl_code[0] = '\0';
	}

	ncbitc_cleanup_str(comments, 255);
#ifdef NCBITC_WITH_COMMENTS
	strncpy(node->comments, comments, 255);
#endif

	return 0;
};

int ncbitc_read_name_entry(char *line, struct names_dmp *name)
{
	sscanf(line, "%d |%[^|]|%[^|]|%[^|]|", &name->tax_id,
	       name->name_txt, name->unique_name, name->name_class);

	ncbitc_cleanup_str(name->name_txt, NCBITC_MAX_LONG_NAME);
	ncbitc_cleanup_str(name->unique_name, NCBITC_MAX_LONG_NAME);
	ncbitc_cleanup_str(name->name_class, NCBITC_MAX_SHORT_NAME);

	return 0;
}

int ncbitc_print_name_entry(struct names_dmp *name)
{
	printf("%d | %s | %s | %s |\n", name->tax_id,
	       name->name_txt, name->unique_name, name->name_class);

	return 0;
}

int ncbitc_search_tax_id(char *filename, int gi, int *tax_id)
{
	FILE *fp;
	long offset;
	struct gi_taxid_nucl_dmp entry;

	fp = fopen(filename, "rb");
	if (fp == NULL) {
		perror("fopen");
		return -1;
	}

	/* seek to line */
	offset = (gi - 1) * sizeof(struct gi_taxid_nucl_dmp);
	if (fseek(fp, offset, SEEK_SET) < 0) {
		perror("fseek");
		fclose(fp);
		return -1;
	}

	fread(&entry, sizeof(struct gi_taxid_nucl_dmp), 1, fp);
	if (verbose_flag) {
#ifdef NCBITC_WITH_GIINDEX
		printf("%d\t%d\n", entry.gi, entry.tax_id);
#else
		printf("%d\t%d\n", gi, entry.tax_id);
#endif
	}
	*tax_id = entry.tax_id;

	fclose(fp);
	return 0;
}

int ncbitc_search_node(char *filename, int tax_id, struct nodes_dmp *node)
{
	FILE *fp;
	long offset;

	fp = fopen(filename, "rb");
	if (fp == NULL) {
		perror("fopen");
		return -1;
	}

	/* seek to line */
	offset = (tax_id - 1) * sizeof(struct nodes_dmp);
	if (fseek(fp, offset, SEEK_SET) < 0) {
		perror("fseek");
		fclose(fp);
		return -1;
	}

	fread(node, sizeof(struct nodes_dmp), 1, fp);
	if (verbose_flag) {
		printf("%d\n", node->tax_id);
	}

	fclose(fp);

	return 0;
}

int ncbitc_seek_name(FILE * fp, int pos, struct names_dmp *name)
{
	long offset;

	/* seek to line */
	offset = (pos - 1) * sizeof(struct names_dmp) + sizeof(int);
	if (fseek(fp, offset, SEEK_SET) < 0) {
		perror("fseek");
		fclose(fp);
		return -1;
	}

	fread(name, sizeof(struct names_dmp), 1, fp);

	return name->tax_id;
}

int ncbitc_search_name(char *filename, int tax_id)
{
	FILE *fp;
	int val = 0;
	int l1, i, j, flag = 0;
	l1 = 0;
	int num;
	struct names_dmp name;

	fp = fopen(filename, "rb");
	if (fp == NULL) {
		perror("fopen");
		return -1;
	}

	fread(&num, sizeof(int), 1, fp);

	/* binary search */
	i = num - 1;
	while (l1 <= i) {
		j = (l1 + i) / 2;
		val = ncbitc_seek_name(fp, j, &name);
		//ncbitc_print_name_entry(&name);
		if (val == tax_id) {
			flag = 1;
			break;
		} else if (val < tax_id)
			l1 = j + 1;
		else
			i = j - 1;
	}
	if (flag == 0) {
		printf("0\n");
		return 0;
	}

	for (i = j - 1; i > 0; i--) {
		val = ncbitc_seek_name(fp, i, &name);
		if (val != tax_id)
			break;
	}

	for (i = i + 1; i < num; i++) {
		val = ncbitc_seek_name(fp, i, &name);
		if (val != tax_id)
			break;
		ncbitc_print_name_entry(&name);
	}

	fclose(fp);

	return 0;
}

int ncbitc_create_gi_taxid(char *filename)
{
	FILE *fp;
	FILE *fpnew;
	char line[NCBITC_LINE_SIZE];
	int i, diff, last;
	char new[128];
	struct gi_taxid_nucl_dmp zero;
	struct gi_taxid_nucl_dmp entry;
	int gi;

	memset(&zero, '\0', sizeof(struct gi_taxid_nucl_dmp));

	sprintf(new, "%s.bin", filename);

	fp = fopen(filename, "r");
	if (fp == NULL) {
		perror("fopen");
		return -1;
	}

	fpnew = fopen(new, "wb");
	if (fp == NULL) {
		perror("fopen");
		return -1;
	}

	last = 0;
	while (fgets(line, NCBITC_LINE_SIZE, fp) != NULL) {
		sscanf(line, "%d\t%d", &gi, &entry.tax_id);
		diff = gi - last;
		if (diff != 1) {
			for (i = 1; i < diff; i++)
				fwrite(&zero, sizeof(struct gi_taxid_nucl_dmp),
				       1, fpnew);
		}
#ifdef NCBITC_WITH_GIINDEX
		entry.gi = gi;
#endif
		fwrite(&entry, sizeof(struct gi_taxid_nucl_dmp), 1, fpnew);
		last = gi;
	}

	fclose(fp);
	fclose(fpnew);

	return 0;
}

int ncbitc_create_nodes(char *filename)
{
	FILE *fp;
	FILE *fpnew;
	char line[NCBITC_LINE_SIZE];
	int i, diff, last;
	char new[128];
	struct nodes_dmp node;
	struct nodes_dmp zero;

	memset(&zero, 0, sizeof(struct nodes_dmp));

	sprintf(new, "%s.bin", filename);

	fp = fopen(filename, "r");
	if (fp == NULL) {
		perror("fopen");
		return -1;
	}

	fpnew = fopen(new, "wb");
	if (fp == NULL) {
		perror("fopen");
		return -1;
	}

	last = 0;
	while (fgets(line, NCBITC_LINE_SIZE, fp) != NULL) {
		ncbitc_read_node_entry(line, &node);
		//ncbitc_print_node_entry(&node);
		diff = node.tax_id - last;
		if (diff != 1) {
			for (i = 1; i < diff; i++)
				fwrite(&zero, sizeof(struct nodes_dmp), 1,
				       fpnew);
		}
		fwrite(&node, sizeof(struct nodes_dmp), 1, fpnew);
		last = node.tax_id;
	}

	fclose(fp);
	fclose(fpnew);

	return 0;
}

int ncbitc_create_names(char *filename)
{
	FILE *fp;
	FILE *fpnew;
	char line[NCBITC_LINE_SIZE];
	//int i, diff, last;
	char new[128];
	struct names_dmp name;
	struct names_dmp zero;
	int num = 0;

	memset(&zero, 0, sizeof(struct names_dmp));

	sprintf(new, "%s.bin", filename);

	fp = fopen(filename, "r");
	if (fp == NULL) {
		perror("fopen");
		return -1;
	}

	fpnew = fopen(new, "wb");
	if (fp == NULL) {
		perror("fopen");
		return -1;
	}

	fwrite(&num, sizeof(int), 1, fpnew);

	while (fgets(line, NCBITC_LINE_SIZE, fp) != NULL) {
		ncbitc_read_name_entry(line, &name);
		//ncbitc_print_name_entry(&name);
		fwrite(&name, sizeof(struct names_dmp), 1, fpnew);
		num++;
	}

	rewind(fpnew);
	fwrite(&num, sizeof(int), 1, fpnew);

	fclose(fp);
	fclose(fpnew);

	return 0;
}

#define NCBITC_OPER_SEARCH		0x0
#define NCBITC_OPER_SEARCH_GI	0x1
#define NCBITC_OPER_SEARCH_NODE	0x2
#define NCBITC_OPER_SEARCH_NAME	0x3
#define NCBITC_OPER_CREATE		0x4

void print_help()
{
	printf("Usage: ncbitc [options]\n");
	printf("Options:\n");
	printf("   -s --search id         search all tree using a gi index\n");
	printf("   -g --search-gi id      search a tax id using a gi index\n");
	printf("   -t --search-node id    search a node entry using a tax id\n");
	printf("   -n --search-name id    search a name entry using a tax id\n");
	printf("   -v --verbose           turn on verbose output\n");
	printf("   -h --help              print this help message\n");
	exit(0);
}

int main(int argc, char *argv[])
{
	int index;
	int tax_id;
	int ret;
	int c;
	int oper = -1;
	int option_index = 0;

	while (1) {
		static struct option long_options[] = {
			{"help", no_argument, 0, 'h'},
			{"verbose", no_argument, 0, 'v'},
			{"create", no_argument, 0, 'c'},
			{"search", required_argument, 0, 's'},
			{"search-gi", required_argument, 0, 'g'},
			{"search-name", required_argument, 0, 'n'},
			{"search-node", required_argument, 0, 't'},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, "hcs:n:vt:g:", long_options,
				&option_index);

		if (c == -1)
			break;

		switch (c) {
		case 'c':
			oper = NCBITC_OPER_CREATE;
			break;

		case 's':
			oper = NCBITC_OPER_SEARCH;
			index = atoi(optarg);
			break;

		case 'g':
			oper = NCBITC_OPER_SEARCH_GI;
			index = atoi(optarg);
			break;

		case 'v':
			verbose_flag = 1;
			break;

		case 'n':
			oper = NCBITC_OPER_SEARCH_NAME;
			tax_id = atoi(optarg);
			break;

		case 't':
			oper = NCBITC_OPER_SEARCH_NODE;
			tax_id = atoi(optarg);
			break;

		case '?':
		case 'h':
			print_help();
			break;

		default:
			print_help();
		}
	}

	if (verbose_flag)
		printf("verbose flag is set\n");

	if (oper == NCBITC_OPER_SEARCH) {
		struct nodes_dmp node;

		ret = ncbitc_search_tax_id(NCBITC_FILE_BIN_GI_TAXID, index, &tax_id);
		if (ret < 0) {
			printf("Error.\n");
			return -1;
		}

		if (tax_id == 0) {
			printf("0\n");
			return 0;
		}
		
		while (tax_id != 1) {
			ret = ncbitc_search_node(NCBITC_FILE_BIN_NODES, tax_id, &node);
			if (ret < 0) {
				printf("Error.\n");
				return -1;
			}

			tax_id = node.parent_tax_id;
			if (tax_id != 1)
				ncbitc_print_node_entry(&node);
		}

	} else if (oper == NCBITC_OPER_SEARCH_GI) {
		struct nodes_dmp node;

		ret = ncbitc_search_tax_id(NCBITC_FILE_BIN_GI_TAXID, index, &tax_id);
		if (ret < 0) {
			printf("Error.\n");
			return -1;
		}

		if (tax_id == 0) {
			printf("0\n");
			return 0;
		}

		ret = ncbitc_search_node(NCBITC_FILE_BIN_NODES, tax_id, &node);
		if (ret < 0) {
			printf("Error.\n");
			return -1;
		}

		ncbitc_print_node_entry(&node);

	} else if (oper == NCBITC_OPER_SEARCH_NODE) {
		struct nodes_dmp node;

		ret = ncbitc_search_node(NCBITC_FILE_BIN_NODES, tax_id, &node);
		if (ret < 0) {
			printf("Error.\n");
			return -1;
		}

		ncbitc_print_node_entry(&node);

	} else if (oper == NCBITC_OPER_SEARCH_NAME) {
		ret = ncbitc_search_name(NCBITC_FILE_BIN_NAMES, tax_id);
		if (ret < 0) {
			printf("Error.\n");
			return -1;
		}

	} else if (oper == NCBITC_OPER_CREATE) {
		ncbitc_create_gi_taxid(NCBITC_FILE_DMP_GI_TAXID);
		ncbitc_create_nodes(NCBITC_FILE_DMP_NODES);
		ncbitc_create_names(NCBITC_FILE_DMP_NAMES);
	} else {
		print_help();
	}

	return 0;
}
