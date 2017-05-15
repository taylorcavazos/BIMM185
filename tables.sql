/* This file consists of the framework to create
** the following tables:
** genomes, replicons, genes, exons,synonyms, external references, and functions
** 
** This file can either be directly ran in mySQL or individual table structures
** can be copied and executed alone.      
*/

/* This table contains genome information such as id, name, size, ext
** genome_id: an automatic key used to link tables together
*/
CREATE TABLE genomes (
  genome_id        INT(10) UNSIGNED NOT NULL,
  name             VARCHAR(256) NOT NULL,
  tax_id           INT(10) UNSIGNED NOT NULL,
  domain           ENUM('bacteria','archaea','eukarya') NOT NULL,
  num_replicons    SMALLINT(5) UNSIGNED NOT NULL,
  num_genes        INT(10) UNSIGNED NOT NULL,
  size_bp          BIGINT(15) UNSIGNED NOT NULL,
  assembly         VARCHAR(25) NOT NULL,
  PRIMARY KEY (genome_id),
  KEY tax_id (tax_id)
) ENGINE=InnoDB;

/* This table contains the different replicons per genome
** replicons are linked to genomes using the genome id
*/
CREATE TABLE replicons (
  replicon_id    INT(10) UNSIGNED NOT NULL,
  genome_id      INT(10) UNSIGNED NOT NULL,
  name           VARCHAR(256) NOT NULL,
  type           ENUM('chromosome','plasmid') NOT NULL,
  shape          ENUM('circular','linear') NOT NULL,
  num_genes      INT(10) UNSIGNED NOT NULL,
  size_bp        BIGINT(15) UNSIGNED NOT NULL,
  accession      VARCHAR(25) NOT NULL,
  release_date   VARCHAR(25) NOT NULL,
  PRIMARY KEY (replicon_id),
  KEY(genome_id)
) ENGINE=InnoDB;

/* The genes table stores information about each gene in a genome*/
CREATE TABLE genes (
  gene_id     INT(10) UNSIGNED NOT NULL,
  genome_id   INT(10) UNSIGNED NOT NULL,
  replicon_id INT(10) UNSIGNED NOT NULL,
  locus_tag   CHAR(25) NOT NULL,
  protein_id  CHAR(25) NOT NULL,
  name        CHAR(10) NOT NULL,
  strand      ENUM('F','R') NOT NULL,
  num_exons   SMALLINT(5) UNSIGNED NOT NULL,
  length      MEDIUMINT(7) UNSIGNED NOT NULL,
  product     VARCHAR(1024) NOT NULL,
  PRIMARY KEY (gene_id),
  KEY (genome_id),
  KEY (replicon_id),
  KEY (locus_tag),
  KEY (protein_id)
) ENGINE=InnoDB;

/* The exon table holds information about each exon in a genome*/
CREATE TABLE exons(
  gene_id INT (10) UNSIGNED NOT NULL,
  exon VARCHAR (100) NOT NULL,
  left_pos INT (10) UNSIGNED NOT NULL,
  right_pos INT (10) UNSIGNED NOT NULL,
  length INT (10) UNSIGNED NOT NULL,
  KEY (gene_id)
) ENGINE=InnoDB;

/* This table holds synonyms for each gene*/
CREATE TABLE gene_synonyms (
  gene_id INT (10) UNSIGNED NOT NULL,
  synonyms VARCHAR (100) NOT NULL,
  KEY (gene_id) 
) ENGINE=InnoDB;

/* This table holds the external references for each gene*/
CREATE TABLE gene_xrefs (
  gene_id INT(10) UNSIGNED NOT NULL,
  xdb VARCHAR(32) NOT NULL,
  xid VARCHAR(24) NOT NULL,
  KEY (gene_id),
  KEY (xid)
) ENGINE=InnoDB;

/* This table holds the functions for each gene*/
CREATE TABLE functions(
  gene_id INT (10) UNSIGNED NOT NULL,
  function VARCHAR (100) NOT NULL,
  KEY (gene_id)
) ENGINE=InnoDB;
