/*
-- Run the following commands in the sqlite3's shell:
CREATE TABLE entrez_string(entrez_id TEXT, string_id TEXT);
.separator ","
.import /Users/igormartire/Projects/tcc-code-py/data/entrez_id-string_id-old.csv entrez_string
-- make sure the file's rows are ending with LF and not with CR
.import /Users/igormartire/Projects/tcc-code-py/data/genage_models2.csv genage_models
.separator "\t"
.import /Users/igormartire/Projects/tcc-code-py/data/entrez_id-go_id.csv entrez_go
.separator " "
.import /Users/igormartire/Projects/tcc-code-py/data/raw_ppi/CE_6239.protein.links.detailed.v10.txt ppi_ce
.import /Users/igormartire/Projects/tcc-code-py/data/raw_ppi/DM_7227.protein.links.detailed.v10.txt ppi_dm
.import /Users/igormartire/Projects/tcc-code-py/data/raw_ppi/MM_10090.protein.links.detailed.v10.txt ppi_mm
.import /Users/igormartire/Projects/tcc-code-py/data/raw_ppi/SC_4932.protein.links.detailed.v10.txt ppi_sc
.read prepare_tables.sql
*/

-- gene_class: gene_id - class - org
-- gene_go: gene_id - go_id - evid - cat - org
-- gene_prot_ce: gene_id - prot_id - score
-- gene_prot_dm: gene_id - prot_id - score
-- gene_prot_mm: gene_id - prot_id - score
-- gene_prot_sc: gene_id - prot_id - score


CREATE TABLE gene_class AS
  SELECT string_id AS gene_id, "longevity influence" AS class, organism AS org
  FROM entrez_string AS a INNER JOIN genage_models AS b
    ON a.entrez_id = b."Entrez gene id";

UPDATE gene_class
  SET org = CASE org
              WHEN 'Caenorhabditis elegans' THEN 'ce'
              WHEN 'Drosophila melanogaster' THEN 'dm'
              WHEN 'Mus musculus' THEN 'mm'
              WHEN 'Saccharomyces cerevisiae' THEN 'sc'
            END,
      class = CASE class
                WHEN 'Anti-Longevity' THEN 'a'
                WHEN 'Pro-Longevity' THEN 'p'
                ELSE NULL
               END;

DELETE FROM gene_class
  WHERE class IS NULL;

CREATE TABLE gene_go AS
  SELECT gene_id, GO_ID AS go_id, Evidence AS evid, Category AS cat, org
  FROM entrez_string AS a INNER JOIN entrez_go AS b
    ON a.entrez_id = b.GeneID
  INNER JOIN gene_class AS c
    ON a.string_id = c.gene_id;

UPDATE gene_go
  SET cat = CASE cat
              WHEN 'Process' THEN 'bp'
              WHEN 'Function' THEN 'mf'
              WHEN 'Component' THEN 'cc'
            END;

CREATE TABLE gene_prot_ce AS
  SELECT protein1 AS gene_id, protein2 AS prot_id, combined_score AS score
  FROM ppi_ce
  WHERE protein1 IN (
    SELECT gene_id
    FROM gene_class
  );

CREATE TABLE gene_prot_dm AS
  SELECT protein1 AS gene_id, protein2 AS prot_id, combined_score AS score
  FROM ppi_dm
  WHERE protein1 IN (
    SELECT gene_id
    FROM gene_class
  );

CREATE TABLE gene_prot_mm AS
  SELECT protein1 AS gene_id, protein2 AS prot_id, combined_score AS score
  FROM ppi_mm
  WHERE protein1 IN (
    SELECT gene_id
    FROM gene_class
  );

CREATE TABLE gene_prot_sc AS
  SELECT protein1 AS gene_id, protein2 AS prot_id, combined_score AS score
  FROM ppi_sc
  WHERE protein1 IN (
    SELECT gene_id
    FROM gene_class
  );
