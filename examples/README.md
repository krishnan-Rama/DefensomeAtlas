# examples

`metadata_35species.tsv` groups the 35-proteome aquatic/terrestrial set used to
develop the cross-lineage path. Columns:

| column | what it is |
|---|---|
| `phylum`, `class`, `order` | taxonomy; `class` is the usual covariate |
| `habitat` | where the animal lives as an adult |
| `life_stage_habitat` | where it lives as a larva, which is where exposure happens |
| `annotation` | BRAKER or Ensembl; a technical covariate, not biology |

**Why `life_stage_habitat` is a separate column.** An adult mayfly flies and an
adult toad sits on land, but both spend their formative life in water, and that
is the exposure that a defensome would have to answer. Grouping on adult
habitat alone would call *Cloeon* terrestrial.

**The contrast this set supports.** Aquatic versus terrestrial, with lineage as
a covariate, because three arthropods (*Chironomus*, *Cloeon*, *Parapoynx*) are
aquatic. Without them "aquatic" would just mean "vertebrate or mollusc" and the
comparison would be phylogeny wearing a different label.

**The contrast it does not support cleanly.** Within Lepidoptera the only
aquatic species, *Parapoynx*, is Ensembl-annotated while seven of the eight
terrestrial moths are BRAKER. *Lacanobia wlatinum* is the single
annotation-matched terrestrial moth, so the clean within-order, within-pipeline
comparison is one species against one. Treat any lepidopteran aquatic effect as
a hypothesis, not a result, and say so.

Nothing in this file is produced by the pipeline. It is a scientific claim you
are making about your species, and you should check it before using it.

## Sourcing proteomes for an endpoint dataset

Reference proteome catalogues do not all cover the same species, and
ecotoxicology endpoint species are mostly laboratory models:

- **Ensembl Rapid Release / Darwin Tree of Life** has broad wild-species
  coverage (1,576 proteomes in the July 2026 listing) but omits the classic
  models: no *Danio rerio*, no *Drosophila melanogaster*, no *Daphnia*.
- **UniProt** has the models. Of 17 species needed for one chlorpyrifos LC50
  set, Rapid Release had 0 and UniProt had 9, including the first crustaceans
  (*Daphnia magna*, *Hyalella azteca*).
- Neither had *Ctenopharyngodon idella*, *Pimephales promelas*, *Pungitius
  pungitius*, *Melanotaenia boesemani*, *Daphnia carinata*, *Eriocheir
  sinensis*, *Penaeus monodon* or *Procambarus clarkii*, for which UniProt
  lists only their viruses.

Check both, and check subspecies naming before concluding a species is absent:
*Cyprinus carpio carpio* and *Cyprinus carpio* are the same animal, and a
join on the exact string silently drops it.
